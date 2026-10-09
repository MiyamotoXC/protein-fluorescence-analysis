# -*- coding: utf-8 -*-
"""
蛋白质荧光图像分析 -- 核/细胞实例分割
=====================================
从多通道荧光 TIFF 直接产出下游定量所需的**实例 mask**，无需训练任何模型：

  1) 核分割：核通道 Otsu 阈值 -> 形态学开运算去噪 -> 距离变换 + 局部极大值做种子
     -> 分水岭，把粘连的核拆开
  2) 细胞分割：以核为种子，在「细胞前景」（蛋白 A 通道 Otsu）内做分水岭，
     得到与核编号一一对应的细胞区域（胞质 = 细胞 - 核）

输出（路径见 protein.yaml）：
  protein_data/masks_nuc/*.tif   核实例 mask（uint16，0=背景，n=第 n 个细胞）
  protein_data/masks/*.tif       细胞实例 mask（编号与核 mask 一一对应）

用法：
  python segment_nuclei.py
  python segment_nuclei.py --min_distance 9 --min_size 60
"""
import argparse
import glob
import os

import numpy as np
import tifffile
import yaml
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu
from skimage.morphology import (binary_dilation, closing, disk, opening,
                                remove_small_objects)
from skimage.segmentation import watershed

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "protein.yaml")


def load_config():
    with open(CONFIG, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def segment_nuclei(nuc, min_distance, min_size):
    """核通道 -> 核实例 mask。"""
    smooth = ndi.gaussian_filter(nuc.astype(np.float32), 1.5)
    binary = smooth > threshold_otsu(smooth)
    binary = opening(binary, disk(2))
    binary = remove_small_objects(binary, max_size=min_size)
    if binary.sum() == 0:
        return None

    dist = ndi.distance_transform_edt(binary)
    coords = peak_local_max(dist, min_distance=min_distance, labels=binary)
    if len(coords) == 0:
        return None
    markers = np.zeros(dist.shape, dtype=np.int32)
    markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    return watershed(-dist, markers, mask=binary).astype(np.uint16)


def segment_cells(foreground, nuclei, min_size):
    """以核为种子，在细胞前景内做分水岭 -> 细胞实例 mask。"""
    smooth = ndi.gaussian_filter(foreground.astype(np.float32), 1.5)
    mask = smooth > threshold_otsu(smooth)
    # 核必然属于细胞：核区并入前景（蛋白 A 偏胞质时核信号低于阈值，会被误排除），
    # 再外扩一圈保证细胞区域严格大于核，胞质环非空
    mask = mask | binary_dilation(nuclei > 0, disk(2))
    mask = binary_dilation(mask, disk(3))
    mask = closing(mask, disk(4))
    # 注意：skimage 0.26 里 min_size 已废弃，新参数名 max_size 的实际行为是
    # 「移除此值以下的小对象」（等价旧 min_size），别被名字误导
    mask = remove_small_objects(mask, max_size=min_size)
    if mask.sum() == 0:
        return None
    dist = ndi.distance_transform_edt(mask)
    labels = watershed(-dist, nuclei, mask=mask).astype(np.uint16)
    # 保证每个细胞都带胞质环：小核的周围像素会被相邻大核在分水岭里抢走，
    # 这类细胞只剩核本身，下游的核/质定位比会除零变成 NaN。
    # 这里把它的核外 2 px 环从邻域里划回来（每个核最多让邻居让出 2 px，影响可忽略）。
    for lab in np.unique(nuclei):
        if lab == 0:
            continue
        nm = nuclei == lab
        if nm.any() and not (labels == lab)[~nm].any():
            labels[binary_dilation(nm, disk(2)) & ~nm] = lab
    return labels


def main():
    ap = argparse.ArgumentParser(description="核/细胞实例分割（无需训练）")
    ap.add_argument("--min_distance", type=int, default=7, help="核种子最小间距（像素）")
    ap.add_argument("--min_size", type=int, default=40, help="核最小面积（像素）")
    ap.add_argument("--data_dir", default=None,
                    help="覆盖 protein.yaml 的数据目录（其下 images/）；"
                         "跑合成演示数据用 protein_data_synth")
    ap.add_argument("--out_dir", default=None,
                    help="产出目录（默认 protein.yaml 的 output_dir）；"
                         "合成演示建议 outputs_synth，避免覆盖真实数据的结果")
    args = ap.parse_args()

    cfg = load_config()
    ch = cfg["channels"]
    image_dir = (os.path.join(args.data_dir, "images") if args.data_dir
                 else os.path.join(HERE, cfg["image_dir"]))
    out_root = os.path.join(HERE, args.out_dir or cfg["output_dir"])
    # 分割产出放在 outputs/segmentation/ 下，与定量结果分开
    mask_dir = os.path.join(out_root, "segmentation", "masks")
    nuc_dir = os.path.join(out_root, "segmentation", "masks_nuc")
    os.makedirs(mask_dir, exist_ok=True)
    os.makedirs(nuc_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(image_dir, "*.tif")) +
                   glob.glob(os.path.join(image_dir, "*.tiff")))
    if not files:
        print(f"未找到图像于 {image_dir}。真实数据：python download_protein.py；"
              f"合成演示：python make_sample_data.py")
        return

    for path in files:
        img = np.asarray(tifffile.imread(path)).astype(np.float32)
        if img.ndim == 2:
            img = img[None, ...]
        nuc_ch = int(ch["nucleus"])
        # 细胞前景优先用蛋白 A（细胞体整体有表达）
        fg_ch = int(ch["protein_a"]) if ch.get("protein_a") is not None else nuc_ch

        nuclei = segment_nuclei(img[nuc_ch], args.min_distance, args.min_size)
        if nuclei is None:
            print(f"[skip] {os.path.basename(path)}：核通道未分割出前景，请检查通道索引或调小 --min_size")
            continue
        cells = segment_cells(img[fg_ch], nuclei, args.min_size)
        if cells is None:
            print(f"[skip] {os.path.basename(path)}：细胞前景为空")
            continue

        name = os.path.splitext(os.path.basename(path))[0] + ".tif"
        tifffile.imwrite(os.path.join(nuc_dir, name), nuclei, photometric="minisblack")
        tifffile.imwrite(os.path.join(mask_dir, name), cells, photometric="minisblack")

        n = int(cells.max())
        print(f"  {os.path.basename(path)}: 分割出 {n} 个细胞 "
              f"(核 {int(nuclei.max())})")

    print(f"\n细胞 mask -> {mask_dir}")
    print(f"核 mask -> {nuc_dir}")
    print("下一步：python analyze_protein.py"
          + (f" --data_dir {args.data_dir}" if args.data_dir else "")
          + (f" --out_dir {args.out_dir}" if args.out_dir else ""))


if __name__ == "__main__":
    main()
