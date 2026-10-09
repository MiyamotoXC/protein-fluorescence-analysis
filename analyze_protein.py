# -*- coding: utf-8 -*-
"""
蛋白质荧光图像分析：表达定量 + 亚细胞定位 + 双蛋白共定位
======================================================
在实例 mask（segment_nuclei.py 产出）的基础上，对荧光通道做逐细胞定量：

  1) 蛋白 A 表达强度：细胞区域内的 mean / std（已扣除背景基线）
  2) 亚细胞定位：核内 vs 胞质 平均强度比 localization_ratio（>1 偏核，<1 偏质）
  3) 双蛋白共定位：Pearson's r 与 Manders' M1（蛋白 A 落在蛋白 B 高表达区的比例）
  4) 形态：面积、偏心率、实心度、核面积

输入（路径与通道顺序见 protein.yaml）：
  protein_data/images/*.tif       多通道荧光
  protein_data/masks/*.tif        细胞实例 mask（0=背景，n=第 n 个细胞）
  protein_data/masks_nuc/*.tif    核实例 mask（编号与细胞 mask 一一对应）

用法：
  python analyze_protein.py
  python analyze_protein.py --out my_features.csv
输出：
  protein_features.csv            每细胞一行，可直接拿去做聚类/UMAP
"""
import argparse
import csv
import glob
import os

import numpy as np
import tifffile
import yaml
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "protein.yaml")


def load_config():
    with open(CONFIG, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_multichannel(path):
    """读取多通道 TIFF，返回 (C, H, W) float 数组。"""
    img = np.asarray(tifffile.imread(path))
    if img.ndim == 2:
        img = img[None, ...]
    elif img.ndim == 3 and img.shape[-1] <= 4:     # 通道在最后一维（如 HPA 的 4 通道图）
        img = np.moveaxis(img, -1, 0)
    return img.astype(np.float32)


def subtract_background(channel, percentile):
    """扣除背景基线（用整图低分位估计），负值截断到 0。"""
    bg = np.percentile(channel, percentile)
    return np.clip(channel.astype(np.float64) - bg, 0.0, None)


def pearson_r(a, b):
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def manders(m, p, thr_p):
    """Manders M1：在蛋白 P 高于阈值处，蛋白 M 的强度占比。"""
    mask_p = p > thr_p
    total = m.sum()
    if mask_p.sum() == 0 or total == 0:
        return 0.0
    return float(m[mask_p].sum() / total)


def analyze_one(image_path, mask_path, nuc_mask_path, cfg):
    img = load_multichannel(image_path)                       # (C, H, W)
    cells = np.asarray(tifffile.imread(mask_path)).astype(np.int32)
    nuclei = np.asarray(tifffile.imread(nuc_mask_path)).astype(np.int32)

    ch = cfg["channels"]
    bg_pct = float(cfg.get("background_percentile", 5))
    m_pct = float(cfg.get("manders_percentile", 50))
    pa_ch = int(ch["protein_a"])
    pb_ch = ch.get("protein_b")

    pa = subtract_background(img[pa_ch], bg_pct)
    pb = subtract_background(img[int(pb_ch)], bg_pct) if pb_ch is not None else None

    rows = []
    for rp in measure.regionprops(cells):
        cell = cells == rp.label
        nuc = (nuclei == rp.label) & cell          # 核 mask 与细胞 mask 编号一一对应
        cyto = cell & ~nuc

        nuc_int = float(pa[nuc].mean())
        cyto_int = float(pa[cyto].mean())
        row = {
            "image": os.path.basename(image_path),
            "label": int(rp.label),
            "area": int(rp.area),
            "eccentricity": float(rp.eccentricity),
            "solidity": float(rp.solidity),
            "nucleus_area": int(nuc.sum()),
            "proteinA_mean": float(pa[cell].mean()),
            "proteinA_std": float(pa[cell].std()),
            "nucleus_int": nuc_int,
            "cyto_int": cyto_int,
            "localization_ratio": float(nuc_int / cyto_int),
        }
        if pb is not None:
            a_vals = pa[cell]
            b_vals = pb[cell]
            row["pearson_AB"] = pearson_r(a_vals, b_vals)
            row["mander_M1"] = manders(a_vals, b_vals, np.percentile(b_vals, m_pct))
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description="蛋白质荧光逐细胞定量")
    ap.add_argument("--out", default=os.path.join(HERE, "protein_features.csv"),
                    help="输出 CSV 路径")
    ap.add_argument("--data_dir", default=None,
                    help="覆盖 protein.yaml 的数据目录（其下 images/ masks/ masks_nuc/）；"
                         "跑合成演示数据用 protein_data_synth")
    args = ap.parse_args()

    cfg = load_config()
    if args.data_dir:
        image_dir = os.path.join(args.data_dir, "images")
        mask_dir = os.path.join(args.data_dir, "masks")
        nuc_dir = os.path.join(args.data_dir, "masks_nuc")
    else:
        image_dir = os.path.join(HERE, cfg["image_dir"])
        mask_dir = os.path.join(HERE, cfg["mask_dir"])
        nuc_dir = os.path.join(HERE, cfg["nucleus_mask_dir"])

    files = sorted(glob.glob(os.path.join(image_dir, "*.tif")) +
                   glob.glob(os.path.join(image_dir, "*.tiff")))
    if not files:
        print(f"未找到图像于 {image_dir}。真实数据：python download_protein.py；"
              f"合成演示：python make_sample_data.py")
        return

    all_rows = []
    for imp in files:
        stem = os.path.splitext(os.path.basename(imp))[0]
        mp = os.path.join(mask_dir, stem + ".tif")
        np_ = os.path.join(nuc_dir, stem + ".tif")
        if not os.path.exists(mp) or not os.path.exists(np_):
            print(f"[skip] 缺少 mask: {stem}（请先运行 python segment_nuclei.py）")
            continue
        all_rows.extend(analyze_one(imp, mp, np_, cfg))

    if not all_rows:
        print("未找到可分析的图像/mask。")
        return

    # 字段取并集，避免个别图缺通道时丢列
    fields = []
    for r in all_rows:
        for k in r.keys():
            if k not in fields:
                fields.append(k)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)

    def avg(key):
        vals = [r[key] for r in all_rows if key in r]
        return sum(vals) / len(vals) if vals else float("nan")

    print(f"已写出 {len(all_rows)} 个细胞的特征 -> {args.out}")
    print(f"  蛋白A 平均表达 : {avg('proteinA_mean'):.1f}")
    print(f"  核/质定位比    : {avg('localization_ratio'):.2f}  (>1 偏核, <1 偏质)")
    if "pearson_AB" in fields:
        print(f"  Pearson(A,B)   : {avg('pearson_AB'):.3f}")
        print(f"  Manders M1     : {avg('mander_M1'):.3f}")


if __name__ == "__main__":
    main()
