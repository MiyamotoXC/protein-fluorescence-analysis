# -*- coding: utf-8 -*-
"""
蛋白质荧光图像分析 -- 合成样本数据生成器
========================================
生成多通道荧光 TIFF，使「分割 -> 定量 -> 定位 -> 共定位」整条链路在没有任何
外部数据集的情况下即可端到端跑通。

每个合成视野包含若干细胞，并刻意混入三类真实场景：
  - 蛋白 A 的亚细胞定位不同：偏核（nuclear）/ 偏质（cytoplasmic）/ 均匀（mixed）
  - 蛋白 B 与 A 的关系不同：共定位（coloc）或非共定位（独立分布）
这样下游算出的 localization_ratio、pearson_AB、mander_M1 会呈现可区分的分布，
聚类才有意义（而不是一堆噪声）。

产出：
  protein_data/images/field_XXX.tif    多通道 uint16 TIFF（通道顺序见 protein.yaml）

真实数据接入：把你自己的多通道免疫荧光 TIFF 直接放进 protein_data/images/ 即可，
后续命令完全一致（无需下载脚本）。

用法：
  python make_sample_data.py                      # 默认 6 个视野
  python make_sample_data.py --n_images 12 --n_cells 30
"""
import argparse
import os

import numpy as np
import tifffile
import yaml
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "protein.yaml")


def load_config():
    with open(CONFIG, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ellipse_mask(shape, center, axes, angle_deg):
    """生成椭圆布尔掩码（用参数方程，避免每个项目都依赖 OpenCV）。"""
    h, w = shape
    cy, cx = center
    a, b = axes
    t = np.deg2rad(angle_deg)
    yy, xx = np.mgrid[0:h, 0:w]
    dy, dx = yy - cy, xx - cx
    cos_t, sin_t = np.cos(t), np.sin(t)
    rx = (dx * cos_t + dy * sin_t) / a
    ry = (-dx * sin_t + dy * cos_t) / b
    return (rx ** 2 + ry ** 2) <= 1.0


def place_cells(shape, n_cells, rng):
    """随机放置互不重叠的细胞，并为每个细胞分配定位模式与共定位属性。"""
    h, w = shape
    placed = np.zeros(shape, dtype=bool)
    cells = []
    tries = 0
    while len(cells) < n_cells and tries < n_cells * 80:
        tries += 1
        a = int(rng.integers(22, 34))
        b = int(rng.integers(16, 28))
        if b > a:
            a, b = b, a
        angle = float(rng.integers(0, 180))
        cx = int(rng.integers(a + 4, w - a - 4))
        cy = int(rng.integers(a + 4, h - a - 4))
        cell_m = ellipse_mask(shape, (cy, cx), (a, b), angle)
        if (cell_m & placed).sum() > 0.15 * cell_m.sum():
            continue
        placed |= cell_m
        cells.append({
            "center": (cy, cx),
            "axes": (a, b),
            "angle": angle,
            "nuc_axes": (max(4, int(a * 0.45)), max(4, int(b * 0.45))),
            "mode": str(rng.choice(["nuclear", "cytoplasmic", "mixed"])),
            "coloc": bool(rng.random() < 0.6),
        })
    return cells


def render(cells, shape, rng):
    """渲染三个通道：0=核, 1=蛋白A, 2=蛋白B。返回 (3, H, W) uint16。"""
    background = 120.0
    nuc = np.full(shape, background)
    pa = np.full(shape, background)
    pb = np.full(shape, background)

    for c in cells:
        cell_m = ellipse_mask(shape, c["center"], c["axes"], c["angle"])
        nuc_m = ellipse_mask(shape, c["center"], c["nuc_axes"], c["angle"])
        cyto_m = cell_m & ~nuc_m

        nuc[nuc_m] += 2600 + rng.normal(0, 200, int(nuc_m.sum()))

        if c["mode"] == "nuclear":
            # 即便偏核的蛋白，胞质也有可检测的基础表达（真实荧光图像均如此）
            pa[nuc_m] += 2000
            pa[cyto_m] += 1000
        elif c["mode"] == "cytoplasmic":
            pa[nuc_m] += 450
            pa[cyto_m] += 1900
        else:
            pa[nuc_m] += 1300
            pa[cyto_m] += 1300

        if c["coloc"]:
            # 与 A 同分布 -> Pearson 高
            if c["mode"] == "nuclear":
                pb[nuc_m] += 2000
                pb[cyto_m] += 1000
            elif c["mode"] == "cytoplasmic":
                pb[nuc_m] += 450
                pb[cyto_m] += 1900
            else:
                pb[nuc_m] += 1300
                pb[cyto_m] += 1300
        else:
            # 只占胞质的一半区域 -> 与 A 部分错开，Pearson 低
            ys, xs = np.nonzero(cyto_m)
            half = xs > c["center"][1]
            sel = np.zeros(shape, dtype=bool)
            sel[ys[half], xs[half]] = True
            pb[sel] += 1800
            pb[nuc_m] += 300

    chans = []
    for ch in (nuc, pa, pb):
        ch = ndi.gaussian_filter(ch, 1.2)                  # 模拟显微 PSF
        ch = ch + rng.normal(0, 60, shape)                 # 读出噪声
        ch = np.clip(ch, 0, 65535).astype(np.uint16)
        chans.append(ch)
    return np.stack(chans, axis=0)


def main():
    ap = argparse.ArgumentParser(description="生成合成多通道荧光 TIFF")
    ap.add_argument("--n_images", type=int, default=6, help="视野数")
    ap.add_argument("--n_cells", type=int, default=24, help="每视野细胞数")
    ap.add_argument("--imgsz", type=int, default=512)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--data_dir", default="protein_data_synth",
                    help="合成数据输出目录（与真实数据 protein_data 隔离）")
    args = ap.parse_args()

    # 合成数据固定写到独立目录，避免和 download_protein.py 下载的真实荧光图混在一起
    image_dir = os.path.join(HERE, args.data_dir, "images")
    os.makedirs(image_dir, exist_ok=True)
    shape = (args.imgsz, args.imgsz)
    rng = np.random.default_rng(args.seed)

    for i in range(args.n_images):
        cells = place_cells(shape, args.n_cells, rng)
        arr = render(cells, shape, rng)
        path = os.path.join(image_dir, f"field_{i:03d}.tif")
        tifffile.imwrite(path, arr, photometric="minisblack")
        modes = {m: sum(1 for c in cells if c["mode"] == m) for m in
                 ("nuclear", "cytoplasmic", "mixed")}
        print(f"  field_{i:03d}.tif: {len(cells)} 个细胞 {modes}, "
              f"共定位 {sum(1 for c in cells if c['coloc'])} 个")

    print(f"\n已生成 {args.n_images} 个视野 -> {image_dir}")
    print("下一步（合成数据的共定位需要第 3 通道，记得把 protein.yaml 的 protein_b 改成 2）：")
    print(f"  python segment_nuclei.py --data_dir {args.data_dir} --out_dir outputs_synth")
    print(f"  python analyze_protein.py --data_dir {args.data_dir} --out_dir outputs_synth")


if __name__ == "__main__":
    main()
