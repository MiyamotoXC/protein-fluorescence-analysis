#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
真实数据下载：人源细胞双通道荧光图像（Hoechst 核 + CellMask 胞质）
==================================================================
数据源（真实公开数据集，CC BY 4.0，非合成）：
  Hugging Face  einarolafsson/cross-channel-cell-from-hoechst
      images/<field>.tif    Hoechst 通道（DNA / 细胞核染料）
  Hugging Face  einarolafsson/cross-channel-nuclei-from-cellmask
      images/<field>.tif    CellMask 通道（细胞质 / 细胞膜染料）
      masks/<field>.tif     uint16 核实例标注（0=背景，1..N=细胞核）

两个仓库来自同一批采集（CSA_screen，3029 个视野，1998x1998 uint16），
**视野名一一对应**，因此可以把 Hoechst 和 CellMask 两个真实荧光通道
配对成一张多通道 TIFF，用来做：
  - 蛋白表达强度（CellMask 通道）
  - 核 / 质定位比（Hoechst 与 CellMask 的相对分布）
  - 共定位（仅当数据提供第二个蛋白通道时才有意义，本数据集只有这两个通道，
    故 protein.yaml 里 protein_b 设为 null，共定位字段不输出）

输出 protein_data/images/<field>.tif，通道顺序与 protein.yaml 一致：
  通道 0 = 核（Hoechst），通道 1 = 蛋白 A（CellMask）

用法：
  python download_protein.py --list                 # 只看规模，不下载
  python download_protein.py                        # 下载 3 个视野
  python download_protein.py --max_fields 10
"""
import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "protein_data", "images")

NUC_REPO = "einarolafsson/cross-channel-cell-from-hoechst"        # Hoechst 通道
CYTO_REPO = "einarolafsson/cross-channel-nuclei-from-cellmask"    # CellMask 通道
FIELDS_CSV = "fields.csv"


def read_fields():
    """读 fields.csv（name, split, n_objects），返回全部视野名。"""
    from huggingface_hub import hf_hub_download
    path = hf_hub_download(CYTO_REPO, FIELDS_CSV, repo_type="dataset")
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append((r["name"], int(r["n_objects"])))
    return rows


def fetch(repo, repo_file, retries=4):
    """HF 走代理时连接会偶发被重置，这里做退避重试。"""
    from huggingface_hub import hf_hub_download
    last = None
    for _ in range(retries):
        try:
            return hf_hub_download(repo, repo_file, repo_type="dataset")
        except Exception as e:      # noqa: BLE001 - 网络错误类型不固定
            last = e
    raise RuntimeError(f"下载失败 {repo}/{repo_file}: {last}")


def download_one(name):
    """下载同一视野的两个通道，合成 (H, W, 2) uint16 TIFF。"""
    nuc = tifffile.imread(fetch(NUC_REPO, f"images/{name}.tif"))
    cyto = tifffile.imread(fetch(CYTO_REPO, f"images/{name}.tif"))
    stack = np.stack([nuc, cyto], axis=0)          # 通道在前，与 protein.yaml 一致
    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, name + ".tif")
    # 用 tifffile 而非 cv2：cv2.imwrite 在中文路径下会静默失败
    tifffile.imwrite(out, stack)
    return stack.shape


def main():
    ap = argparse.ArgumentParser(description="下载真实双通道荧光图像")
    ap.add_argument("--max_fields", type=int, default=3,
                    help="下载多少个视野（每个视野约 10 MB）")
    ap.add_argument("--workers", type=int, default=4, help="并发下载线程数")
    ap.add_argument("--list", action="store_true", help="只打印清单，不下载")
    args = ap.parse_args()

    fields = read_fields()
    print(f"数据源: {NUC_REPO} (Hoechst) + {CYTO_REPO} (CellMask)")
    print(f"可用视野: {len(fields)}，细胞核总数: {sum(n for _, n in fields)}")
    for name, n in fields[:5]:
        print(f"  {name}: {n} 个细胞核")

    if args.list:
        return

    # 优先取细胞核多的视野，样本少时也能得到足够的分析对象
    fields.sort(key=lambda x: -x[1])
    todo = fields[: args.max_fields] if args.max_fields > 0 else fields

    print(f"\n准备下载 {len(todo)} 个视野（并发 {args.workers}）")
    ok = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(download_one, n): n for n, _ in todo}
        for fu in as_completed(futs):
            try:
                shape = fu.result()
                print(f"  完成 {futs[fu]}  {shape}")
                ok += 1
            except Exception as e:      # noqa: BLE001
                print(f"  失败 {futs[fu]}: {e}")

    print(f"\n完成 {ok} 个视野 -> {OUT_DIR}")
    print("下一步：python segment_nuclei.py && python analyze_protein.py")


if __name__ == "__main__":
    import csv
    import numpy as np
    main()
