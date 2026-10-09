# 蛋白质荧光图像分析（单细胞分割项目的下游扩展）

> 单细胞分割的"逐细胞 mask + 裁剪"是公共底座，本项目是挂在其后的**路线 C**：
> 在荧光图里做蛋白质表达定量、亚细胞定位与双蛋白共定位分析。
> 上游分割项目见 `../细胞分析/`（路线 A 聚类、路线 B 时序、分割底座）。

## 它做什么

输入多通道荧光 TIFF（通道顺序在 `protein.yaml` 里配置，默认 0=核, 1=蛋白A, 2=蛋白B），
先做**无需训练的实例分割**（核通道 Otsu + 分水岭），再逐细胞输出：

1. **蛋白 A 表达强度**：细胞区域内 mean / std（已扣除背景基线）
2. **亚细胞定位**：核内 vs 胞质 强度比 `localization_ratio`（>1 偏核，<1 偏质）
3. **双蛋白共定位**：Pearson's r 与 Manders' M1
4. **形态**：面积、偏心率、实心度、核面积

## 快速开始（零下载）

仓库自带合成数据生成器，会刻意混入不同定位模式（偏核/偏质/均匀）与不同共定位关系，
让下游指标呈现可区分的分布：

```bash
pip install -r requirements.txt

python make_sample_data.py     # 生成 protein_data/images/*.tif（多通道荧光）
python segment_nuclei.py       # 生成 protein_data/masks/ 与 masks_nuc/（实例 mask）
python analyze_protein.py      # 逐细胞定量 -> protein_features.csv
```

实测（6 个视野 × 24 细胞 = 144 个细胞）跑出来的分布：

| 指标 | min | p25 | p50 | p75 | max |
|---|---|---|---|---|---|
| `localization_ratio` | 0.41 | 0.51 | 1.29 | 2.12 | 2.22 |
| `pearson_AB` | 0.01 | 0.36 | 0.98 | 0.99 | 0.99 |

定位比呈双峰（偏质 vs 偏核）、Pearson 明显分离（共定位 vs 各自独立），
说明这条链路真的在测量生物学信号，而不是输出一堆噪声。

## 用你自己的数据

不需要任何下载脚本：把多通道荧光 TIFF 直接放进 `protein_data/images/`，
按你的显微镜通道顺序改 `protein.yaml`，然后照常跑后两步：

```bash
python segment_nuclei.py
python analyze_protein.py
```

| 公开数据源 | 说明 |
|---|---|
| Human Protein Atlas 单细胞荧光 | 每张图 4 通道（蓝=核, 红=微管, 黄=ER, 绿=目标蛋白），研究蛋白质亚细胞定位的权威公开数据 https://www.proteinatlas.org/ |
| Kaggle `human-protein-atlas-image-classification` | 相同多通道图，需 Kaggle 凭证 |
| 任意 immunofluorescence 双/三通道 TIFF | 直接放进 `protein_data/images/` 即可 |

只需要共定位以外的分析时，把 `protein.yaml` 里 `channels.protein_b` 设为 `null`，
输出 CSV 会自动去掉 `pearson_AB` / `mander_M1` 两列。

## 分割是怎么做的（为什么不需要训练）

`segment_nuclei.py` 完全无监督，两步分水岭：

1. **核分割**：核通道 Otsu 阈值 → 形态学开运算去噪 → 距离变换 + 局部极大值做种子 → 分水岭，
   把粘连的核拆开
2. **细胞分割**：以核为种子，在"细胞前景"（蛋白 A 通道 Otsu）内做分水岭，
   得到与核编号**一一对应**的细胞区域（胞质 = 细胞 − 核）

两个关键细节（都是实测踩出来的）：

- 细胞前景必须**并入核区并外扩**。蛋白 A 偏胞质时核信号低于 Otsu 阈值，核会被整体排除出前景，
  导致"细胞区域反而比核小"、胞质为空、核/质比算出 NaN。
- 分割结果写入 `protein_data/masks/`（细胞）与 `protein_data/masks_nuc/`（核），编号一一对应，
  `analyze_protein.py` 直接用 `nuclei == label & cell` 取核区。

调参：`python segment_nuclei.py --min_distance 9 --min_size 60`

## 目录结构

```
蛋白质分析/
├── protein.yaml           # 通道索引 + 数据路径 + 分析参数（三个脚本共用）
├── make_sample_data.py    # 合成多通道荧光数据（零下载跑通）
├── segment_nuclei.py      # 核/细胞实例分割（无监督，无需训练）
├── analyze_protein.py     # 表达定量 / 定位比 / 共定位 -> protein_features.csv
└── requirements.txt

../细胞分析/               # 上游：分割底座 + 路线 A 聚类 + 路线 B 时序
```

## 输出字段

`protein_features.csv` 每细胞一行：

| 字段 | 含义 |
|---|---|
| `image`, `label` | 来源视野、细胞编号 |
| `area`, `eccentricity`, `solidity`, `nucleus_area` | 形态 |
| `proteinA_mean`, `proteinA_std` | 蛋白 A 表达强度（已扣背景） |
| `nucleus_int`, `cyto_int`, `localization_ratio` | 核内/胞质强度与定位比 |
| `pearson_AB`, `mander_M1` | 双蛋白共定位（`protein_b` 为 null 时不出现） |

后续可直接拿这份 CSV 做聚类/UMAP（复用 `../细胞分析/聚类分析/` 的思路）。

## 依赖

| 包 | 用途 |
|---|---|
| numpy, scipy | 数组运算、高斯滤波、距离变换 |
| scikit-image | Otsu 阈值、分水岭、形态学、regionprops |
| tifffile, imageio | 多通道 TIFF 读写 / 图像落盘 |
| pyyaml | 读取 protein.yaml |

> 图像落盘统一用 `tifffile` / `imageio`：**`cv2.imwrite` 在中文路径下会返回 True 却不写文件**。
