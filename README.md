# 蛋白质荧光图像分析（单细胞分割项目的下游扩展）

> 单细胞分割的"逐细胞 mask + 裁剪"是公共底座，本项目是挂在其后的**路线 C**：
> 在荧光图里做蛋白质表达定量、亚细胞定位与双蛋白共定位分析。
> 上游分割项目见 `../细胞分析/`（路线 A 聚类、路线 B 时序、分割底座）。

## 它做什么

输入多通道荧光 TIFF（通道顺序在 `protein.yaml` 里配置），先做**无需训练的实例分割**
（核通道 Otsu + 分水岭），再逐细胞输出：

1. **蛋白表达强度**：细胞区域内 mean / std（已扣除背景基线）
2. **亚细胞定位**：核内 vs 胞质 强度比 `localization_ratio`（>1 偏核，<1 偏质）
3. **双蛋白共定位**：Pearson's r 与 Manders' M1（需要第二个蛋白通道）
4. **形态**：面积、偏心率、实心度、核面积

## 快速开始：真实数据

默认走**真实公开数据集**，不需要任何合成数据：

```bash
pip install -r requirements.txt

python download_protein.py        # 下载真实双通道荧光（Hoechst + CellMask）
python segment_nuclei.py          # 核/细胞实例分割 -> protein_data/masks/
python analyze_protein.py         # 逐细胞定量 -> protein_features.csv
```

**真实数据源**（两个 Hugging Face 数据集，CC BY 4.0，来自同一批采集，视野名一一对应）：

| 用途 | 数据集 | 内容 |
|---|---|---|
| 核通道 | `einarolafsson/cross-channel-cell-from-hoechst` | Hoechst（DNA）通道，1998×1998 uint16 |
| 胞质通道 | `einarolafsson/cross-channel-nuclei-from-cellmask` | CellMask 通道 + 核实例标注 |

共 3029 个视野、264,503 个细胞核；`download_protein.py` 默认取 3 个视野
（每个约 10 MB，可用 `--max_fields` 调整，用 `--list` 只看清单）。

本机实测（3 个真实视野）：分割出 1784 个细胞，其中 `A02_1_1` 视野分出 742 个细胞、
数据集自带标注为 683 个细胞核——分割与真值同量级。定量结果：

| 指标 | min | p50 | max |
|---|---|---|---|
| `localization_ratio`（CellMask 核/质比） | 0.38 | 1.55 | 4.12 |

关于**共定位**：这份公开数据只有 Hoechst 与 CellMask 两个通道，没有第二个蛋白通道，
所以 `protein.yaml` 里 `protein_b` 为 `null`，CSV 不输出 `pearson_AB` / `mander_M1`。
要跑共定位，用合成演示数据（3 通道）或换成你自己的双蛋白荧光图。

## 用你自己的数据

不需要任何下载脚本：把多通道荧光 TIFF 放进 `protein_data/images/`，
按显微镜的通道顺序改 `protein.yaml`，然后照常跑后两步：

```bash
python segment_nuclei.py
python analyze_protein.py
```

| 公开数据源 | 说明 |
|---|---|
| Human Protein Atlas 单细胞荧光 | 每张图 4 通道（蓝=核, 红=微管, 黄=ER, 绿=目标蛋白），研究蛋白质亚细胞定位的权威公开数据 https://www.proteinatlas.org/ |
| Kaggle `human-protein-atlas-image-classification` | 与 HPA 相同的多通道图，需 Kaggle 凭证 |
| 任意 immunofluorescence 双/三通道 TIFF | 直接放进 `protein_data/images/` 即可 |

## 合成演示数据（零下载，含共定位）

仓库自带生成器，会刻意混入不同定位模式（偏核/偏质/均匀）与不同共定位关系，
**3 个通道**，因此可以用来验证共定位指标。它写到独立目录 `protein_data_synth/`，
不会和真实数据混在一起：

```bash
python make_sample_data.py                                  # -> protein_data_synth/
# 把 protein.yaml 的 channels.protein_b 改成 2 以启用共定位
python segment_nuclei.py --data_dir protein_data_synth
python analyze_protein.py --data_dir protein_data_synth
```

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
├── download_protein.py    # 下载真实双通道荧光（Hoechst + CellMask）
├── make_sample_data.py    # 合成多通道荧光（零下载，含共定位）-> protein_data_synth/
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
| huggingface_hub | 真实数据集下载 |

> 图像落盘统一用 `tifffile` / `imageio`：**`cv2.imwrite` 在中文路径下会返回 True 却不写文件**。
