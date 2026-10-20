# 肺癌与癌旁组织的差异表达分析

## 1\. 这个仓库回答什么问题

同一份数据里，肺癌组织和癌旁正常组织之间，哪些基因的表达差异最大？差异够不够稳定，值不值得继续做富集分析和生存分析？

## 2. 目录

    .
    ├── README.md
    ├── requirements.txt
    ├── data/                  # 放原始数据（不进版本库）
    │   └── README.md          # 下载地址与文件命名
    ├── src/
    │   ├── geo_io.py          # 读 GEO series matrix / CSV，给样本分组
    │   ├── stats.py           # Welch t 检验 + BH 校正（纯 numpy）
    │   ├── plots.py           # 火山图、热图
    │   └── run_all.py         # 一键跑通
    ├── docs/
    │   └── 说明.md            # 三页说明：目的、方法、结果与局限
    └── results/               # 输出（不进版本库，跑完自动生成）
## 3\. 数据

* 数据集：GSE10072（Landi 等，2008），107 例样本：58 例肺腺癌、49 例正常肺组织。
* 下载页面：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
拉到底部 Series Matrix File(s)，下载 GSE10072\_series\_matrix.txt.gz。
* 命令行下载：

  curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072\_series\_matrix.txt.gz

* 放置：放进 data/，文件名保持原样，即 data/GSE10072\_series\_matrix.txt.gz。
* 大小：约 20 MB。原始数据不进仓库，.gitignore 已排除。
* 可选注释：想按基因符号看结果，再下载 GPL96 注释表（两列：探针 ID、Gene Symbol），跑的时候用 --annot 指过去；不提供也能跑，只是图上的行名是探针号。

## 4\. 怎么跑

需要 Python 3.9 及以上。

&#x20;   pip install -r requirements.txt

&#x20;   # 第一步：先用模拟数据确认环境和流程没问题（不联网、不需要真实数据）
    python src/run\\\_all.py --demo

    # 第二步：跑真实数据
    python src/run\\\_all.py --input data/GSE10072\\\_series\\\_matrix.txt.gz

    # 可选：带基因注释
    python src/run\\\_all.py --input data/GSE10072\\\_series\\\_matrix.txt.gz --annot data/GPL96\\\_annot.txt







参数都能改：

* \--fdr：显著性阈值，默认 0.05
* \--lfc：|log2FC| 阈值，默认 1.0（即 2 倍）
* \--top：热图画多少个基因，默认 30
* \--outdir：输出目录，默认 results

## 5\. 输出

* results/deg\_results.csv：每个探针一行，含 log2FC、t 值、自由度、p、FDR、上下调
* results/volcano.png：火山图，横轴 log2FC，纵轴 -log10 FDR
* results/heatmap.png：FDR 最小的基因，按行 z-score 画热图
* results/summary.txt：样本数、显著基因数、FDR 最小的 10 个基因

## 6\. 方法

1. 分组：从 !Sample\_title 和 !Sample\_source\_name\_ch1 里按关键词判断（tumor / cancer / adenocarcinoma 归肿瘤；normal / adjacent 归正常）。规则在 src/geo\_io.py 顶部的两个常量里，换数据集时改这里。
2. 预处理：剔除缺失超过 20% 的探针；如果数值最大值远大于 50，认为还是线性尺度，先做 log2。
3. 统计：每个探针对两组做 Welch t 检验（不假设方差相等），再对全部 p 值做 Benjamini-Hochberg 校正得到 FDR。
4. 阈值：FDR < 0.05 且 |log2FC| ≥ 1。
5. 不依赖 scipy：t 分布的双侧 p 值用不完全贝塔函数自己实现（src/stats.py），连分式在 \_betacf 里。好处是少装一个大包，代价是代码长一点——公式写在注释里。

## 7\. 局限

* 只有一个数据集，没有独立队列验证，也没做批次校正。
* 差异表达只说明两组之间相关，不说明因果。
* 探针未必全部注释到基因，同一基因可能有多个探针，目前没有合并。
* 没做通路富集和生存分析——这是下一步。

## 8\. 参考

* Landi MT, Dracheva T, Rotunno M, et al. Gene expression signature of cigarette smoking and its role in lung adenocarcinoma development and survival. PLoS One. 2008;3(2):e1651.
* Barrett T, et al. NCBI GEO: archive for functional genomics data sets. Nucleic Acids Research. 2013.

