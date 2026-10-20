# data/

原始数据放这里，但不提交到仓库（见根目录 .gitignore）。

## 需要什么

GSE10072 的 series matrix 文件：

- 页面：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072
- 拉到页面底部 Series Matrix File(s)，下载 GSE10072_series_matrix.txt.gz
- 命令行：

    curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz

放好之后目录应该是这样：

    data/
    └── GSE10072_series_matrix.txt.gz

## 可选：基因注释

想按基因符号而不是探针号出结果，再下一份 GPL96 注释表（制表符分隔，至少两列：探针 ID、Gene Symbol），命名随意，跑的时候用 --annot 指过去即可。

## 为什么不上传

- 文件几十 MB，超过 GitHub 常规建议体积；
- 数据是公共资源，别人按上面地址也能下到；
- 仓库里只留「怎么拿到数据」的说明，比留一份副本更清楚。
