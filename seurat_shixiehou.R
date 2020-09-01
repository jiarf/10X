####将数据整合为表达矩阵####
rm(list=ls())
setwd("/data/jiarf/10X/results/outs/filtered_feature_bc_matrix")
library(Matrix)
cellbarcodes <- read.table("./barcodes.tsv.gz")
genenames <- read.table("./features.tsv.gz")
molecules <- Matrix::readMM("./matrix.mtx.gz")

rownames(molecules) <- genenames[,2]
colnames(molecules) <- cellbarcodes[,1]
#expression_matrix_df 是该样品的表达矩阵,行是gene，列是cell_barcode
expression_matrix_df <- as.data.frame(as.matrix(molecules))
expression_matrix_df_cut <- expression_matrix_df[(rowSums(expression_matrix_df)>0),] 

#############seurat单细胞注释分析###################
library(Seurat)
packageVersion("Seurat")
library(dplyr)
library(patchwork)
##构建seurat对象#####
seurat <- CreateSeuratObject(counts = expression_matrix_df_cut, project = "SSC", min.cells = 3, min.features = 200)
seurat
# An object of class Seurat 
# 17222 features across 8682 samples within 1 assay 
# Active assay: RNA (17222 features, 0 variable features)#22104x3894
# save(seurat,file = "./shixiehou_seurat_out/Rdata/seurat_Rdata.rds")
# load(file=("./shixiehou_seurat_out/Rdata/seurat_Rdata.rds"))#重新读入


expression_matrix_df_cut[c("PGBD2", "ZNF692", "ZNF672"), 1:30]
#           AAACCCAAGGGAGTGG-1 AAACCCAAGGTGCGAT-1 AAACCCACAAAGGCAC-1 AAACCCACACGAGGTA-1 AAACCCAGTCAACATC-1 AAACCCAGTCACCGAC-1
# PGBD2                   0                  0                  0                  0                  0                  0
# ZNF692                  0                  0                  0                  0                  0                  0
# ZNF672                  0                  0                  0                  0                  0                  0
#           AAACCCAGTGAACCGA-1 AAACCCAGTGCACATT-1 AAACCCATCACGGGAA-1 AAACCCATCATTTCCA-1 AAACCCATCGCAGTTA-1 AAACCCATCTAGACAC-1
# PGBD2                   0                  0                  0                  0                  0                  0
# ZNF692                  0                  0                  0                  0                  0                  0
# ZNF672                  0                  0                  0                  0                  0                  0
#           AAACGAAAGATACTGA-1 AAACGAAAGGAACATT-1 AAACGAAAGGTTGACG-1 AAACGAACAATTGTGC-1 AAACGAACAGAGCTAG-1 AAACGAAGTACAGGTG-1
# PGBD2                   0                  1                  0                  2                  1                  0
# ZNF692                  0                  0                  0                  0                  0                  0
# ZNF672                  0                  0                  0                  0                  0                  0
#           AAACGAAGTAGAGATT-1 AAACGCTAGCATGCGA-1 AAACGCTAGGACAGTC-1 AAACGCTAGTTAGTAG-1 AAACGCTCAAACCATC-1 AAACGCTCAATGCTCA-1
# PGBD2                   0                  0                  0                  0                  0                  0
# ZNF692                  0                  0                  0                  0                  0                  0
# ZNF672                  0                  0                  0                  0                  0                  0
#          AAACGCTCAGTGACCC-1 AAACGCTGTGCGGTAA-1 AAACGCTTCCACGGGT-1 AAACGCTTCTGAATCG-1 AAAGAACAGGACAGCT-1 AAAGAACAGGATTCAA-1
# PGBD2                   0                  0                  0                  0                  0                  0
# ZNF692                  0                  0                  0                  0                  0                  0
# ZNF672                  0                  0                  0                  0                  0                  0

#在Seurat中可以使用PercentageFeatureSet函数计算每个细胞中线粒体的含量：在人类参考基因中线粒体基因是以“MT-”开头的，而在小鼠中是以“mt-”开头的。
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")#没意义，因为 All cells have the same value of percent.mt
# Show QC metrics for the first 5 cells
head(seurat@meta.data, 5)
# orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGGGAGTGG-1        SSC       1429          695          0
# AAACCCAAGGTGCGAT-1        SSC        771          471          0
# AAACCCACAAAGGCAC-1        SSC        808          515          0
# AAACCCACACGAGGTA-1        SSC       8775         1663          0
# AAACCCAGTCAACATC-1        SSC       8186          932          0
# Visualize QC metrics as a violin plot
p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf('./shixiehou_seurat_out/SSC.vln1.pdf', width = 16, height = 20)
print(p1)
dev.off()
# All cells have the same value of percent.mt.
###我们还可以使用FeatureScatter函数来对不同特征-特征之间的关系进行可视化#######
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pdf('./shixiehou_seurat_out/SSC.featurescatter.pdf',width = 16, height = 8)
print(plot1 + plot2)
dev.off()

#######根据QC指标进行细胞和基因的过滤###########
# 可视化QC指标，并用它们来过滤细胞：
# 1）将unique基因count数超过2500，或者小于200的细胞过滤掉
# 2）把线粒体含量超过5%以上的细胞过滤掉
SSC_seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
saveRDS(SSC_seurat,file = './shixiehou_seurat_out/Rdata/SSC_seurat_subset.Rdata')
#######数据的归一化#########
# 默认情况下，Seurat使用global-scaling的归一化方法，称为“LogNormalize”，这种方法是利用总的表达量对每个细胞里的基因表达值进行归一化，乘以一个scale factor（默认值是10000），再用log转换一下。归一化后的数据存放在pbmc[["RNA"]]@data里
SSC_seurat_scale <- NormalizeData(SSC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#########鉴定高可变基因（特征选择）###########
# Seurat使用FindVariableFeatures函数鉴定高可变基因，这些基因在PBMC不同细胞之间的表达量差异很大（在一些细胞中高表达，在另一些细胞中低表达）。默认情况下，会返回2,000个高可变基因用于下游的分析，如PCA等
SSC_seurat <- FindVariableFeatures(SSC_seurat_scale, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SSC_seurat), 10)
# plot variable features with and without labels
p1 <- VariableFeaturePlot(SSC_seurat)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
p1 + p2
pdf('./shixiehou_seurat_out/SSC.vln2.pdf', width = 16, height = 8)
print(p1)
dev.off()

pdf('./shixiehou_seurat_out/SSC.labelpoints.pdf', width = 16, height = 8)
print(p2)
dev.off()
#######################数据的标准化##############
# Seurat使用ScaleData函数对归一化后的count矩阵进行一个线性的变换(“scaling”)，将数据进行标准化：
# 1）shifting每个基因的表达，使细胞间的平均表达为0
# 2）scaling每个基因的表达，使细胞间的差异为1
# ScaleData默认对之前鉴定到的2000个高可变基因进行标准化，也可以通过vars.to.regress参数指定其他的变量对数据进行标准化，表达矩阵进行scaling后，其结果存储在pbmc[["RNA"]]@scale.data中。


SSC_seurat <- ScaleData(SSC_seurat)
all.genes <- rownames(SSC_seurat)
SSC_seurat <- ScaleData(SSC_seurat, features = all.genes)
SSC_seurat <- ScaleData(SSC_seurat, vars.to.regress = "percent.mt")


#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
########################进行PCA线性降维#####################
# Seurat使用RunPCA函数对标准化后的表达矩阵进行PCA降维处理，默认情况下，只对前面选出的2000个高可变基因进行线性降维，也可以通过feature参数指定想要降维的数据集。
SSC_seurat <- RunPCA(SSC_seurat, features = VariableFeatures(object = SSC_seurat))
print(SSC_seurat[["pca"]], dims = 1:5, nfeatures = 5)
save(SSC_seurat,file = "./shixiehou_seurat_out/Rdata/SSC.PCA.Rdata")
# PC_ 1 
# Positive:  HSP90AA1, PTGES3, ENSMFAG00000000497, FAM216A, ACTG1 
# Negative:  HMGB4, ENSMFAG00000042551, TNP2, RNF151, ENSMFAG00000035848 
# PC_ 2 
# Positive:  IZUMO4, DYDC1, SPINK2, LDHC, KCNK4 
# Negative:  RPL10, CST3, ITGB1, CD63, APOE 
# PC_ 3 
# Positive:  IGFBP7, DCN, MGP, TIMP3, MYH11 
# Negative:  ENSMFAG00000001879, ENSMFAG00000034362, ENSMFAG00000002462, HERC5, PAGE1 
# PC_ 4 
# Positive:  COTL1, CLDN11, ENSMFAG00000032810, SLC7A5, INHA 
# Negative:  GML, IQCB1, ENSMFAG00000031138, ENSMFAG00000045427, MEIOB 
# PC_ 5 
# Positive:  ENSMFAG00000037049, ENSMFAG00000034362, ENSMFAG00000001879, PLK1, CENPF 
# Negative:  ENSMFAG00000032810, CCK, INHA, CLDN11, ENSMFAG00000032762
########Seurat可以使用VizDimReduction, DimPlot, 和DimHeatmap函数对PCA的结果进行可视化#####
p3 <- VizDimLoadings(SSC_seurat, dims = 1, reduction = "pca")
p4 <- VizDimLoadings(SSC_seurat, dims = 2, reduction = "pca")
p5 <- DimPlot(SSC_seurat, reduction = "pca")#每个点是一个细胞
p3 +p4 
pdf('./shixiehou_seurat_out/SSC.vln3.pdf', width = 16, height = 9)
print(p3+p4)
dev.off()
pdf('./shixiehou_seurat_out/SSC.dimplot1.pdf', width = 16, height = 9)#PCA得分图能将对照组和实验组样本区分开。在PCA图中，如果样本之间聚集在一起，说明这些样本差异性小；反之样本之间距离越远，说明样本之间差异性越大。
print(p5)
dev.off()


pdf('./shixiehou_seurat_out/SSC.dimheatmap1.pdf', width = 16, height = 9)
DimHeatmap(SSC_seurat, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf('./shixiehou_seurat_out/SSC.dimheatmap2.pdf', width = 16, height = 9)
DimHeatmap(SSC_seurat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

####选择PCA降维的维数用于后续的分析######
# Seurat可以使用两种方法确定PCA降维的维数用于后续的聚类分析##########
# 使用JackStrawPlot函数
# 使用JackStraw函数计算每个PC的P值的分布，显著的PC会有较低的p-value
SSC_seurat <- JackStraw(SSC_seurat, num.replicate = 100)##时间很久
SSC_seurat <- ScoreJackStraw(SSC_seurat, dims = 1:20)
# 使用JackStrawPlot函数进行可视化
pdf('./shixiehou_seurat_out/SSC.jackstrawplot.pdf', width = 16, height = 8)
JackStrawPlot(SSC_seurat, dims = 1:15)
dev.off()
#使用ElbowPlot函数
# 使用ElbowPlot函数查看在哪一个PC处出现平滑的挂点：
pdf('./shixiehou_seurat_out/SSC.elbowplot.pdf', width = 16, height = 8)
ElbowPlot(SSC_seurat)
dev.off()
####细胞的聚类分群jiangwei keshihua######
# Seurat使用图聚类的方法对降维后的表达数据进行聚类分群。KNN
SSC_seurat <- FindNeighbors(SSC_seurat, dims = 1:10)
SSC_seurat <- FindClusters(SSC_seurat, resolution = 0.5)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 7575
# Number of edges: 237907
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.9125
#    Number of communities: 19
#    Elapsed time: 0 seconds
head(Idents(SSC_seurat), 5)
# AAACCCAAGGGAGTGG-1 AAACCCACAAAGGCAC-1 AAACCCACACGAGGTA-1 AAACCCAGTCAACATC-1 AAACCCAGTCACCGAC-1 
# 3                  1                 14                  9                  1 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
#非线性降维可视化（UMAP/tSNE）
# Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space.
# UMAP降维可视化
seurat_runUMAP <- RunUMAP(SSC_seurat, dims = 1:15)
pdf('./shixiehou_seurat_out/SSC.UMAP.pdf', width = 16, height = 8)
DimPlot(seurat_runUMAP, reduction = "umap")
dev.off()
#tSNE降维可视化
seurat_RunTSNE <- RunTSNE(SSC_seurat, dims = 1:15)
pdf('./shixiehou_seurat_out/SSC.TSNE.pdf', width = 16, height = 8)
DimPlot(seurat_RunTSNE, reduction = "tsne", label = TRUE)
dev.off()
saveRDS(seurat_RunTSNE,file='./shixiehou_seurat_out/Rdata/seurat_RunTSNE.Rdata')
###########鉴定不同类群之间的差异表达基因###########
# Seurat使用FindMarkers和FindAllMarkers函数进行差异表达基因的筛选，可以通过test.use参数指定不同的差异表达基因筛选的方法。
# find all markers of cluster 1
cluster1.markers <- FindMarkers(seurat_RunTSNE, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# p_val
# PRM2               0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001637327
# ENSMFAG00000040785 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000026294644419117482778978284282645058
# TPPP2              0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000087415613244988002004988769241033615739731492050620049
# C3orf22            0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000016083546723470532063389054955571984616758741145906608266
# RPS8               0.00000000000000000000000000000000000000000000000000000000000000000000000000000000004266991241471333295914784657573187928714365332604707075171860811702
# avg_logFC pct.1 pct.2
# PRM2                0.5828498 1.000 1.000
# ENSMFAG00000040785  0.3989413 1.000 1.000
# TPPP2               0.8850507 0.881 0.886
# C3orf22             0.8942210 0.917 0.925
# RPS8               -0.6145655 0.137 0.552
# p_val_adj
# PRM2               0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002819805
# ENSMFAG00000040785 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000045284636618604125712301246666738419
# TPPP2              0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000150547169130518340400313602924565641393755850589151700
# C3orf22            0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000027699084167160949372704811022673555566983708850985055628
# RPS8               0.000000000000000000000000000000000000000000000000000000000000000000000000000000734861231606193079094524373981648389067994376854040519167516354328
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_RunTSNE, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# > head(cluster5.markers, n = 5)
# p_val avg_logFC pct.1 pct.2 p_val_adj
# ENSMFAG00000042475     0  1.655843 0.994 0.122         0
# ENSMFAG00000031138     0  1.393691 0.968 0.065         0
# ENSMFAG00000037178     0  1.340153 0.960 0.050         0
# ENSMFAG00000045427     0  1.171999 0.945 0.054         0
# TERF2IP                0  1.150245 0.955 0.050         0
# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.markers <- FindAllMarkers(seurat_RunTSNE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#慢
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# # A tibble: 38 x 7
# # Groups:   cluster [19]
# p_val         avg_logFC pct.1 pct.2 p_val_adj cluster gene   
# <dbl>            <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
#   1 0.            1.61  0.89  0.52  0.        0       SPACA1 
# 2 0.            1.58  0.885 0.466 0.        0       TMEM210
# 3 1.60e- 45     0.970 0.659 0.659 2.75e- 41 1       DCUN1D1
# 4 3.24e- 38     1.25  0.653 0.686 5.58e- 34 1       TSPAN6 
# 5 1.16e-303     2.41  0.735 0.257 1.99e-299 2       CCDC168
# 6 6.74e-211     1.67  0.744 0.347 1.16e-206 2       PRSS37 
# 7 1.56e- 71     0.390 1     1     2.69e- 67 3       PRM2   
# 8 1.15e- 45     0.361 0.925 0.882 1.97e- 41 3       TPPP2  
# 9 0.            1.27  0.985 0.243 0.        4       AURKA  
# 10 1.79e-291     1.32  1     0.344 3.08e-287 4       CCNA1  
# … with 28 more rows
table(seurat.markers$cluster)

# 0   1   2   3   4   5   6   7   8   9  10  11  12 
# 9 377 432 114 743 722 441 364 297 495  70 208  10
# write.table(seurat.markers,file='./seurat.markers.table')
#      p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene    
#      <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
# 1 5.81e- 54     0.571 0.982 0.973 1.00e- 49 0       COX3    
# 2 2.28e- 47     0.716 0.858 0.833 3.93e- 43 0       ND3     
# 3 0.            1.57  0.891 0.352 0.        1       SPACA1  
# 4 0.            1.45  0.716 0.165 0.        1       ACRV1   
# 5 0.            2.49  0.753 0.166 0.        2       CCDC168 
# 6 4.25e-223     1.64  0.758 0.313 7.31e-219 2       PRSS37  
# 7 1.03e-114     1.07  0.997 0.992 1.78e-110 3       CRISP2  
# 8 1.41e- 65     1.33  0.722 0.61  2.43e- 61 3       TSPAN6  
# 9 6.38e-149     1.21  0.962 0.536 1.10e-144 4       HSP90AA1
# 10 2.13e-127     1.05  0.905 0.514 3.68e-123 4       LYAR    
# … with 16 more rows

save(seurat,seurat.markers,file = "./shixiehou_seurat_out/Rdata/SSC.marker.Rdata")
#pbmc <- readRDS(file = "../data/pbmc3k_final.rds")加载数据


# 0   1   2   3   4   5   6   7   8   9  10  11  12 
# 0  94 104  13 222 142 111 301 145 437   9  50   6 

top.gene <- seurat.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)#每个簇里排序后的第一个

VlnPlot(seurat, features = top.gene$gene)
#####Marker基因的可视化##########
# Seurat可以使用VlnPlot，FeaturePlot，RidgePlot，DotPlot和DoHeatmap等函数对marker基因的表达进行可视化

p6 <- VlnPlot(seurat_RunTSNE, features = c("LYAR", "SPACA1"))
p6
# you can plot raw counts as well
p7 <- VlnPlot(seurat_RunTSNE, features = c("LYAR", "SPACA1"), slot = "counts", log = TRUE)

p8 <- FeaturePlot(seurat_RunTSNE, features = top.gene$gene)
p9 <- RidgePlot(seurat_RunTSNE, features = top.gene$gene)
p10 <- DotPlot(seurat_RunTSNE, features = c("SPACA1","TSPAN6","CCDC168","PRM2","CCNA1","ENSMFAG00000003237","DCN","COX1","RAD51AP2","TNP2","CETN3","ENSMFAG00000039259","ENSMFAG00000001879","FAM24A","ENSMFAG00000042652","CST9L","CLU"))
top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
p11 <- DoHeatmap(seurat_RunTSNE, features = top10$gene) #+ NoLegend()

pdf('./shixiehou_seurat_out/SSC.vlnplot4.pdf', width = 16, height = 9)
print(p6)
dev.off()

pdf('./shixiehou_seurat_out/SSC.vlnplot5.pdf', width = 16, height = 9)
print(p7)
dev.off()

pdf('./shixiehou_seurat_out/SSC.featureplot.pdf', width = 16, height = 9)
print(p8)
dev.off()

pdf('./shixiehou_seurat_out/SSC.RidgePlot.pdf', width = 16, height = 9)## 峰峦图（RidgePlot）可视化marker基因的表达
print(p9)
dev.off()

pdf('./shixiehou_seurat_out/SSC.DotPlot.pdf', width = 40, height = 9)
print(p10)
dev.off()

pdf('./shixiehou_seurat_out/SSC.DoHeatmap.pdf', width = 16, height = 20)
print(p11)
dev.off()
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
