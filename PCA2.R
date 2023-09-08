library(devtools)
install_github("vqv/ggbiplot")


#绘制PCA图进行分析，导入ggbiplot 包
library(ggbiplot)
library(PCAtools)
library(cowplot)
library(ggsci)
library(tidyverse)
library(xlsx)
library(showtext)
#install.packages("Cairo")
library(Cairo)

setwd("D:\\桌面\\all_for_jia")
#导入基因表达量图谱，数据在原文链接中,密码yaoy
sample_count <- read.xlsx("小红参pca-修.xlsx",1,row.names=1)
#sample_count <- sample_count[, -ncol(sample_count)]
head(sample_count)
#导入组名信息，数据在原文链接中,密码yaoy
sample_info <- read.xlsx(file = '组名信息-pca-修.xlsx',1,row.names=1)
#rownames(sample_info) <- sample_info[, 1]
sample_info
wine.pca <- prcomp(sample_count, scale. = TRUE)
wine.pca
# 获取解释方差贡献
variance_explained <- wine.pca$sdev^2
total_variance <- sum(variance_explained)
percent_variance_explained <- variance_explained / total_variance * 100

# 输出各主成分的解释方差百分比
for (i in seq_along(wine.pca$sdev)) {
  cat(paste0("PC", i, ": ", percent_variance_explained[i], "% variance explained\n"))
}

plot(wine.pca$rotation, main="after PCA")

wine.pca2 <- as.data.frame(wine.pca$rotation)


#加载tidyverse进行处理数据
library(tidyverse)
#首先进行行名的移入，列名用sample_name，然后再把PCA分析结果与组名信息相互结合
pca_rotated_plus <- rownames_to_column(wine.pca2, 
                                       var = 'sample_name') %>%
  left_join(rownames_to_column(sample_info, var = 'sample_name'), 
            by = 'sample_name')

#进行绘图，加载cowplot、ggsci包
#palette.color <- randomColor(count = 18)  #随机生成60种颜色，其实里面有重复的
#palette.color <- distinctColorPalette(18)
library(RColorBrewer)
colourCount = length(unique(sample_info$组名))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

#使用PC1和PC2映射到x、y轴，然后画点图，点大小为8，
#形状用strain进行映射，填充色用stage进行映射，
#映射形状的标度用21和22号点，填充色的标度用cowplot中的set2
#主题用cowplot中的theme_half_open
pca_rotated_plus <- pca_rotated_plus %>%
  mutate(组名 = factor(组名, levels = unique(组名) ))

library(dplyr)

library(ggsci)
ggplot(pca_rotated_plus, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, aes(color = 组名)) +
  ggforce::geom_mark_ellipse(aes(fill= 组名,color = 组名)) +
  labs(x = 'PC1 (34.6% variance explained)',
       y = 'PC2 (25.16% variance explained)') +
  scale_shape_manual(values = 19)  +theme_classic()+
  theme(
    axis.text.x = element_text(size = 10),  # 调整x轴标签的字体大小
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 10))
ggsave("小红参PCA2.png",width = 10,height = 6)
 