setwd("I:\\科研工作进展\\1-灯心草\\1.文章\\重复序列注释\\Repeats")
library(ggplot2)
dir()
# 导入数据
df <- read.table("all.sum.tab", header = TRUE, sep = "\t")

# 获取列名（除去第一列"species"）
var_names <- names(df)[-1]

# 创建存储图表的列表
plots <- list()

# 遍历每一列进行分析和绘图
for (var in var_names) {
  # 创建线性回归图
  p <- ggplot(df, aes(x = genome_size, y = df[[var]])) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    geom_text(aes(label = paste0("R^2 = " , 
                                 round(cor(genome_size, df[[var]])^2, 2), 
                                 " , p = ", 
                                 formatC(summary(lm(df[[var]] ~ genome_size))$coef[2,4], 
                                         format = "e", digits = 2))), 
              x = min(df$genome_size), y = max(df[[var]]), hjust = 0, vjust = 1, 
              color = "blue") +  # 设置标签颜色为蓝色
    labs(title = paste0(var, " vs. Genome Size"),
         x = "Genome Size(bp)",
         y = var) +
    theme(plot.title.position = "plot")
  
  plots[[var]] <- p
  print(p)
  
  # 保存图表为图片文件
  ggsave(filename = paste0(var, "_plot.png"), plot = p,width = 8, height = 6)
}

# 显示图表列表
names(plots) <- var_names
print(plots)
