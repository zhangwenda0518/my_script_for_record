setwd("D://桌面")
dir()
# 读取数据
data <- read.table("repeat_stats_summary.tsv", header = TRUE)

# 根据Total列排序
data <- data[order(data$Size), ]

# 加载所需的包
library(ggplot2)
library(reshape2)
# 将数据转换为长格式，除去"Total"列
data_long <- reshape2::melt(data, id.vars = "File",
                            measure.vars = c("LTR", "TIR", "nonLTR", "nonTIR", "repeat.region"))
#data_long<- reshape2::melt(data, id.vars = "File",
#                            measure.vars = c("Size","LTR", "TIR", "nonLTR", "nonTIR", "repeat.region"))
# 定义自定义配色方案
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# 绘制水平堆叠条形图
p <- ggplot(data_long, aes(x = reorder(File, value), y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Genome Size", x = "Species", y = "Size (Bp)") +
  theme_minimal() +
  scale_fill_manual(values = my_colors,
                    labels = c("LTR", "TIR", "nonLTR", "nonTIR", "repeat_region"))+ 
  theme(legend.direction = "vertical", legend.position = "right") +
  coord_flip()


ggsave(file="repeats.bar.pdf",plot = p)
ggsave("repeats.bar.png", plot = p,width = 8,height = 8)

# 获取Size值
total_values <- subset(data_long2, variable == "Size")
total_values <- total_values[order(total_values$value), ]
# 添加Size的文本标签
p <- p + geom_text(data = total_values, aes(x = reorder(File, -value), y = value, label = paste0(round(value), " Mb")),
                   hjust = -0.2, size = 3,
                   nudge_x = 0.1)

# 显示绘图
print(p)

