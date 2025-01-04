#合并多组GO，生成气泡图
#在进行气泡图分析之前，可以手动修改上一步输出的富集分析的csv文件中的GO/KEGG term的词条，把长词条变成短词条，以方便在气泡图上展示。
#从csv文件读取富集的信息到数据框（count_num可以指定获取的GO/KEGG term的数量，默认是15;enrich_type指定富集类型GO/KEGG）
#原始富集信息已经保存到本地，可根据需要自定义绘图
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(cowplot))

getdata_GO <- function(Prefix,enrich_type="GO",count_num=10){
trd_GO <- read.csv(paste(Prefix,"_",enrich_type,".csv",sep = ""),header = T)
pic_trd_GO <- head(trd_GO,count_num)
pic_trd_GO$adjust <- -log10(pic_trd_GO$p.adjust)
pic_trd_GO$type <- Prefix
return(pic_trd_GO)
}

# 检查命令行参数数量 
args <- commandArgs(trailingOnly = TRUE) 
if (length(args) < 2) { 
    stop("Usage: Rscript this_script.R deg_up_prefix deg_down_prefix") 
    } 
# 从命令行参数获取前缀 
deg_up <- args[1] 
deg_down <- args[2]

#读取GO
deg_up_GO <- getdata_GO(deg_up)
deg_down_GO <- getdata_GO(deg_down)
#合并5组GO数据
data_GO_all <- rbind(deg_up_GO,deg_down_GO)
#气泡图可视化富集结果
data_GO_all$Description_short <- strtrim(data_GO_all$Description, width = 50)
#截断注释40个字符
ggplot(data_GO_all,aes(type,Description_short)) +
geom_point(aes(fill=adjust,size=Count),alpha=0.9,pch=21) + #fill对应点的填充色，colour对应点的边框色
scale_fill_gradient(low='SpringGreen',high='DeepPink') + #设定颜色的变化范围
scale_size_area() + #设定点的大小比例和图例上显示的间隔
labs(y='GO term',fill='-log10(P.adjust)',size='Metabolites number')+
theme_bw()+
theme(axis.text.y = element_text(size = 8))
ggsave("diff_GO_all.pdf",dpi=300)
write.csv("data_GO_all",file="diff_GO_all.csv")

#合并多组KEGG数据,kegg富集

getdata_KEGG <- function(Prefix,enrich_type="KEGG",count_num=50){
trd_KEGG <- read.csv(paste(Prefix,"_",enrich_type,".csv",sep = ""),header = T)
pic_trd_KEGG <- head(trd_KEGG,count_num)
pic_trd_KEGG$adjust <- -log10(pic_trd_KEGG$p.adjust)
pic_trd_KEGG$type <- Prefix
return(pic_trd_KEGG)
}
#读取KEGG
deg_up_KEGG <- getdata_KEGG(deg_up)
deg_down_KEGG <- getdata_KEGG(deg_down)

#合并5组KEGG数据
data_kegg_all <- rbind(deg_up_KEGG,deg_down_KEGG)
#气泡图可视化富集结果
ggplot(data_kegg_all,aes(type,Description)) +
geom_point(aes(fill=adjust,size=Count),alpha=0.9,pch=21) + #fill对应点的填充色，colour对应点的边框色
scale_fill_gradient(low='SpringGreen',high='DeepPink') + #设定颜色的变化范围
scale_size_area() + #设定点的大小比例和图例上显示的间隔
labs(y='KEGG term',fill='-log10(P.adjust)',size='Metabolites number')+
theme_bw()+ theme(axis.text = element_text(size = 12))+ labs(x = "")
ggsave("diff_KEGG_all.pdf",dpi=1200)
ggsave("diff_KEGG_all.png",dpi=1200)
write.csv("data_kegg_all",file="diff_KEGG_all.csv")
