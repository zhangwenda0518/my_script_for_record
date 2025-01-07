# 获取命令行参数，如果没有提供参数，则使用默认文件名
args <- commandArgs(trailingOnly = TRUE)
# 参数说明
if (length(args) == 0) {
  cat("没有提供参数\n")
  cat("使用帮助\n")
  cat("  第一个参数（filename）\n")
  cat("  第二个参数（sample_name）\n")
  cat("  第三个参数（KEGG_list）:kegg注释文件 ，默认为 'KEGG.info'\n")
  q()
}
#提供参数，否则使用默认文件名
filename <- args[1]
sample <- args[2]
KEGG_list <- ifelse(length(args) > 2, args[3], "KEGG.info")

# 富集分析
suppressMessages(library(clusterProfiler))
suppressMessages(library(tidyverse))
# 可视化
suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(ggridges))

print("##############加载数据完成####################")

##获取基因的注释文件，需要输入基因列表和genetype(用于输出文件名的前缀)
get_gesa_kegg <- function(filename,sample,KEGG_list){
#定义需要分析的基因列表
data <- read.table(filename,header=TRUE,sep="\t")
data_sort <- data %>% arrange(desc(logFC))
gene_list <- data_sort$logFC
names(gene_list) <- rownames(data_sort)
sample_name <- sample
#读取KEGG信息
kegglist <- fread(KEGG_list,header = FALSE)
#GSEA-KEGG富集
kegg2gene <- data.frame(kegglist$V1,kegglist$V3) #keggid,Geneid
kegg2name <- data.frame(kegglist$V1,kegglist$V2) #keggid,kegg描述
enrichment_gsea <- GSEA(geneList = gene_list,
                        TERM2GENE = kegg2gene,
                        TERM2NAME = kegg2name,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",by = "fgsea")
#保存结果
write.csv(file = paste0(sample_name,".GSEA-kegg_results.csv"), x = enrichment_gsea@result)

#绘图
if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    #绘制山脊图
    p_ridgeplot <- ridgeplot(enrichment_gsea, core_enrichment= FALSE, fill="p.adjust", orderBy = "NES", showCategory=10) + ggtitle("Ridge plot for GSEA") 

    ggsave(filename = paste0(sample_name,".GSEA_kegg-ridgeplot.pdf"),                # EDIT THIS
           plot =  p_ridgeplot,  dpi = 1200)
    #绘制气泡图
    p_dotplot <- dotplot(enrichment_gsea,showCategory=10, label_format = 40,,font.size=8,title = "KEGG Enrich" )
    ggsave(filename = paste0(sample_name,".GSEA_kegg-dotplot.pdf"),                # EDIT THIS
           plot =  p_dotplot,  dpi = 1200)
}
}
# 调用函数并传参
get_gesa_kegg(filename, sample,KEGG_list)
print(运行结束)
