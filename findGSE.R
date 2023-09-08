library("findGSE");
#纯合基因组
findGSE(histo="../jellyfish/fastq.counts.histo", sizek=21, outdir="findGSE_21mer");
findGSE(histo="../jellyfish/fastq.counts.histo", sizek=31, outdir="findGSE_31mer");
#杂合基因组
findGSE(histo="../jellyfish/fastq.counts.histo", sizek=21, outdir="findGSE_21mer_exp_90", exp_hom=90);
findGSE(histo="../jellyfish/fastq.counts.histo", sizek=31, outdir="findGSE_31mer_exp_90", exp_hom=90);