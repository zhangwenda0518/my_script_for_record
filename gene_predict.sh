#!/usr/bin/bash

#文件准备
#1.基因组文件：genome.fasta
#2.软屏蔽文件：genome.softmask.fasta
#3.同源蛋白文件：homolog.fasta
#4.转录组文件：rna.1.fq
#5.转录组文件：rna.2.fq
#6.线程数：
$genome
$soft_masked_genome
$homolog
$rna1
$rna2
$cpu

mkdir result 
#**********************************************************************************************************************#
####################################  近缘物种蛋白预测     ###########################
                            mkdir homolog_predict;cd homolog_predict
####################################                      ###########################
#**********************************************************************************************************************#

                    ################     miniprot       #################
mkdir miniprot ;cd miniprot

miniprot -I -u --gff -t 60  $soft_masked_genome $homolog |sed '/##PAF/d'  >miniprot.gff3 

cp miniprot.gff3 ../../result
cd ..

                    ################      genewise      #################
mkdir genewise ;cd genewise
#同源蛋白和基因组序列比对注释
homolog_genewise --cpu 64 --coverage_ratio 0.4 --evalue 1e-9 --max_intron_length 20000 $homolog $soft_masked_genome 
#将结果转化为gff3格式        输入homolog_genewiseGFF2GFF3 查看帮助信息
homolog_genewiseGFF2GFF3 --input_genewise_start_info genewise.start_info.txt --genome $soft_masked_genome  --min_score 15 \
                         --gene_prefix genewise genewise.gff > genewise.gff3 2> gene_id.with_stop_codon.list
cp genewise.gff3 ../../result

cd ..
cd ..
#**********************************************************************************************************************#
#####################################         转录本预测         ##################################
                                    mkdir rna_predict;cd rna_predict
#####################################                           ##################################       
#**********************************************************************************************************************#                            
#making_hints
mkdir -p making_hints ; cd making_hints

hisat2-build -p 60 $soft_masked_genome genome

hisat2 -x genome -p 80 --dta --dta-cufflinks -1 $rna1 -2 $rna2 \
            -S rna.sam --new-summary --summary-file rna.hisat2.summary --very-sensitive -t 

samtools sort -@ 80 -O bam -o rna.sort.bam rna.sam
#生成转录本hints文件
bam2hints --intronsonly --in=rna.sort.bam --out=hints.gff
cd ..

#转录组组装
mkdir RNA_assembly;cd RNA_assembly
#######无参组装
Trinity --seqType fq --max_memory 80G --left $rna1  --right $rna2  --CPU 120  --normalize_reads  \
        --bflyCalculateCPU --output trinity_denovo --full_cleanup  &> trinity_denovo.log
#######有参组装
Trinity --max_memory 50G --CPU 40 --jaccard_clip --normalize_reads --genome_guided_bam merged.sort.bam \
        --genome_guided_max_intron 20000 --output trinity_genomeGuided --bflyCalculateCPU &> trinity_genomeGuided.log

cd ..
#**********************************************************************************************************************#
#####################################      1.先比对-组装-预测     #################################
#**********************************************************************************************************************#
mkdir transdecoder ;cd transdecoder
#转录本组装 
stringtie -p 10 -o stringtie.gtf ../making_hints/rna.sort.bam
#第一步: 从GTF文件中提取FASTA序列
perl ~/biosoft/gene_prediction/transDecoder/util/gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fasta > stringtie.transcripts.fasta
#第二步: 将GTF文件转成GFF3格式
perl ~/biosoft/gene_prediction/transDecoder/util/gtf_to_alignment_gff3.pl stringtie.gtf > stringtie.gff3
#第三步: 根据6个读码框翻译得到蛋白序列，预测转录本中长的开放阅读框, 默认最小蛋白长度是100个氨基酸，可以用-m修改
perl ~/biosoft/gene_prediction/transDecoder/TransDecoder.LongOrfs -t stringtie.transcripts.fasta
###对序列名，替换操作，仅保留序列名，以免程序运行失败，只适用Trinity
###perl -pe 's/^>(\S+).*/>$1/' longest_orfs.pep > longest_orfs.pep
#第四步:将ORF蛋白序列比对到公共数据库，能比对上的则表示存在对应转录本序列
# 比对到SwissProt数据库
ln -s ~/database/swiss_port/Plants.fa ./uniprot_sprot.fasta

diamond makedb --threads 80 --db uniprot_sprot --in uniprot_sprot.fasta

diamond blastp --db uniprot_sprot --query stringtie.transcripts.fasta.transdecoder_dir/longest_orfs.pep --out blast.xml \
        --outfmt 5 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 20 --index-chunks 1 --threads 80

parsing_blast_result.pl --no-header --max-hit-num 20 --query-coverage 0.2 --subject-coverage 0.1 blast.xml > blastp.outfmt6

# 比对到PFAM数据库
para_hmmscan.pl --cpu 80 stringtie.transcripts.fasta.transdecoder_dir/longest_orfs.pep  --hmm_db /opt/biosoft/hmmer-3.3.2/Pfam-A.hmm >pfam.domtbl

#第五步：去除假阳性ORF，保留正确的ORF，使用上面两个比对结果，进行假阳性去除
TransDecoder.Predict -t stringtie.transcripts.fasta --retain_pfam_hits pfam.domtbl --retain_blastp_hits blastp.outfmt6

#第六步: 生成基于参考基因组的编码区注释文件
perl ~/biosoft/gene_prediction/transDecoder/util/cdna_alignment_orf_to_genome_orf.pl stringtie.transcripts.fasta.transdecoder.gff3 \
                                         stringtie.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

GFF3Clear --genome $soft_masked_genome --gene_prefix transDecoder --GFF3_source transDecoder *.transdecoder.gff3 > out; mv out transDecoder.gff3

cd ..
#**********************************************************************************************************************#
######################################      2.先组装-pasa比对-预测      #############################################33
#**********************************************************************************************************************#
ln -s ../making_hints/stringtie.gtf ./
# 将 RNA-Seq de novo 组装序列和参考组装序列合并到一个文件中
cat Trinity-denovo.fasta  Trinity-GG.fasta >transcripts.fasta
#获得Denovo组装结果的序列名，每行一条序列名，用于构建综合的转录组数据库,不能用有参转录本结果
perl -e 'while (<>) { print "$1\n" if />(\S+)/ }' Trinity-denovo.fasta > tdn.accs
# 对 transcripts，去除污染序列，截断低质量 end-trimming (vector, adaptor, primer, polyA/T tails)
seqclean transcripts.fasta -v ~/biosoft/gene_prediction/pasa-2.5/UniVec/UniVec -c 10 
#-c 线程数，不能太大会报错 20就报错了
rm -rf cleaning_*
#####################################
#修改sqlite数据库配置
#vi sqlite.confs/alignAssembly.config
#DATABASE=/tmp/dxc_mydb_pasa.sqlite
#或者删掉
rm /tmp/sample_mydb_pasa.sqlite
## 运行 PASA 主程序，将 transcripts 序列比对到基因组上，得到去冗余的转录子序列、转录子和基因组的比对结果和可变剪接信息
perl ~/biosoft/gene_prediction/pasa-2.5/Launch_PASA_pipeline.pl -c sqlite.confs/alignAssembly.config -R -C -g genome.fasta \
    -t transcripts.fasta.clean -T -u transcripts.fasta  --ALIGNERS gmap,blat,minimap2 --CPU 40 --stringent_alignment_overlap 30.0 \
    --TDN tdn.accs --MAX_INTRON_LENGTH 20000 --TRANSDECODER  &> pasa.log

rm tmp-*

DBName=`ls *.assemblies.fasta | perl -pe 's/.assemblies.fasta\s*//'`

perl ~/biosoft/gene_prediction/pasa-2.5/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta $DBName.assemblies.fasta \
    --pasa_transcripts_gff3 $DBName.pasa_assemblies.gff3
### 10 min
ln -fs $DBName.assemblies.fasta.transdecoder.genome.gff3 pasa.gff3
##对GFF3文件格式进行修改
GFF3Clear --genome genome.fasta --gene_prefix pasa --GFF3_source PASA pasa.gff3 > out; mv out pasa.gff3
# 进一步提取完好的基因模型（此脚本在geta里），用于重头预测的训练
geneModels2AugusutsTrainingInput --cpu 180 --min_evalue 1e-9 --min_identity 0.7 --min_coverage_ratio 0.7 --min_cds_num 1 \
            --min_cds_length 450 --min_cds_exon_ratio 0.4 --keep_ratio_for_excluding_too_long_gene 0.99 pasa.gff3 genome.fasta


#**********************************************************************************************************************#
#####################################          重头预测          ##################################
                                    mkdir ab_predict ; cd ab_predict
#####################################                           ##################################
#**********************************************************************************************************************#


##使用 GeneMark-ES 进行基因预测  (http://topaz.gatech.edu/GeneMark/)
mkdir -p genemark_es_et/es ; cd genemark_es_et/es
#数据准备
gmes_petap.pl --sequence $soft_masked_genome --ES --cores 60  
#--cores 60 过大报错，160
#将gtf转换成gff3格式（由pasa提供）
perl /opt/biosoft/PASApipeline.v2.4.1/misc_utilities/gtf_to_gff3_format.pl genemark.gtf  $soft_masked_genome > genemark.gff3
#格式整理
GFF3Clear --genome $soft_masked_genome --no_attr_add --GFF3_source GeneMarkES --gene_prefix gmes genemark.gff3 > ../genemarkES.gff3

cd ..

#使用 GeneMark-ET 进行基因预测
mkdir -p et;cd et
#使用ET算法时，需要输入含有内含子信息的GFF文件，可以由Tophat的结果文件转化而来，或者使用geta中augustus的中间文件（推荐）
ln -s ../making_hints/hints.gff ./

hints2genemarkETintron.pl $soft_masked_genome hints.gff > genemartET.intron.gff

gmes_petap.pl --sequence $soft_masked_genome --ET genemartET.intron.gff --et_score 10 --cores 60 

#将gtf转换成gff3格式（由pasa提供）
/opt/biosoft/PASApipeline.v2.4.1/misc_utilities/gtf_to_gff3_format.pl genemark.gtf $soft_masked_genome > genemark.gff3
#格式整理
GFF3Clear --genome $soft_masked_genome --no_attr_add --GFF3_source GeneMarkET --gene_prefix gmet genemark.gff3 > ../genemarkET.gff3

cd ..

cp genemarkET.gff3 genemarkES.gff3 ../../result/

cd ..

#############   snap预测    ################

mkdir snap ; cd snap
ln -s ../pasa/out.filter2.gff3 ./
export PERL5LIB=$PERL5LIB:/opt/biosoft/EVM_r2012-06-25/PerlLib/
##1.将gff3文件转化成zff文件，并准备对应的参考序列文件，推荐pasa，august，GeneMark-ET结合RNA预测的结果进行预测。
/opt/biosoft/EVM_r2012-06-25/OtherGeneFinderTrainingGuide/SNAP/gff3_to_SNAP_train.pl out.filter2.gff3 $soft_masked_genome
#生成genome.ann genome.dna
###2.检测基因注释错误，然后对其进行修正和丢弃
#第一步是查看基因的一些特征：
fathom genome.ann genome.dna -gene-stats &> gene-stats.log
#接下来，验证基因没有明显的错误：
fathom genome.ann genome.dna -validate &> validate.log
perl -ne 'print "$1\n" if /.*:\s+(\S+)\s+OK/' validate.log > zff2keep.txt
perl -e 'open IN, "zff2keep.txt"; while (<IN>) { chomp; $keep{$_} = 1; } while (<>) { if (m/>/) { print; } \
        else { chomp; @_ = split /\t/; print "$_\n" if exists $keep{$_[-1]}; } }' \
        genome.ann > out;mv out genome.ann
###3.1 将序列分解为每个序列一个基因的片段：
fathom genome.ann genome.dna -categorize 1000
rm  err.* olp.* wrn.*

###3.2 将 uni 基因转换为正链：
fathom uni.ann uni.dna -export 1000 -plus

###4.运行参数计算程序
mkdir params; cd params
forge ../export.ann ../export.dna
cd ..
###5.构建HMM文件
hmm-assembler.pl species params > species.hmm

#2，使用snap命令进行最后的基因预测

snap species.hmm $soft_masked_genome > snap_out.zff

###   20min
##不能用nohup 会引入错误字符串 scoring....decoding.10.20.30.40.50.60.70.80.90.100
#3，将预测得到的zff结果转换成GFF3结果。使用evm的脚本，同时添加ID和Parent值，以避免下一步报错

export PERL5LIB=$PERL5LIB:/opt/biosoft/EVM_r2012-06-25/PerlLib/
/opt/biosoft/EVM_r2012-06-25/OtherGeneFinderTrainingGuide/SNAP/SNAP_output_to_gff3.pl snap_out.zff $soft_masked_genome > snap_out.gff3

perl -i -pe 's/scoring....decoding.10.20.30.40.50.60.70.80.90.100 done//' snap_out.zff
perl -i -ne 'print unless /^$/' snap_out.zff ##去空行

###最后格式化gff3文件
GFF3Clear --genome $soft_masked_genome --no_attr_add --GFF3_source SNAP --gene_prefix snap snap_out.gff3 > snap.gff3
cp snap.gff3 ../../result/
cd ..


#augustus








提取蛋白
gff3_to_protein.pl $soft_masked_genome miniprot.gff3 > miniprot.pep.fa

busco.sh -i miniprot.pep.fa -m 3

cd ..