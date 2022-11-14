#!/usr/bin/bash
thread=60
threads=$3
##0.数据qc
mkdir 0.fastqc
fastqc -t 20 -o 0.fastqc $1  $2 &
###1.1数据矫正，使用rcorrector软件对数据进行矫正
perl ~/biosoft/RNA-seq/Rcorrector/run_rcorrector.pl -1 $1  -2 $2 -t $threads -od 1.rcorrector

###2.1 trimmomatic过滤
mkdir 2.trimmomatic&&cd 2.trimmomatic

java -jar ~/biosoft/trimmomatic-0.40/trimmomatic-0.40-rc1.jar PE -threads $threads \
../1.rcorrector/*1.cor.fq ../1.rcorrector/*2.cor.fq   \
trim.1.fq fragment.1.unpaired.fastq trim.2.fq fragment.2.unpaired.fastq \
ILLUMINACLIP:/home/zhangwenda/biosoft/trimmomatic-0.40/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:25 TOPHRED33  \
-trimlog trim.log -summary summary.log

cd ..
### 3.数据qc
mkdir 3.fastqc
fastqc -t 40 -o 3.fastqc 2.trimmomatic/trim.1.fq  2.trimmomatic/trim.2.fq &

### 4.trinity脚本归一化去除过多重复过表达序列，减小组装消耗
perl ~/biosoft/RNA-seq/trinityrnaseq-v2.14.0/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
--left 2.trimmomatic/trim.1.fq --right 2.trimmomatic/trim.2.fq \
--pairs_together --PARALLEL_STATS --CPU $threads --output 4.normalized_data

###5.使用sortmerna去除rrna
ln -s ~/database/sortmerna-db/idx/ ./

sortmerna --index 0 --threads $threads --ref ~/biosoft/RNA-seq/sortmerna-4.3.4-Linux/database/smr_v4.3_sensitive_db_rfam_seeds.fasta \
--reads 4.normalized_data/*1*fq --reads 4.normalized_data/*2*fq  \
--workdir ./ --aligned "5.sortmerna/rRNA" --other "5.sortmerna/clean" --paired_in --fastx --out2

###6.结果软链接
fq1=`basename $1`
fq2=`basename $2`
mkdir 6.result
ln -s 5.sortmerna/clean_fwd.fq  6.result\/$fq1
ln -s 5.sortmerna/clean_rev.fq  6.result\/$fq2
