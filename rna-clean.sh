#!/usr/bin/bash
#1   数据qc(后台)
#2   fastp过滤
#3   使用rcorrector软件对数据进行矫正
#4   数据qc
#5   使用sortmerna去除rrna

threads=$3

#1   数据qc(后台)
mkdir 1.fastqc
fastqc -t $threads -o 1.fastqc $1  $2 &
#2   fastp过滤
mkdir 2.fastp
fastp -i $1 -o 2.fastp/1.fq -I $2 -O 2.fastp/2.fq --thread=10 -c -x
mv fastp* 2.fastp
#3   使用rcorrector软件对数据进行矫正
perl ~/biosoft/RNA-seq/Rcorrector/run_rcorrector.pl -1 2.fastp/1.fq  -2 2.fastp/2.fq -t $threads -od 3.rcorrector
#4   数据qc
mkdir 4.fastqc
fastqc -t $threads -o 4.fastqc 3.rcorrector/1.cor.fq  3.rcorrector/2.cor.fq &
#5   使用sortmerna去除rrna
mkdir 5.sortmerna
cd 5.sortmerna
ln -s ~/database/sortmerna-db/idx/ ./

sortmerna --index 0  --ref ~/database/sortmerna-db/smr_v4.3_sensitive_db.fasta --reads ../3.rcorrector/1.cor.fq --reads ../3.rcorrector/2.cor.fq  --workdir ./ --aligned rRNA --other non_rRNA --paired_in --fastx --out2 --threads $threads 
cd ..
cp 5.sortmerna/non_rRNA_fwd.fq clean_R1.fastq
cp 5.sortmerna/non_rRNA_rev.fq clean_R2.fastq
