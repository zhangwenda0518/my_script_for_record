#!/bin/bash
set -e
#2022/04/12
#1.fastqc质检-2.fastauniq去重-3.bless修正-4.musket修正-5.trimmomatic质控-6.fastp质控
#bless kmerlength 25/31

#目录准备
#ln -s $1 sample_R1.fastq
#ln -s $2 sample_R2.fastq
#1.fastq
#mkdir fastqc
#nohup fastqc -t 80 -l 23 -o fastqc sample_R1.fastq  sample_R2.fastq &

#2.fastuniq
#mkdir fastuniq && cd fastuniq
#ls ../*.fastq > fragment.list
#fastuniq -i fragment.list -o uniq.1.fastq -p uniq.2.fastq >fastuniq.log
#cd ..

#3.1 bless-25
mkdir bless-25&&cd bless-25
bless
/opt/biosoft/BLESS/bless -read1 ../fastuniq/uniq.1.fastq -read2 ../fastuniq/uniq.2.fastq -kmerlength 21 -smpthread 30 -prefix out -max_mem 80 &>bless.log
mkdir fastqc
fastqc -t 80 -o fastqc out.1.corrected.fastq out.2.corrected.fastq &
cd ..
#3.2 bless-31
mkdir bless-31&&cd bless-31
bless
/opt/biosoft/BLESS/bless -read1 ../bless-25/out.1.corrected.fastq -read2 ../bless-25/out.2.corrected.fastq -kmerlength 31 -smpthread 30 -prefix out -max_mem 80 &>bless.log
mkdir fastqc
fastqc -t 80 -o fastqc out.1.corrected.fastq out.2.corrected.fastq &
cd ..
#4.musket
mkdir  musket && cd musket
musket ../bless-31/out.1.corrected.fastq ../bless-31/out.2.corrected.fastq -omulti musket -inorder -p 30 &>musket.log
mv musket.0 out.1.fastq;
mv musket.1 out.2.fastq;
mkdir fastqc
fastqc -t 80 -o fastqc out.1.fastq out.2.fastq &
cd ..
#5.trimmomatic
mkdir trimmomatic && cd trimmomatic
java -jar ~/biosoft/trimmomatic-0.40/trimmomatic-0.40-rc1.jar PE -threads 120 ../musket/out.1.fastq ../musket/out.2.fastq trimmomatic.1.fastq trimmomatic.1.unpaired.fastq trimmomatic.2.fastq trimmomatic.2.unpaired.fastq ILLUMINACLIP:/opt/biosoft/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50 TOPHRED33 &>trimmomatic.log
mkdir fastqc
fastqc -t 80 -o fastqc trimmomatic.1.fastq trimmomatic.2.fastq &
cd ..
#6.fastp
mkdir fastp && cd fastp
fastp -i ../musket/out.1.fastq -o fastp_1.fastq -I ../musket/out.2.fastq -O fastp_2.fastq --thread=16 &>fastp.log
mkdir fastqc
fastqc -t 80 -o fastqc fastp_1.fastq fastp_2.fastq &
cd ..
echo -e "sucessful ! all done/n"
