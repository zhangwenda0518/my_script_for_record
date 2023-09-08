#!/usr/bin/bash
:<<EOF
为简化GATK文件准备，写的流程脚本，使用前，请注意安装好bowtie2,gatk,samtools,Picard(并修改脚本picard位置)
EOF

#参数解析
if [ $# -lt 1 ]; then
    echo "useage: gatk-bam.sh genome.fa illumina.1.fastq illumina.2.fastq name cpu"
    exit 1
fi

genome=`realpath $1`
fq1=`realpath $2`
fq2=`realpath $3`
name=$4
cpu=60
cpu=$5

#z工作目录
mkdir $name ;cd $name

#0.文件准备
echo "#0.文件准备"

mkdir 0.data ;cd 0.data
    ln -s $fq1 sample_R1.fastq
    ln -s $fq2 sample_R2.fastq
    ln -s $genome genome.fasta
    bowtie2-build --threads $cpu genome.fasta genome &>bowtie2-build.log 
cd ..

echo "#0.文件准备 ,完成"

#1.trim序列修剪
echo "#1.trim序列修剪"
mkdir 1.trim-fastp;cd 1.trim-fastp
    fastp -i ../0.data/sample_R1.fastq -o fastp_1.fastq -I ../0.data/sample_R2.fastq -O fastp_2.fastq  --thread=16 1>$name.fastp.log 2>$name.fastp.errors.log
cd ..
echo "#1.trim序列修剪 ，完成"

#2.利用bowtie2进行不严格比对（默认参数）
echo "bowtie2不严格比对"

mkdir 2.bowtie2 ;cd 2.bowtie2
    bowtie2 -p $cpu -x ../0.data/genome -1 ../1.trim-fastp/fastp_1.fastq -2 ../1.trim-fastp/fastp_2.fastq -S $name.sam --rg-id $name --rg "PL:Illumina" --rg "SM:$name" 1>$name.bowtie2.log 2> $name.bowtie2.errors.log
    samtools sort -@ $cpu -O BAM -o $name.bam $name.sam
cd ..
echo "bowtie2不严格比对 ，完成"

#3.利用Picard软件去除PCR重复
echo "#3.利用Picard软件去除PCR重复 "

mkdir 3.picard ;cd 3.picard
    java -jar ~/biosoft/binary/picard.jar MarkDuplicates --REMOVE_DUPLICATES true -I ../2.bowtie2/$name.bam -O $name.MD.bam -M $name.metrics 1>picard.md1.log 2>picard.md1.errors.log
cd ..
echo "#3.利用Picard软件去除PCR重复 ,完成"

#4.将去除了PCR重复的SAM、BAM文件转换得到比对的fastq文件
echo "#4.将去除了PCR重复的SAM、BAM文件转换得到比对的fastq文件"
mkdir 4.sam2fq;cd 4.sam2fq
    java -jar ~/biosoft/binary/picard.jar SamToFastq -I ../3.picard/$name.MD.bam -F $name.MD.1.fastq -F2 $name.MD.2.fastq 1>picard.sam2fq.log 2>picard.sam2fq.errors.log
cd ..
echo "#4.将去除了PCR重复的SAM、BAM文件转换得到比对的fastq文件 ,完成"

#5.利用BLESS软件对去除了PCR重复的数据进行reads修正

echo "#5.利用BLESS软件对去除了PCR重复的数据进行reads修正,运行"
mkdir 5.bless;cd 5.bless
	/opt/biosoft/BLESS/bless -read1 ../4.sam2fq/$name.MD.1.fastq -read2 ../4.sam2fq/$name.MD.2.fastq -kmerlength 21 -smpthread $cpu -prefix $name.bless -max_mem 80  1>$name.bless.log 2>$name.bless.errors.log
cd ..
echo "#5.利用BLESS软件对去除了PCR重复的数据进行reads修正,运行完成"

#6.严格的阈值进行bowtie2比对
echo "#6.严格的阈值进行bowtie2比对,运行"
mkdir 6.bowtie2-strict ;cd 6.bowtie2-strict
    bowtie2 -p $cpu -x ../0.data/genome --score-min L,-0.3,-0.3 --1 ../5.bless/$name.bless.1.corrected.fastq -2 ../5.bless/$name.bless.2.corrected.fastq -S $name.cor.sam --rg-id $name --rg "PL:Illumina" --rg "SM:$name" 1>$name.cor.bowtie2.log 2> $name.cor.bowtie2.error.log
    samtools sort -@ $cpu -O BAM -o $name.cor.bam $name.cor.sam
cd ..
echo "#6.严格的阈值进行bowtie2比对,运行完成"

#7.再次利用Picard软件去除PCR重复，最终得到用于下游分析的GATK的BAM文件
echo "#7.再次利用Picard软件去除PCR重复，最终得到用于下游分析的GATK的BAM文件,运行"
mkdir 7.picard-cor ;cd 7.picard-cor
    java -jar ~/biosoft/binary/picard.jar MarkDuplicates --REMOVE_DUPLICATES true -I ../6.bowtie2-strict/$name.cor.bam -O $name.cor.MD.bam -M $name.cor.metrics 1>picard.md2.log 2>picard.md2.errors.log
    bam=`realpath $name.cor.MD.bam`
cd ..
echo "#7.再次利用Picard软件去除PCR重复，最终得到用于下游分析的GATK的BAM文件,运行完成"

echo "----------------------------------------------------------------------------------------------------------"
echo "                                           DATA Preparation ALL DONE SUCCESSFULLY                                          "
echo "----------------------------------------------------------------------------------------------------------"

#8.单个样本的变异位点分析
#.文件准备
echo "#8.gatk文件准备"

mkdir 8.gatk.data ;cd 8.gatk.data
    ln -s $genome genome.fasta
    ln -s $bam $name.bam
    #构建参考基因组的索引文件
    samtools faidx genome.fasta
    #构建bam文件的索引文件
    samtools index $name.bam
    #构建参考基因组的字典：
    java -jar ~/biosoft/binary/picard.jar CreateSequenceDictionary -R genome.fasta -O genome.dict
cd ..

echo "#8.gatk文件准备 ,完成"

#9.单个样本变异分析-HaplotypeCaller
mkdir 9.HaplotypeCaller;cd 9.HaplotypeCaller
    gatk HaplotypeCaller -R ../8.gatk.data/genome.fasta -I ../8.gatk.data/$name.bam -ERC GVCF -O $name.g.vcf --pcr-indel-model CONSERVATIVE \
         --sample-ploidy 2 --min-base-quality-score 10 --kmer-size 10 --kmer-size 25 1>gatk.log 2>gatk.error.log
cd ..

echo "----------------------------------------------------------------------------------------------------------"
echo "                                           GATK HaplotypeCaller DONE SUCCESSFULLY                                          "
echo "----------------------------------------------------------------------------------------------------------"