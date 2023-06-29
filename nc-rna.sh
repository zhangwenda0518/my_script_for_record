#/usr/bin/bash
if [ $# -lt 1 ] ; then 
    echo "nc-rna.sh genome.fa threads"
    exit 1
fi

genome=`realpath $1`
threads=30
threads=$2

mkdir Non-coding-RNA ;cd Non-coding-RNA

#1.1 使用rnammer 进行rRAN分析,包括5S/8S、16S/18S、23S/28S rRNA。依赖hmmer2.3

mkdir rRNAmmer ; cd rRNAmmer
    rnammer -S euk -multi -m tsu,lsu,ssu -f rRNA.fasta -h rRNA.hmmreport -xml rRNA.xml -gff rRNA.gff2 $genome &
    ## 30 min ##
        rnammer_pid=$!
cd ..

####1.2  使用barrnap进行rRNA预测分析####

mkdir barrnap ;cd barrnap
    barrnap --kingdom euk --threads $threads  $genome --outseq rRNA.fasta --quiet > rRNA.gff3 &
    barrnap_pid=$!
cd ..

#2.使用tRNAscan-SE进行tRNA分析（支持原核、真核、线粒体等），包括所在基因组位置、碱基组成、反密码子和携带氨基酸类型、以及二级结构信息。

mkdir tRNAscan ;cd tRNAscan
    tRNAscan-SE -E -o tRNA.tblout -f tRNA.ss -m tRNA.stats -a tRNA.fasta --thread $threads $genome 2>tRNAscan.log &
    tRNAscan_pid=$!
cd ..

#3.使用Rfam数据库对非编码RNA（rRNA基因、snRNA基因和miRNA基因）进行分析

mkdir rfam ; cd rfam
    ln -s  $genome ./
    ln -s ~/biosoft/infernal-1.1.4/Rfam.cm* ./
    #cmsearch
        cmsearch  --tblout cmsearch_rfam_out.tab  --rfam --cut_ga --nohmmonly --noali --cpu $threads Rfam.cm $genome > cmsearch_rfam_out.txt &>cmsearch.log  &
        cmsearch_pid=$!
    #cmpress
        cmscan  --tblout cmscan_rfam_out.tab  --rfam --cut_ga --nohmmonly --noali --cpu $threads Rfam.cm $genome > cmscan_rfam_out.txt &>cmscan.log &
        cmscan_pid=$!
cd ..

#结果解析
wait  $rnammer_pid
    cd rRNAmmer
        rRNAmmer_gff2gff3.pl rRNA.gff2 > rRNA.gff3
    cd ..
echo "rnammer run "

wait $barrnap_pid

echo "barrnap run"

wait $tRNAscan_pid
    cd tRNAscan
        tRNAscanSE2GFF3.pl tRNA.tblout tRNA.ss > tRNA.gff3 2>tRNAscan.stat
    cd ..
echo "tRNAscan run "

#陈老师脚本，提取与统计相关信息（rRNA miRNA snRNA）
wait  $cmsearch_pid
    cd rfam
        #统计rRNA
       perl /home/zhangwenda/bin/Rfam_rRNA_stats.pl cmsearch_rfam_out.tab
        #统计miRNA
       perl /home/zhangwenda/bin/Rfam_miRNA_stats.pl cmsearch_rfam_out.tab ~/biosoft/infernal-1.1.4/Rfam.cm
        #统计snRNA
       perl /home/zhangwenda/bin/Rfam_snRNA_stats.pl cmsearch_rfam_out.tab ~/biosoft/infernal-1.1.4/snoRNA_type.txt
    cd ..
echo "rnafam_cmsearch run "

wait  $cmscan_pid
    cd rfam
        #统计rRNA
        perl /home/zhangwenda/bin/Rfam_rRNA_stats.pl cmscan_rfam_out.tab
        #统计miRNA
        perl /home/zhangwenda/bin/Rfam_miRNA_stats.pl cmscan_rfam_out.tab ~/biosoft/infernal-1.1.4/Rfam.cm
        #统计snRNA
        perl /home/zhangwenda/bin/Rfam_snRNA_stats.pl cmscan_rfam_out.tab ~/biosoft/infernal-1.1.4/snoRNA_type.txt
    cd ..
echo "rnafam_cmscan run "

echo "##############################################################################"
echo "###############################ALL DONE!######################################"
echo "##############################################################################"
