#!/usr/bin/bash

function usage
{
echo '''
############################################################
haplotigs.sh 这个脚本是为了简化我常使用的haplotigs命令，它常用于基因组的去杂合处理。
#2023.04.25
############################################################
Usage: haplotigs.sh -g genome.fa -g genome.fa [-l hifi.fa |-i illumina.1.fq -I illumina.2.fq] | [-lam long.aligned.bam |-sbam short.aligned.bam ]-m mode -t threads 
-g 输入基因组
-l 长读长数据
-s1 短读前端
-s2 短读后端
-lbam 长读比对bam
-sbam 短读比对bam
-m 运行模式
   1 长读模式
   2 短读模式
   3 混合运行模式
-t 线程数，默认60
'''
}
############################################################
function run_long_map
{	
	genome=$1
	long=$2
	threads=$3
	minimap2 -t $threads -x map-hifi  $genome $long -a \
		| samtools view -@ $threads -hF 256 - \
		| samtools sort -@ $threads  -o long.aligned.bam -T tmp.ali 

	purge_haplotigs readhist -b long.aligned.bam -g $genome -t $threads
}

############################################################
function run_short_map
{
	genome=$1
	s1=$2
	s2=$3
	threads=$4
	~/biosoft/bwa-mem2-2.2.1/bwa-mem2 index  $genome -p genome
	~/biosoft/bwa-mem2-2.2.1/bwa-mem2 mem -t $threads genome $s1 $s2 \
		| samtools view -@ $threads -hF 256 - \
		| samtools sort -@ $threads  -o short.aligned.bam -T tmp.ali

	purge_haplotigs readhist -b short.aligned.bam -g $genome -t $threads
}
############################################################
####默认参数设置####
mode=short
threads=60
###########################################################
intRe='^[0-9]+$'
###########################################################
while getopts ":g:l::i::I::L::S::m::t::" opt
do
  case $opt in
    g)
      genome=$OPTARG
      if [ ! -f "$genome" ]
      then
          echo "ERROR: $genome is not input"
          exit 1
      fi
      ;;
    l)
      long=$OPTARG
      ;;
	i)
	  short1=$OPTARG
	  ;;
    I)
      short2=$OPTARG
      ;;
	L)
	  lbam=$OPTARG
	  ;;
    S)
      sbam=$OPTARG
      ;;
	m)
	  mode=$OPTARG
	  ;;
    t)
      threads=$OPTARG
      if ! [[ $threads =~ $intRe ]]
      then
          echo "ERROR: threads should be an integer, $THREADS is not an integer"
          exit 1
      fi
      ;;
	\?)
		usage
		exit 1
		;;
	esac
done
######################################################################
if [ $# -lt 1 ]; then
	usage
	exit
fi
######################################################################
if [ $mode = 1 ] ||[ $mode = long ] 
	then mode=long
elif [ $mode = 2 ] || [ $mode = short ] 
	then mode=short
elif [ $mode = 3 ] || [ $mode = both ]
	then mode=both
else 
	usage
	exit 1
fi


source /home/zhangwenda/sysoft/anaconda3/bin/activate purge_haplotigs

echo "###################################################################"
echo "                        Beigin Run                             "
echo "                        Wrate  ...                             "
echo "###################################################################"

#长读数据模式或混合模式
if [ $mode = long ] || [ $mode = both ]
	then 
		if  [ -n "$lbam" ]
			then 
				echo "long_bam_mode running ....."
				echo "########################"
				purge_haplotigs readhist -b $lbam -g $genome -t $threads &>long.bam.log
			else
				echo "long_map_mode running ....."
				echo "########################"
				run_long_map $genome $long $threads &>long.map.log
		fi
fi
#短读数据模式或混合模式
if [ $mode = short ] || [ $mode = both ]
	then 
		if  [ -n "$sbam" ]
			then
			   echo "short_bam_mode running ....."	
			   echo "########################"
				purge_haplotigs readhist -b $sbam -g $genome -t $threads &>short.bam.log
			else
				echo "short_map_mode running ....."
				echo "########################"
				run_short_map $genome $short1 $short2 $threads &>short.map.log
		fi
fi
echo "###################################################################"
echo "                        ALL succesffully Done                      "
echo "###################################################################"

source /home/zhangwenda/sysoft/anaconda3/bin/deactivate
