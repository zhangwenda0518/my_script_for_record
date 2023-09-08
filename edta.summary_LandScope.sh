#!/usr/bin/bash
mkdir TE-REPORT;cd TE-REPORT
ln -s ../genome.fasta.mod.EDTA.anno/genome.fasta.mod.cat.gz .
mkdir Repeats-complete;cd Repeats-complete
#########################生成完整报告####################################
perl ~/biosoft/genome-feature/Plant_Annotation_TEs/ProcessRepeats/ProcessRepeats-complete.pl -species viridiplantae \
                        -nolow -noint ../genome.fasta.mod.cat.gz &
complete=$!
cd ..
#########################生成简单报告####################################
mkdir Repeats-lite ;cd Repeats-lite
perl ~/biosoft/genome-feature/Plant_Annotation_TEs/ProcessRepeats/ProcessRepeats-lite.pl -species viridiplantae \
                        -nolow -noint -a ../genome.fasta.mod.cat.gz&
lite=$!
cd ..

#
wait $lite
cd Repeats-lite
    # Rename the result file and move it to the main EDTA folder
    cp genome.fasta.mod.tbl TEs-Report-lite.txt
    cp TEs-Report-lite.txt ../../
    #ProcessRepeats-lite.pl脚本生成一个名为：TEs-Report-lite.txt 的结果文件。
    #绘制
    bash /home/zhangwenda/bin/RepeatLandScape.sh
    mv Rplots.pdf RepeatLandScape-lite.pdf
    cp RepeatLandscape.html RepeatLandscape-lite.html
    #rm align2.txt
    #rm tmp.txt
cd ..
wait $complete
cd Repeats-complete
    # Rename the result file and move it to the main EDTA folder
    cp genome.fasta.mod.tbl TEs-Report-Complete.txt
    cp TEs-Report-lite.txt ../../
    #ProcessRepeats-complete.pl脚本生成一个名为：TEs-Report-Complete.txt 的结果文件。
    #在本报告中，部分元素将以“-like”后缀命名（例如，Angela-like）
    #绘制
    bash /home/zhangwenda/bin/RepeatLandScape.sh
    mv Rplots.pdf RepeatLandScape-lite.pdf
    cp RepeatLandscape.html RepeatLandscape-lite.html
    #移动到EDTA主目录
    #rm align2.txt
    #rm tmp.txt
cd ..
    cp Repeats-*/RepeatLandscape* ../
cd ..
###################################   LTR-AGE绘制   ######################################
mkdir TE-LTR-AGE;cd TE-LTR-AGE
ln -s ../genome.fasta.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list .
#
ln -s ~/biosoft/genome-feature/Plant_Annotation_TEs/Rscripts/plot-AGE-Gypsy.R .
ln -s ~/biosoft/genome-feature/Plant_Annotation_TEs/Rscripts/plot-AGE-Copia.R .
#
# Preparing the file
cat -n genome.fasta.mod.pass.list | grep Gypsy | cut -f 1,13 | sed 's# ##g' | sed 's#^#Cluster_#g' | awk '{if ($2 > 0) print $n}' > AGE-Gypsy.txt
cat -n genome.fasta.mod.pass.list | grep Copia | cut -f 1,13 | sed 's# ##g' | sed 's#^#Cluster_#g' | awk '{if ($2 > 0) print $n}' > AGE-Copia.txt
#
# Generating the plots
Rscript plot-AGE-Gypsy.R
Rscript plot-AGE-Copia.R
cp AGE-Copia.pdf ../RepeatAGE-Copia.pdf
cp AGE-Gypsy.pdf ../RepeatAGE-Gypsy.pdf
cd ..