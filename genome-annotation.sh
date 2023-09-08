#!/usr/bin/bash
if [ $# -lt 1 ]; then
        echo "nohup genome-annotation.sh proteins.fasta species cpu &"
        exit 1
fi
proteins=$1
species=$2
threads=30
threads=$3
##去除非法字符 * 
sed  's/\*//g' $proteins >proteins.fasta
## 01 对基因组蛋白质序列进行Nr注释
mkdir 01.Nr  ; cd 01.Nr
ln -s ~/database/nr-db/diamond-nr/diamond-nr.dmnd ./
diamond blastp --db diamond-nr.dmnd --query ../proteins.fasta --out Nr.xml --outfmt 5 \
               --sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1 --threads 60 &>diamond.nr.log &
Nr_annotation=$!
echo "01.Nr_annotation is running"
cd ..

#02 Swiss-Prot注释

mkdir -p 02.Swiss-Prot  ; cd 02.Swiss-Prot
ln -s  ~/database/funannotate_db/uniprot.dmnd ./

diamond blastp --db uniprot.dmnd --query ../proteins.fasta --out uniprot_sprot.xml --outfmt 5 \
               --sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1 --threads 10 &>diamond.kog.log &
swissport_annotation=$!
echo "02.Swissport_annotation is running"
cd ..

#03  KOG 注释
mkdir -p 03.KOG ; cd 03.KOG
ln -sf ~/database/kog/kog.dmnd ./
diamond blastp --db kog --query ../proteins.fasta --out kog.xml --outfmt 5 \
               --sensitive --max-target-seqs 200 --evalue 1e-5 --id 10 --tmpdir /dev/shm --index-chunks 1 --threads 10 &>diamond.kog.log &
kog_annotation=$!
echo "03.Kog_annotation is running"
cd ..

#04 eggNOG 注释
mkdir -p 04.eggNOG  ; cd 04.eggNOG
ln -sf /home/zhangwenda/database/eggnog/eggnog_proteins.dmnd ./
diamond blastp --db eggnog_proteins  --query ../proteins.fasta --out eggNOG.xml --outfmt 5 \
               --sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --tmpdir /dev/shm \
               --index-chunks 1 --threads 30 &>diamond.eggNOG.log &
eggNOG_annotation=$!
echo "04.EggNOG_annotation is running"
cd ..

#05 Interpro 注释
mkdir -p 05.InterPro  ; cd 05.InterPro
export PATH=/home/zhangwenda/sysoft/jdk-11.0.1/bin:$PATH
python3 ~/database/interproscan-5.55-88.0/interproscan.sh --output-file-base interpro -dp --cpu 30 \
            --formats TSV,XML,GFF3 --goterms --input ../proteins.fasta &>interproscan.log &
interproscan_annotation=$!
echo "05.Interproscan_annotation is running"
cd ..

#6. Pfam 注释
mkdir -p 06.Pfam  ; cd 06.Pfam
para_hmmscan.pl --outformat --cpu 20 --hmm_db ~/database/funannotate_db/Pfam-A.hmm ../proteins.fasta > Pfam.tab 2>pfam.log &
pfam_annotation=$!
echo "06.pfam_annotation is running"

cd ..

#07.kegg注释
mkdir  07.KEGG;cd 07.KEGG
kofamscan -o kofam.out --ko-list /home/zhangwenda/database/kofam_scan/db/ko_list \
          --profile /home/zhangwenda/database/kofam_scan/db/profiles --cpu 10 \
          --format mapper -e 1e-5 ../proteins.fasta &>kofam.log &
kegg_annotation=$!
echo "07.kegg_annotation is running"
cd ..


####################################        结果处理           ###########################################
##########
wait $Nr_annotation
cd 01.Nr 
    parsing_blast_result.pl --out-hit-confidence --suject-annotation Nr.xml > Nr.tab
    nr_species_distribution.pl Nr.tab > Nr_species_distribution.txt
    gene_annotation_from_Nr.pl Nr.tab > Nr.txt
    cp Nr.tab ../functional_annotation.Nr.tab
    cp Nr.txt ../functional_annotation.Nr.txt
    cp Nr_species_distribution.txt ../functional_annotation.Nr.species_distribution
echo "01.Nr_annotation is done!"
cd ..
#########
cd 02.Swiss-Prot
wait $swissport_annotation
    parsing_blast_result.pl --out-hit-confidence --suject-annotation uniprot_sprot.xml > uniprot_sprot.tab
    #对BLAST的xml或tab格式的结果进行解析和过滤，得到更准确的BLAST结果。结果为表格形式（BLAST outfmt6），结果按query序列的ID排序，每个query序列的比对结果按得分排序。
    gene_annotation_from_SwissProt.pl uniprot_sprot.tab > SwissProt.txt
    #程序用于提取基因的SwissPort注释结果。Blast分析后的SwissProt注释结果中含有多种结果，本程序提取得到基因最可能的注释结果。
    #程序输入blast outfmt 6表格形式的注释结果，分析输入文件最后一列表示SubjectAnnotation的描述结果，获得基因的注释结果。
    #其结果信息分两列，第一列基因ID，第二列注释的蛋白名称。每个基因一行，若基因有多个可能的注释结果，以分号分割。
    cp uniprot_sprot.tab ../functional_annotation.Swiss-Prot.tab
    cp SwissProt.txt ../functional_annotation.Swiss-Prot.txt
echo "02.Swissport_annotation is done!"
cd ..
########
wait $kog_annotation
cd 03.KOG
    cog_from_xml.pl --coverage 0.2 --evalue 1e-5  --db-fasta ~/database/cog/kyva \
                --db-class ~/database/cog/kog --fun-txt ~/database/cog/kog.fun.txt kog.xml

    cut -f 1,3,4 out.annot | gene_annotation_from_table.pl - > KOG.txt

    cog_R.pl --title "KOG Function Classification of Whole Genome Genes of $species" --y-name "Number of Genes" out.class

    cp out.annot ../functional_annotation.KOG.tab
    cp KOG.txt ../functional_annotation.KOG.txt
    cp out.pdf ../functional_annotation.KOG.pdf
    #cp out.png ../functional_annotation.KOG.png
echo "03.Kog_annotation is done!"
cd ..
##########
wait $eggNOG_annotation
cd 04.eggNOG
    #解析
    parsing_blast_result.pl --no-header --max-hit-num 1 --HSP-num 1 eggNOG.xml | cut -f 1,2,11,12 > eggNOG.emapper.seed_orthologs
    python3 ~/biosoft/eggnog-mapper/emapper.py -m no_search --annotate_hits_table eggNOG.emapper.seed_orthologs \
                                                -o eggNOG --dbmem --override --cpu 60
    #结果处理
    ln -sf eggNOG.emapper.annotations eggNOG.annot
    cp eggNOG.annot ../functional_annotation.eggNOG.tab
    grep -v -P "^#" eggNOG.emapper.annotations | cut -f 1,8 > eggNOG.txt
    cp eggNOG.txt ../functional_annotation.eggNOG.txt
    mkdir tbtools-go_kegg
    java -cp /home/zhangwenda/biosoft/tools/TBtools-1.098769/TBtools_JRE1.6.jar biocjava.bioIO.BioSoftPipeServer.eggNogMapperResult \
         --inEggNOGAnno eggNOG.emapper.annotations --outDir tbtools-go_kegg
    mv tbtools-go_kegg ../
echo "04.EggNOG_annotation is done!"
cd ..
##########
wait $interproscan_annotation
cd 05.InterPro
    grep IPR interpro.tsv | cut -f 1,12,13 | gene_annotation_from_table.pl - > InterPro.txt
    cp interpro.tsv ../functional_annotation.InterPro.tab
    cp InterPro.txt ../functional_annotation.InterPro.txt
    cp interpro.gff3 ../functional_annotation.InterPro.gff3
    #cp interpro.html.tar.gz ../functional_annotation.InterPro.html.tar.gz
    #cp interpro.svg.tar.gz ../functional_annotation.InterPro.svg.tar.gz
echo "05.Interproscan_annotation is done"
cd ..
###########
wait $pfam_annotation
cd 06.Pfam
    cut -f 1,2,7 Pfam.tab | perl -e '<>; print <>' | gene_annotation_from_table.pl - > Pfam.txt
    cp Pfam.tab ../functional_annotation.Pfam.tab
    cp Pfam.txt ../functional_annotation.Pfam.txt
echo "06.pfam_annotation is done!"
cd ..
############
###结果解析
wait $kegg_annotation
cd 07.KEGG
    gene_annotation_from_kaas.pl kofam.out >kofam.KEGG.out
    python3 ~/myscripts/kofamscan_plus/kofamscan_plus.py  -K ~/database/kegg/ko00001.keg -i kofam.out -o kegg.anno.xls
    cp kofam.KEGG.out ../functional_annotation.KEGG.txt
echo "07.kegg_annotation is done!"
cd ..

#GO注释
mkdir 08.GO  ; cd 08.GO
#整合
go_from_eggNOG_and_interpro.pl ../04.eggNOG/eggNOG.emapper.annotations ../05.InterPro/interpro.tsv > go.annot
#去除上层的GO编号，保留注释更精细的子级编号
go_reducing_go_number_para.pl /opt/biosoft/go_class/bin/go.obo go.annot 80 > go_reduced.annot
sort go_reduced.annot > go.annot; rm go_reduced.annot
#go_reducing_go_number.pl程序通过gene_ontology_edit.obo文件解析GO term之间的is_a关系；然后对单个基因的GO注释编号之间的关系进行解析，去除上层的GO编号，保留注释更精细的子级编号。
gene_annotation_from_table.pl go.annot > GO.txt



# GO wego分类图（按照goslim_agr.obo方法分了53类）
perl /opt/biosoft/go_class/bin/annot2wego.pl go.annot > go.wego

get_Genes_From_GO.pl /opt/biosoft/go_class/bin/go.obo go.wego > go_class.tab

go_svg.pl --outdir ./ --name out --color "green" --mark "Whole Genome Genes" --note "GO Class of whole genome genes of $species " go.wego

perl -p -i -e 's/MovePer:0.125/MovePer:0.5/; s/FontSize:.*/FontSize:24/;' out.lst

perl /opt/biosoft/go_class/svg/distributing_svg.pl out.lst out.svg

perl /opt/biosoft/go_class/bin/changsvgsize.pl out.svg 150 -100

convert out.svg out.png

# GO分类（按照goslim_agr.obo方法分了53类）
#程序输入全基因组的GO注释结果和OBO文件，对GO注释进行分类统计。
go_class.pl /opt/biosoft/go_class/bin/go.obo go.annot --go_class_method  /opt/biosoft/go_class/bin/goslim_plant.obo > go_class.tab
#输出的结果文件有5列：第1列是GO三大类的名称；第2列是GO分类编号；第3列是对GO编号的描述；第4列是该GO编号下包含的基因数目；第5列是由逗号分割的所有基因IDs。
cp go.annot ../functional_annotation.GO.tab
cp GO.txt ../functional_annotation.GO.txt
cp out.svg ../functional_annotation.GO_class.svg
cp out.png ../functional_annotation.GO_class.png
cp go_class.tab ../functional_annotation.GO_class.tab
echo "08.go_annotation is done!"
cd ..


