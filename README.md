# my_script_for_record
some scripts for myself !
any question issue me!

有朋自远方来，不亦乐乎！

WeChat：15651721689


#GFAP-linux2.py 使用

GFAP 是一款新的植物基因功能注释的软件，我常用于KEGG和GO的分析工作。

但是原脚本限制只能在程序的目录进行使用，当我进行多个基因组的注释时，会遇到很多问题。因此，我尝试进行了代码的补充。

1.可以在任意工作目录使用

2.可以指定输出文件目录

3.运行日志，进行了存放，以保持页面整洁。

示例：

python3 ~/biosoft/GFAP-linux/GFAP-linux2.py -go  -pfam -kegg -awd nr -cpu 20 -qp test.fa -o gfap-nr
