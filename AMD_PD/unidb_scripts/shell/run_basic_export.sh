#make sure you are runing the basic_export in python2
cat gene_list|xargs -I% echo "python basic_export.py "
