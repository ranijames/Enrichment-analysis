#!/bin/bash#
echo "The input path where you can find resuts from basic export $1"
echo "The mode used $2"
echo "The path for saving the outputs $3"
if [ -e lancet_res_ctrl_$2.sh ] 
then
   rm -- lancet_res_ctrl_$2.sh
fi
filenames=$(ls -1 $1/*csv)
for file in $filenames
do  B="$(cut -d'.' -f1 <<<"$file")"
    gene="$(cut -d'_' -f5 <<<"$B")"
        echo "python ../pd_pats.py -genename=$gene -mode=$2 $1/$file $3 > err_pats_$2_$gene &" >>lancet_pats_$2.sh
done 

#sh /data/ajames/GD_PD/scripts/lancet_pats_class.sh
               
