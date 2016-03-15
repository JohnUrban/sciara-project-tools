#!/bin/bash


if [ $# -eq 0 ];then echo; echo "characterize-qv-PilonChanges.sh qv_pilon.changes"; echo; echo "qv_pilon.changes is output of qv-of-pilon-changes.py" ; echo; exit 0; fi

#echo Total Number of Changes
n_changes=`cat $1 | wc -l`
#echo n_changes $n_changes
#echo

#echo Number of Substitution Changes \(with range of median QV and mean of median QVs\)
n_sub=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub"' | wc -l`; #echo n_sub $n_sub
n_sub_1_for_1=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == 1 && $3 == 1' | wc -l`; #echo n_sub_1_for_1 $n_sub_1_for_1
sub_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub"' | cut -f 5 | awkRange`; #echo sub_med_Q_range $sub_med_Q_range
sub_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub"' | cut -f 5 | awkMean`; #echo sub_med_Q_mean $sub_med_Q_mean
n_sub_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $5 > 40' | wc -l`; #echo n_sub_med_Q_gt_40 $n_sub_med_Q_gt_40
sub_size_range_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $5 > 40' | cut -f 2 | awkRange`; #echo sub_size_range_med_Q_gt_40 $sub_size_range_med_Q_gt_40
#echo

#echo Number of Substitution Changes where the sub is the same length as the old seq -- and size range -- \(and range of median QV and mean of median QVs\)
n_sub_e=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3' | wc -l`; #echo n_sub_e $n_sub_e
sub_e_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3' | cut -f 2 | awkRange`; #echo sub_e_size_range $sub_e_size_range
sub_e_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3' | cut -f 5 | awkRange`; #echo sub_e_med_Q_range $sub_e_med_Q_range
sub_e_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3' | cut -f 5 | awkMean`; #echo sub_e_med_Q_mean $sub_e_med_Q_mean
n_sub_e_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3 && $5 > 40' | wc -l`; #echo n_sub_e_med_Q_gt_40 $n_sub_e_med_Q_gt_40
sub_e_Q_gt_40_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 == $3 && $5 > 40' | cut -f 2 | awkRange`; #echo sub_e_Q_gt_40_size_range $sub_e_Q_gt_40_size_range
#echo

#echo Number of RevComp subs --same length and Revcomp-- and size range
n_sub_rc=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1' | wc -l`; #echo n_sub_rc $n_sub_rc
sub_rc_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1' | cut -f 2 | awkRange`; #echo sub_rc_size_range $sub_rc_size_range
sub_rc_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1' | cut -f 5 | awkRange`; #echo sub_rc_med_Q_range $sub_rc_med_Q_range
sub_rc_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1' | cut -f 5 | awkMean`; #echo sub_rc_med_Q_mean $sub_rc_med_Q_mean
if [ -z $sub_rc_med_Q_mean ]; then sub_rc_med_Q_mean=-; fi
n_sub_rc_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1 && $5 > 40' | wc -l`; #echo n_sub_rc_med_Q_gt_40 $n_sub_rc_med_Q_gt_40
sub_rc_med_Q_gt_40_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$8 == "True"' | awk '$2 > 1 && $5 > 40' | cut -f 2 | awkRange`; #echo sub_rc_med_Q_gt_40_size_range $sub_rc_med_Q_gt_40_size_range
#echo

#echo Number of Substitution Changes where the sub is the different length than old seq -- and size range of difference between them \(new-old\)
n_sub_ne=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3' | wc -l`; #echo n_sub_ne $n_sub_ne
sub_ne_original_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3' | awk '{print $2}' | awkRange`; #echo sub_ne_original_size_range $sub_ne_original_size_range
sub_ne_size_diff_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3' | awk '{print $3-$2}' | awkRange`; #echo sub_ne_size_diff_range $sub_ne_size_diff_range
sub_ne_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3' | cut -f 5 | awkRange`; #echo sub_ne_med_Q_range $sub_ne_med_Q_range
sub_ne_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3' | cut -f 5 | awkMean`; #echo sub_ne_med_Q_mean $sub_ne_med_Q_mean
if [ -z $sub_ne_med_Q_mean ]; then sub_ne_med_Q_mean=-; fi
n_sub_ne_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3 && $5 > 40' | wc -l`; #echo n_sub_ne_med_Q_gt_40 $n_sub_ne_med_Q_gt_40
sub_ne_med_Q_gt_40_original_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3 && $5 > 40' | awk '{print $2}' | awkRange`; #echo sub_ne_med_Q_gt_40_original_size_range $sub_ne_med_Q_gt_40_original_size_range
sub_ne_med_Q_gt_40_diff_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "sub" && $2 != $3 && $5 > 40' | awk '{print $3-$2}' | awkRange`; #echo sub_ne_med_Q_gt_40_diff_range $sub_ne_med_Q_gt_40_diff_range
#echo

#echo Number of Deletion Changes \(old=seq,new=.\)
n_del=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del"' | wc -l`; #echo n_del $n_del
n_del_single_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 == 1' | wc -l`; #echo n_del_single_bp $n_del_single_bp
n_del_le_18_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 <= 18' | wc -l`; #echo n_del_le_18_bp $n_del_le_18_bp
n_del_gt_18_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 > 18' | wc -l`; #echo n_del_gt_18_bp $n_del_gt_18_bp
n_del_ge_100_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 >= 100' | wc -l`; #echo n_del_ge_100_bp $n_del_ge_100_bp
n_del_ge_500_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 >= 500' | wc -l`; #echo n_del_ge_500_bp $n_del_ge_500_bp
n_del_ge_1000_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 >= 1000' | wc -l`; #echo n_del_ge_1000_bp $n_del_ge_1000_bp
n_del_ge_2000_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 >= 2000' | wc -l`; #echo n_del_ge_2000_bp $n_del_ge_2000_bp
n_del_ge_5000_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $2 >= 5000' | wc -l`; #echo n_del_ge_5000_bp $n_del_ge_5000_bp
del_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del"' | cut -f 2 | awkRange`; #echo del_size_range $del_size_range
del_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del"' | cut -f 5 | awkRange`; #echo del_med_Q_range $del_med_Q_range
del_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del"' | cut -f 5 | awkMean`; #echo del_med_Q_mean $del_med_Q_mean
n_del_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $5 > 40' | wc -l`; #echo n_del_med_Q_gt_40 $n_del_med_Q_gt_40
del_med_Q_gt_40_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "del" && $5 > 40' | cut -f 2 | awkRange`; #echo del_med_Q_gt_40_size_range $del_med_Q_gt_40_size_range
#echo

#echo Number of Insertion Changes \(old=.,new=seq\)
n_ins=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins"' | wc -l`; #echo n_ins $n_ins
n_ins_single_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 == 1' | wc -l`; #echo n_ins_single_bp $n_ins_single_bp
n_ins_le_18_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 <= 18' | wc -l`; #echo n_ins_le_18_bp $n_ins_le_18_bp
n_ins_gt_18_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 > 18' | wc -l`; #echo n_ins_gt_18_bp $n_ins_gt_18_bp
n_ins_ge_100_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 >= 100' | wc -l`; #echo n_ins_gt_100_bp $n_ins_gt_100_bp
n_ins_ge_500_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 >= 500' | wc -l`; #echo n_ins_gt_500_bp $n_ins_gt_500_bp
n_ins_ge_1000_bp=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $3 >= 1000' | wc -l`; #echo n_ins_gt_1000_bp $n_ins_gt_1000_bp
ins_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins"' | cut -f 3 | awkRange`; #echo ins_size_range $ins_size_range
ins_med_Q_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins"' | cut -f 5 | awkRange`; #echo ins_med_Q_range $ins_med_Q_range
ins_med_Q_mean=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins"' | cut -f 5 | awkMean`; #echo ins_med_Q_mean $ins_med_Q_mean
n_ins_med_Q_gt_40=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $5 > 40' | wc -l`; #echo n_ins_med_Q_gt_40 $n_ins_med_Q_gt_40
ins_med_Q_gt_40_size_range=`awk 'OFS="\t" {print $5,$6,$7,$8,$9,$10,$11,$12}' $1 | awk '$1 == "ins" && $5 > 40' | cut -f 3 | awkRange`; #echo ins_med_Q_gt_40_size_range $ins_med_Q_gt_40_size_range
#echo


NAMES="n_changes n_sub n_sub_1_for_1 sub_med_Q_min sub_med_Q_max sub_med_Q_mean n_sub_med_Q_gt_40 sub_size_min_med_Q_gt_40 sub_size_max_med_Q_gt_40 n_sub_e sub_e_size_min sub_e_size_max sub_e_med_Q_min sub_e_med_Q_max sub_e_med_Q_mean n_sub_e_med_Q_gt_40 sub_e_Q_gt_40_size_min sub_e_Q_gt_40_size_max n_sub_rc sub_rc_size_min sub_rc_size_max sub_rc_med_Q_min sub_rc_med_Q_max sub_rc_med_Q_mean n_sub_rc_med_Q_gt_40 sub_rc_med_Q_gt_40_size_min sub_rc_med_Q_gt_40_size_max n_sub_ne sub_ne_original_size_min sub_ne_original_size_max sub_ne_size_diff_min sub_ne_size_diff_max sub_ne_med_Q_min sub_ne_med_Q_max sub_ne_med_Q_mean n_sub_ne_med_Q_gt_40 sub_ne_med_Q_gt_40_original_size_min sub_ne_med_Q_gt_40_original_size_max sub_ne_med_Q_gt_40_diff_min sub_ne_med_Q_gt_40_diff_max n_del n_del_single_bp n_del_le_18_bp n_del_gt_18_bp n_del_ge_100_bp n_del_ge_500_bp n_del_ge_1000_bp n_del_ge_2000_bp n_del_ge_5000_bp del_size_min del_size_max del_med_Q_min del_med_Q_max del_med_Q_mean n_del_med_Q_gt_40 del_med_Q_gt_40_size_min del_med_Q_gt_40_size_max n_ins n_ins_single_bp n_ins_le_18_bp n_ins_gt_18_bp n_ins_ge_100_bp n_ins_ge_500_bp n_ins_ge_1000_bp ins_size_min ins_size_max ins_med_Q_min ins_med_Q_max ins_med_Q_mean n_ins_med_Q_gt_40 ins_med_Q_gt_40_size_min ins_med_Q_gt_40_size_max" 
VALUES="$n_changes $n_sub $n_sub_1_for_1 $sub_med_Q_range $sub_med_Q_mean $n_sub_med_Q_gt_40 $sub_size_range_med_Q_gt_40 $n_sub_e $sub_e_size_range $sub_e_med_Q_range $sub_e_med_Q_mean $n_sub_e_med_Q_gt_40 $sub_e_Q_gt_40_size_range $n_sub_rc $sub_rc_size_range $sub_rc_med_Q_range $sub_rc_med_Q_mean $n_sub_rc_med_Q_gt_40 $sub_rc_med_Q_gt_40_size_range $n_sub_ne $sub_ne_original_size_range $sub_ne_size_diff_range $sub_ne_med_Q_range $sub_ne_med_Q_mean $n_sub_ne_med_Q_gt_40 $sub_ne_med_Q_gt_40_original_size_range $sub_ne_med_Q_gt_40_diff_range $n_del $n_del_single_bp $n_del_le_18_bp $n_del_gt_18_bp $n_del_ge_100_bp $n_del_ge_500_bp $n_del_ge_1000_bp $n_del_ge_2000_bp $n_del_ge_5000_bp $del_size_range $del_med_Q_range $del_med_Q_mean $n_del_med_Q_gt_40 $del_med_Q_gt_40_size_range $n_ins $n_ins_single_bp $n_ins_le_18_bp $n_ins_gt_18_bp $n_ins_ge_100_bp $n_ins_ge_500_bp $n_ins_ge_1000_bp $ins_size_range $ins_med_Q_range $ins_med_Q_mean $n_ins_med_Q_gt_40 $ins_med_Q_gt_40_size_range"


#echo $NAMES | awk '{gsub(/\ /,"\t"); print}'
#echo $VALUES | awk '{gsub(/\ /,"\t"); print}'
#echo
#echo $NAMES | awk '{gsub(/\ /,"\n"); print}'
#echo $VALUES | awk '{gsub(/\ /,"\n"); print}'
#echo
paste <(echo $NAMES | awk '{gsub(/\ /,"\n"); print}') <(echo $VALUES | awk '{gsub(/\ /,"\n"); print}')



#for var in n_changes n_sub n_sub_1_for_1 sub_med_Q_min sub_med_Q_max sub_med_Q_mean n_sub_med_Q_gt_40 sub_size_min_med_Q_gt_40 sub_size_max_med_Q_gt_40 n_sub_e sub_e_size_min sub_e_size_max sub_e_med_Q_min sub_e_med_Q_max sub_e_med_Q_mean n_sub_e_med_Q_gt_40 sub_e_Q_gt_40_size_min sub_e_Q_gt_40_size_max n_sub_rc sub_rc_size_min sub_rc_size_max sub_rc_med_Q_min sub_rc_med_Q_max sub_rc_med_Q_mean n_sub_rc_med_Q_gt_40 sub_rc_med_Q_gt_40_size_min sub_rc_med_Q_gt_40_size_max n_sub_ne sub_ne_original_size_min sub_ne_original_size_max sub_ne_size_diff_min sub_ne_size_diff_max sub_ne_med_Q_min sub_ne_med_Q_max sub_ne_med_Q_mean n_sub_ne_med_Q_gt_40 sub_ne_med_Q_gt_40_original_size_min sub_ne_med_Q_gt_40_original_size_max sub_ne_med_Q_gt_40_diff_min sub_ne_med_Q_gt_40_diff_max n_del n_del_single_bp n_del_le_18_bp n_del_gt_18_bp n_del_ge_100_bp n_del_ge_500_bp n_del_ge_1000_bp n_del_ge_2000_bp n_del_ge_5000_bp del_size_min del_size_max del_med_Q_min del_med_Q_max del_med_Q_mean n_del_med_Q_gt_40 del_med_Q_gt_40_size_min del_med_Q_gt_40_size_max n_ins n_ins_single_bp n_ins_le_18_bp ins_size_min ins_size_max ins_med_Q_min ins_med_Q_max ins_med_Q_mean n_ins_med_Q_gt_40 ins_med_Q_gt_40_size_min ins_med_Q_gt_40_size_max ; do 
#eval "echo -e $var '\t' \$$var"
#done
