mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS

cd 05_TE_POSITIONS

mkdir 00_cp

cd 00_cp

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

#move over full output files to appropriate directory
cd /projects/phillipslab/shared/TE_SLiM/full_output

cp -v *_full_DOM_neutral.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_NEUTRAL
cp -v *_full_DOM_sel_all.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL
cp -v *_full_DOM_sel_all_WEAK.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL_WEAK
cp -v *_full_DOM_sel_all_WEAK_POS_MUT.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL_WEAK_POS_MUT
cp -v *_full_DOM_sel_arm_WEAK_center_Strong.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ARM_WEAK_CENTER_STRONG
cp -v *_full_DOM_sel_arm_WEAK_center_Strong_LOSS.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
cp -v *_full_DOM_sel_center_weak.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_CENTR_WEAK
cp -v *_full_NO_DOM_neutral.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_NEUTRAL
cp -v *_full_no_dom_sel_all.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL
cp -v *_full_no_dom_sel_all_weak.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL_WEAK
cp -v *_full_no_dom_sel_all_weak_pos_Mut.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL_WEAK_POS_MUT
cp -v *_full_no_dom_sel_arm_weak_Center_Strong.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
cp -v *_full_Nno_dom_sel_arm_weak_Center_Strong-LOSS.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
cp -v *_full_no_dom_sel_center_weak.txt /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_CENTR_WEAK


#get only the "mutations" part of the outout

mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_NEUTRAL
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL_WEAK
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/DOM_SEL_CENTR_WEAK
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_NEUTRAL
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL_WEAK
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/00_cp/NO_DOM_SEL_CENTR_WEAK
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_CENTR_WEAK/$i; done


#remove first and last lines

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_NEUTRAL
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_CENTR_WEAK
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_NEUTRAL
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_CENTR_WEAK
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done


#ok, awk the shit

mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK


cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_NEUTRAL",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_ALL",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_ALL_WEAK",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_ALL_WEAK_POS_MUT",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_ARM_WEAK_CENTER_STRONG",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"DOM_SEL_CENTR_WEAK",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_NEUTRAL",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_ALL",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_ALL_WEAK",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_ALL_WEAK_POS_MUT",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_ARM_WEAK_CENTR_STRONG",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,"NO_DOM_SEL_CENTR_WEAK",FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_CENTR_WEAK/$i; done


#make the windows file for the simulated genome for bedtools

sed -n '/I\t0\t10000/,/I\t2990000\t3000000/p' /projects/phillipslab/gavincw/repeats_12-18-18/29_prep_plots/repeatmasker_ii/11_bedtools/windows/inopinata.10kb.windows > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows

grep -v 'II' /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows.tmp

mv simulated_genome.10kb.windows.tmp simulated_genome.10kb.windows

#bedtools sort

module load bedtools/2.25.0


mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_NEUTRAL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/DOM_SEL_CENTR_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_CENTR_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_NEUTRAL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/02_awk/NO_DOM_SEL_CENTR_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_CENTR_WEAK/$i;done



#ok, get number of sites that have a  TE ; this I think is the thing that Anastasia is measuring.


mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_CENTR_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/03_bedtools_sort/NO_DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,"1"}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_CENTR_WEAK/$i;done




#bedtools map


mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_NEUTRAL
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL_WEAK
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/DOM_SEL_CENTR_WEAK
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_NEUTRAL
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL_WEAK
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/07_sites_awk/NO_DOM_SEL_CENTR_WEAK
for i in *; do bedtools map -o sum -c 8 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_CENTR_WEAK/$i; done

#awk

mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK



cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_NEUTRAL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/DOM_SEL_CENTR_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_NEUTRAL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/08_sites_bedtools_map/NO_DOM_SEL_CENTR_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_CENTR_WEAK/$i; done

#sed to replace . with zero and make tab delimited

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_NEUTRAL
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_CENTR_WEAK
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_NEUTRAL
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_CENTR_WEAK
for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done



#ok, cool maybe we have something now. cat!



mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/


cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_NEUTRAL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_NEUTRAL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_ALL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_ALL_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ALL_WEAK_POS_MUT
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_ALL_WEAK_POS_MUT
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_ARM_WEAK_CENTER_STRONG
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/DOM_SEL_CENTR_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/DOM_SEL_CENTR_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_NEUTRAL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_NEUTRAL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_ALL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_ALL_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_ALL_WEAK_POS_MUT
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/09_sites_awk/NO_DOM_SEL_CENTR_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/NO_DOM_SEL_CENTR_WEAK




cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/10_sites_cat/

for i in *; do echo -e "Chr\tBP\ttot_TE\tnorm_dist_center\tsim_id" | cat - $i > $i.tmp && mv $i.tmp $i; done
	#these files are in the folder "data", files combined to make early_dec_simulations_TE_sites.tsv

#newer evolutionary scenarios


wkdir="/projects/phillipslab/gavincw/repeats_12-18-18/"

mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019

mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING

cd /projects/phillipslab/shared/TE_SLiM/full_output

cp -v *_full_DOM_sel_all_WEAK_SEL_-0.0002.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0002
cp -v *_full_DOM_sel_all_WEAK_SEL_-0.0005.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0005
cp -v *_full_DOM_sel_all_WEAK_SEL_-0.0015.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0015
cp -v *_full_DOM_sel_all_WEAK_SEL_-0.001.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.001
cp -v *_full_DOM_sel_all_WEAK_SEL_-0.002.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.002
cp -v *_full_DOM_sel_all_WEAK_SELFING.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SELFING
cp -v *_full_DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
cp -v *_full_DOM_sel_arm_WEAK_center_Strong_SELFING.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_arm_WEAK_center_Strong_SELFING
cp -v *_full_no_dom_sel_all_weak_SELFING.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_all_weak_SELFING
cp -v *_full_no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
cp -v *_full_no_dom_sel_arm_weak_Center_Strong_SELFING.txt $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_arm_weak_Center_Strong_SELFING


mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/

cd  $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING


cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_all_WEAK_SELFING
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_all_weak_SELFING
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/00_cp/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do sed -n '/Mutations/,/Individuals/p' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done




cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SELFING
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_all_weak_SELFING
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do sed -i '1d'  $i; done
for i in *; do sed -i '$ d' $i; done

mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/

cd  $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING


cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SEL_-0.0002",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SEL_-0.0005",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SEL_-0.0015",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SEL_-0.001",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SEL_-0.002",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_all_WEAK_SELFING
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_all_WEAK_SELFING",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_sel_arm_WEAK_center_Strong_SELFING",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_all_weak_SELFING
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"no_dom_sel_all_weak_SELFING",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/01_sed/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"no_dom_sel_arm_weak_Center_Strong_SELFING",FILENAME,1}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done


#bedtools sort

module load bedtools/2.25.0


mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/


mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING


cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_all_WEAK_SELFING
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_all_weak_SELFING
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/02_awk/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do bedtools sort -i $i  > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done



mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SELFING
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_all_weak_SELFING
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do bedtools map -o sum -c 10 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done



mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING


cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_all_WEAK_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_all_weak_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_TE_sites/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done



mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0002
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SEL_-0.0002
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0005
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SEL_-0.0005
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.0015
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SEL_-0.0015
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.001
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SEL_-0.001
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SEL_-0.002
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SEL_-0.002
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_all_WEAK_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_all_WEAK_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/DOM_sel_arm_WEAK_center_Strong_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/DOM_sel_arm_WEAK_center_Strong_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_all_weak_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/no_dom_sel_all_weak_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/08_awk_TE_sites/no_dom_sel_arm_weak_Center_Strong_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/no_dom_sel_arm_weak_Center_Strong_SELFING


cd  $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/09_cat_TE_sites/


for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done

for i in *; do echo -e "Chr\tBP\ttot_TE\tnorm_dist_center\tsim_id" | cat - $i > $i.tmp && mv $i.tmp $i; done
	#files combined to make all_dec_27_simulations_TE_sites.tsv


#ages

mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_all_WEAK_SELFING
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_all_weak_SELFING
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/03_bedtools_sort/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done



mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/

mkdir DOM_sel_all_WEAK_SEL_-0.0002
mkdir DOM_sel_all_WEAK_SEL_-0.0005
mkdir DOM_sel_all_WEAK_SEL_-0.0015
mkdir DOM_sel_all_WEAK_SEL_-0.001
mkdir DOM_sel_all_WEAK_SEL_-0.002
mkdir DOM_sel_all_WEAK_SELFING
mkdir DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
mkdir DOM_sel_arm_WEAK_center_Strong_SELFING
mkdir no_dom_sel_all_weak_SELFING
mkdir no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
mkdir no_dom_sel_arm_weak_Center_Strong_SELFING


cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0002
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0005
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0005/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.0015
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0015/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.001
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.001/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SEL_-0.002
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.002/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_all_WEAK_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/DOM_sel_arm_WEAK_center_Strong_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_arm_WEAK_center_Strong_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_all_weak_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_all_weak_SELFING/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN/$i;done
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/07_bedtools_map_ages/no_dom_sel_arm_weak_Center_Strong_SELFING
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_arm_weak_Center_Strong_SELFING/$i;done



mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0002
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SEL_-0.0002
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0005
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SEL_-0.0005
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.0015
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SEL_-0.0015
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.001
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SEL_-0.001
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SEL_-0.002
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SEL_-0.002
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_all_WEAK_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_all_WEAK_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/DOM_sel_arm_WEAK_center_Strong_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/DOM_sel_arm_WEAK_center_Strong_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_all_weak_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/no_dom_sel_all_weak_SELFING
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN
cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/11_awk_ages/no_dom_sel_arm_weak_Center_Strong_SELFING
cat * > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/no_dom_sel_arm_weak_Center_Strong_SELFING


mkdir $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/13_awk_ages/

cd  $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/12_cat_ages/

for i in *; do awk '$3 != "."'  $i > $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/13_awk_ages/$i; done

cd $wkdir/69_SLiM_simulations_12-20-19/06_SIMULATIONS_DEC_27_2019/13_awk_ages/

for i in *; do sed -i -e 's/ /\t/g' $i; done

for i in *; do echo -e "Chr\tBP\tage\tnorm_dist_center\tsim_id" | cat - $i > $i.tmp && mv $i.tmp $i; done
	#these are in folder "data/ages"



#ages


mkdir  /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/
cd  /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK


cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_NEUTRAL",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_ALL",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_ALL_WEAK",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_ALL_WEAK_POS_MUT",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_ARM_WEAK_CENTER_STRONG",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"DOM_SEL_CENTR_WEAK",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_NEUTRAL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_NEUTRAL",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_ALL",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_ALL_WEAK",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_ALL_WEAK_POS_MUT",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_ARM_WEAK_CENTR_STRONG",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/01_sed/NO_DOM_SEL_CENTR_WEAK
for i in *; do awk 'BEGIN  {OFS="\t"} {print "I",$4+1,$4+1,$9,$3,50000-$8,$5,"NO_DOM_SEL_CENTR_WEAK",FILENAME}' $i >/projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_CENTR_WEAK/$i; done


#bedtools sort



mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_NEUTRAL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/DOM_SEL_CENTR_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_CENTR_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_NEUTRAL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_NEUTRAL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL_WEAK/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i;done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/00_awk/NO_DOM_SEL_CENTR_WEAK
for i in *; do bedtools sort -i $i  > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_CENTR_WEAK/$i;done


#bedtools map



mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_NEUTRAL
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL_WEAK
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/DOM_SEL_CENTR_WEAK
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_NEUTRAL
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL_WEAK
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/01_bedtools_sort/NO_DOM_SEL_CENTR_WEAK
for i in *; do bedtools map -o mean -c 6 -a /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/simulated_genome.10kb.windows -b $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_CENTR_WEAK/$i; done


#awk


mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/

cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/

mkdir DOM_NEUTRAL
mkdir DOM_SEL_ALL
mkdir DOM_SEL_ALL_WEAK
mkdir DOM_SEL_ALL_WEAK_POS_MUT
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG
mkdir DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
mkdir DOM_SEL_CENTR_WEAK
mkdir NO_DOM_NEUTRAL
mkdir NO_DOM_SEL_ALL
mkdir NO_DOM_SEL_ALL_WEAK
mkdir NO_DOM_SEL_ALL_WEAK_POS_MUT
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
mkdir NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
mkdir NO_DOM_SEL_CENTR_WEAK


cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_NEUTRAL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/DOM_SEL_CENTR_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_CENTR_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_NEUTRAL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_NEUTRAL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL_WEAK/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ALL_WEAK_POS_MUT
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS/$i; done
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/02_bedtools_map/NO_DOM_SEL_CENTR_WEAK
for i in *; do awk 'function abs(v) {return v < 0 ? -v : v} {print $1,$2+1,$4,(abs(1500000-$2)/1500000)/2,FILENAME}' $i > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_CENTR_WEAK/$i; done



mkdir /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/


cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_NEUTRAL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_NEUTRAL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_ALL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_ALL_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ALL_WEAK_POS_MUT
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_ALL_WEAK_POS_MUT
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_ARM_WEAK_CENTER_STRONG
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/DOM_SEL_CENTR_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/DOM_SEL_CENTR_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_NEUTRAL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_NEUTRAL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_ALL
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_ALL_WEAK
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ALL_WEAK_POS_MUT
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_ALL_WEAK_POS_MUT
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS
cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/03_awk/NO_DOM_SEL_CENTR_WEAK
cat * > /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat/NO_DOM_SEL_CENTR_WEAK




cd /projects/phillipslab/gavincw/repeats_12-18-18/69_SLiM_simulations_12-20-19/05_TE_POSITIONS/12_AGES_NEW_SIMULATIONS/04_cat

for i in *; do sed -i -e 's/ \. / 0 /g' $i; done
for i in *; do sed -i -e 's/ /\t/g' $i; done

for i in *; do echo -e "Chr\tBP\tage\tnorm_dist_center\tsim_id" | cat - $i > $i.tmp && mv $i.tmp $i; done

	#in "data/ages"


#on local
cd 09_cat_TE_sites 

for i in *; do awk 'BEGIN  {OFS="\t"} {print $0,FILENAME}' $i > /home/gavin/genome/genome/repeats_12-18-18/revisions/simulations_12-20-19/work_on_macbook_12-22-19/06_SIMULATIONS_DEC_27_2019/local_00_awk_TE_sites/$i;done

cd local_00_awk_TE_sites

for i in *; do sed -i '1d' $i; done

cat * > all_dec_27_simulations_TE_sites.tsv

echo -e "Chr\tBP\ttot_TE\tnorm_dist_center\tsim_id\tscenario" | cat - all_dec_27_simulations_TE_sites.tsv > all_dec_27_simulations_TE_sites.tsv.tmp && mv all_dec_27_simulations_TE_sites.tsv.tmp all_dec_27_simulations_TE_sites.tsv

cat all_dec_27_simulations_TE_sites.tsv early_dec_simulations_TE_sites.tsv > all_slim_te_site_data.tsv