#############
#############
#############
#############
#Revision: Extract Kimura distances (as a proxy for age) of repetitive elements
#############
#############
#############
#############

mkdir 63_RepeatMasker_align_file_output

wkdir="/projects/phillipslab/gavincw/repeats_12-18-18/"

#get ".align" files from repeatmasker with option -a this time.
cd 63_RepeatMasker_align_file_output

mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

RepeatMasker -a -s -lib $wkdir/24_fasta_filter/briggsae_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/briggsae_genome.fa
RepeatMasker -a -s -lib $wkdir/24_fasta_filter/elegans_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/elegans_genome.fa
RepeatMasker -a -s -lib $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/inopinata_genome.fa
RepeatMasker -a -s -lib $wkdir/24_fasta_filter/nigoni_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/nigoni_genome.fa
RepeatMasker -a -s -lib $wkdir/24_fasta_filter/remanei_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/remanei_genome.fa


cd $wkdir/
mkdir 65_custom_parse_RM_align
cd 65_custom_parse_RM_align

mkdir 00_links

ln -s  /projects/phillipslab/gavincw/repeats_12-18-18/63_RepeatMasker_align_file_output/inopinata/inopinata_genome.align /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/00_links/inopinata_genome.align

ln -s  /projects/phillipslab/gavincw/repeats_12-18-18/63_RepeatMasker_align_file_output/briggsae/briggsae_genome.fa.align /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/00_links/briggsae_genome.align

ln -s  /projects/phillipslab/gavincw/repeats_12-18-18/63_RepeatMasker_align_file_output/elegans/elegans_genome.fa.align /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/00_links/elegans_genome.align

ln -s  /projects/phillipslab/gavincw/repeats_12-18-18/63_RepeatMasker_align_file_output/nigoni/nigoni_genome.fa.align /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/00_links/nigoni_genome.align

ln -s  /projects/phillipslab/gavincw/repeats_12-18-18/63_RepeatMasker_align_file_output/remanei/remanei_genome.fa.align /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/00_links/remanei_genome.align



#split each alignment into its own file for extracting the information i want (genomic coordinates, cluster id, and kimura distance)
cd $wkdir/65_custom_parse_RM_align/
mkdir 01_csplit
cd 01_csplit
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei


cd $wkdir/65_custom_parse_RM_align/01_csplit/inopinata
csplit -f RMalign_ -z $wkdir/65_custom_parse_RM_align/00_links/inopinata_genome.align  /^[0-9]/ {*}

cd $wkdir/65_custom_parse_RM_align/01_csplit/elegans
csplit -f RMalign_ -z $wkdir/65_custom_parse_RM_align/00_links/elegans_genome.align  /^[0-9]/ {*}

cd $wkdir/65_custom_parse_RM_align/01_csplit/briggsae
csplit -f RMalign_ -z $wkdir/65_custom_parse_RM_align/00_links/briggsae_genome.align  /^[0-9]/ {*}

cd $wkdir/65_custom_parse_RM_align/01_csplit/nigoni
csplit -f RMalign_ -z $wkdir/65_custom_parse_RM_align/00_links/nigoni_genome.align  /^[0-9]/ {*}

cd $wkdir/65_custom_parse_RM_align/01_csplit/remanei
csplit -f RMalign_ -z $wkdir/65_custom_parse_RM_align/00_links/remanei_genome.align  /^[0-9]/ {*}



#copy only insertions that have kimura distances (ie, not low complexity or simple repeats or satellite DNA)
cd $wkdir/65_custom_parse_RM_align/
mkdir 02_cp_grep
cd 02_cp_grep
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/65_custom_parse_RM_align/01_csplit/inopinata/

for i in *; do 
 num_det=$(grep -c 'Kimura' $i)
 if (($num_det > 0)) ; then
 cp $i $wkdir/65_custom_parse_RM_align/02_cp_grep/inopinata/$i
 fi
done

cd $wkdir/65_custom_parse_RM_align/01_csplit/briggsae/

for i in *; do 
 num_det=$(grep -c 'Kimura' $i)
 if (($num_det > 0)) ; then
 cp $i $wkdir/65_custom_parse_RM_align/02_cp_grep/briggsae/$i
 fi
done

cd $wkdir/65_custom_parse_RM_align/01_csplit/elegans/

for i in *; do 
 num_det=$(grep -c 'Kimura' $i)
 if (($num_det > 0)) ; then
 cp $i $wkdir/65_custom_parse_RM_align/02_cp_grep/elegans/$i
 fi
done

cd $wkdir/65_custom_parse_RM_align/01_csplit/nigoni/

for i in *; do 
 num_det=$(grep -c 'Kimura' $i)
 if (($num_det > 0)) ; then
 cp $i $wkdir/65_custom_parse_RM_align/02_cp_grep/nigoni/$i
 fi
done

cd $wkdir/65_custom_parse_RM_align/01_csplit/remanei/

for i in *; do 
 num_det=$(grep -c 'Kimura' $i)
 if (($num_det > 0)) ; then
 cp $i $wkdir/65_custom_parse_RM_align/02_cp_grep/remanei/$i
 fi
done


#ok, now get first line and the cluster id and the coordinates for each alignment
cd $wkdir/65_custom_parse_RM_align/
mkdir 03_grep_ids_coordinates
cd 03_grep_ids_coordinates
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/inopinata/
for i in *; do grep -E '^[0-9]' $i | awk 'BEGIN {FS="[ ]"} {OFS="\t"} {print $5,$6,$7}' > $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/inopinata/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/briggsae/
for i in *; do grep -E '^[0-9]' $i | awk 'BEGIN {FS="[ ]"} {OFS="\t"} {print $5,$6,$7}' > $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/briggsae/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/elegans/
for i in *; do grep -E '^[0-9]' $i | awk 'BEGIN {FS="[ ]"} {OFS="\t"} {print $5,$6,$7}' > $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/elegans/$i; done


cd $wkdir/65_custom_parse_RM_align/02_cp_grep/remanei/
for i in *; do grep -E '^[0-9]' $i | awk 'BEGIN {FS="[ ]"} {OFS="\t"} {print $5,$6,$7}' > $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/remanei/$i; done


cd $wkdir/65_custom_parse_RM_align/02_cp_grep/nigoni/
for i in *; do grep -E '^[0-9]' $i | awk 'BEGIN {FS="[ ]"} {OFS="\t"} {print $5,$6,$7}' > $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/nigoni/$i; done


#now get the kimura distance for each alignment
cd $wkdir/65_custom_parse_RM_align/
mkdir 04_grep_kimura
cd 04_grep_kimura
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/inopinata/
for i in *; do grep 'Kimura' $i | sed -e 's/Kimura.*= //g' > $wkdir/65_custom_parse_RM_align/04_grep_kimura/inopinata/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/briggsae/
for i in *; do grep 'Kimura' $i | sed -e 's/Kimura.*= //g' > $wkdir/65_custom_parse_RM_align/04_grep_kimura/briggsae/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/elegans/
for i in *; do grep 'Kimura' $i | sed -e 's/Kimura.*= //g' > $wkdir/65_custom_parse_RM_align/04_grep_kimura/elegans/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/nigoni/
for i in *; do grep 'Kimura' $i | sed -e 's/Kimura.*= //g' > $wkdir/65_custom_parse_RM_align/04_grep_kimura/nigoni/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/remanei/
for i in *; do grep 'Kimura' $i | sed -e 's/Kimura.*= //g' > $wkdir/65_custom_parse_RM_align/04_grep_kimura/remanei/$i; done

#get the cluster id
cd $wkdir/65_custom_parse_RM_align/
mkdir 05_grep_id
cd 05_grep_id
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/inopinata/
for i in *; do grep -E '^[0-9]' $i | sed -e 's/.*Cluster/Cluster/g' | sed -e 's/ .*//g' > $wkdir/65_custom_parse_RM_align/05_grep_id/inopinata/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/briggsae/
for i in *; do grep -E '^[0-9]' $i | sed -e 's/.*Cluster/Cluster/g' | sed -e 's/ .*//g' > $wkdir/65_custom_parse_RM_align/05_grep_id/briggsae/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/elegans/
for i in *; do grep -E '^[0-9]' $i | sed -e 's/.*Cluster/Cluster/g' | sed -e 's/ .*//g' > $wkdir/65_custom_parse_RM_align/05_grep_id/elegans/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/nigoni/
for i in *; do grep -E '^[0-9]' $i | sed -e 's/.*Cluster/Cluster/g' | sed -e 's/ .*//g' > $wkdir/65_custom_parse_RM_align/05_grep_id/nigoni/$i; done

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/remanei/
for i in *; do grep -E '^[0-9]' $i | sed -e 's/.*Cluster/Cluster/g' | sed -e 's/ .*//g' > $wkdir/65_custom_parse_RM_align/05_grep_id/remanei/$i; done

#confirm all have same number of files....

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/inopinata/
ls -lh | wc -l
#250239

cd $wkdir/65_custom_parse_RM_align/04_grep_kimura/inopinata/
ls -lh | wc -l
#250239

cd $wkdir/65_custom_parse_RM_align/05_grep_id/inopinata/
ls -lh | wc -l
# 250239

cd $wkdir/65_custom_parse_RM_align/02_cp_grep/inopinata/
ls -lh | wc -l
#250239

#ok, good, now paste

#paste!
cd $wkdir/65_custom_parse_RM_align/
mkdir 06_paste
cd 06_paste
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/inopinata/
for i in *; do paste $i $wkdir/65_custom_parse_RM_align/04_grep_kimura/inopinata/$i $wkdir/65_custom_parse_RM_align/05_grep_id/inopinata/$i > $wkdir/65_custom_parse_RM_align/06_paste/inopinata/$i; done

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/briggsae/
for i in *; do paste $i $wkdir/65_custom_parse_RM_align/04_grep_kimura/briggsae/$i $wkdir/65_custom_parse_RM_align/05_grep_id/briggsae/$i > $wkdir/65_custom_parse_RM_align/06_paste/briggsae/$i; done

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/elegans/
for i in *; do paste $i $wkdir/65_custom_parse_RM_align/04_grep_kimura/elegans/$i $wkdir/65_custom_parse_RM_align/05_grep_id/elegans/$i > $wkdir/65_custom_parse_RM_align/06_paste/elegans/$i; done

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/nigoni/
for i in *; do paste $i $wkdir/65_custom_parse_RM_align/04_grep_kimura/nigoni/$i $wkdir/65_custom_parse_RM_align/05_grep_id/nigoni/$i > $wkdir/65_custom_parse_RM_align/06_paste/nigoni/$i; done

cd $wkdir/65_custom_parse_RM_align/03_grep_ids_coordinates/remanei/
for i in *; do paste $i $wkdir/65_custom_parse_RM_align/04_grep_kimura/remanei/$i $wkdir/65_custom_parse_RM_align/05_grep_id/remanei/$i > $wkdir/65_custom_parse_RM_align/06_paste/remanei/$i; done

#cat!

cd $wkdir/65_custom_parse_RM_align/
mkdir 07_cat

cd $wkdir/65_custom_parse_RM_align/06_paste/inopinata/
cat * > $wkdir/65_custom_parse_RM_align/07_cat/inopinata

cd $wkdir/65_custom_parse_RM_align/06_paste/briggsae/
cat * > $wkdir/65_custom_parse_RM_align/07_cat/briggsae

cd $wkdir/65_custom_parse_RM_align/06_paste/elegans/
cat * > $wkdir/65_custom_parse_RM_align/07_cat/elegans

cd $wkdir/65_custom_parse_RM_align/06_paste/nigoni/
cat * > $wkdir/65_custom_parse_RM_align/07_cat/nigoni

cd $wkdir/65_custom_parse_RM_align/06_paste/remanei/
cat * > $wkdir/65_custom_parse_RM_align/07_cat/remanei


#replace # with tab

cd $wkdir/65_custom_parse_RM_align/07_cat/
for i in *; do sed -i -e 's/#/\t/g' $i; done

#sort

mkdir $wkdir/65_custom_parse_RM_align/08_bedtools_sort/

module load bedtools/2.25.0

cd $wkdir/65_custom_parse_RM_align/07_cat/

bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i inopinata > $wkdir/65_custom_parse_RM_align/08_bedtools_sort/inopinata

bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i briggsae > $wkdir/65_custom_parse_RM_align/08_bedtools_sort/briggsae
bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i elegans > $wkdir/65_custom_parse_RM_align/08_bedtools_sort/elegans
bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i nigoni > $wkdir/65_custom_parse_RM_align/08_bedtools_sort/nigoni
bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i remanei > $wkdir/65_custom_parse_RM_align/08_bedtools_sort/remanei

#just _try_ bedtools intersect to see what happens....
#with the disjoined bed file

mkdir $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/

cd $wkdir/65_custom_parse_RM_align/08_bedtools_sort/

bedtools intersect -wao -a $wkdir/60_disjoined_gff/14_awk_bed/inopinata -b inopinata > $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/inopinata
bedtools intersect -wao -a $wkdir/60_disjoined_gff/14_awk_bed/briggsae -b briggsae > $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/briggsae
bedtools intersect -wao -a $wkdir/60_disjoined_gff/14_awk_bed/elegans -b elegans > $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/elegans
bedtools intersect -wao -a $wkdir/60_disjoined_gff/14_awk_bed/nigoni -b nigoni > $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/nigoni
bedtools intersect -wao -a $wkdir/60_disjoined_gff/14_awk_bed/remanei -b remanei > $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/remanei




#oh yeah... this is awesome

#get just the clusters that are the same in the final disjoined gff/bed and my kimura distances....

cd $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' inopinata > inopinata_col_4
sed -i -e 's/_.*//g' inopinata_col_4
paste inopinata inopinata_col_4 > inopinata.tmp

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' briggsae > briggsae_col_4
sed -i -e 's/_.*//g' briggsae_col_4
paste briggsae briggsae_col_4 > briggsae.tmp

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' elegans > elegans_col_4
sed -i -e 's/_.*//g' elegans_col_4
paste elegans elegans_col_4 > elegans.tmp

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' nigoni > nigoni_col_4
sed -i -e 's/_.*//g' nigoni_col_4
paste nigoni nigoni_col_4 > nigoni.tmp

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' remanei > remanei_col_4
sed -i -e 's/_.*//g' remanei_col_4
paste remanei remanei_col_4 > remanei.tmp

#keep if col 15 = col 18 !!!!

mkdir $wkdir/65_custom_parse_RM_align/10_awk/

cd $wkdir/65_custom_parse_RM_align/09_bedtools_intersect/

awk ' BEGIN {FS="\t"} {OFS="\t"} $15==$18 {print $0}' inopinata.tmp > $wkdir/65_custom_parse_RM_align/10_awk/inopinata
awk ' BEGIN {FS="\t"} {OFS="\t"} $15==$18 {print $0}' briggsae.tmp > $wkdir/65_custom_parse_RM_align/10_awk/briggsae
awk ' BEGIN {FS="\t"} {OFS="\t"} $15==$18 {print $0}' elegans.tmp > $wkdir/65_custom_parse_RM_align/10_awk/elegans
awk ' BEGIN {FS="\t"} {OFS="\t"} $15==$18 {print $0}' nigoni.tmp > $wkdir/65_custom_parse_RM_align/10_awk/nigoni
awk ' BEGIN {FS="\t"} {OFS="\t"} $15==$18 {print $0}' remanei.tmp > $wkdir/65_custom_parse_RM_align/10_awk/remanei


#i think we got it!!!

#Chr_gff	BP_start_gff	BP_end_gff	insertion_id	RM_classification_gff	repeat_class	repeat_subclass	repeat_order	repeat_superfamily	repeat_family	Chr_kimura_bed	BP_start_kimura_bed	BP_end_kimura_bed	kimura_distance	cluster_id_kimura_bed	RM_classification_kimura_bed	Num_bp_overlap	cluster_id_gff

#the above is the header for the files in $wkdir/65_custom_parse_RM_align/10_awk/ (what the columns mean)

#cat header

mkdir $wkdir/65_custom_parse_RM_align/11_header/

cd $wkdir/65_custom_parse_RM_align/10_awk/

for i in *; do echo -e "Chr\tBP\tbp_rep\tspecies\trep_class" | cat - $i > $i.tmp; done

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff" | cat - inopinata > $wkdir/65_custom_parse_RM_align/11_header/inopinata

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff" | cat - elegans > $wkdir/65_custom_parse_RM_align/11_header/elegans

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff" | cat - briggsae > $wkdir/65_custom_parse_RM_align/11_header/briggsae

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff" | cat - nigoni > $wkdir/65_custom_parse_RM_align/11_header/nigoni

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff" | cat - remanei > $wkdir/65_custom_parse_RM_align/11_header/remanei

#make a big tsv with everything...


mkdir $wkdir/65_custom_parse_RM_align/12_awk_species/

cd $wkdir/65_custom_parse_RM_align/10_awk/

rm header

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,FILENAME}' $i > $wkdir/65_custom_parse_RM_align/12_awk_species/$i; done

cd $wkdir/65_custom_parse_RM_align/12_awk_species/

cat * > all

#add header

echo -e "Chr_gff\tBP_start_gff\tBP_end_gff\tinsertion_id\tRM_classification_gff\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tChr_kimura_bed\tBP_start_kimura_bed\tBP_end_kimura_bed\tkimura_distance\tcluster_id_kimura_bed\tRM_classification_kimura_bed\tNum_bp_overlap\tcluster_id_gff\tspecies" | cat - all > all.tmp

mv all.tmp all



#ok, genomic landscapes.....


mkdir $wkdir/66_kimura_distance_genomic_landscapes/
mkdir $wkdir/66_kimura_distance_genomic_landscapes/00_links

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/all /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/all

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/briggsae /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/briggsae

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/elegans /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/elegans

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/inopinata /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/inopinata

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/nigoni /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/nigoni

ln -s /projects/phillipslab/gavincw/repeats_12-18-18/65_custom_parse_RM_align/12_awk_species/remanei /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/00_links/remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/
mkdir 01_class
mkdir 02_order
mkdir 03_superfamily
mkdir 04_family

cd 01_class

mkdir 00_awk
cd 00_awk
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#class $6, order $8, superfamily $9, family $10

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/briggsae/
awk '{print > $6}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/briggsae
sed -i -e 's/NA/Unknown/' NA
mv I class_I_retrotransposon
mv II class_II_DNA_transposon
mv NA Unknown

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/elegans/
awk '{print > $6}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/elegans
sed -i -e 's/NA/Unknown/' NA
mv I class_I_retrotransposon
mv II class_II_DNA_transposon
mv NA Unknown

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/inopinata/
awk '{print > $6}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/inopinata
sed -i -e 's/NA/Unknown/' NA
mv I class_I_retrotransposon
mv II class_II_DNA_transposon
mv NA Unknown

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/nigoni/
awk '{print > $6}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/nigoni
sed -i -e 's/NA/Unknown/' NA
mv I class_I_retrotransposon
mv II class_II_DNA_transposon
mv NA Unknown

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/remanei/
awk '{print > $6}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/remanei
sed -i -e 's/NA/Unknown/' NA
mv I class_I_retrotransposon
mv II class_II_DNA_transposon
mv NA Unknown

#ok, bedtools 

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei



cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/briggsae
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/elegans
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/inopinata
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/nigoni
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/00_awk/remanei
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/remanei/$i; done


mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"briggsae"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"elegans"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"inopinata"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"nigoni"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/01_bedtools_map/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"remanei"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk/remanei/$i; done






#ok cool..... just get the columns we want

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk//briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"class"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk//elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"class"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk//inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"class"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk//nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"class"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/02_sed_awk//remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"class"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk/remanei/$i; done

#ok cat species

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk//briggsae
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/briggsae
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk//elegans
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/elegans
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk//inopinata
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/inopinata
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk//nigoni
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/nigoni
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/04_awk//remanei
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/remanei

#cat them all

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/06_cat_all/
cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/05_cat/
cat * > $wkdir/66_kimura_distance_genomic_landscapes/01_class/06_cat_all/all


#awk by repeat type

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/07_awk_repeat_type

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/07_awk_repeat_type



awk '{print > $4}' $wkdir/66_kimura_distance_genomic_landscapes/01_class/06_cat_all/all

#add header

for i in *; do echo -e "Chr\tBP\tkimura_distance\trep_type\tspecies\trep_taxononmic_rank" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done

#make pdfs

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/08_genomic_landscape_rep_type_pdf/

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/01_class/distance_genomic_landscape_pdf.R $i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/08_genomic_landscape_rep_type_pdf/

pdfunite *.pdf genomic_landscape_distances_class.pdf



#orders

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order

mkdir 00_awk
cd 00_awk
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#class $6, order $8, superfamily $9, family $10

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/briggsae/
awk '{print > $8}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/briggsae
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/elegans/
awk '{print > $8}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/elegans
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/inopinata/
awk '{print > $8}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/inopinata
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/nigoni/
awk '{print > $8}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/nigoni
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/remanei/
awk '{print > $8}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/remanei
rm NA

#ok, bedtools 

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei



cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/briggsae
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/elegans
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/inopinata
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/nigoni
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/00_awk/remanei
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/remanei/$i; done


mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"briggsae"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"elegans"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"inopinata"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"nigoni"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/01_bedtools_map/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"remanei"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk/remanei/$i; done






#ok cool..... just get the columns we want

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk//briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"order"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk//elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"order"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk//inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"order"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk//nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"order"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/02_sed_awk//remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"order"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk/remanei/$i; done

#ok cat species

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk//briggsae
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/briggsae
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk//elegans
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/elegans
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk//inopinata
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/inopinata
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk//nigoni
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/nigoni
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/04_awk//remanei
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/remanei

#cat them all

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/06_cat_all/
cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/05_cat/
cat * > $wkdir/66_kimura_distance_genomic_landscapes/02_order/06_cat_all/all


#awk by repeat type

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/07_awk_repeat_type

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/07_awk_repeat_type



awk '{print > $4}' $wkdir/66_kimura_distance_genomic_landscapes/02_order/06_cat_all/all

#add header

for i in *; do echo -e "Chr\tBP\tkimura_distance\trep_type\tspecies\trep_taxononmic_rank" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/distance_genomic_landscape_pdf.R $wkdir/66_kimura_distance_genomic_landscapes/02_order/distance_genomic_landscape_pdf.R

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/08_genomic_landscape_rep_type_pdf/

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/02_order/distance_genomic_landscape_pdf.R $i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/08_genomic_landscape_rep_type_pdf/

pdfunite *.pdf genomic_landscape_distances_order.pdf


#superfamilies


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily

mkdir 00_awk
cd 00_awk
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#class $6, order $8, superfamily $9, family $10

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/briggsae/
awk '{print > $9}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/briggsae
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/elegans/
awk '{print > $9}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/elegans
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/inopinata/
awk '{print > $9}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/inopinata
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/nigoni/
awk '{print > $9}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/nigoni
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/remanei/
awk '{print > $9}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/remanei
rm NA

#ok, bedtools 

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei



cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/briggsae
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/elegans
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/inopinata
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/nigoni
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/00_awk/remanei
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/remanei/$i; done


mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"briggsae"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"elegans"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"inopinata"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"nigoni"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/01_bedtools_map/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"remanei"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk/remanei/$i; done






#ok cool..... just get the columns we want

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk//briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"superfamily"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk//elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"superfamily"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk//inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"superfamily"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk//nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"superfamily"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/02_sed_awk//remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"superfamily"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk/remanei/$i; done

#ok cat species

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk//briggsae
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/briggsae
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk//elegans
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/elegans
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk//inopinata
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/inopinata
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk//nigoni
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/nigoni
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/04_awk//remanei
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/remanei

#cat them all

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/06_cat_all/
cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/
cat * > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/06_cat_all/all


#awk by repeat type

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/07_awk_repeat_type

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/07_awk_repeat_type



awk '{print > $4}' $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/06_cat_all/all

#add header

for i in *; do echo -e "Chr\tBP\tkimura_distance\trep_type\tspecies\trep_taxononmic_rank" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/distance_genomic_landscape_pdf.R $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/distance_genomic_landscape_pdf.R

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/08_genomic_landscape_rep_type_pdf/

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/distance_genomic_landscape_pdf.R $i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/08_genomic_landscape_rep_type_pdf/

pdfunite *.pdf genomic_landscape_distances_superfamily.pdf



#families


cd $wkdir/66_kimura_distance_genomic_landscapes/04_family

mkdir 00_awk
cd 00_awk
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#class $6, order $8, superfamily $9, family $10

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/briggsae/
awk '{print > $10}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/briggsae
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/elegans/
awk '{print > $10}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/elegans
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/inopinata/
awk '{print > $10}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/inopinata
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/nigoni/
awk '{print > $10}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/nigoni
rm NA

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/remanei/
awk '{print > $10}' $wkdir/66_kimura_distance_genomic_landscapes/00_links/remanei
rm NA

#ok, bedtools 

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei



cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/briggsae
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/elegans
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/inopinata
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/nigoni
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/00_awk/remanei
for i in *; do bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/remanei/$i; done


mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"briggsae"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"elegans"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"inopinata"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"nigoni"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/01_bedtools_map/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $0,FILENAME,"remanei"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk/remanei/$i; done






#ok cool..... just get the columns we want

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk//briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"family"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/briggsae/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk//elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"family"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/elegans/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk//inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"family"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/inopinata/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk//nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"family"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/nigoni/$i; done
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/02_sed_awk//remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4,$5,$6,"family"}' $i > $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk/remanei/$i; done

#ok cat species

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk//briggsae
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/briggsae
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk//elegans
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/elegans
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk//inopinata
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/inopinata
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk//nigoni
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/nigoni
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/04_awk//remanei
cat * >  $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/remanei

#cat them all

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/06_cat_all/
cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/05_cat/
cat * > $wkdir/66_kimura_distance_genomic_landscapes/04_family/06_cat_all/all


#awk by repeat type

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/07_awk_repeat_type

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/07_awk_repeat_type



awk '{print > $4}' $wkdir/66_kimura_distance_genomic_landscapes/04_family/06_cat_all/all

#add header

for i in *; do echo -e "Chr\tBP\tkimura_distance\trep_type\tspecies\trep_taxononmic_rank" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/distance_genomic_landscape_pdf.R $wkdir/66_kimura_distance_genomic_landscapes/04_family/distance_genomic_landscape_pdf.R

mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/08_genomic_landscape_rep_type_pdf/

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/04_family/distance_genomic_landscape_pdf.R $i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/08_genomic_landscape_rep_type_pdf/

pdfunite *.pdf genomic_landscape_distances_family.pdf




#sina plots with x=species, y=kimura distances for all repeat taxa



mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/09_species_sina/


cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/01_class/sina_species.R  $i; done &

cd  $wkdir/66_kimura_distance_genomic_landscapes/01_class/09_species_sina/

pdfunite *.pdf class_species_sina.pdf






mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/09_species_sina/

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/sina_species.R  $wkdir/66_kimura_distance_genomic_landscapes/02_order/sina_species.R 

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/02_order/sina_species.R  $i; done &

cd  $wkdir/66_kimura_distance_genomic_landscapes/02_order/09_species_sina/

pdfunite *.pdf order_species_sina.pdf




mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/09_species_sina/

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/sina_species.R  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/sina_species.R 


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/sina_species.R  $i; done &


cd  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/09_species_sina/

pdfunite *.pdf superfamily_species_sina.pdf







mkdir $wkdir/66_kimura_distance_genomic_landscapes/04_family/09_species_sina/

cp $wkdir/66_kimura_distance_genomic_landscapes/01_class/sina_species.R  $wkdir/66_kimura_distance_genomic_landscapes/04_family/sina_species.R 


cd $wkdir/66_kimura_distance_genomic_landscapes/04_family/07_awk_repeat_type


for i in *; do Rscript $wkdir/66_kimura_distance_genomic_landscapes/04_family/sina_species.R  $i; done &


cd  $wkdir/66_kimura_distance_genomic_landscapes/04_family/09_species_sina/

pdfunite *.pdf family_species_sina.pdf






#norm dist chr center class



mkdir  $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/
cd  $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/briggsae

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/05_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/elegans

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/05_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/inopinata

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/05_cat/inopinata

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/nigoni

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/05_cat/nigoni

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/remanei

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/05_cat/remanei

mkdir $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/nigoni


cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/11_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/



cat * > all.tsv


echo -e "Chr\tBP\tkimura_distance\trep_type\tspecies\trep_taxononmic_rank\tnorm_dist_center" | cat - all.tsv > all.tsv.tmp && mv all.tsv.tmp all.tsv

sed -i -e 's/ /\t/g' all.tsv



#order



mkdir  $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/
cd  $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/briggsae

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/02_order/05_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/elegans

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/02_order/05_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/inopinata

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/02_order/05_cat/inopinata

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/nigoni

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/02_order/05_cat/nigoni

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/remanei

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/02_order/05_cat/remanei

mkdir $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/nigoni


cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/11_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/



cat * > all.tsv


#superfamily


mkdir  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/
cd  $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/briggsae

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/elegans

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/inopinata

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/inopinata

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/nigoni

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/nigoni

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/remanei

awk '{print > $1}' /projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/03_superfamily/05_cat/remanei

mkdir $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/nigoni


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/11_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/



cat * > all.tsv


#put all together....

cat $wkdir/66_kimura_distance_genomic_landscapes/01_class/12_cat/all.tsv $wkdir/66_kimura_distance_genomic_landscapes/02_order/12_cat/all.tsv $wkdir/66_kimura_distance_genomic_landscapes/03_superfamily/12_cat/all.tsv $wkdir/66_kimura_distance_genomic_landscapes/04_family/12_cat/all.tsv > $wkdir/66_kimura_distance_genomic_landscapes/kimura_distances.tsv

	#THIS IS THE FILE USED FOR MOST ANALYSES OF REPEAT AGE. Here, averages of Kimura distances of elements of the same type across 10kb windows.

#global landscape



mkdir $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape
mkdir $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/
cd $wkdir/66_kimura_distance_genomic_landscapes/00_links
bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b briggsae > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/briggsae
bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b elegans > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/elegans
bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b inopinata > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/inopinata
bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b nigoni > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/nigoni
bedtools map -o mean -c 14 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b remanei > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/remanei


mkdir $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/00_bedtools/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} $4!="."  {print $1,$2+1,$4,FILENAME}' $i > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/$i; done

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk

cat * > all


echo -e "Chr\tBP\tkimura_distance\tspecies" | cat - all > all.tmp

mv all.tmp kimura_distances_global_genomic_landscape.tsv

	#This is the global genomic landscape including all repetitive elements with Kimura distances in RepeatMasker output.



#global landscape with norm dist center

mkdir  $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/
cd  $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/briggsae

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/elegans

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/inopinata

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/inopinata

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/nigoni

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/nigoni

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/remanei

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/01_awk/remanei

mkdir $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/elegans

cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/nigoni


cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/02_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/05_global_distance_genomic_landscape/03_cat/



cat * > all.tsv


echo -e "Chr\tBP\tkimura_distance\tspecies\tnorm_dist_center" | cat - all.tsv > all.tmp

mv all.tmp kimura_distances_global_landscape_norm_dist_center.tsv


#clusters

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/
cd  $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir//54_cluster_distributions/04_size_distributions/01_awk_cluster/briggsae

for i in *; do awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/briggsae/$i; done &


cd $wkdir//54_cluster_distributions/04_size_distributions/01_awk_cluster/elegans

for i in *; do awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/elegans/$i; done &


cd $wkdir//54_cluster_distributions/04_size_distributions/01_awk_cluster/inopinata

for i in *; do awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/inopinata/$i; done &


cd $wkdir//54_cluster_distributions/04_size_distributions/01_awk_cluster/nigoni

for i in *; do awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/nigoni/$i; done &


cd $wkdir//54_cluster_distributions/04_size_distributions/01_awk_cluster/remanei

for i in *; do awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/remanei/$i; done &

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/
cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/briggsae/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/briggsae/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/elegans/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/elegans/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/inopinata/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/inopinata/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/nigoni/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/nigoni/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/remanei/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/remanei/$i; done &


mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/briggsae/

cat * > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/briggsae



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/elegans/

cat * > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/elegans


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/inopinata/

cat * > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/nigoni/

cat * > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/nigoni




cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/01_awk_filename/remanei/

cat * > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/remanei



#just grab number of clusters...

cd $wkdir/54_cluster_distributions/03_count_general_clusters/

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/02_grep_cluster_insertion_count/

for i in *; do grep Cluster $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/02_grep_cluster_insertion_count/$i; done &



#total bp covered....

#get the disjoined bed file ready

mkdir $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster
mkdir $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/00_awk

cd $wkdir/60_disjoined_gff/14_awk_bed/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' $i > $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/00_awk/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/00_awk/

for i in *; do sed -i -e  's/_.*//g' $i; done &

mkdir $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/01_paste/

cd $wkdir/60_disjoined_gff/14_awk_bed/

for i in *; do paste $i $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/00_awk/$i > $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/01_paste/$i; done &

#just the clusters

mkdir $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/02_grep/

cd $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/01_paste/

for i in *; do grep Cluster $i > $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/02_grep/$i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/00_links/00_disjoined_gff_w_cluster/02_grep/

for i in *; do sed -i -e  's/\t\t/\t/g' $i; done &



#total bp covered....

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk

cd $wkdir/66_kimura_distance_genomic_landscapes/00_links//00_disjoined_gff_w_cluster/02_grep/

#this gets length of each insertion in disjoined gff

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} function abs(v) {return v < 0 ? -v : v} {print $11, (abs($3-$2))}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/$i; done &

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk


#this breaks them up by cluster

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/

cd  $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/briggsae

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/briggsae


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/elegans

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/elegans


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/inopinata

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/inopinata



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/nigoni

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/nigoni



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/remanei

awk '{print > $1}' $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/00_awk/remanei




#get the sum


mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/

cd  $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/briggsae

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {sum+=$2;} END{print FILENAME,sum;}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/briggsae/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/elegans

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {sum+=$2;} END{print FILENAME,sum;}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/elegans/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/inopinata

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {sum+=$2;} END{print FILENAME,sum;}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/inopinata/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/nigoni

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {sum+=$2;} END{print FILENAME,sum;}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/nigoni/$i; done &


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/01_awk/remanei

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {sum+=$2;} END{print FILENAME,sum;}' $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/remanei/$i; done &



#put together and get percentage of genome...

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/briggsae/

cat * | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, (($2/105183150)*100)}' > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/briggsae



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/elegans/

cat * | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, (($2/100272607)*100)}' > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/elegans




cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/inopinata/

cat * | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, (($2/122994328)*100)}' > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/inopinata




cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/nigoni/

cat * | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, (($2/117831421)*100)}' > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/nigoni



cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/02_awk/remanei/

cat * | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, (($2/124856471)*100)}' > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/remanei



#get the repeat taxonomy for the clusters

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/04_taxonomy

cd  $wkdir/66_kimura_distance_genomic_landscapes/00_links//00_disjoined_gff_w_cluster/02_grep/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $11,$5,$6,$7,$8,$9,$10}' $i | sort | uniq > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/04_taxonomy/$i; done &

#header for summary age stats "rep_type	average	median	minimum	maximum	standard_deviation	iqr	num_alignments	species"

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/00_cluster_kimura_distance_means_local/cluster_summary_stats/

for i in *; do sed -i '1d' $i; done &

#they all have the same line count, woop!


#what will the header _actually_ be? "cluster_id	rm_classification	repeat_class	repeat_subclass	repeat_order	repeat_superfamily	repeat_family	average	median	minimum	maximum	standard_deviation	iqr	num_alignments	species	mean_insertion_length	insertion_count	total_genome_bp_cluster	total_genome_percent_cluster"

mkdir  $wkdir//66_kimura_distance_genomic_landscapes/06_cluster/05_merge


cd  $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/04_taxonomy


for i in *; do perl $wkdir/merge.pl -k -e "no_key" $i $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/00_cluster_kimura_distance_means_local/cluster_summary_stats/$i $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/01_mean_insertion_length_cluster/02_cat/$i $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/02_grep_cluster_insertion_count/$i $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/03_cluster_insertion_total_bp/03_cat/$i 2> $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/$i.merge_pl.error > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/05_merge/$i; done &


#add header 


cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/05_merge/

for i in *; do echo -e "cluster_id\trm_classification\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\taverage\tmedian\tminimum\tmaximum\tstandard_deviation\tiqr\tnum_alignments\tspecies\tmean_insertion_length\tinsertion_count\ttotal_genome_bp_cluster\ttotal_genome_percent_cluster" | cat - $i > $i.tmp && mv $i.tmp $i; done


#so... I'll just exclude those one-insertion clusters. I'm somehwat confused as to how they'd exist at all considering they are supposed to be repetitive. I guess it's cause they have homology to pre-existing repeats in the db.

mkdir $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/06_awk/

mv $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/06_awk/ $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/06_grep/

cd $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/05_merge/

for i in *; do grep -v "no_key" $i > $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/06_grep/$i; done 

cd  $wkdir/66_kimura_distance_genomic_landscapes/06_cluster/06_grep/

for i in *; do grep "no_key" $i; done 

#cool

#ok, put the  together

cat * > all

#remove the headers

grep -v "cluster_id" all > all.tmp

#put the header back on

echo -e "cluster_id\trm_classification\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\taverage\tmedian\tminimum\tmaximum\tstandard_deviation\tiqr\tnum_alignments\tspecies\tmean_insertion_length\tinsertion_count\ttotal_genome_bp_cluster\ttotal_genome_percent_cluster" | cat - all.tmp > all_cluster_insertion_age_abundance.tsv

	#This file is used for supplemental figures 25-27.


#discrete repeat count landscape

mkdir $wkdir/70_global_repeat_count/

cd $wkdir/70_global_repeat_count

mkdir 00_awk

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed

for i in *; do awk 'BEGIN {OFS="\t"} {print $0,1}' $i > $wkdir/70_global_repeat_count/00_awk/$i; done

cd $wkdir/70_global_repeat_count/00_awk

mv briggsae_mask_bed briggsae
mv elegans_mask_bed elegans
mv inopinata_mask_bed inopinata
mv nigoni_mask_bed nigoni
mv remanei_mask_bed remanei

cd $wkdir/70_global_repeat_count/
mkdir 01_bedtools_map

cd  $wkdir/70_global_repeat_count/00_awk/

bedtools map -o sum -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b briggsae > $wkdir/70_global_repeat_count/01_bedtools_map/briggsae

bedtools map -o sum -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b elegans > $wkdir/70_global_repeat_count/01_bedtools_map/elegans

bedtools map -o sum -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b inopinata > $wkdir/70_global_repeat_count/01_bedtools_map/inopinata

bedtools map -o sum -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b nigoni > $wkdir/70_global_repeat_count/01_bedtools_map/nigoni

bedtools map -o sum -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b remanei > $wkdir/70_global_repeat_count/01_bedtools_map/remanei


cd $wkdir/70_global_repeat_count/01_bedtools_map/

for i in *; do sed -i -e 's/\./0/g' $i; done

cd $wkdir/70_global_repeat_count/
mkdir 02_awk


cd $wkdir/70_global_repeat_count/01_bedtools_map/

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,FILENAME}' $i > $wkdir/70_global_repeat_count/02_awk/$i; done 

cd wkdir/70_global_repeat_count/02_awk/

cat * > all_te_count.tsv


echo -e "Chr\tBP\telement_count\tspecies" | cat - all_te_count.tsv > all_te_count.tsv.tmp
mv all_te_count.tsv.tmp all_te_count.tsv


#TE length landscape

mkdir $wkdir/71_global_repeat_length/

cd $wkdir/71_global_repeat_length

mkdir 00_awk

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed

for i in *; do awk 'BEGIN {OFS="\t"} {print $0,$3-$2}' $i > $wkdir/71_global_repeat_length/00_awk/$i; done

cd $wkdir/71_global_repeat_length/00_awk

mv briggsae_mask_bed briggsae
mv elegans_mask_bed elegans
mv inopinata_mask_bed inopinata
mv nigoni_mask_bed nigoni
mv remanei_mask_bed remanei

cd $wkdir/71_global_repeat_length/
mkdir 01_bedtools_map

cd  $wkdir/71_global_repeat_length/00_awk/

bedtools map -o mean -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/briggsae.10kb.windows -b briggsae > $wkdir/71_global_repeat_length/01_bedtools_map/briggsae

bedtools map -o mean -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/elegans.10kb.windows -b elegans > $wkdir/71_global_repeat_length/01_bedtools_map/elegans

bedtools map -o mean -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/inopinata.10kb.windows -b inopinata > $wkdir/71_global_repeat_length/01_bedtools_map/inopinata

bedtools map -o mean -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/nigoni.10kb.windows -b nigoni > $wkdir/71_global_repeat_length/01_bedtools_map/nigoni

bedtools map -o mean -c 5 -a $wkdir/33_remanei_prep_plots/11_bedtools/windows/remanei.10kb.windows -b remanei > $wkdir/71_global_repeat_length/01_bedtools_map/remanei

cd $wkdir/71_global_repeat_length/
mkdir 02_awk

cd $wkdir/71_global_repeat_length/01_bedtools_map/

for i in *; do awk '$4 != "."' $i > $wkdir/71_global_repeat_length/02_awk/$i ; done

cd $wkdir/71_global_repeat_length/
mkdir 03_awk


cd $wkdir/71_global_repeat_length/02_awk/

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,FILENAME}' $i > $wkdir/71_global_repeat_length/03_awk/$i; done 

cd wkdir/71_global_repeat_length/03_awk/

cat * > all_te_length.tsv


echo -e "Chr\tBP\telement_length\tspecies" | cat - all_te_length.tsv > all_te_length.tsv.tmp
mv all_te_length.tsv.tmp all_te_length.tsv

