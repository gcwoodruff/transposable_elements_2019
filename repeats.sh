#This is the code used to generate repeat masks and data for Woodruff and Teterina 2019.

#set working directory

wkdir="/projects/phillipslab/gavincw/repeats_12-18-18/"

cd $wkdir

#make script directory

mkdir $wkdir/scripts
	#copy scripts in "~/scripts" to this directory

#############
#############
#############
#############
#Find repetitive sequences in five Caenorhabditis genome assemblies
#############
#############
#############
#############

#make folder for assemblies, retrieve assemblies
mkdir 00_wget
cd 00_wget

wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_briggsae_CB4.scaffolds.fa.gz
wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_elegans_WBcel235.scaffolds.fa.gz
wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp34_NK74SC_v710.scaffolds.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/caenorhabditis_nigoni/PRJNA384657/caenorhabditis_nigoni.PRJNA384657.WBPS12.genomic.fa.gz

#C. remanei assembly can be found here: XXX

#unzip

gunzip *

#rename genomes

mv Caenorhabditis_briggsae_CB4.scaffolds.fa briggsae_genome.fa
mv Caenorhabditis_elegans_WBcel235.scaffolds.fa elegans_genome.fa
mv Caenorhabditis_sp34_NK74SC_v710.scaffolds.fa inopinata_genome.fa
mv caenorhabditis_nigoni.PRJNA384657.WBPS12.genomic.fa nigoni_genome.fa

#rename chromosomes to I,II,III,IV,V,X

sed -i -e 's/>I caenorhabditis_elegans_core_32_85_250 scaffold/>I/g' elegans_genome.fa
sed -i -e 's/>II caenorhabditis_elegans_core_32_85_250 scaffold/>II/g' elegans_genome.fa
sed -i -e 's/>III caenorhabditis_elegans_core_32_85_250 scaffold/>III/g' elegans_genome.fa
sed -i -e 's/>IV caenorhabditis_elegans_core_32_85_250 scaffold/>IV/g' elegans_genome.fa
sed -i -e 's/>V caenorhabditis_elegans_core_32_85_250 scaffold/>V/g' elegans_genome.fa
sed -i -e 's/>X caenorhabditis_elegans_core_32_85_250 scaffold/>X/g' elegans_genome.fa
sed -i -e 's/>I_caenorhabditis_briggsae_core_32_85_250_scaffold/>I/g' briggsae_genome.fa
sed -i -e 's/>II_caenorhabditis_briggsae_core_32_85_250_scaffold/>II/g' briggsae_genome.fa
sed -i -e 's/>III_caenorhabditis_briggsae_core_32_85_250_scaffold/>III/g' briggsae_genome.fa
sed -i -e 's/>IV_caenorhabditis_briggsae_core_32_85_250_scaffold/>IV/g' briggsae_genome.fa
sed -i -e 's/>V_caenorhabditis_briggsae_core_32_85_250_scaffold/>V/g' briggsae_genome.fa
sed -i -e 's/>X_caenorhabditis_briggsae_core_32_85_250_scaffold/>X/g' briggsae_genome.fa
sed -i -e 's/>CSP34.Sp34_Chr1 caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>I/g' inopinata_genome.fa
sed -i -e 's/>CSP34.Sp34_Chr2 caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>II/g' inopinata_genome.fa
sed -i -e 's/>CSP34.Sp34_Chr3 caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>III/g' inopinata_genome.fa
sed -i -e 's/>CSP34.Sp34_Chr4 caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>IV/g' inopinata_genome.fa
sed -i -e 's/>CSP34.Sp34_Chr5 caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>V/g' inopinata_genome.fa
sed -i -e 's/>CSP34.Sp34_ChrX caenorhabditis_sp34_NK74SC_v710_core_32_85_1 scaffold/>X/g' inopinata_genome.fa
sed -i -e 's/>CM008509.1 length=16740257/>I/g' nigoni_genome.fa
sed -i -e 's/>CM008510.1 length=19229515/>II/g' nigoni_genome.fa
sed -i -e 's/>CM008511.1 length=15535478/>III/g' nigoni_genome.fa
sed -i -e 's/>CM008512.1 length=20390332/>IV/g' nigoni_genome.fa
sed -i -e 's/>CM008513.1 length=22287381/>V/g' nigoni_genome.fa
sed -i -e 's/>CM008514.1 length=23648458/>X/g' nigoni_genome.fa

#extract just chromosomal sequences
	#copy file "chromosome_sequences.txt" and "fasta_filter.pl" to folder $wkdir/scripts/

#make folder for chromosome-only fasta files

mkdir $wkdir/01_fasta_filter

cd $wkdir/00_wget

#extract the chromosomal sequences for each assembly with fasta_filter.pl
	#ensure only genome fastas with renamed chromsomes are in folder
for i in *; do perl $wkdir/scripts/fasta_filter.pl $wkdir/scripts/chromosome_sequences.txt $i $wkdir/01_fasta_filter/$i; done


#building de novo repeat libraries with RepeatModeler. This requires RepeatModeler installation. We used version 1.0.11.

#This pipeline is inspired by Berriman M, Coghlan A, Tsai IJ. "Creation of a comprehensive repeat library for a newly sequenced parasitic worm genome" https://doi.org/10.1038/protex.2018.054

#build the repeatmodeler databases

mkdir $wkdir/02_repeatmodeler_builddatabase
cd $wkdir/02_repeatmodeler_builddatabase

BuildDatabase -name "elegans" $wkdir/01_fasta_filter/elegans_genome.fa
BuildDatabase -name "briggsae" $wkdir/01_fasta_filter/briggsae_genome.fa
BuildDatabase -name "inopinata" $wkdir/01_fasta_filter/inopinata_genome.fa
BuildDatabase -name "nigoni" $wkdir/01_fasta_filter/nigoni_genome.fa
BuildDatabase -name "remanei" $wkdir/01_fasta_filter/remanei_genome.fa

#create folders for repeatmodeler output

mkdir $wkdir/03_RepeatModeler
cd $wkdir/03_RepeatModeler
mkdir elegans
mkdir inopinata
mkdir briggsae
mkdir nigoni
mkdir remanei

#run repeatmodeler. This can take a long time. In practice this was done in parallel for each assembly.

cd $wkdir/03_RepeatModeler/briggsae
RepeatModeler -engine ncbi -database $wkdir/02_repeatmodeler_builddatabase/briggsae > briggsae.out

cd $wkdir/03_RepeatModeler/elegans
RepeatModeler -engine ncbi -database $wkdir/02_repeatmodeler_builddatabase/elegans > elegans.out

cd $wkdir/03_RepeatModeler/inopinata
RepeatModeler -engine ncbi -database $wkdir/02_repeatmodeler_builddatabase/inopinata > inopinata.out

cd $wkdir/03_RepeatModeler/nigoni
RepeatModeler -engine ncbi -database $wkdir/02_repeatmodeler_builddatabase/nigoni > nigoni.out

cd $wkdir/03_RepeatModeler/remanei
RepeatModeler -engine ncbi -database $wkdir/02_repeatmodeler_builddatabase/remanei > remanei.out


	#for those wondering, directory 04 was where the repeatmasker Rhabditida library was stored. It is now included in directory "additional_repeat_libraries"

#run TransposonPSI for DNA transposons. This requires TransposonPSI installation. We used version TransposonPSI_08222010

mkdir $wkdir/05_transposonPSI
cd $wkdir/05_transposonPSI
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/05_transposonPSI/briggsae
transposonPSI.pl $wkdir/01_fasta_filter/briggsae_genome.fa nuc

cd $wkdir/05_transposonPSI/elegans
transposonPSI.pl $wkdir/01_fasta_filter/elegans_genome.fa nuc

cd $wkdir/05_transposonPSI/inopinata
transposonPSI.pl $wkdir/01_fasta_filter/inopinata_genome.fa nuc

cd $wkdir/05_transposonPSI/nigoni
transposonPSI.pl $wkdir/01_fasta_filter/nigoni_genome.fa nuc

cd $wkdir/05_transposonPSI/remanei
transposonPSI.pl $wkdir/01_fasta_filter/remanei_genome.fa nuc



#now turn the results into a fasta and remove seqs less than 50 bp. This requires bedtools installation. We used version 2.25.0.

cd $wkdir/05_transposonPSI/briggsae

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' briggsae_genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 > briggsae_genome.fa.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/01_fasta_filter/briggsae_genome.fa -bed briggsae_genome.fa.TPSI.allHits.chains.bestPerLocus.bed -fo  briggsae.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' briggsae.TPSI.allHits.chains.bestPerLocus.fa  > briggsae.TPSI.allHits.chains.bestPerLocus_50bp.fa 

cd $wkdir/05_transposonPSI/elegans

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' elegans_genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 > elegans_genome.fa.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/01_fasta_filter/elegans_genome.fa -bed elegans_genome.fa.TPSI.allHits.chains.bestPerLocus.bed -fo  elegans.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' elegans.TPSI.allHits.chains.bestPerLocus.fa  > elegans.TPSI.allHits.chains.bestPerLocus_50bp.fa 

cd $wkdir/05_transposonPSI/inopinata

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' inopinata_genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 > inopinata_genome.fa.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/01_fasta_filter/inopinata_genome.fa -bed inopinata_genome.fa.TPSI.allHits.chains.bestPerLocus.bed -fo  inopinata.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' inopinata.TPSI.allHits.chains.bestPerLocus.fa  > inopinata.TPSI.allHits.chains.bestPerLocus_50bp.fa 

cd $wkdir/05_transposonPSI/nigoni

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' nigoni_genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 > nigoni_genome.fa.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/01_fasta_filter/nigoni_genome.fa -bed nigoni_genome.fa.TPSI.allHits.chains.bestPerLocus.bed -fo  nigoni.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' nigoni.TPSI.allHits.chains.bestPerLocus.fa  > nigoni.TPSI.allHits.chains.bestPerLocus_50bp.fa 

cd $wkdir/05_transposonPSI/remanei

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' remanei_genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 > remanei_genome.fa.TPSI.allHits.chains.bestPerLocus.bed
bedtools getfasta -name -fi $wkdir/01_fasta_filter/remanei_genome.fa -bed remanei_genome.fa.TPSI.allHits.chains.bestPerLocus.bed -fo  remanei.TPSI.allHits.chains.bestPerLocus.fa 
awk '!/^>/ { next } { getline seq } length(seq) >= 50 { print $0 "\n" seq }' remanei.TPSI.allHits.chains.bestPerLocus.fa  > remanei.TPSI.allHits.chains.bestPerLocus_50bp.fa 

#run genometools suffixerator for downstream ltrharvest and ltrdigest. This requires genometools installation. We used version 1.5.9.

mkdir $wkdir/06_gt_suffixerator

cd $wkdir/06_gt_suffixerator

gt suffixerator -db $wkdir/01_fasta_filter/briggsae_genome.fa -indexname briggsae_genome.fa -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db $wkdir/01_fasta_filter/elegans_genome.fa -indexname elegans_genome.fa -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db $wkdir/01_fasta_filter/inopinata_genome.fa -indexname inopinata_genome.fa -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db $wkdir/01_fasta_filter/nigoni_genome.fa -indexname nigoni_genome.fa -tis -suf -lcp -des -ssp -sds -dna
gt suffixerator -db $wkdir/01_fasta_filter/remanei_genome.fa -indexname remanei_genome.fa -tis -suf -lcp -des -ssp -sds -dna

#run ltrharvest for LTR retrotransposons

mkdir $wkdir/07_gt_ltrharvest
cd $wkdir/07_gt_ltrharvest

gt ltrharvest -index $wkdir/06_gt_suffixerator/briggsae_genome.fa -seqids yes -tabout no > briggsae.ltrharvest.out
gt ltrharvest -index $wkdir/06_gt_suffixerator/briggsae_genome.fa -gff3 briggsae_genome.fa.ltr.gff -out briggsae_genome.fa.ltr.fa
gt ltrharvest -index $wkdir/06_gt_suffixerator/elegans_genome.fa -seqids yes -tabout no > elegans.ltrharvest.out
gt ltrharvest -index $wkdir/06_gt_suffixerator/elegans_genome.fa -gff3 elegans_genome.fa.ltr.gff -out elegans_genome.fa.ltr.fa
gt ltrharvest -index $wkdir/06_gt_suffixerator/inopinata_genome.fa -seqids yes -tabout no > inopinata.ltrharvest.out
gt ltrharvest -index $wkdir/06_gt_suffixerator/inopinata_genome.fa -gff3 inopinata_genome.fa.ltr.gff -out inopinata_genome.fa.ltr.fa
gt ltrharvest -index $wkdir/06_gt_suffixerator/nigoni_genome.fa -seqids yes -tabout no > nigoni.ltrharvest.out
gt ltrharvest -index $wkdir/06_gt_suffixerator/nigoni_genome.fa -gff3 nigoni_genome.fa.ltr.gff -out nigoni_genome.fa.ltr.fa
gt ltrharvest -index $wkdir/06_gt_suffixerator/remanei_genome.fa -seqids yes -tabout no > remanei.ltrharvest.out
gt ltrharvest -index $wkdir/06_gt_suffixerator/remanei_genome.fa -gff3 remanei_genome.fa.ltr.gff -out remanei_genome.fa.ltr.fa

#sort gff ltrharvest output

mkdir $wkdir/08_gt_gff3_sort

cd $wkdir/07_gt_ltrharvest

for i in *.fa.ltr.gff; do gt gff3 -sort $i > $wkdir//08_gt_gff3_sort/${i%.gff}.s.gff; done

#get the pfam hmm's from Steinbiss, et al. 2009 for ltrdigest


#make directory for pfam domains
mkdir $wkdir/09_ltr_pfam_domains
cd $wkdir/09_ltr_pfam_domains
mkdir 00_wget
mkdir 01_hmmconvert

#download and name the domain hmm's
cd 00_wget

wget http://pfam.xfam.org/family/PF00075/hmm
mv hmm PF00075.hmm
wget http://pfam.xfam.org/family/PF00077/hmm
mv hmm PF00077.hmm
wget http://pfam.xfam.org/family/PF00078/hmm
mv hmm PF00078.hmm
wget http://pfam.xfam.org/family/PF00098/hmm
mv hmm PF00098.hmm
wget http://pfam.xfam.org/family/PF00385/hmm
mv hmm PF00385.hmm
wget http://pfam.xfam.org/family/PF00429/hmm
mv hmm PF00429.hmm
wget http://pfam.xfam.org/family/PF00516/hmm
mv hmm PF00516.hmm
wget http://pfam.xfam.org/family/PF00517/hmm
mv hmm PF00517.hmm
wget http://pfam.xfam.org/family/PF00552/hmm
mv hmm PF00552.hmm
wget http://pfam.xfam.org/family/PF00607/hmm
mv hmm PF00607.hmm
wget http://pfam.xfam.org/family/PF00665/hmm
mv hmm PF00665.hmm
wget http://pfam.xfam.org/family/PF00692/hmm
mv hmm PF00692.hmm
wget http://pfam.xfam.org/family/PF01021/hmm
mv hmm PF01021.hmm
wget http://pfam.xfam.org/family/PF01140/hmm
mv hmm PF01140.hmm
wget http://pfam.xfam.org/family/PF01141/hmm
mv hmm PF01141.hmm
wget http://pfam.xfam.org/family/PF01393/hmm
mv hmm PF01393.hmm
wget http://pfam.xfam.org/family/PF02022/hmm
mv hmm PF02022.hmm
wget http://pfam.xfam.org/family/PF02093/hmm
mv hmm PF02093.hmm
wget http://pfam.xfam.org/family/PF02337/hmm
mv hmm PF02337.hmm
wget http://pfam.xfam.org/family/PF03056/hmm
mv hmm PF03056.hmm
wget http://pfam.xfam.org/family/PF03078/hmm
mv hmm PF03078.hmm
wget http://pfam.xfam.org/family/PF03408/hmm
mv hmm PF03408.hmm
wget http://pfam.xfam.org/family/PF03732/hmm
mv hmm PF03732.hmm
wget http://pfam.xfam.org/family/PF04094/hmm
mv hmm PF04094.hmm
wget http://pfam.xfam.org/family/PF04195/hmm
mv hmm PF04195.hmm
wget http://pfam.xfam.org/family/PF05380/hmm
mv hmm PF05380.hmm
wget http://pfam.xfam.org/family/PF06815/hmm
mv hmm PF06815.hmm
wget http://pfam.xfam.org/family/PF06817/hmm
mv hmm PF06817.hmm
wget http://pfam.xfam.org/family/PF07253/hmm
mv hmm PF07253.hmm
wget http://pfam.xfam.org/family/PF07727/hmm
mv hmm PF07727.hmm
wget http://pfam.xfam.org/family/PF08284/hmm
mv hmm PF08284.hmm
wget http://pfam.xfam.org/family/PF08330/hmm
mv hmm PF08330.hmm
wget http://pfam.xfam.org/family/PF08791/hmm
mv hmm PF08791.hmm
wget http://pfam.xfam.org/family/PF09590/hmm
mv hmm PF09590.hmm

#use HMMER to convert the domains from HMMER3/f to HMMER2.0. requires hmmer installation. used hmmer version 3.1b2.

for i in *; do hmmconvert -2 $i > $wkdir/09_ltr_pfam_domains/01_hmmconvert/$i; done

#run ltrharvest for LTR retrotransposons. this requires hmmer.

cd $wkdir/10_gt_ltrdigest/briggsae
gt ltrdigest -hmms $wkdir/09_ltr_pfam_domains/01_hmmconvert/*.hmm -outfileprefix briggsae_ltrdigest $wkdir/08_gt_gff3_sort/briggsae_genome.fa.ltr.s.gff $wkdir/06_gt_suffixerator/briggsae_genome.fa >briggsae_ltrdigest_output_gff

cd $wkdir/10_gt_ltrdigest/elegans
gt ltrdigest -hmms $wkdir/09_ltr_pfam_domains/01_hmmconvert/*.hmm -outfileprefix elegans_ltrdigest $wkdir/08_gt_gff3_sort/elegans_genome.fa.ltr.s.gff $wkdir/06_gt_suffixerator/elegans_genome.fa >elegans_ltrdigest_output_gff

cd $wkdir/10_gt_ltrdigest/inopinata
gt ltrdigest -hmms $wkdir/09_ltr_pfam_domains/01_hmmconvert/*.hmm -outfileprefix inopinata_ltrdigest $wkdir/08_gt_gff3_sort/inopinata_genome.fa.ltr.s.gff $wkdir/06_gt_suffixerator/inopinata_genome.fa >inopinata_ltrdigest_output_gff

cd $wkdir/10_gt_ltrdigest/nigoni
gt ltrdigest -hmms $wkdir/09_ltr_pfam_domains/01_hmmconvert/*.hmm -outfileprefix nigoni_ltrdigest $wkdir/08_gt_gff3_sort/nigoni_genome.fa.ltr.s.gff $wkdir/06_gt_suffixerator/nigoni_genome.fa >nigoni_ltrdigest_output_gff

cd $wkdir/10_gt_ltrdigest/remanei
gt ltrdigest -hmms $wkdir/09_ltr_pfam_domains/01_hmmconvert/*.hmm -outfileprefix remanei_ltrdigest $wkdir/08_gt_gff3_sort/remanei_genome.fa.ltr.s.gff $wkdir/06_gt_suffixerator/remanei_genome.fa > remanei_ltrdigest_output_gff

#now, gt select to remove all LTR retrotransposon candidates that don't have any domain hit at all (to help get rid of things that might not be LTR retrotransposon insertions)


mkdir $wkdir/11_gt_select
cd $wkdir/11_gt_select

#this requires script filter_protein_match.lua ; this is from http://avrilomics.blogspot.com/2015/09/ltrharvest.html ; originally by Sascha Kastens; can be found in folder $wkdir/scripts/

gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < $wkdir/10_gt_ltrdigest/briggsae/briggsae_ltrdigest_output_gff > briggsae_ltrdigest_output_gff2
gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < $wkdir/10_gt_ltrdigest/elegans/elegans_ltrdigest_output_gff > elegans_ltrdigest_output_gff2
gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < $wkdir/10_gt_ltrdigest/inopinata/inopinata_ltrdigest_output_gff > inopinata_ltrdigest_output_gff2
gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < $wkdir/10_gt_ltrdigest/nigoni/nigoni_ltrdigest_output_gff > nigoni_ltrdigest_output_gff2
gt select -rule_files $wkdir/scripts/filter_protein_match.lua -- < $wkdir/10_gt_ltrdigest/remanei/remanei_ltrdigest_output_gff > remanei_ltrdigest_output_gff2


#get the ltr retrotransposon seqs. requires script get_ltr_retrotransposon_seqs.pl ; originally by Avril Coghlan ; this is from http://avrilomics.blogspot.com/2015/09/ltrharvest.html ; can be found in folder $wkdir/scripts/

mkdir $wkdir/12_get_ltr_retrotransposon_seqs
cd $wkdir/12_get_ltr_retrotransposon_seqs

perl -w $wkdir/scripts/get_ltr_retrotransposon_seqs.pl $wkdir/07_gt_ltrharvest/briggsae_genome.fa.ltr.fa $wkdir/11_gt_select/briggsae_ltrdigest_output_gff2 > briggsae_ltrdigest_output_gff2.fa
perl -w $wkdir/scripts/get_ltr_retrotransposon_seqs.pl $wkdir/07_gt_ltrharvest/elegans_genome.fa.ltr.fa $wkdir/11_gt_select/elegans_ltrdigest_output_gff2 > elegans_ltrdigest_output_gff2.fa
perl -w $wkdir/scripts/get_ltr_retrotransposon_seqs.pl $wkdir/07_gt_ltrharvest/inopinata_genome.fa.ltr.fa $wkdir/11_gt_select/inopinata_ltrdigest_output_gff2 > inopinata_ltrdigest_output_gff2.fa
perl -w $wkdir/scripts/get_ltr_retrotransposon_seqs.pl $wkdir/07_gt_ltrharvest/nigoni_genome.fa.ltr.fa $wkdir/11_gt_select/nigoni_ltrdigest_output_gff2 > nigoni_ltrdigest_output_gff2.fa
perl -w $wkdir/scripts/get_ltr_retrotransposon_seqs.pl $wkdir/07_gt_ltrharvest/remanei_genome.fa.ltr.fa $wkdir/11_gt_select/remanei_ltrdigest_output_gff2 > remanei_ltrdigest_output_gff2.fa

#now use usearch to generate non-redunant repeat libraries

#repeat library "Rhabditida.repeatmasker" retrieved from the Rhabditida-specific repeats in RepeatMasker.
	#Used this command to generate: $EBROOTREPEATMASKER/util/queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

#repeat library "ce_cb_repbase.fa" retreived from RepBase (23.03) by concatenating these files: cbrapp.ref celapp.ref cbrrep.ref celrep.ref
	#used this command to generate: cat cbrapp.ref celapp.ref cbrrep.ref celrep.ref > ce_cb_repbase.fa

#these additional repeat families are in folder "additional_repeat_families"

#prepare directories and make links

mkdir $wkdir/15_usearch
mkdir $wkdir/15_usearch/00_links
cd $wkdir/15_usearch/00_links
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

	#repeatmodeler

ln -s $wkdir/02_repeatmodeler_builddatabase/briggsae-families.fa $wkdir/15_usearch/00_links/briggsae/briggsae-families.fa
ln -s $wkdir/02_repeatmodeler_builddatabase/elegans-families.fa $wkdir/15_usearch/00_links/elegans/elegans-families.fa
ln -s $wkdir/02_repeatmodeler_builddatabase/inopinata-families.fa $wkdir/15_usearch/00_links/inopinata/inopinata-families.fa
ln -s $wkdir/02_repeatmodeler_builddatabase/nigoni-families.fa $wkdir/15_usearch/00_links/nigoni/nigoni-families.fa
ln -s $wkdir/02_repeatmodeler_builddatabase/remanei-families.fa $wkdir/15_usearch/00_links/remanei/remanei-families.fa

	#04_Rhabditida_repeatmasker

ln -s $wkdir/additional_repeat_families/Rhabditida.repeatmasker $wkdir/15_usearch/00_links/briggsae/Rhabditida.repeatmasker.fa
ln -s $wkdir/additional_repeat_families/Rhabditida.repeatmasker $wkdir/15_usearch/00_links/elegans/Rhabditida.repeatmasker.fa
ln -s $wkdir/additional_repeat_families/Rhabditida.repeatmasker $wkdir/15_usearch/00_links/inopinata/Rhabditida.repeatmasker.fa
ln -s $wkdir/additional_repeat_families/Rhabditida.repeatmasker $wkdir/15_usearch/00_links/nigoni/Rhabditida.repeatmasker.fa
ln -s $wkdir/additional_repeat_families/Rhabditida.repeatmasker $wkdir/15_usearch/00_links/remanei/Rhabditida.repeatmasker.fa

	#05_transposonPSI

ln -s $wkdir/05_transposonPSI/briggsae/briggsae.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/15_usearch/00_links/briggsae/briggsae.TPSI.allHits.chains.bestPerLocus_50bp.fa
ln -s $wkdir/05_transposonPSI/elegans/elegans.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/15_usearch/00_links/elegans/elegans.TPSI.allHits.chains.bestPerLocus_50bp.fa
ln -s $wkdir/05_transposonPSI/inopinata/inopinata.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/15_usearch/00_links/inopinata/inopinata.TPSI.allHits.chains.bestPerLocus_50bp.fa
ln -s $wkdir/05_transposonPSI/nigoni/nigoni.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/15_usearch/00_links/nigoni/nigoni.TPSI.allHits.chains.bestPerLocus_50bp.fa
ln -s $wkdir/05_transposonPSI/remanei/remanei.TPSI.allHits.chains.bestPerLocus_50bp.fa $wkdir/15_usearch/00_links/remanei/remanei.TPSI.allHits.chains.bestPerLocus_50bp.fa

	#12_get_ltr_retrotransposon_seqs (ie, ltrdigest and ltrharvest)

ln -s $wkdir/12_get_ltr_retrotransposon_seqs/briggsae_ltrdigest_output_gff2.fa $wkdir/15_usearch/00_links/briggsae/briggsae_ltrdigest_output_gff2.fa
ln -s $wkdir/12_get_ltr_retrotransposon_seqs/elegans_ltrdigest_output_gff2.fa $wkdir/15_usearch/00_links/elegans/elegans_ltrdigest_output_gff2.fa
ln -s $wkdir/12_get_ltr_retrotransposon_seqs/inopinata_ltrdigest_output_gff2.fa $wkdir/15_usearch/00_links/inopinata/inopinata_ltrdigest_output_gff2.fa
ln -s $wkdir/12_get_ltr_retrotransposon_seqs/nigoni_ltrdigest_output_gff2.fa $wkdir/15_usearch/00_links/nigoni/nigoni_ltrdigest_output_gff2.fa
ln -s $wkdir/12_get_ltr_retrotransposon_seqs/remanei_ltrdigest_output_gff2.fa $wkdir/15_usearch/00_links/remanei/remanei_ltrdigest_output_gff2.fa

	# 14_repbase

ln -s $wkdir/additional_repeat_families/ce_cb_repbase.fa $wkdir/15_usearch/00_links/briggsae/ce_cb_repbase.fa
ln -s $wkdir/additional_repeat_families/ce_cb_repbase.fa $wkdir/15_usearch/00_links/elegans/ce_cb_repbase.fa
ln -s $wkdir/additional_repeat_families/ce_cb_repbase.fa $wkdir/15_usearch/00_links/inopinata/ce_cb_repbase.fa
ln -s $wkdir/additional_repeat_families/ce_cb_repbase.fa $wkdir/15_usearch/00_links/nigoni/ce_cb_repbase.fa
ln -s $wkdir/additional_repeat_families/ce_cb_repbase.fa $wkdir/15_usearch/00_links/remanei/ce_cb_repbase.fa

#combine the libraries

mkdir $wkdir/15_usearch/01_cat

cd $wkdir/15_usearch/00_links/briggsae
cat * > $wkdir/15_usearch/01_cat/briggsae_COMBO_repeats_unclass.fasta

cd $wkdir/15_usearch/elegans
cat * > $wkdir/15_usearch/01_cat/elegans_COMBO_repeats_unclass.fasta

cd $wkdir/15_usearch/inopinata
cat * > $wkdir/15_usearch/01_cat/inopinata_COMBO_repeats_unclass.fasta

cd $wkdir/15_usearch/nigoni
cat * > $wkdir/15_usearch/01_cat/nigoni_COMBO_repeats_unclass.fasta

cd $wkdir/15_usearch/remanei
cat * > $wkdir/15_usearch/01_cat/remanei_COMBO_repeats_unclass.fasta

#run usearch to get non-redunant repeat libraries. requires usearch installation. we used version 8.0 .

mkdir $wkdir/15_usearch/02_usearch
cd $wkdir/15_usearch/02_usearch
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/15_usearch/02_usearch/briggsae
usearch -cluster_fast $wkdir/15_usearch/01_cat/briggsae_COMBO_repeats_unclass.fasta -id 0.8 -centroids briggsae_COMBO_rep_centroids.fasta -uc briggsae_clusters.uc -consout briggsae_repeats_clust_80_unclass.fasta -msaout briggsae_clust_aligned_80.fasta

cd $wkdir/15_usearch/02_usearch/elegans
usearch -cluster_fast $wkdir/15_usearch/01_cat/elegans_COMBO_repeats_unclass.fasta -id 0.8 -centroids elegans_COMBO_rep_centroids.fasta -uc elegans_clusters.uc -consout elegans_repeats_clust_80_unclass.fasta -msaout elegans_clust_aligned_80.fasta

cd $wkdir/15_usearch/02_usearch/inopinata
usearch -cluster_fast $wkdir/15_usearch/01_cat/inopinata_COMBO_repeats_unclass.fasta -id 0.8 -centroids inopinata_COMBO_rep_centroids.fasta -uc inopinata_clusters.uc -consout inopinata_repeats_clust_80_unclass.fasta -msaout inopinata_clust_aligned_80.fasta

cd $wkdir/15_usearch/02_usearch/nigoni
usearch -cluster_fast $wkdir/15_usearch/01_cat/nigoni_COMBO_repeats_unclass.fasta -id 0.8 -centroids nigoni_COMBO_rep_centroids.fasta -uc nigoni_clusters.uc -consout nigoni_repeats_clust_80_unclass.fasta -msaout nigoni_clust_aligned_80.fasta

cd $wkdir/15_usearch/02_usearch/remanei
usearch -cluster_fast $wkdir/15_usearch/01_cat/remanei_COMBO_repeats_unclass.fasta -id 0.8 -centroids remanei_COMBO_rep_centroids.fasta -uc remanei_clusters.uc -consout remanei_repeats_clust_80_unclass.fasta -msaout remanei_clust_aligned_80.fasta


#next run repeatclassifier, but first, split the files to run in parallel. requires script fasta-splitter.pl by Kirill Kryukov ; can be found in folder $wkdir/scripts/

mkdir $wkdir/16_fasta-splitter
cd $wkdir/16_fasta-splitter
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --out-dir $wkdir/16_fasta-splitter/briggsae/ --nopad --measure count $wkdir/15_usearch/02_usearch/briggsae/briggsae_repeats_clust_80_unclass.fasta

perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --out-dir $wkdir/16_fasta-splitter/elegans/ --nopad --measure count $wkdir/15_usearch/02_usearch/elegans/elegans_repeats_clust_80_unclass.fasta

perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --out-dir $wkdir/16_fasta-splitter/inopinata/ --nopad --measure count $wkdir/15_usearch/02_usearch/inopinata/inopinata_repeats_clust_80_unclass.fasta

perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --out-dir $wkdir/16_fasta-splitter/nigoni/ --nopad --measure count $wkdir/15_usearch/02_usearch/nigoni/nigoni_repeats_clust_80_unclass.fasta

perl $wkdir/scripts/fasta-splitter.pl --n-parts 250 --line-length 0 --out-dir $wkdir/16_fasta-splitter/remanei/ --nopad --measure count $wkdir/15_usearch/02_usearch/remanei/remanei_repeats_clust_80_unclass.fasta

#run repeatclassifier. requires repeatmodeler. we used version 1.0.11

cd $wkdir/16_fasta-splitter/briggsae/
for i in *; do RepeatClassifier -consensi $i -engine ncbi; done


cd $wkdir/16_fasta-splitter/elegans/
for i in *; do RepeatClassifier -consensi $i -engine ncbi; done


cd $wkdir/16_fasta-splitter/inopinata/
for i in *; do RepeatClassifier -consensi $i -engine ncbi; done

cd $wkdir/16_fasta-splitter/nigoni/
for i in *; do RepeatClassifier -consensi $i -engine ncbi; done


cd $wkdir/16_fasta-splitter/remanei/
for i in *; do RepeatClassifier -consensi $i -engine ncbi; done

#combine classified repeats

mkdir $wkdir/17_cat_RepeatClassifer

cd $wkdir/16_fasta-splitter/briggsae
cat *.fasta.classified >  $wkdir/17_cat_RepeatClassifer/briggsae_combo_repeats_classified.fa

cd $wkdir/16_fasta-splitter/elegans
cat *.fasta.classified >  $wkdir/17_cat_RepeatClassifer/elegans_combo_repeats_classified.fa

cd $wkdir/16_fasta-splitter/inopinata
cat *.fasta.classified >  $wkdir/17_cat_RepeatClassifer/inopinata_combo_repeats_classified.fa

cd $wkdir/16_fasta-splitter/nigoni
cat *.fasta.classified >  $wkdir/17_cat_RepeatClassifer/nigoni_combo_repeats_classified.fa

cd $wkdir/16_fasta-splitter/remanei/
cat *.fasta.classified >  $wkdir/17_cat_RepeatClassifer/remanei_combo_repeats_classified.fa

#next, will align unclassified repeats to predicted proteins to exclude such sequences from being called "unknown" repetitive elements. this requires ncbi toolkit. we used version 18.0.0. Proteins are in directory $wkdir/protein

#make blast db's

mkdir $wkdir/18_makeblastdb

cd $wkdir/18_makeblastdb


makeblastdb -in $wkdir/protein/sp34.fa -parse_seqids -out $wkdir/18_makeblastdb/inopinata -title inopinata_prot -dbtype prot

makeblastdb -in $wkdir/protein/elegans.fa -parse_seqids -out $wkdir/18_makeblastdb/elegans -title elegans_prot -dbtype prot

makeblastdb -in $wkdir/protein/briggsae.fa -parse_seqids -out $wkdir/18_makeblastdb/briggsae -title briggsae_prot -dbtype prot

makeblastdb -in $wkdir/protein/nigoni.fa -parse_seqids -out $wkdir/18_makeblastdb/nigoni -title nigoni_prot -dbtype prot

makeblastdb -in $wkdir/protein/remanei.fa -parse_seqids -out $wkdir/18_makeblastdb/remanei -title remanei_prot -dbtype prot

#blast repeats to prot database to find overlap. this requires ncbi toolkit. we used version 18.0.0.

mkdir $wkdir/19_blastx

cd $wkdir/19_blastx

blastx -db $wkdir/18_makeblastdb/briggsae -num_threads 24 -outfmt 6 -query $wkdir/17_cat_RepeatClassifer/briggsae_combo_repeats_classified.fa -out BLAST_clust80_briggsae.res.txt -evalue 0.001

blastx -db $wkdir/18_makeblastdb/elegans -num_threads 24 -outfmt 6 -query $wkdir/17_cat_RepeatClassifer/elegans_combo_repeats_classified.fa -out BLAST_clust80_elegans.res.txt -evalue 0.001

blastx -db $wkdir/18_makeblastdb/inopinata -num_threads 24 -outfmt 6 -query $wkdir/17_cat_RepeatClassifer/inopinata_combo_repeats_classified.fa -out BLAST_clust80_inopinata.res.txt -evalue 0.001

blastx -db $wkdir/18_makeblastdb/nigoni -num_threads 24 -outfmt 6 -query $wkdir/17_cat_RepeatClassifer/nigoni_combo_repeats_classified.fa -out BLAST_clust80_nigoni.res.txt -evalue 0.001

blastx -db $wkdir/18_makeblastdb/remanei -num_threads 24 -outfmt 6 -query $wkdir/17_cat_RepeatClassifer/remanei_combo_repeats_classified.fa -out BLAST_clust80_remanei.res.txt -evalue 0.001


#get list of proteins to exclude

mkdir $wkdir/20_cut

cd $wkdir/19_blastx

#these are the repeats that blast to proteins

for i in *; do cut -f1 $i |uniq > $wkdir/20_cut/${i%.res.txt}_list_of_prot_match; done

#these are the repeats that blast to proteins and are unknown

mkdir $wkdir/21_grep

cd $wkdir/20_cut

for i in *; do grep "Unknown" $i | sort > $wkdir/21_grep/${i%_list_of_prot_match}_match_prot_and_unknown; done

#these are lists of all the repeats

mkdir $wkdir/22_grep

cd $wkdir/17_cat_RepeatClassifer

for i in *; do grep ">" $i | sed -e 's/>//g' | sort > $wkdir/22_grep/${i%.fa}_list_all; done

#these are repeats that exclude those that match proteins and are unknown

mkdir $wkdir/23_comm

cd $wkdir/22_grep

comm -23 <(grep -Po '\S+' briggsae_combo_repeats_classified_list_all | sort) <(grep -Po '\S+' $wkdir/21_grep/BLAST_clust80_briggsae_match_prot_and_unknown | sort)  > $wkdir/23_comm/briggsae_repeat_list

comm -23 <(grep -Po '\S+' elegans_combo_repeats_classified_list_all | sort) <(grep -Po '\S+' $wkdir/21_grep/BLAST_clust80_elegans_match_prot_and_unknown | sort)  > $wkdir/23_comm/elegans_repeat_list

comm -23 <(grep -Po '\S+' inopinata_combo_repeats_classified_list_all | sort) <(grep -Po '\S+' $wkdir/21_grep/BLAST_clust80_inopinata_match_prot_and_unknown | sort)  > $wkdir/23_comm/inopinata_repeat_list

comm -23 <(grep -Po '\S+' nigoni_combo_repeats_classified_list_all | sort) <(grep -Po '\S+' $wkdir/21_grep/BLAST_clust80_nigoni_match_prot_and_unknown | sort)  > $wkdir/23_comm/nigoni_repeat_list

comm -23 <(grep -Po '\S+' remanei_combo_repeats_classified_list_all | sort) <(grep -Po '\S+' $wkdir/21_grep/BLAST_clust80_remanei_match_prot_and_unknown | sort)  > $wkdir/23_comm/remanei_repeat_list

#extract those sequences that exclude those that match proteins and are unknown. thanks Kevin Nyberg for fasta_filter.pl

mkdir $wkdir/24_fasta_filter

perl $wkdir/scripts/fasta_filter.pl $wkdir/23_comm/briggsae_repeat_list $wkdir/17_cat_RepeatClassifer/briggsae_combo_repeats_classified.fa $wkdir/24_fasta_filter/briggsae_combo_repeats_classified_no_prot.fa

perl $wkdir/scripts/fasta_filter.pl $wkdir/23_comm/elegans_repeat_list $wkdir/17_cat_RepeatClassifer/elegans_combo_repeats_classified.fa $wkdir/24_fasta_filter/elegans_combo_repeats_classified_no_prot.fa

perl $wkdir/scripts/fasta_filter.pl $wkdir/23_comm/inopinata_repeat_list $wkdir/17_cat_RepeatClassifer/inopinata_combo_repeats_classified.fa $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa

perl $wkdir/scripts/fasta_filter.pl $wkdir/23_comm/nigoni_repeat_list $wkdir/17_cat_RepeatClassifer/nigoni_combo_repeats_classified.fa $wkdir/24_fasta_filter/nigoni_combo_repeats_classified_no_prot.fa

perl $wkdir/scripts/fasta_filter.pl $wkdir/23_comm/remanei_repeat_list $wkdir/17_cat_RepeatClassifer/remanei_combo_repeats_classified.fa $wkdir/24_fasta_filter/remanei_combo_repeats_classified_no_prot.fa

#not all sequences were extracted faithfully ; here are the sequences that were added by hand

#this sequence was added to $wkdir/24_fasta_filter/briggsae_combo_repeats_classified_no_prot.fa :

#>Cluster570#ARTEFACT 
#GGTGATGCTGCCAACTTACTGATTTAGTGTATGATGGTGTTTTTGAGGTG
#CTCCAGTGGCTTCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGAC
#GGGGTGGTGCGTAACGGCAAAAGCACCGCCGGACATCAGCGCTATCTCTG
#CTCTCACTGCCGTAAAACATGGCAACTGCAGTTCACTTACACCGCTTCTC
#AACCCGGTACGCACCAGAAAATCATTGATATGGCCATGAATGGCGTTGGA
#TGCCGGGCAACAGCCCGCATTATGGGCGTTGGCCTCAACACGATTTTACG
#TCACTTAAAAAACTCAGGCCGCAGTCGGTAACCTCGCGCATACAGCCGGG
#CAGTGACGTCATCGTCTGCGCGGAAATGGACGAACAGTGGGGCTATGTCG
#GGGCTAAATCGCGCCAGCGCTGGCTGTTTTACGCGTATGACAGTCTCCGG
#AAGACGGTTGTTGCGCACGTATTCGGTGAACGCACTATGGCGACGCTGGG
#GCGTCTTATGAGCCTGCTGTCACCCTTTGACGTGGTGATATGGATGACGG
#ATGGCTGGCCGCTGTATGAATCCCGCCTGAAGGGAAAGCTGCACGTAATC
#AGCAAGCGATATACGCAGCGAATTGAGCGGCATAACCTGAATCTGAGGCA
#GCACCTGGCACGGCTGGGACGGAAGTCGCTGTCGTTCTCAAAATCGGTGG
#AGCTGCATGACAAAGTCATCGGGCATTATCTGAACATAAAACACTATCAA
#TAAGTTGGAGTCATTACC

#this sequence was added to $wkdir/24_fasta_filter/elegans_combo_repeats_classified_no_prot.fa

#>Cluster550#DNA 
#TGGGGTTATTCATGACATAAAAAATGTATATTTTTTTGGAATTTTTAAAC
#ACAAAAAATTAATGCATGCATGCGAAAAAATGTATTTCATGTATTTTTTC
#GTATTTTTTTGTATTTAAAAAACAAATCATGCATAACCCCA

#this sequence was added to $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa

#>Cluster1170#DNA/MULE-MuDR 
#TAATTTCAAAAAAACATAGAAAGTAAAATTATTTCAGGCAGAAATTATTG
#TGTTATTACTAACGAAAACCAAACAAAACAAATAATGCGAAACGTAAAAT
#CTAGCCTTGGTATTTTTGTATGAACAATAAGTGTTCAACAGATCGAGCAA
#TTATAAACTATTTTTCAAGAGAATTTGAGTCTTTTTGGAGGTGCAACGTT
#TCTCATACTCACAGAGTACACAGGGACTTCTTTTGATTCTTTCATTTAAA
#ATCTATTCAATGATTGATTCTGAATCAGTGGGTGTTCTATCTAAATTCTC
#ATCCTGGTCAGTTTTGTATTCAGGCTAAACGGAAAAAACACAGAATTAAT
#ATGATAAAAACATATTAAACGAATCATATTCAATCACGAATACTATTTCT
#GTTTCAAAGGAGAATACTTCTTAAAATTTATATTTTCTGCTATTCTTCGC
#AGTAACGATTGATAAACGTAAGAATAAAAAACTCTTGTTTGTTCATATTA
#ATTGAAACCCTAAGGCCTACTTGTGGTACCCCCGAATCATAATTGTGGTA
#CCCCCAAATAATATCGTTGGTACCCCCTGAATCATAGCGTCATGGTACCC
#CTCTGTCATTGTACCCTAAAAATATAGTGCCCCCGAATAATAGTACCCCC
#CGAATCATGAACTTGGTACCCCCGAATCGTTTGTACCCCGTTGTACCCCC
#GAATCATTGTACCCCCCGAATTAAATGCC

#this sequence was added to $wkdir/24_fasta_filter/remanei_combo_repeats_classified_no_prot.fa

#>Cluster870#DNA/hAT 
#GGTACTTATGGGTTTCGTTCCCCCCAAAATGTTCATTCAATTATTTAATA
#CTGAATTTTTAATTTTAATCCACACGTGAAAGTTTATTTTAATACTGTTT
#TCATTTTCAGGCTTAGGAAACCATCTTCCTAAGCCTGAGAATGAAAAAAA
#AGTTCACGTCTCTATTAAAAACTATGGAAAACATACTGGGGGGAACGAAA
#CCCATAAGTACC

#run of RepeatMasker with the previously generated repeat libraries. repeatmodeler version 1.0.11.

mkdir $wkdir/27_RepeatMasker_II
cd $wkdir/27_RepeatMasker_II
mkdir hard
mkdir soft

cd $wkdir/27_RepeatMasker_II/hard/

RepeatMasker -s -lib $wkdir/24_fasta_filter/briggsae_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/briggsae_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/elegans_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/elegans_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/inopinata_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/nigoni_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/nigoni_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/remanei_combo_repeats_classified_no_prot.fa  -gff -pa 16  $wkdir/01_fasta_filter/remanei_genome.fa

cd $wkdir/27_RepeatMasker_II/soft/

RepeatMasker -s -lib $wkdir/24_fasta_filter/briggsae_combo_repeats_classified_no_prot.fa -xsmall -gff -pa 16  $wkdir/01_fasta_filter/briggsae_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/elegans_combo_repeats_classified_no_prot.fa -xsmall -gff -pa 16  $wkdir/01_fasta_filter/elegans_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa -xsmall -gff -pa 16  $wkdir/01_fasta_filter/inopinata_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/nigoni_combo_repeats_classified_no_prot.fa -xsmall -gff -pa 16  $wkdir/01_fasta_filter/nigoni_genome.fa
RepeatMasker -s -lib $wkdir/24_fasta_filter/remanei_combo_repeats_classified_no_prot.fa -xsmall -gff -pa 16  $wkdir/01_fasta_filter/remanei_genome.fa

#############
#############
#############
#############
#Now, generating data needed to plot the global genomic landscape of repeats in all species
#############
#############
#############
#############

#make data for Figure 2, ie the global genomic landscapes of repetitive elements
mkdir $wkdir/29_prep_plots
cd $wkdir/29_prep_plots
mkdir repeatmasker_ii
cd repeatmasker_ii
mkdir 00_grep

#First, we need to format the gff output from RepeatMasker to give us usable columns of cluster ids and RepeatClassifier classifications

#get repeat gff without header


cd $wkdir/27_RepeatMasker_II/hard

grep -v "#" briggsae_genome.fa.out.gff > $wkdir//29_prep_plots/repeatmasker_ii/00_grep/briggsae
grep -v "#" elegans_genome.fa.out.gff > $wkdir//29_prep_plots/repeatmasker_ii/00_grep/elegans
grep -v "#" inopinata_genome.fa.out.gff > $wkdir//29_prep_plots/repeatmasker_ii/00_grep/inopinata
grep -v "#" nigoni_genome.fa.out.gff > $wkdir//29_prep_plots/repeatmasker_ii/00_grep/nigoni
grep -v "#" remanei_genome.fa.out.gff > $wkdir//29_prep_plots/repeatmasker_ii/00_grep/remanei

#clean up list of repeat id's for grep and awk

mkdir $wkdir//29_prep_plots/repeatmasker_ii/01_awk

cd $wkdir//29_prep_plots/repeatmasker_ii/00_grep

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $9}' $i > $wkdir/29_prep_plots/repeatmasker_ii/01_awk/$i; done

cd ..
cd 01_awk

for i in *; do sed -i -e 's/ [0-9].*//g' $i; done
for i in *; do sed -i -e 's/Target "Motif://g' $i; done
for i in *; do sed -i -e 's/"//g' $i; done

#get just the simple repeats (tandem repeats of kmers)

cd ..

mkdir $wkdir//29_prep_plots/repeatmasker_ii/02_grep_simple_repeats

cd $wkdir//29_prep_plots/repeatmasker_ii/01_awk

for i in *; do grep 'n' $i | sort | uniq > $wkdir/29_prep_plots/repeatmasker_ii/02_grep_simple_repeats/$i; done

#add "Simple_repeat"

cd  $wkdir/29_prep_plots/repeatmasker_ii/02_grep_simple_repeats

for i in *; do sed -i -e 's/$/#Simple_repeat/g' $i; done


#get the cluster ids

mkdir $wkdir/29_prep_plots/repeatmasker_ii/03_grep_clusters


cd $wkdir/17_cat_RepeatClassifer/

for i in *; do grep ">" $i > $wkdir/29_prep_plots/repeatmasker_ii/03_grep_clusters/${i%_combo_repeats_classified.fa}; done

#clean up cluster id list for grep and awk

cd $wkdir/29_prep_plots/repeatmasker_ii/03_grep_clusters

for i in *; do sed -i -e 's/>//g' $i; done



#need to grep rich! for things like "AT-rich", "A-rich" 

mkdir 04_grep_rich

cd $wkdir/29_prep_plots/repeatmasker_ii/01_awk

for i in *; do grep 'rich' $i | sort | uniq > $wkdir/29_prep_plots/repeatmasker_ii/04_grep_rich/$i; done


cd $wkdir/29_prep_plots/repeatmasker_ii/04_grep_rich
mkdir awk

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"#",$0}' $i > $i.tmp && mv $i.tmp $i; done
for i in *; do sed -i -e 's/\t//g' $i; done


#put together the cluster ids

mkdir $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids
cd $wkdir/29_prep_plots/repeatmasker_ii/02_grep_simple_repeats

for i in *; do cat $i $wkdir/29_prep_plots/repeatmasker_ii/03_grep_clusters/$i $wkdir/29_prep_plots/repeatmasker_ii/04_grep_rich/$i > $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/$i; done


#remove "#"


cd $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids
for i in *; do sed -i -e 's/#/\t/g' $i; done




#get the repeat classifications for each record in the gff!

mkdir $wkdir/29_prep_plots/repeatmasker_ii/06_grep
cd $wkdir/29_prep_plots/repeatmasker_ii/01_awk

while read line; do grep -w  $line $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/briggsae; done < briggsae > $wkdir/29_prep_plots/repeatmasker_ii/06_grep/briggsae &

while read line; do grep -w  $line $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/elegans; done < elegans > $wkdir/29_prep_plots/repeatmasker_ii/06_grep/elegans &

while read line; do grep -w $line $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/inopinata; done < inopinata > $wkdir/29_prep_plots/repeatmasker_ii/06_grep/inopinata &

while read line; do grep -w $line $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/nigoni; done < nigoni > $wkdir/29_prep_plots/repeatmasker_ii/06_grep/nigoni &

while read line; do grep -w  $line $wkdir/29_prep_plots/repeatmasker_ii/05_cat_repeat_ids/remanei; done < remanei > $wkdir/29_prep_plots/repeatmasker_ii/06_grep/remanei &

#woo!

#sed on the things to remove painful characters....

cd $wkdir/29_prep_plots/repeatmasker_ii/06_grep/

for i in *; do sed -i -e 's/(//g' $i; done
for i in *; do sed -i -e 's/)//g' $i; done
for i in *; do sed -i -e 's/?/_q/g' $i; done
for i in *; do sed -i -e 's/\//-/g' $i; done

#now paste!

mkdir $wkdir/29_prep_plots/repeatmasker_ii/07_paste

cd $wkdir/29_prep_plots/repeatmasker_ii/00_grep

for i in *; do paste $i $wkdir/29_prep_plots/repeatmasker_ii/06_grep/$i > $wkdir/29_prep_plots/repeatmasker_ii/07_paste/$i; done &

#now, we have gff's with cluster and classification columns.

#next, we turn this into a bed file for use with bedtools

cd $wkdir/29_prep_plots/repeatmasker_ii/

mkdir $wkdir/29_prep_plots/repeatmasker_ii/17_awk

cd $wkdir/29_prep_plots/repeatmasker_ii/07_paste

#extract the columns we need in the right order for bedtools ; add length of insertion and species column

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$4,$5,$14,$5-$4,FILENAME}' $i > $wkdir/29_prep_plots/repeatmasker_ii/17_awk/$i; done

cd $wkdir/29_prep_plots/repeatmasker_ii/17_awk

cat * > $wkdir/29_prep_plots/repeatmasker_ii/18_cat/all


#change unknown bases in assembly (N) to (X) for bedtools nuc such that only repetitive bases are N and are counted as such.

mkdir 43_bedtools_nuc_perc_N_3-11-19

cd 43_bedtools_nuc_perc_N_3-11-19

mkdir 00_genome_links

ln -s $wkdir/01_fasta_filter/briggsae_genome.fa $wkdir/43_bedtools_nuc_perc_N_3-11-19/00_genome_links/briggsae_genome.fa
ln -s $wkdir/01_fasta_filter/elegans_genome.fa $wkdir/43_bedtools_nuc_perc_N_3-11-19/00_genome_links/elegans_genome.fa
ln -s $wkdir/01_fasta_filter/nigoni_genome.fa $wkdir/43_bedtools_nuc_perc_N_3-11-19/00_genome_links/nigoni_genome.fa
ln -s $wkdir/01_fasta_filter/remanei_genome.fa $wkdir/43_bedtools_nuc_perc_N_3-11-19/00_genome_links/remanei_genome.fa
ln -s $wkdir/01_fasta_filter/inopinata_genome.fa $wkdir/43_bedtools_nuc_perc_N_3-11-19/00_genome_links/inopinata_genome.fa

mkdir 01_sed_genomes_X

cd 00_genome_links

for i in *; do sed -e 's/N/X/g' $i > $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/$i; done &

#get simplified bed files for masking with bedtools maskfasta

mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed

cd $wkdir/29_prep_plots/repeatmasker_ii/17_awk/

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' briggsae | bedtools sort -i  > $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/briggsae_mask_bed

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' elegans | bedtools sort -i  > $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/elegans_mask_bed

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' inopinata | bedtools sort -i  > $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/inopinata_mask_bed

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' nigoni | bedtools sort -i  > $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/nigoni_mask_bed

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' briggsae | bedtools sort -i  > $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/nigoni_mask_bed

#mask genome assemblies with prepped bed files and bedtools maskfasta (bedtools version 2.25.0)

mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all

bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/briggsae_genome.fa -bed $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/briggsae_mask_bed -fo $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/briggsae &

bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/elegans_genome.fa -bed $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/elegans_mask_bed -fo $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/elegans &

bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/inopinata_genome.fa -bed $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/inopinata_mask_bed -fo $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/inopinata  &

bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/nigoni_genome.fa -bed $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/nigoni_mask_bed -fo $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/nigoni  &


bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/remanei_genome.fa -bed $wkdir/43_bedtools_nuc_perc_N_3-11-19/08_prep_mask_bed/remanei_mask_bed -fo $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/remanei  &

#generate 10kb and 100kb genomic windows with bedtools 

mkdir $wkdir/windows

cd $wkdir/01_fasta_filter/

#get length of all chromosomes per species
	#.fai index files should exist because upstream processing ; otherwise can be made with samtools faidx (samtools faidx <ref.fasta>).

for i in *.fa.fai; do cut -f1,2 $i > $wkdir/windows/${i%_genome.fa.fai}.chrom.sizes; done

#make the windows for bedtools with bedtools makewindows
cd $wkdir/windows/

	#10kb
for i in *.chrom.sizes; do bedtools makewindows -w 10000 -g $i > ${i%.chrom.sizes}.10kb.windows; done

	#100 kb
for i in *.chrom.sizes; do bedtools makewindows -w 100000 -g $i > ${i%.chrom.sizes}.100kb.windows; done


#bedtools nuc to calculate the number of repetitive bases in 10kb windows in all species

mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all

bedtools nuc -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/briggsae -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all/briggsae &

bedtools nuc -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/nigoni -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all/nigoni &

bedtools nuc -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/remanei -bed $wkdir/windows/remanei.10kb.windows > $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all/remanei &

bedtools nuc -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/elegans -bed $wkdir/windows/elegans.10kb.windows > $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all/elegans &

bedtools nuc -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/09_bedtools_maskfasta_all/inopinata -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all/inopinata &

#prep for plotting and stats (here, get just chromosome, position, number of N in window, and species)

mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/10_bedtools_nuc_all

awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"briggsae"}' briggsae > $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/briggsae
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"elegans"}' elegans > $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/elegans
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"inopinata"}' inopinata > $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/inopinata
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"nigoni"}' nigoni > $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/nigoni
awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"remaeni"}' remanei > $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/remanei

#remove header and clean up spelling mistake
cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk

for i in *; do sed -i '1d' $i; done

for i in *; do sed -i -e 's/remaeni/remanei/g' $i; done

#combine into one file
cat * > all

#add header
for i in *; do echo -e "Chr\tBP\tbp_rep\tspecies" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done

#now, normalize by chromosomal position.....

mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos
mkdir $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat
cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#split up data by chromosome for each species
cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/briggsae

awk '{print > $1}' $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/briggsae


cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/elegans

awk '{print > $1}' $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/elegans

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/inopinata

awk '{print > $1}' $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/inopinata

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/nigoni

awk '{print > $1}' $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/nigoni

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/remanei

awk '{print > $1}' $wkdir/43_bedtools_nuc_perc_N_3-11-19/11_awk/remanei

#forgot to remove header so get rid of file "Chr"
cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/briggsae

rm Chr


cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/elegans

rm Chr

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/inopinata

rm Chr

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/nigoni

rm Chr

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/remanei

rm Chr


##normalize by chromosome position for each chromosome in each species by ((|(chromosome length-position)|/chromosome length)/2). centers will be 0, ends will be 0.5.

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/briggsae


cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/elegans

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/inopinata



cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/nigoni


cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/12_awk_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/remanei

#combine all species files
cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/

cat * > global_repeat_density_10kb_win_norm_dist_cent.tsv
######The above is the data for Figure 2, Figure 4a and Supplemental Figure 1-2! [formerly "all" , "all_N_account_for_unassembled_bases_4-16-19"]


#add headers
for i in *; do echo -e "Chr\tBP\tnum_rep\tspecies\tnorm_dist_center" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done
#fix tabs
for i in *; do sed -i -e 's/ /\t/g' $i; done


#############
#############
#############
#############
#make a disjoined TE gff
#############
#############
#############
#############

#For the purposes of estimating the perentage of genomic windows beloning to specific repeat taxa (ie, TE superfamilies, classes, orders, etc.), I set out to make a disjoined TE gff3 that accounts for TE's nested within each other so that bases are not counted multiple times. This is inspired by Stitzer et al. 2019 "The Genomic Ecosystem of Transposable Elements in Maize" (https://www.biorxiv.org/content/10.1101/559922v1). See https://mcstitzer.github.io/maize_TEs/ for nice graphical explanation and https://github.com/mcstitzer/w22_te_annotation/blob/master/combine_all_TEs.R for code.

mkdir $wkdir/60_disjoined_gff

cd $wkdir/60_disjoined_gff

mkdir $wkdir/60_disjoined_gff/00_links

ln -s $wkdir/29_prep_plots/repeatmasker_ii/07_paste/briggsae $wkdir/60_disjoined_gff/00_links/briggsae
ln -s $wkdir/29_prep_plots/repeatmasker_ii/07_paste/elegans $wkdir/60_disjoined_gff/00_links/elegans
ln -s $wkdir/29_prep_plots/repeatmasker_ii/07_paste/inopinata $wkdir/60_disjoined_gff/00_links/inopinata
ln -s $wkdir/29_prep_plots/repeatmasker_ii/07_paste/nigoni $wkdir/60_disjoined_gff/00_links/nigoni
ln -s $wkdir/29_prep_plots/repeatmasker_ii/07_paste/remanei $wkdir/60_disjoined_gff/00_links/remanei

mkdir $wkdir/60_disjoined_gff/01_awk

cd $wkdir/60_disjoined_gff/00_links

#convert to format for gff

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$11,$10}' $i > $wkdir/60_disjoined_gff/01_awk/$i; done

#ok, awk by cluster... here the goal is to give each insertion unique id

mkdir $wkdir/60_disjoined_gff/02_awk/
mkdir $wkdir/60_disjoined_gff/02_awk/briggsae
mkdir $wkdir/60_disjoined_gff/02_awk/elegans
mkdir $wkdir/60_disjoined_gff/02_awk/inopinata
mkdir $wkdir/60_disjoined_gff/02_awk/nigoni
mkdir $wkdir/60_disjoined_gff/02_awk/remanei

cd $wkdir/60_disjoined_gff/02_awk/briggsae

awk '{print > $10}' $wkdir/60_disjoined_gff/01_awk/briggsae

cd $wkdir/60_disjoined_gff/02_awk/elegans
awk '{print > $10}' $wkdir/60_disjoined_gff/01_awk/elegans

cd $wkdir/60_disjoined_gff/02_awk/inopinata
awk '{print > $10}' $wkdir/60_disjoined_gff/01_awk/inopinata

cd $wkdir/60_disjoined_gff/02_awk/nigoni
awk '{print > $10}' $wkdir/60_disjoined_gff/01_awk/nigoni

cd $wkdir/60_disjoined_gff/02_awk/remanei
awk '{print > $10}' $wkdir/60_disjoined_gff/01_awk/remanei

#add unique insertion ids (numbers)

mkdir $wkdir/60_disjoined_gff/03_nl/
mkdir $wkdir/60_disjoined_gff/03_nl/briggsae
mkdir $wkdir/60_disjoined_gff/03_nl/elegans
mkdir $wkdir/60_disjoined_gff/03_nl/inopinata
mkdir $wkdir/60_disjoined_gff/03_nl/nigoni
mkdir $wkdir/60_disjoined_gff/03_nl/remanei

cd $wkdir/60_disjoined_gff/02_awk/briggsae
for i in *; do nl -n ln $i > $wkdir/60_disjoined_gff/03_nl/briggsae/$i; done

cd $wkdir/60_disjoined_gff/02_awk/elegans
for i in *; do nl -n ln $i > $wkdir/60_disjoined_gff/03_nl/elegans/$i; done

cd $wkdir/60_disjoined_gff/02_awk/inopinata
for i in *; do nl -n ln $i > $wkdir/60_disjoined_gff/03_nl/inopinata/$i; done

cd $wkdir/60_disjoined_gff/02_awk/nigoni
for i in *; do nl -n ln $i > $wkdir/60_disjoined_gff/03_nl/nigoni/$i; done

cd $wkdir/60_disjoined_gff/02_awk/remanei
for i in *; do nl -n ln $i > $wkdir/60_disjoined_gff/03_nl/remanei/$i; done

#move id numbers

mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/
mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/briggsae
mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/elegans
mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/inopinata
mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/nigoni
mkdir $wkdir/60_disjoined_gff/04_awk_move_id_numbers/remanei


cd $wkdir/60_disjoined_gff/03_nl/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' $i > $wkdir/60_disjoined_gff/04_awk_move_id_numbers/briggsae/$i; done

cd $wkdir/60_disjoined_gff/03_nl/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' $i > $wkdir/60_disjoined_gff/04_awk_move_id_numbers/elegans/$i; done

cd $wkdir/60_disjoined_gff/03_nl/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' $i > $wkdir/60_disjoined_gff/04_awk_move_id_numbers/inopinata/$i; done

cd $wkdir/60_disjoined_gff/03_nl/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' $i > $wkdir/60_disjoined_gff/04_awk_move_id_numbers/nigoni/$i; done

cd $wkdir/60_disjoined_gff/03_nl/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$1}' $i > $wkdir/60_disjoined_gff/04_awk_move_id_numbers/remanei/$i; done



#replace last tab with underscore

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/briggsae
for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/elegans
for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/inopinata
for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/nigoni
for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/remanei
for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done


#ok, now cat
mkdir $wkdir/60_disjoined_gff/05_cat_gff/

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/briggsae
cat * > $wkdir/60_disjoined_gff/05_cat_gff/briggsae

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/elegans
cat * > $wkdir/60_disjoined_gff/05_cat_gff/elegans

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/inopinata
cat * > $wkdir/60_disjoined_gff/05_cat_gff/inopinata

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/nigoni
cat * > $wkdir/60_disjoined_gff/05_cat_gff/nigoni

cd $wkdir/60_disjoined_gff/04_awk_move_id_numbers/remanei
cat * > $wkdir/60_disjoined_gff/05_cat_gff/remanei

#get right for gff sorting

mkdir $wkdir/60_disjoined_gff/06_awk/

cd $wkdir/60_disjoined_gff/05_cat_gff/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,"ID="$10";Name="$9}' $i |sed -e 's/     /\g/' > $wkdir/60_disjoined_gff/06_awk/$i; done



#get a maize control up there... /B73.structuralTEv2.fulllength.2018-09-19.gff3 maize_control


#add "ID" and "Name" for disjoined TE gff3 and disjoin_gff_args.R

mkdir $wkdir/60_disjoined_gff/11_awk/

cd $wkdir/60_disjoined_gff/05_cat_gff/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,"ID="$10";Name="$9";Name="$9}' $i | sed 's/ \+//g'  >$wkdir/60_disjoined_gff/11_awk/$i; done

cd $wkdir/60_disjoined_gff/11_awk/

mkdir $wkdir/60_disjoined_gff/12_gt_gff3

#add XXXX for downstream processing and ##gff-version so R likes the format

for i in *; do sed -i -e 's/;Name=/XXXX/1' $i; done

for i in *; do echo -e "##gff-version 3" | cat - $i > $i.gff3; done


#sort gff with genometools
for i in *gff3; do gt gff3 -sort -retainids $i > $wkdir/60_disjoined_gff/12_gt_gff3/$i; done



#clean up weirdness...

cd $wkdir/60_disjoined_gff/12_gt_gff3/

for i in *; do sed -i -e 's/###//g' $i; done

for i in *; do sed -i -e '/^$/d' $i; done


#disjoin

module load gcc/7.3
mkdir 13_disjoin_gff_args_R

cd $wkdir/60_disjoined_gff/12_gt_gff3/

for i in *; do Rscript $wkdir/60_disjoined_gff/disjoin_gff_args.R $i; done &


#ok, great now we got something. get in bed format.


mkdir $wkdir/60_disjoined_gff/14_awk_bed/

cd  $wkdir/60_disjoined_gff/13_disjoin_gff_args_R/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$9}' $i | grep -v "#" > $wkdir/60_disjoined_gff/14_awk_bed/$i; done

cd $wkdir/60_disjoined_gff/14_awk_bed/

#clean up

rm maize_control
mv briggsae.gff3 briggsae
mv elegans.gff3 elegans
mv inopinata.gff3 inopinata
mv nigoni.gff3 nigoni
mv remanei.gff3 remanei

for i in *; do sed -i -e 's/X=.*;ID/ID/g' $i; done
for i in *; do sed -i -e 's/;intact.*//g' $i; done
for i in *; do sed -i -e 's/ID=//g' $i; done
for i in *; do sed -i -e 's/XXXX/\t/g' $i; done


#add the repeat taxonomy (repeat_taxonomy.tsv) for analysis of various repeat taxa
cd $wkdir/60_disjoined_gff/14_awk_bed/


for i in *; do sed -i -e 's/GA-rich/GA-rich_only\t/2' $i; done
for i in *; do sed -i -e 's/A-rich$/A-rich_only\t/g' $i; done
for i in *; do sed -i -e 's/G-rich$/G-rich_only\t/g' $i; done
for i in *; do sed -i -e 's/LINE$/LINE_only\t/g' $i; done
for i in *; do sed -i -e 's/LTR$/LTR_only\t/g' $i; done
for i in *; do sed -i -e 's/Satellite$/Satellite_only\t/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar$/DNA-TcMar_only\t/g' $i; done
for i in *; do sed -i -e 's/DNA$/DNA_only\t/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT$/DNA-hAT_only\t/g' $i; done
for i in *; do sed -i -e 's/DNA-Sola$/DNA-Sola_only\t/g' $i; done
for i in *; do sed -i -e 's/LTR-Gypsy$/LTR-Gypsy_only\t/g' $i; done
for i in *; do sed -i -e 's/LINE-L1$/LINE-L1_only\t/g' $i; done
for i in *; do sed -i -e 's/LINE-R1$/LINE-R1_only\t/g' $i; done
for i in *; do sed -i -e 's/LINE-R2$/LINE-R2_only\t/g' $i; done
for i in *; do sed -i -e 's/SINE-tRNA$/SINE-tRNA_only\t/g' $i; done
for i in *; do sed -i -e 's/A-rich_only/A-rich\tLow_complexity\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/ARTEFACT/ARTEFACT\tARTEFACT\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-CMC-Chapaev/DNA-CMC-Chapaev\tII\t1\tTIR\tCACTA\tChapaev/g' $i; done
for i in *; do sed -i -e 's/DNA-CMC-EnSpm/DNA-CMC-EnSpm\tII\t1\tTIR\tCACTA\tEnSpm/g' $i; done
for i in *; do sed -i -e 's/DNA-CMC-Mirage/DNA-CMC-Mirage\tII\t1\tTIR\tCACTA\tMirage/g' $i; done
for i in *; do sed -i -e 's/DNA-Crypton-S/DNA-Crypton-S\tII\t1\tCrypton\tCrypton\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-Dada/DNA-Dada\tII\t1\tTIR\tDada\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Ac/DNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Blackjack/DNA-hAT-Blackjack\tII\t1\tTIR\thAT\tBlackjack/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Charlie/DNA-hAT-Charlie\tII\t1\tTIR\thAT\tCharlie/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-hATx/DNA-hAT-hATx\tII\t1\tTIR\thAT\thATx/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-hobo/DNA-hAT-hobo\tII\t1\tTIR\thAT\thobo/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT_only/DNA-hAT\tII\t1\tTIR\thAT\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Restless/DNA-hAT-Restless\tII\t1\tTIR\thAT\tRestless/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Tag1/DNA-hAT-Tag1\tII\t1\tTIR\thAT\tTag1/g' $i; done
for i in *; do sed -i -e 's/DNA-hAT-Tip100/DNA-hAT-Tip100\tII\t1\tTIR\thAT\tTip100/g' $i; done
for i in *; do sed -i -e 's/DNA_only/DNA\tII\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-Kolobok-T2/DNA-Kolobok-T2\tII\t1\tTIR\tKolobok\tT2/g' $i; done
for i in *; do sed -i -e 's/DNA-Maverick/DNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-Merlin/DNA-Merlin\tII\t1\tTIR\tMerlin\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-MuLE-MuDR/DNA-MuLE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' $i; done
for i in *; do sed -i -e 's/DNA-MULE-MuDR/DNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' $i; done
for i in *; do sed -i -e 's/DNA-MULE-NOF/DNA-MULE-NOF\tII\t1\tTIR\tMutator\tNOF/g' $i; done
for i in *; do sed -i -e 's/DNA-P-Fungi/DNA-P-Fungi\tII\t1\tTIR\tP\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-PIF-Harbinger/DNA-PIF-Harbinger\tII\t1\tTIR\tPIF-Harbinger\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-PIF-ISL2EU/DNA-PIF-ISL2EU\tII\t1\tTIR\tPIF-Harbinger\tISL2EU/g' $i; done
for i in *; do sed -i -e 's/DNA-PiggyBac/DNA-PiggyBac\tII\t1\tTIR\tPiggyBac\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-Sola-3/DNA-Sola-3\tII\t1\tTIR\tSola\tSola3/g' $i; done
for i in *; do sed -i -e 's/DNA-Sola_only/DNA-Sola\tII\t1\tTIR\tSola\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Fot1/DNA-TcMar-Fot1\tII\t1\tTIR\tTc1-Mariner\tFot1/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar_only/DNA-TcMar\tII\t1\tTIR\tTc1-Mariner\tNA/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-m44/DNA-TcMar-m44\tII\t1\tTIR\tTc1-Mariner\tm44/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Mariner/DNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Pogo/DNA-TcMar-Pogo\tII\t1\tTIR\tTc1-Mariner\tPogo/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Tc1/DNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Tc2/DNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/g' $i; done
for i in *; do sed -i -e 's/DNA-TcMar-Tc4/DNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' $i; done
for i in *; do sed -i -e 's/G-rich_only/G-rich\tLow_complexity\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-CR1/LINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' $i; done
for i in *; do sed -i -e 's/LINE-CRE/LINE-CRE\tI\tNA\tLINE\tR2\tCRE/g' $i; done
for i in *; do sed -i -e 's/LINE-Dong-R4/LINE-Dong-R4\tI\tNA\tLINE\tR2\tDong/g' $i; done
for i in *; do sed -i -e 's/LINE-I-Jockey/LINE-I-Jockey\tI\tNA\tLINE\tI\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE_only/LINE\tI\tNA\tLINE\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-Jockey/LINE-Jockey\tI\tNA\tLINE\tJockey\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-L1_only/LINE-L1\tI\tNA\tLINE\tL1\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-L1-Tx1/LINE-L1-Tx1\tI\tNA\tLINE\tL1\tTx1/g' $i; done
for i in *; do sed -i -e 's/LINE-L2/LINE-L2\tI\tNA\tLINE\tJockey\tL2/g' $i; done
for i in *; do sed -i -e 's/LINE-Penelope/LINE-Penelope\tI\tNA\tPLE\tPenelope\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-R1_only/LINE-R1\tI\tNA\tLINE\tI\tR1/g' $i; done
for i in *; do sed -i -e 's/LINE-R1-LOA/LINE-R1-LOA\tI\tNA\tLINE\tI\tR1/g' $i; done
for i in *; do sed -i -e 's/LINE-R2_only/LINE-R2\tI\tNA\tLINE\tR2\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-R2-NeSL/LINE-R2-NeSL\tI\tNA\tLINE\tR2\tNeSL/g' $i; done
for i in *; do sed -i -e 's/LINE-RTE-BovB/LINE-RTE-BovB\tI\tNA\tLINE\tRTE\tBovB/g' $i; done
for i in *; do sed -i -e 's/LINE-RTE-RTE/LINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' $i; done
for i in *; do sed -i -e 's/LINE-RTE-X/LINE-RTE-X\tI\tNA\tLINE\tRTE\tX/g' $i; done
for i in *; do sed -i -e 's/LTR-Copia/LTR-Copia\tI\tNA\tLTR\tCopia\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR-DIRS/LTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR-ERV1/LTR-ERV1\tI\tNA\tLTR\tERV\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR-ERVL/LTR-ERVL\tI\tNA\tLTR\tERV\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR-Foamy/LTR-Foamy\tI\tNA\tLTR\tERV\tFoamy/g' $i; done
for i in *; do sed -i -e 's/LTR-Gypsy-Cigr/LTR-Gypsy-Cigr\tI\tNA\tLTR\tGypsy\tCigr/g' $i; done
for i in *; do sed -i -e 's/LTR-Gypsy_only/LTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR_only/LTR\tI\tNA\tLTR\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/LTR-Pao/LTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' $i; done
for i in *; do sed -i -e 's/RC-Helitron/RC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' $i; done
for i in *; do sed -i -e 's/rRNA/rRNA\trRNA\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/Satellite_only/Satellite\tSatellite\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/Satellite-W-chromosome/Satellite-W-chromosome\tSatellite\tNA\tNA\tNA\tW-chromosome/g' $i; done
for i in *; do sed -i -e 's/Simple_repeat/Simple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/SINE-ID/SINE-ID\tI\tNA\tSINE\ttRNA\tID/g' $i; done
for i in *; do sed -i -e 's/SINE_q/SINE_q\tI\tNA\tSINE\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/SINE-R2/SINE-R2\tI\tNA\tSINE\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/SINE-tRNA_only/SINE-tRNA\tI\tNA\tSINE\ttRNA\tNA/g' $i; done
for i in *; do sed -i -e 's/SINE-tRNA-RTE/SINE-tRNA-RTE\tI\tNA\tSINE\ttRNA\tNA/g' $i; done
for i in *; do sed -i -e 's/snRNA/snRNA\tsnRNA\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/Unknown/Unknown\tNA\tNA\tNA\tNA\tNA/g' $i; done
for i in *; do sed -i -e 's/\t_only//g' $i; done


#############
#############
#############
#############
#now, generate landscapes for all repeat taxa (classes, orders, superfamilies, families)
#############
#############
#############
#############


#superfamilies first




mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links

cp -r $wkdir/60_disjoined_gff/14_awk_bed $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/

cp -r $wkdir/46_new_taxonomy_num_bp/00_links/X_genomes $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/

#superfamily first. this is where the major discrepancies are anyway


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/remanei

#awk by superfamily
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/briggsae
awk '{print > $9}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/briggsae
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/elegans
awk '{print > $9}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/elegans
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/inopinata
awk '{print > $9}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/inopinata
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/nigoni
awk '{print > $9}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/nigoni
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/remanei
awk '{print > $9}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/remanei

#remove NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/briggsae
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/elegans
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/inopinata
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/nigoni
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/remanei
rm NA


#now, sort bed

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort//elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/remanei/$i; done &



#now, maskfasta

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/briggsae

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/briggsae.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/elegans

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/elegans.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/inopinata

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/inopinata.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/inopinata/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/nigoni

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/nigoni.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/nigoni/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/01_bedtools_sort/remanei

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/remanei.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/remanei/$i; done &



#bedtools nuc

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/briggsae

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/elegans

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/elegans.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/inopinata

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/inopinata/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/nigoni

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/nigoni/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/02_bedtools_maskfasta/remanei

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/remanei.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/remanei/$i; done &




#now, awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/briggsae

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"briggsae", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/elegans

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"elegans", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/elegans/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/inopinata

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"inopinata", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/inopinata/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/nigoni

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"nigoni", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/nigoni/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/03_bedtools_nuc/remanei

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"remanei", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/remanei/$i; done &




#remove headers

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/briggsae/

for i in *; do sed -i '1d' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/elegans/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/inopinata/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/nigoni/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/remanei/

for i in *; do sed -i '1d' $i; done



#now cat species

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/briggsae

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/elegans


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/inopinata


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/

cat * > all

#awk by repeat type
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/06_awk_class


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/06_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species/all


#make genomic landscape plots for all repeat superfamilies ;  make sure to change output directory in repeat_landscapes_plots.R
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/07_perc_N_rep_class_pdf
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/06_awk_class

for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/07_perc_N_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf


#now, set up for normalized by chromosome position stuff; this time with sed for appropriate values

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/05_cat_species

for i in *; do cp $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/$i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center

rm all

#these should probably be run in parallel, will take a while
cd $wkdir/scripts

./briggsae_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/briggsae
./elegans_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/elegans
./inopinata_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/inopinata
./nigoni_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/nigoni
./remanei_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/08_norm_dist_center/remanei


#awk by class

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/briggsae
awk '{print > $6}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/briggsae
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/elegans
awk '{print > $6}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/elegans
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/inopinata
awk '{print > $6}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/inopinata
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/nigoni
awk '{print > $6}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/nigoni
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/remanei
awk '{print > $6}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/remanei

#now, sort bed

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort//elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/remanei/$i; done &


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/briggsae

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/briggsae.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/elegans

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/elegans.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/inopinata

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/inopinata.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/inopinata/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/nigoni

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/nigoni.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/nigoni/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/01_bedtools_sort/remanei

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/remanei.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/remanei/$i; done &





#bedtools nuc

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/briggsae

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/elegans

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/elegans.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/inopinata

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/inopinata/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/nigoni

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/nigoni/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/remanei

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/remanei.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/remanei/$i; done &



#bedtools nuc

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/briggsae

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/elegans

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/elegans.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/inopinata

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/inopinata/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/nigoni

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/nigoni/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/02_bedtools_maskfasta/remanei

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/remanei.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/remanei/$i; done &



#now, awk to get only columns i need

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/briggsae

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"briggsae", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/elegans

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"elegans", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/elegans/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/inopinata

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"inopinata", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/inopinata/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/nigoni

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"nigoni", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/nigoni/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/03_bedtools_nuc/remanei

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"remanei", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/remanei/$i; done &


#remove headers

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/briggsae/

for i in *; do sed -i '1d' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/elegans/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/inopinata/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/nigoni/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/remanei/

for i in *; do sed -i '1d' $i; done



#now cat species

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/briggsae

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/elegans


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/inopinata


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/

cat * > all



#awk by repeat type
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/06_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/07_perc_N_rep_class_pdf

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/06_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species/all


#add header
for i in *; do echo -e "Chr\tBP\tbp_rep\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done

	#remanei is spelled wrong...

for i in *; do sed -i -e 's/remaeni/remanei/g' $i; done


#change some things...

awk 'BEGIN {FS="\t"} {OFS="\t"}$5 == "I" { $5="class_I_retrotransposon" }1' I > class_I_retrotransposon
awk 'BEGIN {FS="\t"} {OFS="\t"}$5 == "II" { $5="class_II_DNA_transposon" }1' II > class_II_DNA_transposon
awk 'BEGIN {FS="\t"} {OFS="\t"}$5 == "NA" { $5="Unclassified" }1' NA > Unclassified

rm I
rm II
rm NA

#make the pdf's ; make sure to change output directory in repeat_landscapes_plots.R
for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/07_perc_N_rep_class_pdf

pdfunite *.pdf all_repeat_classes.pdf


#now, set up for normalized by chromosome position stuff; this time with sed for appropriate values

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/05_cat_species

for i in *; do cp $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/$i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center

rm all


#these should probably be run in parallel
cd $wkdir/scripts

./briggsae_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/briggsae
./elegans_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/elegans
./inopinata_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/inopinata
./nigoni_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/nigoni
./remanei_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/08_norm_dist_center/remanei


#now, repeat orders

#awk by order

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/briggsae
awk '{print > $8}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/briggsae
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/elegans
awk '{print > $8}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/elegans
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/inopinata
awk '{print > $8}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/inopinata
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/nigoni
awk '{print > $8}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/nigoni
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/remanei
awk '{print > $8}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/remanei

#now, sort bed

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort//elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/remanei/$i; done &

#now, maskfasta

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/briggsae

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/briggsae.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/elegans

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/elegans.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/inopinata

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/inopinata.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/inopinata/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/nigoni

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/nigoni.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/nigoni/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/01_bedtools_sort/remanei

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/remanei.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/remanei/$i; done &

#bedtools nuc

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/briggsae

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/elegans

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/elegans.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/inopinata

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/inopinata/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/nigoni

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/nigoni/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/02_bedtools_maskfasta/remanei

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/remanei.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/remanei/$i; done &

#now, awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/briggsae

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"briggsae", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/elegans

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"elegans", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/elegans/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/inopinata

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"inopinata", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/inopinata/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/nigoni

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"nigoni", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/nigoni/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/03_bedtools_nuc/remanei

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"remanei", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/remanei/$i; done &


#remove headers

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/briggsae/

for i in *; do sed -i '1d' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/elegans/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/inopinata/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/nigoni/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/remanei/

for i in *; do sed -i '1d' $i; done



#now cat species

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/briggsae

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/elegans


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/inopinata


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/

cat * > all



#split up by repeat type for making plots
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/06_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/07_perc_N_rep_class_pdf

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/06_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species/all



for i in *; do echo -e "Chr\tBP\tbp_rep\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done

	#remanei is spelled wrong...

for i in *; do sed -i -e 's/remaeni/remanei/g' $i; done

#make the pdf's ; make sure to change output directory in repeat_landscapes_plots.R
for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/07_perc_N_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf
#for normalize by distance to chromosome center


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/05_cat_species

for i in *; do cp $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/$i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center

rm all



#these should probably be run in parallel
cd $wkdir/scripts

./briggsae_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/briggsae
./elegans_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/elegans
./inopinata_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/inopinata
./nigoni_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/nigoni
./remanei_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/08_norm_dist_center/remanei





#ok, now family

#awk by order
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/briggsae
awk '{print > $10}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/briggsae
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/elegans
awk '{print > $10}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/elegans
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/inopinata
awk '{print > $10}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/inopinata
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/nigoni
awk '{print > $10}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/nigoni
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/remanei
awk '{print > $10}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/disjoined_bed/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/briggsae
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/elegans
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/inopinata
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/nigoni
rm NA
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/remanei
rm NA


#now, sort bed

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort//elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/remanei/$i; done &

#now, maskfasta

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/briggsae

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/briggsae.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/elegans

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/elegans.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/inopinata

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/inopinata.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/inopinata/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/nigoni

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/nigoni.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/nigoni/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/01_bedtools_sort/remanei

for i in *; do bedtools maskfasta -fi $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/X_genomes/remanei.fa -bed $i -fo $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/remanei/$i; done &

#bedtools nuc

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/briggsae

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/elegans

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/elegans.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/elegans/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/inopinata

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/inopinata/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/nigoni

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/nigoni/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/02_bedtools_maskfasta/remanei

for i in *; do bedtools nuc -fi $i -bed $wkdir/windows/remanei.10kb.windows > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/remanei/$i; done &

#now, awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/briggsae
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/briggsae

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"briggsae", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/elegans

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"elegans", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/elegans/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/inopinata

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"inopinata", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/inopinata/$i; done &




cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/nigoni

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"nigoni", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/nigoni/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/03_bedtools_nuc/remanei

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,"remanei", FILENAME}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/remanei/$i; done &


#remove headers

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/briggsae/

for i in *; do sed -i '1d' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/elegans/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/inopinata/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/nigoni/

for i in *; do sed -i '1d' $i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/remanei/

for i in *; do sed -i '1d' $i; done



#now cat species

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/briggsae

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/elegans


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/inopinata


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/

cat * > all


#split up by repeat type for making plots
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/06_awk_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/07_perc_N_rep_class_pdf

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/06_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species/all



for i in *; do echo -e "Chr\tBP\tbp_rep\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done

	#remanei is spelled wrong...

for i in *; do sed -i -e 's/remaeni/remanei/g' $i; done


for i in *; do Rscript $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/repeat_classes_perc_n_plots_4-24-19.R $i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/07_perc_N_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf

#for normalize by distance to chromosome center


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/00_links/norm_dist_center_sed

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/05_cat_species

for i in *; do cp $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/$i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center

rm all


#these should probably be run in parallel
cd $wkdir/scripts

./briggsae_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/briggsae
./elegans_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/elegans
./inopinata_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/inopinata
./nigoni_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/nigoni
./remanei_norm_dist_cen_sed.sh $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/08_norm_dist_center/remanei



#next, repeat counts


#############
#############
#############
#############
#repeat counts (ie, instead of number of repetitive base pairs as above, it is number of insertions)
#############
#############
#############
#############

#here, using the non-disjoined bedfiles to generate counts


mkdir 48_new_taxonomy_count


mkdir $wkdir/48_new_taxonomy_count/01_awk_species
mkdir $wkdir/48_new_taxonomy_count/02_class
mkdir $wkdir/48_new_taxonomy_count/03_order
mkdir $wkdir/48_new_taxonomy_count/04_superfamily
mkdir $wkdir/48_new_taxonomy_count/05_family


mkdir $wkdir/48_new_taxonomy_count/00_links

cp $wkdir/33_remanei_prep_plots/18_cat/all $wkdir/00_links/all_rep_all_species.bed

#get rid of columns I do not need

cd $wkdir/48_new_taxonomy_count/00_links

awk 'BEGIN {OFS="\t"} {print $4, $1,$2,$3,$6}' all_rep_all_species.bed > all_rep_all_species.bed.tmp
mv all_rep_all_species.bed.tmp all_rep_all_species.bed


#then, used sed.srun to replace the RepeatMasker classification with the new taxonomy.
#these are the commands

sed -i -e 's/GA-rich/GA-rich_only\t/g' all_rep_all_species.bed
sed -i -e 's/A-rich\t/A-rich_only\t/g' all_rep_all_species.bed
sed -i -e 's/LINE\t/LINE_only\t/g' all_rep_all_species.bed
sed -i -e 's/LTR\t/LTR_only\t/g' all_rep_all_species.bed
sed -i -e 's/Satellite\t/Satellite_only\t/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar\t/DNA-TcMar_only\t/g' all_rep_all_species.bed
sed -i -e 's/DNA\t/DNA_only\t/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT\t/DNA-hAT_only\t/g' all_rep_all_species.bed
sed -i -e 's/DNA-Sola\t/DNA-Sola_only\t/g' all_rep_all_species.bed
sed -i -e 's/LTR-Gypsy\t/LTR-Gypsy_only\t/g' all_rep_all_species.bed
sed -i -e 's/LINE-L1\t/LINE-L1_only\t/g' all_rep_all_species.bed
sed -i -e 's/LINE-R1\t/LINE-R1_only\t/g' all_rep_all_species.bed
sed -i -e 's/LINE-R2\t/LINE-R2_only\t/g' all_rep_all_species.bed
sed -i -e 's/SINE-tRNA\t/SINE-tRNA_only\t/g' all_rep_all_species.bed

sed -i -e 's/A-rich_only/A-rich\tLow_complexity\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/ARTEFACT/ARTEFACT\tARTEFACT\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-CMC-Chapaev/DNA-CMC-Chapaev\tII\t1\tTIR\tCACTA\tChapaev/g' all_rep_all_species.bed
sed -i -e 's/DNA-CMC-EnSpm/DNA-CMC-EnSpm\tII\t1\tTIR\tCACTA\tEnSpm/g' all_rep_all_species.bed
sed -i -e 's/DNA-CMC-Mirage/DNA-CMC-Mirage\tII\t1\tTIR\tCACTA\tMirage/g' all_rep_all_species.bed
sed -i -e 's/DNA-Crypton-S/DNA-Crypton-S\tII\t1\tCrypton\tCrypton\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-Dada/DNA-Dada\tII\t1\tTIR\tDada\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Ac/DNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Blackjack/DNA-hAT-Blackjack\tII\t1\tTIR\thAT\tBlackjack/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Charlie/DNA-hAT-Charlie\tII\t1\tTIR\thAT\tCharlie/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-hATx/DNA-hAT-hATx\tII\t1\tTIR\thAT\thATx/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-hobo/DNA-hAT-hobo\tII\t1\tTIR\thAT\thobo/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT_only/DNA-hAT\tII\t1\tTIR\thAT\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Restless/DNA-hAT-Restless\tII\t1\tTIR\thAT\tRestless/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Tag1/DNA-hAT-Tag1\tII\t1\tTIR\thAT\tTag1/g' all_rep_all_species.bed
sed -i -e 's/DNA-hAT-Tip100/DNA-hAT-Tip100\tII\t1\tTIR\thAT\tTip100/g' all_rep_all_species.bed
sed -i -e 's/DNA_only/DNA\tII\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-Kolobok-T2/DNA-Kolobok-T2\tII\t1\tTIR\tKolobok\tT2/g' all_rep_all_species.bed
sed -i -e 's/DNA-Maverick/DNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-Merlin/DNA-Merlin\tII\t1\tTIR\tMerlin\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-MuLE-MuDR/DNA-MuLE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' all_rep_all_species.bed
sed -i -e 's/DNA-MULE-MuDR/DNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' all_rep_all_species.bed
sed -i -e 's/DNA-MULE-NOF/DNA-MULE-NOF\tII\t1\tTIR\tMutator\tNOF/g' all_rep_all_species.bed
sed -i -e 's/DNA-P-Fungi/DNA-P-Fungi\tII\t1\tTIR\tP\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-PIF-Harbinger/DNA-PIF-Harbinger\tII\t1\tTIR\tPIF-Harbinger\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-PIF-ISL2EU/DNA-PIF-ISL2EU\tII\t1\tTIR\tPIF-Harbinger\tISL2EU/g' all_rep_all_species.bed
sed -i -e 's/DNA-PiggyBac/DNA-PiggyBac\tII\t1\tTIR\tPiggyBac\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-Sola-3/DNA-Sola-3\tII\t1\tTIR\tSola\tSola3/g' all_rep_all_species.bed
sed -i -e 's/DNA-Sola_only/DNA-Sola\tII\t1\tTIR\tSola\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Fot1/DNA-TcMar-Fot1\tII\t1\tTIR\tTc1-Mariner\tFot1/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar_only/DNA-TcMar\tII\t1\tTIR\tTc1-Mariner\tNA/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-m44/DNA-TcMar-m44\tII\t1\tTIR\tTc1-Mariner\tm44/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Mariner/DNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Pogo/DNA-TcMar-Pogo\tII\t1\tTIR\tTc1-Mariner\tPogo/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Tc1/DNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Tc2/DNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/g' all_rep_all_species.bed
sed -i -e 's/DNA-TcMar-Tc4/DNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' all_rep_all_species.bed
sed -i -e 's/GA-rich_only/GA-rich\tLow_complexity\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/G-rich/G-rich\tLow_complexity\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-CR1/LINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' all_rep_all_species.bed
sed -i -e 's/LINE-CRE/LINE-CRE\tI\tNA\tLINE\tR2\tCRE/g' all_rep_all_species.bed
sed -i -e 's/LINE-Dong-R4/LINE-Dong-R4\tI\tNA\tLINE\tR2\tDong/g' all_rep_all_species.bed
sed -i -e 's/LINE-I-Jockey/LINE-I-Jockey\tI\tNA\tLINE\tI\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE_only/LINE\tI\tNA\tLINE\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-Jockey/LINE-Jockey\tI\tNA\tLINE\tJockey\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-L1_only/LINE-L1\tI\tNA\tLINE\tL1\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-L1-Tx1/LINE-L1-Tx1\tI\tNA\tLINE\tL1\tTx1/g' all_rep_all_species.bed
sed -i -e 's/LINE-L2/LINE-L2\tI\tNA\tLINE\tJockey\tL2/g' all_rep_all_species.bed
sed -i -e 's/LINE-Penelope/LINE-Penelope\tI\tNA\tPLE\tPenelope\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-R1_only/LINE-R1\tI\tNA\tLINE\tI\tR1/g' all_rep_all_species.bed
sed -i -e 's/LINE-R1-LOA/LINE-R1-LOA\tI\tNA\tLINE\tI\tR1/g' all_rep_all_species.bed
sed -i -e 's/LINE-R2_only/LINE-R2\tI\tNA\tLINE\tR2\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-R2-NeSL/LINE-R2-NeSL\tI\tNA\tLINE\tR2\tNeSL/g' all_rep_all_species.bed
sed -i -e 's/LINE-RTE-BovB/LINE-RTE-BovB\tI\tNA\tLINE\tRTE\tBovB/g' all_rep_all_species.bed
sed -i -e 's/LINE-RTE-RTE/LINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' all_rep_all_species.bed
sed -i -e 's/LINE-RTE-X/LINE-RTE-X\tI\tNA\tLINE\tRTE\tX/g' all_rep_all_species.bed
sed -i -e 's/LTR-Copia/LTR-Copia\tI\tNA\tLTR\tCopia\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR-DIRS/LTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR-ERV1/LTR-ERV1\tI\tNA\tLTR\tERV\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR-ERVL/LTR-ERVL\tI\tNA\tLTR\tERV\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR-Foamy/LTR-Foamy\tI\tNA\tLTR\tERV\tFoamy/g' all_rep_all_species.bed
sed -i -e 's/LTR-Gypsy-Cigr/LTR-Gypsy-Cigr\tI\tNA\tLTR\tGypsy\tCigr/g' all_rep_all_species.bed
sed -i -e 's/LTR-Gypsy_only/LTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR_only/LTR\tI\tNA\tLTR\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/LTR-Pao/LTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' all_rep_all_species.bed
sed -i -e 's/RC-Helitron/RC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' all_rep_all_species.bed
sed -i -e 's/rRNA/rRNA\trRNA\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/Satellite_only/Satellite\tSatellite\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/Satellite-W-chromosome/Satellite-W-chromosome\tSatellite\tNA\tNA\tNA\tW-chromosome/g' all_rep_all_species.bed
sed -i -e 's/Simple_repeat/Simple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/SINE-ID/SINE-ID\tI\tNA\tSINE\ttRNA\tID/g' all_rep_all_species.bed
sed -i -e 's/SINE_q/SINE_q\tI\tNA\tSINE\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/SINE-R2/SINE-R2\tI\tNA\tSINE\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/SINE-tRNA_only/SINE-tRNA\tI\tNA\tSINE\ttRNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/SINE-tRNA-RTE/SINE-tRNA-RTE\tI\tNA\tSINE\ttRNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/snRNA/snRNA\tsnRNA\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed
sed -i -e 's/Unknown/Unknown\tNA\tNA\tNA\tNA\tNA/g' all_rep_all_species.bed

#flip the columns around again

awk 'BEGIN {OFS="\t"} {print $7,$8,$9,$10,$1,$2,$3,$4,$5,$6}' all_rep_all_species.bed > all_rep_all_species.bed.tmp
mv all_rep_all_species.bed.tmp all_rep_all_species.bed

mkdir $wkdir/48_new_taxonomy_count/01_awk_species
mkdir $wkdir/48_new_taxonomy_count/02_class
mkdir $wkdir/48_new_taxonomy_count/03_order
mkdir $wkdir/48_new_taxonomy_count/04_superfamily
mkdir $wkdir/48_new_taxonomy_count/05_family


#break up new repeat bed annotation file by species 

cd $wkdir/48_new_taxonomy_count/01_awk_species

awk '{print > $4}' $wkdir/48_new_taxonomy_count/00_links/all_rep_all_species.bed

#now, class

#awk by class

mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class
mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/00_awk_class/remanei


cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/briggsae
awk '{print > $6}' $wkdir/48_new_taxonomy_count/01_awk_species/briggsae
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/elegans
awk '{print > $6}' $wkdir/48_new_taxonomy_count/01_awk_species/elegans
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/inopinata
awk '{print > $6}' $wkdir/48_new_taxonomy_count/01_awk_species/inopinata
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/nigoni
awk '{print > $6}' $wkdir/48_new_taxonomy_count/01_awk_species/nigoni
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/remanei
awk '{print > $6}' $wkdir/48_new_taxonomy_count/01_awk_species/remanei

#now, sort bed

mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort
mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/remanei

cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort//elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/remanei/$i; done &


#awk one for sum with bedtools map

mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk
mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/02_class/02_awk/briggsae/$i; done
cd $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/02_class/02_awk/elegans/$i; done
cd $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/02_class/02_awk/inopinata/$i; done
cd $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/02_class/02_awk/nigoni/$i; done
cd $wkdir/48_new_taxonomy_count/02_class/01_bedtools_sort/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/02_class/02_awk/remanei/$i; done

#ok, cool now bedtools map

mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map
mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/remanei



cd $wkdir/48_new_taxonomy_count/02_class/02_awk/briggsae
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/briggsae.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/02_awk/elegans
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/elegans.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/02_awk/inopinata
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/inopinata.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/02_awk/nigoni
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/nigoni.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/02_class/02_awk/remanei
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/remanei.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/remanei/$i; done &

#grep to remove lines with no data

mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep
mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/04_grep/remanei



cd $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/briggsae
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/02_class/04_grep/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/elegans
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/02_class/04_grep/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/inopinata
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/02_class/04_grep/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/nigoni
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/02_class/04_grep/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map/remanei
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/02_class/04_grep/remanei/$i; done &

#awk to add in species and class data

mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk
mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/05_awk/remanei


cd $wkdir/48_new_taxonomy_count/02_class/04_grep/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/05_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/04_grep/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/05_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/04_grep/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/05_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/04_grep/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/05_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/04_grep/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/05_awk/remanei/$i; done &



#now cat species

mkdir $wkdir/48_new_taxonomy_count/02_class/06_cat_species

cd $wkdir/48_new_taxonomy_count/02_class/05_awk/briggsae/

cat * > $wkdir/48_new_taxonomy_count/02_class/06_cat_species/briggsae

cd $wkdir/48_new_taxonomy_count/02_class/05_awk/elegans/

cat * > $wkdir/48_new_taxonomy_count/02_class/06_cat_species/elegans


cd $wkdir/48_new_taxonomy_count/02_class/05_awk/inopinata/

cat * > $wkdir/48_new_taxonomy_count/02_class/06_cat_species/inopinata


cd $wkdir/48_new_taxonomy_count/02_class/05_awk/nigoni/

cat * > $wkdir/48_new_taxonomy_count/02_class/06_cat_species/nigoni

cd $wkdir/48_new_taxonomy_count/02_class/05_awk/remanei/

cat * > $wkdir/48_new_taxonomy_count/02_class/06_cat_species/remanei

cd $wkdir/48_new_taxonomy_count/02_class/06_cat_species/

cat * > all





mkdir $wkdir/48_new_taxonomy_count/02_class/07_awk_class
mkdir $wkdir/48_new_taxonomy_count/02_class/08_count_rep_class_pdf

cd $wkdir/48_new_taxonomy_count/02_class/07_awk_class

#awk by class again
awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/48_new_taxonomy_count/02_class/06_cat_species/all


#add headers
for i in *; do echo -e "Chr\tBP\telement_count\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs [not included in the manuscript  ;  make sure to change output directory in repeat_landscapes_plots.R]
for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &

mv $wkdir/48_new_taxonomy_count/02_class/08_perc_N_rep_class_pdf $wkdir/48_new_taxonomy_count/02_class/08_count_rep_class_pdf

cd $wkdir/48_new_taxonomy_count/02_class/08_count_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf

#join to normalized distance to chr center and num_bp...


mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp


#copy info re count

cp -r $wkdir/48_new_taxonomy_count/02_class/03_bedtools_map $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/

cd $wkdir/48_new_taxonomy_count/09_join_norm_dist_center_num_bp

mv 03_bedtools_map 01_cp_count

#clean it up -- replace "." with 0

cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do sed -i -e 's/\./0/g' $i; done


mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk
mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/remanei/$i; done &








#now, order counts

mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class
mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/00_awk_class/remanei


cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/briggsae
awk '{print > $8}' $wkdir/48_new_taxonomy_count/01_awk_species/briggsae
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/elegans
awk '{print > $8}' $wkdir/48_new_taxonomy_count/01_awk_species/elegans
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/inopinata
awk '{print > $8}' $wkdir/48_new_taxonomy_count/01_awk_species/inopinata
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/nigoni
awk '{print > $8}' $wkdir/48_new_taxonomy_count/01_awk_species/nigoni
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/remanei
awk '{print > $8}' $wkdir/48_new_taxonomy_count/01_awk_species/remanei

#now, sort bed


mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort
mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/remanei

cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort//elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/remanei/$i; done &


#ok, now order

#awk one for sum with bedtools map

mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk
mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/03_order/02_awk/briggsae/$i; done
cd $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/03_order/02_awk/elegans/$i; done
cd $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/03_order/02_awk/inopinata/$i; done
cd $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/03_order/02_awk/nigoni/$i; done
cd $wkdir/48_new_taxonomy_count/03_order/01_bedtools_sort/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/03_order/02_awk/remanei/$i; done

#ok, cool now bedtools map

mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map
mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/remanei



cd $wkdir/48_new_taxonomy_count/03_order/02_awk/briggsae
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/briggsae.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/02_awk/elegans
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/elegans.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/02_awk/inopinata
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/inopinata.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/02_awk/nigoni
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/nigoni.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/03_order/02_awk/remanei
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/remanei.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/remanei/$i; done &

#grep to remove lines with no data

mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep
mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/04_grep/remanei



cd $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/briggsae
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/03_order/04_grep/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/elegans
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/03_order/04_grep/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/inopinata
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/03_order/04_grep/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/nigoni
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/03_order/04_grep/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map/remanei
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/03_order/04_grep/remanei/$i; done &

#awk to add in species and class data

mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk
mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/05_awk/remanei


cd $wkdir/48_new_taxonomy_count/03_order/04_grep/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/05_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/04_grep/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/05_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/04_grep/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/05_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/04_grep/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/05_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/04_grep/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/05_awk/remanei/$i; done &



#now cat species

mkdir $wkdir/48_new_taxonomy_count/03_order/06_cat_species

cd $wkdir/48_new_taxonomy_count/03_order/05_awk/briggsae/

cat * > $wkdir/48_new_taxonomy_count/03_order/06_cat_species/briggsae

cd $wkdir/48_new_taxonomy_count/03_order/05_awk/elegans/

cat * > $wkdir/48_new_taxonomy_count/03_order/06_cat_species/elegans


cd $wkdir/48_new_taxonomy_count/03_order/05_awk/inopinata/

cat * > $wkdir/48_new_taxonomy_count/03_order/06_cat_species/inopinata


cd $wkdir/48_new_taxonomy_count/03_order/05_awk/nigoni/

cat * > $wkdir/48_new_taxonomy_count/03_order/06_cat_species/nigoni

cd $wkdir/48_new_taxonomy_count/03_order/05_awk/remanei/

cat * > $wkdir/48_new_taxonomy_count/03_order/06_cat_species/remanei

cd $wkdir/48_new_taxonomy_count/03_order/06_cat_species/

cat * > all





mkdir $wkdir/48_new_taxonomy_count/03_order/07_awk_class
mkdir $wkdir/48_new_taxonomy_count/03_order/08_count_rep_class_pdf



cp $wkdir/48_new_taxonomy_count/03_order/repeat_classes_count_plots_4-24-19.R $wkdir/48_new_taxonomy_count/03_order/repeat_classes_count_plots_4-24-19.R

#awk by order again

cd $wkdir/48_new_taxonomy_count/03_order/07_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/48_new_taxonomy_count/03_order/06_cat_species/all


#add headers

for i in *; do echo -e "Chr\tBP\telement_count\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs [not included in the manuscript  ;  make sure to change output directory in repeat_landscapes_plots.R]
for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &

mv $wkdir/48_new_taxonomy_count/03_order/08_perc_N_rep_class_pdf $wkdir/48_new_taxonomy_count/03_order/08_count_rep_class_pdf

cd $wkdir/48_new_taxonomy_count/03_order/08_count_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf

#join to normalized distance to chr center and num_bp...


mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp

#copy info re count

cp -r $wkdir/48_new_taxonomy_count/03_order/03_bedtools_map $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp

mv 03_bedtools_map 01_cp_count

#clean that shit up -- replace "." with 0

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do sed -i -e 's/\./0/g' $i; done


mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk
mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/remanei/$i; done &



#now, superfamily

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/remanei


cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/briggsae
awk '{print > $9}' $wkdir/48_new_taxonomy_count/01_awk_species/briggsae
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/elegans
awk '{print > $9}' $wkdir/48_new_taxonomy_count/01_awk_species/elegans
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/inopinata
awk '{print > $9}' $wkdir/48_new_taxonomy_count/01_awk_species/inopinata
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/nigoni
awk '{print > $9}' $wkdir/48_new_taxonomy_count/01_awk_species/nigoni
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/remanei
awk '{print > $9}' $wkdir/48_new_taxonomy_count/01_awk_species/remanei


cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/briggsae
rm NA
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/elegans
rm NA
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/inopinata
rm NA
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/nigoni
rm NA
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/remanei
rm NA


#now, sort bed


mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/remanei

cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort//elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/remanei/$i; done &


#ok, now superfamily

#awk one for sum with bedtools map

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/briggsae/$i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/elegans/$i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/inopinata/$i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/nigoni/$i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/01_bedtools_sort/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/remanei/$i; done

#ok, cool now bedtools map

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/remanei



cd $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/briggsae
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/briggsae.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/elegans
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/elegans.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/inopinata
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/inopinata.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/nigoni
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/nigoni.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/04_superfamily/02_awk/remanei
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/remanei.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/remanei/$i; done &

#grep to remove lines with no data

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/remanei



cd $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/briggsae
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/elegans
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/inopinata
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/nigoni
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map/remanei
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/remanei/$i; done &

#awk to add in species and class data

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/remanei


cd $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/04_grep/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/remanei/$i; done &



#now cat species

mkdir $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species

cd $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/briggsae/

cat * > $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/briggsae

cd $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/elegans/

cat * > $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/elegans


cd $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/inopinata/

cat * > $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/inopinata


cd $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/nigoni/

cat * > $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/nigoni

cd $wkdir/48_new_taxonomy_count/04_superfamily/05_awk/remanei/

cat * > $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/remanei

cd $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/

cat * > all





mkdir $wkdir/48_new_taxonomy_count/04_superfamily/07_awk_class
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/08_count_rep_class_pdf



cp $wkdir/48_new_taxonomy_count/02_order/repeat_classes_count_plots_4-24-19.R $wkdir/48_new_taxonomy_count/04_superfamily/repeat_classes_count_plots_4-24-19.R


#awk by superfamily again

cd $wkdir/48_new_taxonomy_count/04_superfamily/07_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/48_new_taxonomy_count/04_superfamily/06_cat_species/all


#add headers

for i in *; do echo -e "Chr\tBP\telement_count\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs [not included in the manuscript  ;  make sure to change output directory in repeat_landscapes_plots.R]

for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &

mv $wkdir/48_new_taxonomy_count/04_superfamily/08_perc_N_rep_class_pdf $wkdir/48_new_taxonomy_count/04_superfamily/08_count_rep_class_pdf

cd $wkdir/48_new_taxonomy_count/04_superfamily/08_count_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf
#join to normalized distance to chr center and num_bp...


mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp

#copy info re norm dist center and class info

#copy info re count

cp -r $wkdir/48_new_taxonomy_count/04_superfamily/03_bedtools_map $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp

mv 03_bedtools_map 01_cp_count

#clean that shit up -- replace "." with 0

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do sed -i -e 's/\./0/g' $i; done


mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/remanei/$i; done &





#ok, now family

mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class
mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/00_awk_class/remanei


cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/briggsae
awk '{print > $10}' $wkdir/48_new_taxonomy_count/01_awk_species/briggsae
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/elegans
awk '{print > $10}' $wkdir/48_new_taxonomy_count/01_awk_species/elegans
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/inopinata
awk '{print > $10}' $wkdir/48_new_taxonomy_count/01_awk_species/inopinata
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/nigoni
awk '{print > $10}' $wkdir/48_new_taxonomy_count/01_awk_species/nigoni
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/remanei
awk '{print > $10}' $wkdir/48_new_taxonomy_count/01_awk_species/remanei


cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/briggsae
rm NA
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/elegans
rm NA
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/inopinata
rm NA
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/nigoni
rm NA
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/remanei
rm NA


#now, sort bed


mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort
mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/remanei

cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/briggsae
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/briggsae_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/elegans
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/elegans_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort//elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/inopinata
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/inopinata_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/nigoni
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/nigoni_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/00_awk_class/remanei
for i in *; do bedtools sort -faidx $wkdir/01_fasta_filter/remanei_genome.fa.fai -i $i > $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/remanei/$i; done &


#ok, now superfamily

#awk one for sum with bedtools map

mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk
mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/briggsae
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/05_family/02_awk/briggsae/$i; done
cd $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/elegans
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/05_family/02_awk/elegans/$i; done
cd $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/inopinata
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/05_family/02_awk/inopinata/$i; done
cd $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/nigoni
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/05_family/02_awk/nigoni/$i; done
cd $wkdir/48_new_taxonomy_count/05_family/01_bedtools_sort/remanei
for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"1"}' $i > $wkdir/48_new_taxonomy_count/05_family/02_awk/remanei/$i; done

#ok, cool now bedtools map

mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map
mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/remanei



cd $wkdir/48_new_taxonomy_count/05_family/02_awk/briggsae
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/briggsae.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/briggsae/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/02_awk/elegans
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/elegans.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/elegans/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/02_awk/inopinata
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/inopinata.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/inopinata/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/02_awk/nigoni
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/nigoni.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/nigoni/$i; done &
cd $wkdir/48_new_taxonomy_count/05_family/02_awk/remanei
for i in *; do bedtools map -o sum -c 11 -a $wkdir/windows/remanei.10kb.windows -b $i > $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/remanei/$i; done &

#grep to remove lines with no data

mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep
mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/04_grep/remanei



cd $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/briggsae
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/05_family/04_grep/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/elegans
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/05_family/04_grep/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/inopinata
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/05_family/04_grep/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/nigoni
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/05_family/04_grep/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map/remanei
for i in *; do grep -v '\.' $i > $wkdir/48_new_taxonomy_count/05_family/04_grep/remanei/$i; done &

#awk to add in species and class data

mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk
mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/05_awk/remanei


cd $wkdir/48_new_taxonomy_count/05_family/04_grep/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/05_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/04_grep/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/05_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/04_grep/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/05_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/04_grep/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/05_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/04_grep/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/05_awk/remanei/$i; done &



#now cat species

mkdir $wkdir/48_new_taxonomy_count/05_family/06_cat_species

cd $wkdir/48_new_taxonomy_count/05_family/05_awk/briggsae/

cat * > $wkdir/48_new_taxonomy_count/05_family/06_cat_species/briggsae

cd $wkdir/48_new_taxonomy_count/05_family/05_awk/elegans/

cat * > $wkdir/48_new_taxonomy_count/05_family/06_cat_species/elegans


cd $wkdir/48_new_taxonomy_count/05_family/05_awk/inopinata/

cat * > $wkdir/48_new_taxonomy_count/05_family/06_cat_species/inopinata


cd $wkdir/48_new_taxonomy_count/05_family/05_awk/nigoni/

cat * > $wkdir/48_new_taxonomy_count/05_family/06_cat_species/nigoni

cd $wkdir/48_new_taxonomy_count/05_family/05_awk/remanei/

cat * > $wkdir/48_new_taxonomy_count/05_family/06_cat_species/remanei

cd $wkdir/48_new_taxonomy_count/05_family/06_cat_species/

cat * > all





mkdir $wkdir/48_new_taxonomy_count/05_family/07_awk_class
mkdir $wkdir/48_new_taxonomy_count/05_family/08_count_rep_class_pdf

#awk by superfamily again

cd $wkdir/48_new_taxonomy_count/05_family/07_awk_class

awk 'BEGIN {FS="\t"} {OFS="\t"} {print > $5}' $wkdir/48_new_taxonomy_count/05_family/06_cat_species/all


#add headers

for i in *; do echo -e "Chr\tBP\telement_count\tspecies\trep_class" | cat - $i > $i.tmp; done

for i in *; do mv $i.tmp $i; done


#make pdfs [not included in the manuscript  ;  make sure to change output directory in repeat_landscapes_plots.R]

for i in *; do Rscript $wkdir/scripts/repeat_landscapes_plots.R $i; done &

mv $wkdir/48_new_taxonomy_count/05_family/08_perc_N_rep_class_pdf $wkdir/48_new_taxonomy_count/05_family/08_count_rep_class_pdf

cd $wkdir/48_new_taxonomy_count/05_family/08_count_rep_class_pdf

rm NA.pdf

pdfunite *.pdf all_repeat_classes.pdf

#join to normalized distance to chr center and num_bp...


mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp

#copy info re norm dist center and class info

#copy info re count

cp -r $wkdir/48_new_taxonomy_count/05_family/03_bedtools_map $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp

mv 03_bedtools_map 01_cp_count

#clean that shit up -- replace "." with 0

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do sed -i -e 's/\./0/g' $i; done
cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do sed -i -e 's/\./0/g' $i; done


mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk
mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/briggsae
mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/elegans
mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/inopinata
mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/nigoni
mkdir $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/remanei


cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/briggsae
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"briggsae", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/briggsae/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/elegans
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"elegans", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/elegans/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/inopinata
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/inopinata/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/nigoni
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"nigoni", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/nigoni/$i; done &

cd $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/01_cp_count/remanei
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"remanei", FILENAME}' $i > $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/remanei/$i; done &


#############
#############
#############
#############
#merging repeat counts, bp repeat, and normalized chromosome positions for data analysis (on original file "pipeline_continued_1-2-19.sh", merge steps start at line 3540)
#############
#############
#############
#############


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp

#get normalized chr positions for all species
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center

cd $wkdir/43_bedtools_nuc_perc_N_3-11-19/13_cat/

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2,$5}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/$i; done

#remove headers

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center

rm all

for i in *; do sed -i '1d' $i; done

#add ~ to create key for merge

for i in *; do sed -i -e 's/\t/~/' $i; done

#now merge counts, bp of repeat, and normalized distances to chromsome centers for all repeat classses


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp

cp -r $wkdir/48_new_taxonomy_count/02_class/09_join_norm_dist_center_num_bp/02_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/


cp -r $wkdir/61_repeat_taxon_num_bp_disjoined_gff/02_class/04_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/








#replace first tab with ~ to create key for merge
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/briggsae
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/elegans
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/inopinata
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/nigoni
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/remanei
for i in *; do sed -i -e 's/\t/~/' $i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/remanei

#merge counts with bp repeats

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/02_awk/briggsae
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/briggsae/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/briggsae/$i; done &


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/02_awk/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/elegans/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/elegans/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/02_awk/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/inopinata/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/inopinata/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/02_awk/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/nigoni/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/nigoni/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/00_cp_counts/02_awk/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/01_cp_num_bp/04_awk/remanei/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/remanei/$i; done &


#connect that data with normalized distance from chromosome center...

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/briggsae/
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/briggsae 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/elegans 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/inopinata 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/nigoni 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/02_merge_count_num_rep/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/remanei 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/remanei/$i; done &

#ok, start to get data right


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/briggsae/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"class"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/briggsae/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/elegans/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"class"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/elegans/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/inopinata/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"class"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/inopinata/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/nigoni/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"class"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/nigoni/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/03_merge_norm_dist_center/remanei/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"class"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/remanei/$i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/briggsae


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/elegans

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/inopinata

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/remanei

#fix that ~

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/

for i in *; do sed -i -e 's/~/\t/' $i; done

cat * > all_classes

#replace cryptic class labels
awk 'BEGIN {FS="\t"} {OFS="\t"}$7 == "I" { $7="class_I_retrotransposon" }1' all_classes > all_classes.tmp1
awk 'BEGIN {FS="\t"} {OFS="\t"}$7 == "II" { $7="class_II_DNA_transposon" }1' all_classes.tmp1 > all_classes.tmp2
awk 'BEGIN {FS="\t"} {OFS="\t"}$7 == "NA" { $7="Unclassified" }1' all_classes.tmp2 > all_classes

#add header
echo -e "Chr\tBP\tnorm_dist_chr_center\tnum_bp_rep\tnum_rep_count\tspecies\trep_class\ttaxonomic_rank" | cat - all_classes > all_classes.tsv

rm all_classes.tmp1
rm all_classes.tmp2


#orders

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp

cp -r $wkdir/48_new_taxonomy_count/03_order/09_join_norm_dist_center_num_bp/02_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/


cp -r $wkdir/61_repeat_taxon_num_bp_disjoined_gff/03_order/04_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/briggsae
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/elegans
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/inopinata
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/nigoni
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/remanei
for i in *; do sed -i -e 's/\t/~/' $i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/02_awk/briggsae
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/briggsae/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/02_awk/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/elegans/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/elegans/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/02_awk/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/inopinata/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/inopinata/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/02_awk/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/nigoni/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/nigoni/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/00_cp_counts/02_awk/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/01_cp_num_bp/04_awk/remanei/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/remanei/$i; done &


#connect with norm dist center...

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/briggsae/
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/briggsae 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/elegans 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/inopinata 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/nigoni 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/02_merge_count_num_rep/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/remanei 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/remanei/$i; done &

#ok, start to get data right


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/briggsae/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/briggsae/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/elegans/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/elegans/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/inopinata/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/inopinata/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/nigoni/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/nigoni/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/03_merge_norm_dist_center/remanei/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/remanei/$i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/briggsae


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/elegans

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/inopinata

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/remanei

#fix the shit

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/

for i in *; do sed -i -e 's/~/\t/' $i; done

cat * > all_orders


echo -e "Chr\tBP\tnorm_dist_chr_center\tnum_bp_rep\tnum_rep_count\tspecies\trep_class\ttaxonomic_rank" | cat -  all_orders >  all_orders.tsv



#superfamilies

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp

cp -r $wkdir/48_new_taxonomy_count/04_superfamily/09_join_norm_dist_center_num_bp/02_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/


cp -r $wkdir/61_repeat_taxon_num_bp_disjoined_gff/04_superfamily/04_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/briggsae
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/elegans
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/inopinata
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/nigoni
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/remanei
for i in *; do sed -i -e 's/\t/~/' $i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/02_awk/briggsae
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/briggsae/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/02_awk/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/elegans/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/elegans/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/02_awk/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/inopinata/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/inopinata/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/02_awk/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/nigoni/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/nigoni/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/00_cp_counts/02_awk/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/01_cp_num_bp/04_awk/remanei/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/remanei/$i; done &


#connect with norm dist center...

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/briggsae/
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/briggsae 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/elegans 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/inopinata 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/nigoni 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/02_merge_count_num_rep/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/remanei 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/remanei/$i; done &

#ok, start to get data right


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/briggsae/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/briggsae/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/elegans/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/elegans/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/inopinata/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/inopinata/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/nigoni/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/nigoni/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/03_merge_norm_dist_center/remanei/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/remanei/$i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/briggsae


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/elegans

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/inopinata

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/remanei

#fix the shit

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/

for i in *; do sed -i -e 's/~/\t/' $i; done

cat * > all_orders


echo -e "Chr\tBP\tnorm_dist_chr_center\tnum_bp_rep\tnum_rep_count\tspecies\trep_class\ttaxonomic_rank" | cat -  all_orders >  all_orders.tsv



#families

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp

cp -r $wkdir/48_new_taxonomy_count/05_family/09_join_norm_dist_center_num_bp/02_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/


cp -r $wkdir/61_repeat_taxon_num_bp_disjoined_gff/05_family/04_awk/ $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/briggsae
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/elegans
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/inopinata
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/nigoni
for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/remanei
for i in *; do sed -i -e 's/\t/~/' $i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/02_awk/briggsae
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/briggsae/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/briggsae/$i; done &



cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/02_awk/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/elegans/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/elegans/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/02_awk/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/inopinata/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/inopinata/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/02_awk/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/nigoni/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/nigoni/$i; done &

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/00_cp_counts/02_awk/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/01_cp_num_bp/04_awk/remanei/$i 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/remanei/$i; done &


#connect with norm dist center...

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/remanei


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/briggsae/
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/briggsae 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/briggsae/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/elegans
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/elegans 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/elegans/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/inopinata
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/inopinata 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/inopinata/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/nigoni
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/nigoni 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/nigoni/$i; done &
cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/02_merge_count_num_rep/remanei
for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/remanei 2> $wkdir/$i.merge_pl.error > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/remanei/$i; done &

#ok, start to get data right


mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/briggsae/
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/elegans
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/inopinata
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/nigoni
mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/remanei

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/briggsae/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/briggsae/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/elegans/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/elegans/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/inopinata/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/inopinata/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/nigoni/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/nigoni/$i; done


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/03_merge_norm_dist_center/remanei/
for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$8,$5,$2,$3,$4,"order"}' $i > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/remanei/$i; done

mkdir $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/briggsae/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/briggsae


cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/elegans/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/elegans

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/inopinata/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/inopinata

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/nigoni/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/nigoni

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/04_awk/remanei/

cat * > $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/remanei

#fix the shit

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/

for i in *; do sed -i -e 's/~/\t/' $i; done

cat * > all_orders


echo -e "Chr\tBP\tnorm_dist_chr_center\tnum_bp_rep\tnum_rep_count\tspecies\trep_class\ttaxonomic_rank" | cat -  all_orders >  all_orders.tsv



#combine all data

#make header-less versions

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/

sed '1d' all_classes > all_classes_no_header

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/

sed '1d' all_orders > all_orders_no_header

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/

sed '1d' all_superfamilies > all_superfamilies_no_header

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/

sed '1d' all_families > all_families_no_header

#combine

cd $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/

cat $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/00_class/05_cat_species/all_classes_no_header $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/01_order/05_cat_species/all_orders_no_header $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/02_superfamily/05_cat_species/all_superfamilies_no_header $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/03_family/05_cat_species/all_families_no_header > all_repeat_taxa.tmp

#add header
echo -e "Chr\tBP\tnorm_dist_chr_center\tnum_bp_rep\tnum_rep_count\tspecies\trep_class\ttaxonomic_rank" | cat -  all_repeat_taxa.tmp >  all_repeat_taxa_density.tsv #[formerly all_repeat_taxa_count_num_bp_disjoined_10kb_win_8-16-19.tsv]
	#THIS is the data used for Figure 3 and 4; all statistics regarding repeat taxa and their chromosomal distributions

#make a version of this data to include all "0"'s for every genomic window in briggsae for repeat superfamily RTE (it does not exist in briggsae) for constructing figure 3. it is called "all_repeat_taxa_count_num_bp_disjoined_10kb_win_plus_RTE_briggsae_zeros_8-17-19.tsv"
	#THIS is the modified data used for Figure 3 specifically

#############
#############
#############
#############
#removing four superfamilies to see their impact on the repetitive landscape in C. inopinata 
#############
#############
#############
#############


mkdir $wkdir/49_remove_repeat_families_4-30-19

mkdir $wkdir/49_remove_repeat_families_4-30-19//00_links

cd  $wkdir/29_prep_plots/repeatmasker_ii//07_paste

for i in *; do cp $i $wkdir/49_remove_repeat_families_4-30-19//00_links/$i; done

mkdir $wkdir/49_remove_repeat_families_4-30-19//01_awk

cd $wkdir/49_remove_repeat_families_4-30-19//00_links

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $11,$1,$4,$5,"briggsae"}' briggsae > $wkdir/49_remove_repeat_families_4-30-19//01_awk/briggsae
awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $11,$1,$4,$5,"elegans"}' elegans > $wkdir/49_remove_repeat_families_4-30-19//01_awk/elegans
awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $11,$1,$4,$5,"inopinata"}' inopinata > $wkdir/49_remove_repeat_families_4-30-19//01_awk/inopinata
awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $11,$1,$4,$5,"nigoni"}' nigoni > $wkdir/49_remove_repeat_families_4-30-19//01_awk/nigoni
awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $11,$1,$4,$5,"remanei"}' remanei > $wkdir/49_remove_repeat_families_4-30-19//01_awk/remanei

#the goal here is to remove Tc1-Mariner, RTE, and Bel-Pao superfamilies


mkdir $wkdir/49_remove_repeat_families_4-30-19//02_remove_Tc1-Mariner
mkdir $wkdir/49_remove_repeat_families_4-30-19//03_remove_RTE
mkdir $wkdir/49_remove_repeat_families_4-30-19//04_Bel-Pao


cd $wkdir/49_remove_repeat_families_4-30-19//01_awk

for i in *; do grep -v "DNA-TcMar" $i > $wkdir/49_remove_repeat_families_4-30-19//02_remove_Tc1-Mariner/$i; done

cd $wkdir/49_remove_repeat_families_4-30-19//02_remove_Tc1-Mariner

for i in *; do grep -v "LINE-RTE" $i > $wkdir/49_remove_repeat_families_4-30-19//03_remove_RTE/$i; done

cd $wkdir/49_remove_repeat_families_4-30-19//03_remove_RTE

for i in *; do grep -v "LTR-Pao" $i > $wkdir/49_remove_repeat_families_4-30-19//04_Bel-Pao/$i; done

mkdir $wkdir/49_remove_repeat_families_4-30-19//05b_remove_Gypsy

cd $wkdir/49_remove_repeat_families_4-30-19//04_Bel-Pao
for i in *; do grep -v "LTR-Gypsy" $i > $wkdir/49_remove_repeat_families_4-30-19//05b_remove_Gypsy/$i; done

mkdir $wkdir/49_remove_repeat_families_4-30-19//05_awk

cd $wkdir/49_remove_repeat_families_4-30-19//05b_remove_Gypsy


for i in *; do awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $2,$3,$4,$5,$1}' $i > $wkdir/49_remove_repeat_families_4-30-19//05_awk/$i; done

#maskfasta

mkdir $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta

cd $wkdir/49_remove_repeat_families_4-30-19//05_awk
bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/briggsae_genome.fa -bed briggsae -fo $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta/briggsae &


bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/elegans_genome.fa -bed elegans -fo $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta/elegans & 


bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/inopinata_genome.fa -bed inopinata -fo $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta/inopinata &


bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/nigoni_genome.fa -bed nigoni -fo $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta/nigoni &

bedtools maskfasta -fi $wkdir/43_bedtools_nuc_perc_N_3-11-19/01_sed_genomes_X/remanei_genome.fa -bed remanei -fo $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta/remanei &

#bedtools nuc

mkdir $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc

cd $wkdir/49_remove_repeat_families_4-30-19//06_maskfasta

bedtools nuc -fi briggsae -bed $wkdir/windows/briggsae.10kb.windows > $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc/briggsae &



bedtools nuc -fi elegans -bed $wkdir/windows/elegans.10kb.windows > $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc/elegans &



bedtools nuc -fi inopinata -bed $wkdir/windows/inopinata.10kb.windows > $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc/inopinata &


bedtools nuc -fi nigoni -bed $wkdir/windows/nigoni.10kb.windows > $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc/nigoni &

bedtools nuc -fi remanei -bed $wkdir/windows/remanei.10kb.windows > $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc/remanei &

#get columns needed

mkdir $wkdir/49_remove_repeat_families_4-30-19//08_awk

cd $wkdir/49_remove_repeat_families_4-30-19//07_bedtools_nuc

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$10,"briggsae"}' briggsae > $wkdir/49_remove_repeat_families_4-30-19//08_awk/briggsae

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$10,"elegans"}' elegans > $wkdir/49_remove_repeat_families_4-30-19//08_awk/elegans

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$10,"inopinata"}' inopinata > $wkdir/49_remove_repeat_families_4-30-19//08_awk/inopinata

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$10,"nigoni"}' nigoni > $wkdir/49_remove_repeat_families_4-30-19//08_awk/nigoni

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$10,"remanei"}' remanei > $wkdir/49_remove_repeat_families_4-30-19//08_awk/remanei

#remove header

cd $wkdir/49_remove_repeat_families_4-30-19//08_awk
for i in *; do sed -i '1d' $i; done

#join for norm dist chr center

mkdir $wkdir/49_remove_repeat_families_4-30-19//09_sed_tilde

for i in *; do sed 's/\t/~/' $i >  $wkdir/49_remove_repeat_families_4-30-19//09_sed_tilde/$i; done

#join...

mkdir $wkdir/49_remove_repeat_families_4-30-19//10_merge

cd $wkdir/49_remove_repeat_families_4-30-19//09_sed_tilde/

for i in *; do perl $wkdir/scripts/merge.pl -k -e "no_key" $i $wkdir/61_repeat_taxon_num_bp_disjoined_gff/06_merge_count_num_bp/norm_dist_chrom_center/$i 2> $wkdir/49_remove_repeat_families_4-30-19//$i.merge_pl.error > $wkdir/49_remove_repeat_families_4-30-19//10_merge/$i; done &

##clean up, cat, add header

cd $wkdir/49_remove_repeat_families_4-30-19//10_merge

for i in *; do grep -v "no_key" $i >$i.tmp; done

for i in *.tmp; do sed -i -e 's/~/\t/' $i; done

cat *.tmp > all.tmp

mv all.tmp all
mv briggsae.tmp briggsae
mv elegans.tmp elegans
mv inopinata.tmp inopinata
mv nigoni.tmp nigoni
mv remanei.tmp remanei

for i in *; do echo -e "Chr\tBP\tnum_rep_bp\tspecies\tnorm_dist_center" | cat - $i > $i.tmp && mv $i.tmp $i; done

mv all global_repeat_density_remove_four_superfamilies.tsv # formerly "tot_bp_rep_remove_four_superfamilies_4-30-19.tsv"
	#THIS is the data used for Figure 5!
}


#############
#############
#############
#############
#Aligning protein-coding genes to transposons
#############
#############
#############
#############

#for versions of protein sets retrieved see file "genome_gff_protein_software_versions.txt"

#protein sets were filtered to remove alternative splice variants to create "canonical" protein sets. here, the longest isoform per predicted gene was retained. these were prepared as in Woodruff 2019 "Patterns of putative gene loss suggest rampant developmental system drift in nematodes" (https://www.biorxiv.org/content/10.1101/627620v1.abstract) see (https://github.com/gcwoodruff/gene_loss/blob/master/prepare_protein_sets.sh) for how these protein sets were generated.



mkdir $wkdir/30_blastp

cd $wkdir/30_blastp

mkdir $wkdir/30_blastp/00_links
	#put the canonical protein sets in this directory!

#make a folder for blast databases
mkdir 01_makeblastdb
mkdir transposon_prot_db

#copy the transposonPSI protein database to directory transposon_prot_db
#in my installation of transposonPSI, it was in ".../transposonpsi/transposon_ORF_lib/transposon_db.pep" This is in included in folder "additional_repeat_libraries"

cp /packages/transposonpsi/transposon_ORF_lib/transposon_db.pep $wkdir/30_blastp/transposon_prot_db/transposon_db.fa

#get the filenames right for makeblastdb... there were a number of sequence id's it did not like

sed -i -e 's/\[.*//g' transposon_db.fa

sed -i -e 's/X02417 ORF2/X02417_b ORF2/g' transposon_db.fa

sed -i 's/Pi_model.supercont1.18.5.1 1572/Pi_model.supercont1.18.5.1_b 1572/g' transposon_db.fa

sed -i 's/Pi_model.supercont1.7.16.1 1499/Pi_model.supercont1.7.16.1_b 1499/g' transposon_db.fa

sed -i 's/Pi_model.supercont1.6.8.1 1637/Pi_model.supercont1.6.8.1_b 1637/g' transposon_db.fa

sed -i '0,/56547708_bc/s/56547708_bc/56547708_bcd/1' transposon_db.fa

sed -i 's/gi|56547708_bcd|gb|AAV92925.1| /gi_56547708_bcd_gb_AAV92925.1_/g' transposon_db.fa

sed -i 's/gi|56547708_bc|gb|AAV92925.1| /gi_56547708_bc_gb_AAV92925.1_/g' transposon_db.fa

#make the blast db for the transposonPSI transposon protein set ; this includes 1537 protein sequences

makeblastdb -in $wkdir/30_blastp/transposon_prot_db/transposon_db.fa -parse_seqids -out $wkdir/30_blastp/01_makeblastdb/transposon_db -title transposon_db -dbtype prot

#blastp for all canonical Caenorhabditis (plus Diploscapter coronatus) protein sets to the transposonPSI database
mkdir $wkdir/30_blastp/02_blastp

cd $wkdir/30_blastp/00_links

for i in *; do blastp -db $wkdir/30_blastp/01_makeblastdb/transposon_db -outfmt 6 -query $i -out $wkdir/30_blastp/02_blastp/${i%.fa}.blast_out -evalue 0.001; done &


#get unique proteins that hit

mkdir $wkdir/30_blastp/03_cut_uniq

cd $wkdir/30_blastp/02_blastp

for i in *; do cut -f1 $i | uniq >  $wkdir/30_blastp/03_cut_uniq/${i%.blast_out}; done

#get number of unique proteins that align to transposons per species
cd $wkdir/30_blastp/ 03_cut_uniq

for i in *; do wc -l $i; done


#number of unique proteins aligning to transposons in transposonPSI library 

#162 afra
#158 angaria
#416 brenneri
#290 briggsae
#164 castelli
#445 D_coronatus
#396 doughertyi
#100 elegans
#3619 japonica
#197 kamaaina
#269 latens
#167 monodelphis
#776 nigoni
#117 plicata
#315 remanei
#653 sinica
#165 sp21
#254 sp26
#134 sp28
#187 sp29
#206 sp31
#176 sp32
#3349 sp34
#261 sp38
#161 sp39
#312 sp40
#211 tropicalis
#217 virilis
#368 wallacei
	#this is in the data used for supplemental figure 18 ; see file proteins_hit_transposons.tsv

#get number of proteins in all data sets

cd $wkdir/30_blastp/00_links

for i in *; do grep -c "^>" $i /dev/null; done

#total number of proteins per species

#afra.fa:19834
#angaria.fa:27970
#brenneri.fa:30660
#briggsae.fa:22405
#castelli.fa:19694
#D_coronatus.fa:34421
#doughertyi.fa:30860
#elegans.fa:20236
#japonica.fa:29931
#kamaaina.fa:19448
#latens.fa:24855
#monodelphis.fa:20262
#nigoni.fa:29167
#plicata.fa:14843
#remanei.fa:24867
#sinica.fa:34696
#sp21.fa:16412
#sp26.fa:22198
#sp28.fa:17300
#sp29.fa:21924
#sp31.fa:27614
#sp32.fa:18192
#sp34.fa:21609
#sp38.fa:22278
#sp39.fa:30124
#sp40.fa:24787
#tropicali:22326
#virilis.fa:19447
#wallacei.fa:20628
	#this is in the data used for supplemental figure 18 ; it is in file "proteins_hit_transposons.tsv"

#get the sp34 top hits and their respective ids


mkdir 04_awk_sp34
cd 04_awk_sp34

awk '{print > $1}' $wkdir/30_blastp/02_blastp/sp34.blast_out

cd ..

mkdir 05_head_sp34
cd 04_awk_sp34
for i in *; do head -1 $i > $wkdir/30_blastp/05_head_sp34/$i; done

cd ..

mkdir 06_cat_sp34
cd 05_head_sp34
cat * > $wkdir/30_blastp/06_cat_sp34/sp34_top_hits

awk '{print $2}' sp34_top_hits | sort | uniq > sp34_top_db_hits_ids

grep "GB:" sp34_top_db_hits_ids | sed 's/GB://g'


#############
#############
#############
#############
#Gene density
#############
#############
#############
#############


#gene densities

#gff versions are in genome_gff_protein_software_versions.txt
#gff's were largely previously cleaned up to rename chromosome scaffolds, remove mitochondrial dna, and just include "genes" as defined by third column in gff. gff's that were not previously cleaned up in this way are done as below.

mkdir $wkdir/39_gene_density

mkdir $wkdir/39_gene_density/00_gff
	#put the annotations here



mkdir $wkdir/39_gene_density/01_prep_gff

mkdir $wkdir/39_gene_density/01_prep_gff/briggsae

mkdir $wkdir/39_gene_density/01_prep_gff/nigoni

#

#split by scaffold

cd $wkdir/39_gene_density/01_prep_gff/briggsae

awk '{print > $1}' $wkdir/39_gene_density/00_gff/caenorhabditis_briggsae.PRJNA10731.WBPS12.annotations.gff3 &

cd $wkdir/39_gene_density/01_prep_gff/nigoni

awk '{print > $1}' $wkdir/39_gene_density/00_gff/caenorhabditis_nigoni.PRJNA384657.WBPS12.annotations.gff3 &

#cat only the relevant scaffolds


cd $wkdir/39_gene_density/01_prep_gff/briggsae/

cat I II III IV V X > $wkdir/39_gene_density/01_prep_gff/briggsae_cat.gff

cd $wkdir/39_gene_density/01_prep_gff/nigoni

cat CM008509.1 CM008510.1 CM008511.1 CM008512.1 CM008513.1 CM008514.1  > $wkdir/39_gene_density/01_prep_gff/nigoni_cat.gff

#remove mitochondrial DNA from remanei gff

grep -v "MtDNA" $wkdir/39_gene_density/00_gff/remanei.gff > $wkdir/39_gene_density/01_prep_gff/remanei_no_mito.gff

#get just the genes...

cd $wkdir/39_gene_density/01_prep_gff/

awk '$3 == "gene"' briggsae_cat.gff > briggsae.gff

awk '$3 == "gene"' nigoni_cat.gff > nigoni.gff

awk '$3 == "gene"' remanei_no_mito.gff > remanei.gff

#replace nigoni scaffold name with chr roman numerals

sed -i -e 's/CM008509.1/I/g' nigoni.gff
sed -i -e 's/CM008510.1/II/g' nigoni.gff
sed -i -e 's/CM008511.1/III/g' nigoni.gff
sed -i -e 's/CM008512.1/IV/g' nigoni.gff
sed -i -e 's/CM008513.1/V/g' nigoni.gff
sed -i -e 's/CM008514.1/X/g' nigoni.gff

#remove other things

rm -r briggsae
rm -r nigoni
rm briggsae_cat.gff
rm nigoni_cat.gff
rm remanei_no_mito.gff

#link other gff

ln -s $wkdir/39_gene_density/00_gff/elegans.gff $wkdir/39_gene_density/01_prep_gff/elegans.gff
ln -s $wkdir/39_gene_density/00_gff/inopinata.gff $wkdir/39_gene_density/01_prep_gff/inopinata.gff

#ok, get cols of interest to make it a bed file; add "1" for every gene so it can be summed...

mkdir $wkdir/39_gene_density/02_awk/

cd $wkdir/39_gene_density/01_prep_gff/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,"1"}'  $i > $wkdir/39_gene_density/02_awk/${i%.gff}.bed; done

#bedtools sort

mkdir $wkdir/39_gene_density/03_bedtools_sort/

cd $wkdir/39_gene_density/02_awk/

for i in *; do bedtools sort -i $i > $wkdir/39_gene_density/03_bedtools_sort/$i; done


#bedtools map to sum up number of genes per window

mkdir $wkdir/39_gene_density/04_bedtools_map/

mkdir $wkdir/39_gene_density/04_bedtools_map/10_kb_windows

cd $wkdir/39_gene_density/03_bedtools_sort/

	#didn't change the chromosome names for inopinata
sed -i -e 's/CSP34.Sp34_Chr1/I/g' inopinata.bed
sed -i -e 's/CSP34.Sp34_Chr2/II/g' inopinata.bed
sed -i -e 's/CSP34.Sp34_Chr3/III/g' inopinata.bed
sed -i -e 's/CSP34.Sp34_Chr4/IV/g' inopinata.bed
sed -i -e 's/CSP34.Sp34_Chr5/V/g' inopinata.bed
sed -i -e 's/CSP34.Sp34_ChrX/X/g' inopinata.bed


bedtools map -b briggsae.bed -a $wkdir/windows/briggsae.10kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/briggsae


bedtools map -b elegans.bed -a $wkdir/windows/elegans.10kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/elegans



bedtools map -b inopinata.bed -a $wkdir/windows/inopinata.10kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/inopinata



bedtools map -b nigoni.bed -a $wkdir/windows/nigoni.10kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/nigoni



bedtools map -b remanei.bed -a $wkdir/windows/remanei.10kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/remanei

#replace "." with 0

cd $wkdir/39_gene_density/04_bedtools_map/

for i in *; do sed -i -e 's/\./0/g' $i; done

################

#briefly, all_N repeats with 100kb windows... ; this is for paired data with gene densities in 100kb windows

mkdir $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows


bedtools nuc -fi $wkdir/27_RepeatMasker_II/hard/briggsae_genome.fa.masked -bed $wkdir/windows/briggsae.100kb.windows > $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/briggsae

bedtools nuc -fi $wkdir/27_RepeatMasker_II/hard/elegans_genome.fa.masked -bed $wkdir/windows/elegans.100kb.windows > $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/elegans

bedtools nuc -fi $wkdir/27_RepeatMasker_II/hard/inopinata_genome.fa.masked -bed $wkdir/windows/inopinata.100kb.windows > $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/inopinata

bedtools nuc -fi $wkdir/27_RepeatMasker_II/hard/nigoni_genome.fa.masked -bed $wkdir/windows/nigoni.100kb.windows > $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/nigoni

bedtools nuc -fi $wkdir/27_RepeatMasker_II/hard/remanei_genome.fa.masked -bed $wkdir/windows/remanei.100kb.windows > $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/remanei

#awk sed and cat to make new 100kb window %N file

mkdir $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows

cd $wkdir/33_remanei_prep_plots/15_bedtools_nuc/100_kb_windows/

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$10,FILENAME}' $i > $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/$i; done

cd $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows

#remove and add header
for i in *; do sed -i '1d' $i; done

cat * > all_N

echo -e "Chr\tBP\tnum_rep\tspecies" | cat - all_N > all_N_tmp && mv all_N_tmp all_N

##########

#now bedtools map the gene density with 100bp windows

mkdir $wkdir/39_gene_density/04_bedtools_map/100_kb_windows

cd $wkdir/39_gene_density/03_bedtools_sort/


bedtools map -b briggsae.bed -a $wkdir/windows/briggsae.100kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/100_kb_windows/briggsae


bedtools map -b elegans.bed -a $wkdir/windows/elegans.100kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/100_kb_windows/elegans



bedtools map -b inopinata.bed -a $wkdir/windows/inopinata.100kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/100_kb_windows/inopinata



bedtools map -b nigoni.bed -a $wkdir/windows/nigoni.100kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/100_kb_windows/nigoni



bedtools map -b remanei.bed -a $wkdir/windows/remanei.100kb.windows -c 4 -o sum > $wkdir/39_gene_density/04_bedtools_map/100_kb_windows/remanei

mkdir $wkdir/39_gene_density/04_bedtools_map/10_kb_windows

mv $wkdir/39_gene_density/04_bedtools_map/briggsae  $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/briggsae
mv $wkdir/39_gene_density/04_bedtools_map/elegans $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/elegans
mv $wkdir/39_gene_density/04_bedtools_map/inopinata $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/inopinata
mv $wkdir/39_gene_density/04_bedtools_map/nigoni $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/nigoni
mv $wkdir/39_gene_density/04_bedtools_map/remanei $wkdir/39_gene_density/04_bedtools_map/10_kb_windows/remanei

cd $wkdir/39_gene_density/04_bedtools_map/100_kb_windows
for i in *; do sed -i -e 's/\./0/g' $i; done

#awk

mkdir $wkdir/39_gene_density/05_awk/
mkdir $wkdir/39_gene_density/05_awk/10_kb_windows
mkdir $wkdir/39_gene_density/05_awk/100_kb_windows

cd $wkdir/39_gene_density/04_bedtools_map/10_kb_windows

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,FILENAME}' $i > $wkdir/39_gene_density/05_awk/10_kb_windows/$i; done

cd $wkdir/39_gene_density/04_bedtools_map/100_kb_windows

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,FILENAME}' $i > $wkdir/39_gene_density/05_awk/100_kb_windows/$i; done

#alright, prep merging gene density and repeat density data sets...

mkdir $wkdir/39_gene_density/06_prep_merge/
mkdir $wkdir/39_gene_density/06_prep_merge/10_kb_windows
mkdir $wkdir/39_gene_density/06_prep_merge/100_kb_windows

ln -s $wkdir/39_gene_density/05_awk/10_kb_windows/briggsae $wkdir/39_gene_density/06_prep_merge/10_kb_windows/briggsae_gene_density
ln -s $wkdir/39_gene_density/05_awk/10_kb_windows/elegans $wkdir/39_gene_density/06_prep_merge/10_kb_windows/elegans_gene_density
ln -s $wkdir/39_gene_density/05_awk/10_kb_windows/inopinata $wkdir/39_gene_density/06_prep_merge/10_kb_windows/inopinata_gene_density
ln -s $wkdir/39_gene_density/05_awk/10_kb_windows/nigoni $wkdir/39_gene_density/06_prep_merge/10_kb_windows/nigoni_gene_density
ln -s $wkdir/39_gene_density/05_awk/10_kb_windows/remanei $wkdir/39_gene_density/06_prep_merge/10_kb_windows/remanei_gene_density

ln -s $wkdir/39_gene_density/05_awk/100_kb_windows/briggsae $wkdir/39_gene_density/06_prep_merge/100_kb_windows/briggsae_gene_density
ln -s $wkdir/39_gene_density/05_awk/100_kb_windows/elegans $wkdir/39_gene_density/06_prep_merge/100_kb_windows/elegans_gene_density
ln -s $wkdir/39_gene_density/05_awk/100_kb_windows/inopinata $wkdir/39_gene_density/06_prep_merge/100_kb_windows/inopinata_gene_density
ln -s $wkdir/39_gene_density/05_awk/100_kb_windows/nigoni $wkdir/39_gene_density/06_prep_merge/100_kb_windows/nigoni_gene_density
ln -s $wkdir/39_gene_density/05_awk/100_kb_windows/remanei $wkdir/39_gene_density/06_prep_merge/100_kb_windows/remanei_gene_density

ln -s $wkdir/33_remanei_prep_plots/16_awk/briggsae $wkdir/39_gene_density/06_prep_merge/10_kb_windows/briggsae_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/elegans $wkdir/39_gene_density/06_prep_merge/10_kb_windows/elegans_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/inopinata $wkdir/39_gene_density/06_prep_merge/10_kb_windows/inopinata_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/nigoni $wkdir/39_gene_density/06_prep_merge/10_kb_windows/nigoni_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/remanei $wkdir/39_gene_density/06_prep_merge/10_kb_windows/remanei_perc_n

ln -s $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/briggsae $wkdir/39_gene_density/06_prep_merge/100_kb_windows/briggsae_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/elegans $wkdir/39_gene_density/06_prep_merge/100_kb_windows/elegans_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/inopinata $wkdir/39_gene_density/06_prep_merge/100_kb_windows/inopinata_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/nigoni $wkdir/39_gene_density/06_prep_merge/100_kb_windows/nigoni_perc_n
ln -s $wkdir/33_remanei_prep_plots/16_awk/100_kb_windows/remanei $wkdir/39_gene_density/06_prep_merge/100_kb_windows/remanei_perc_n


#replace first tab with ~ for joining

cd $wkdir/39_gene_density/06_prep_merge/10_kb_windows

for i in *; do sed -i -e 's/\t/~/' $i; done

cd $wkdir/39_gene_density/06_prep_merge/100_kb_windows

for i in *; do sed -i -e 's/\t/~/' $i; done


##merge ; cp /projects/phillipslab/gavincw/reproducing_lost_genes_1-22-19/scripts/merge.pl $wkdir/scripts/merge.pl 

mkdir $wkdir/39_gene_density/07_merge
mkdir $wkdir/39_gene_density/07_merge/10_kb_windows
mkdir $wkdir/39_gene_density/07_merge/100_kb_windows

cd $wkdir/39_gene_density/06_prep_merge/10_kb_windows/

perl $wkdir/scripts/merge.pl -k briggsae_gene_density briggsae_perc_n > $wkdir/39_gene_density/07_merge/10_kb_windows/briggsae 2>$wkdir/39_gene_density/07_merge/10_kb_windows/briggsae_error
perl $wkdir/scripts/merge.pl -k elegans_gene_density elegans_perc_n > $wkdir/39_gene_density/07_merge/10_kb_windows/elegans 2>$wkdir/39_gene_density/07_merge/10_kb_windows/elegans_error
perl $wkdir/scripts/merge.pl -k inopinata_gene_density inopinata_perc_n > $wkdir/39_gene_density/07_merge/10_kb_windows/inopinata 2>$wkdir/39_gene_density/07_merge/10_kb_windows/inopinata_error
perl $wkdir/scripts/merge.pl -k nigoni_gene_density nigoni_perc_n > $wkdir/39_gene_density/07_merge/10_kb_windows/nigoni 2>$wkdir/39_gene_density/07_merge/10_kb_windows/nigoni_error
perl $wkdir/scripts/merge.pl -k remanei_gene_density remanei_perc_n > $wkdir/39_gene_density/07_merge/10_kb_windows/remanei 2>$wkdir/39_gene_density/07_merge/10_kb_windows/remanei_error

cd $wkdir/39_gene_density/06_prep_merge/100_kb_windows/

perl $wkdir/scripts/merge.pl -k briggsae_gene_density briggsae_perc_n > $wkdir/39_gene_density/07_merge/100_kb_windows/briggsae 2>$wkdir/39_gene_density/07_merge/100_kb_windows/briggsae_error
perl $wkdir/scripts/merge.pl -k elegans_gene_density elegans_perc_n > $wkdir/39_gene_density/07_merge/100_kb_windows/elegans 2>$wkdir/39_gene_density/07_merge/100_kb_windows/elegans_error
perl $wkdir/scripts/merge.pl -k inopinata_gene_density inopinata_perc_n > $wkdir/39_gene_density/07_merge/100_kb_windows/inopinata 2>$wkdir/39_gene_density/07_merge/100_kb_windows/inopinata_error
perl $wkdir/scripts/merge.pl -k nigoni_gene_density nigoni_perc_n > $wkdir/39_gene_density/07_merge/100_kb_windows/nigoni 2>$wkdir/39_gene_density/07_merge/100_kb_windows/nigoni_error
perl $wkdir/scripts/merge.pl -k remanei_gene_density remanei_perc_n > $wkdir/39_gene_density/07_merge/100_kb_windows/remanei 2>$wkdir/39_gene_density/07_merge/100_kb_windows/remanei_error

cd $wkdir/39_gene_density/07_merge/100_kb_windows/

rm *_error

cd $wkdir/39_gene_density/07_merge/10_kb_windows/

rm *_error

####put it all together

mkdir $wkdir/39_gene_density/08_cat

cd $wkdir/39_gene_density/07_merge/10_kb_windows/

cat * | awk 'BEGIN {OFS="\t"} {print $1,$2,(($4/10000)*100),$3}' | sed -e 's/~/\t/' > $wkdir/39_gene_density/08_cat/all_gene_num_repeat_perc_N_10_kb_windows


cd $wkdir/39_gene_density/07_merge/100_kb_windows/

cat * | awk 'BEGIN {OFS="\t"} {print $1,$2,(($4/100000)*100),$3}' | sed -e 's/~/\t/' > $wkdir/39_gene_density/08_cat/all_gene_num_repeat_perc_N_100_kb_windows


echo -e "Chr\tBP\tgene_num\tperc_N\tspecies" | cat - all_gene_num_repeat_perc_N_10_kb_windows > all_gene_num_repeat_perc_N_10_kb_windows_tmp && mv all_gene_num_repeat_perc_N_10_kb_windows_tmp all_gene_num_repeat_perc_N_10_kb_windows


echo -e "Chr\tBP\tgene_num\tperc_N\tspecies" | cat - all_gene_num_repeat_perc_N_100_kb_windows > all_gene_num_repeat_perc_N_100_kb_windows_tmp && mv all_gene_num_repeat_perc_N_100_kb_windows_tmp all_gene_num_repeat_perc_N_100_kb_windows


###for the 100bp windows, do normalized distances from chromosome centers

mkdir $wkdir/39_gene_density/09_awk_norm_chr_pos
mkdir $wkdir/39_gene_density/10_cat
cd $wkdir/39_gene_density/09_awk_norm_chr_pos
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei


cd $wkdir/39_gene_density/09_awk_norm_chr_pos/briggsae

awk '{print > $1}' $wkdir/39_gene_density/05_awk/100_kb_windows/briggsae


cd $wkdir/39_gene_density/09_awk_norm_chr_pos/elegans

awk '{print > $1}' $wkdir/39_gene_density/05_awk/100_kb_windows/elegans

cd $wkdir/39_gene_density/09_awk_norm_chr_pos/inopinata

awk '{print > $1}' $wkdir/39_gene_density/05_awk/100_kb_windows/inopinata

cd $wkdir/39_gene_density/09_awk_norm_chr_pos/nigoni

awk '{print > $1}' $wkdir/39_gene_density/05_awk/100_kb_windows/nigoni

cd $wkdir/39_gene_density/09_awk_norm_chr_pos/remanei

awk '{print > $1}' $wkdir/39_gene_density/05_awk/100_kb_windows/remanei



cd $wkdir/39_gene_density/09_awk_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/10_cat/briggsae


cd $wkdir/39_gene_density/09_awk_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/10_cat/elegans

cd $wkdir/39_gene_density/09_awk_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/10_cat/inopinata



cd $wkdir/39_gene_density/09_awk_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/10_cat/nigoni


cd $wkdir/39_gene_density/09_awk_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/10_cat/remanei


cd $wkdir/39_gene_density/10_cat/

cat * > all

for i in *; do echo -e "Chr\tBP\tgene_num\tspecies\tnorm_dist_center" | cat - $i > $i.tmp && mv $i.tmp $i; done

for i in *; do sed -i -e 's/ /\t/g' $i; done

	#file "all" will be used later on to make data for gene density figures and analysis



#############
#############
#############
#############
#Exclude exclusively transposon-aligning genes from C. inopinata gene set and re-measure gene density
#############
#############
#############
#############


#align the C. inopinata transposon-aligning genes to the set of C. elegans proteins that DO NOT align to transposons

#get elegans non-transposon-aligning genes
mkdir $wkdir/38_blast_fused_proteins

#get all the elegans protein ids
grep ">" $wkdir/30_blastp/00_links/elegans.fa | sed -e 's/>//g' | sort | uniq > $wkdir/38_blast_fused_proteins/00_tot_elegans_prot_id

#get unique elegans protein ids that align to transposons
sort $wkdir/30_blastp/03_cut_uniq/elegans | uniq > $wkdir/38_blast_fused_proteins/01_elegans_align_transposon_prot_id

#get the elegans proteins that DO NOT align to transposons
comm -23 $wkdir/38_blast_fused_proteins/00_tot_elegans_prot_id $wkdir/38_blast_fused_proteins/01_elegans_align_transposon_prot_id > $wkdir/38_blast_fused_proteins/02_elegans_non-transposon_prot_id

#extract just the those sequences
cd  $wkdir/38_blast_fused_proteins

perl $wkdir/scripts/fasta_filter.pl 02_elegans_non-transposon_prot_id $wkdir/30_blastp/00_links/elegans.fa 03_elegans_non-transposon_prot.fa &

#make a blast db with the non-transposon elegans proteins
makeblastdb -in 03_elegans_non-transposon_prot.fa -parse_seqids -out elegans_non-transposon_prot_db -title elegans_non-transposon_prot_db -dbtype prot

#extract the inopinata transposon-aligning proteins
perl $wkdir/scripts/fasta_filter.pl $wkdir/30_blastp/03_cut_uniq/sp34 $wkdir/30_blastp/00_links/sp34.fa 04_sp34_transposon_prot.fa &

#do the blast...

blastp -db elegans_non-transposon_prot_db -outfmt 6 -query 04_sp34_transposon_prot.fa -out 05_sp34_transposon_prot_to_elegans_blast_out -evalue 0.001 &

#get the ids of inopinata-transposon-aligning genes that align to elegans-non-transposons

cut -f1 05_sp34_transposon_prot_to_elegans_blast_out | sort | uniq > 06_34_transposon_prot_align_to_elegans_list

#how many?

wc -l 06_34_transposon_prot_align_to_elegans_list


#this is the list of transposon-aligning proteins that also align to elegans non-transposons
	# $wkdir/38_blast_fused_proteins/06_34_transposon_prot_align_to_elegans_list
	#860

#let's get just the transposons that do not align to elegans proteins; this will be used to exclude the exclusively-transposon-aligning genes from the C. inopinata gene set

cd $wkdir/39_gene_density

mkdir $wkdir/39_gene_density/11_inopinata_remove_transposon_cds

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds

comm -23 <(grep -Po '\S+' $wkdir/30_blastp/03_cut_uniq/sp34 | sort) <(grep -Po '\S+' $wkdir/38_blast_fused_proteins/06_34_transposon_prot_align_to_elegans_list | sort)  > 00_inopinata_transposon_protein-coding_genes

#2489 genes ; 3349-860 = 2489 great
	#there are 2489 C. inopinata genes that align to transposons but DO NOT align to any elegans non-transposon-aligning protein-coding genes!! I want to exclude these transposons from the gene set!

#add "ID=" to the front for grep

sed -e 's/^/ID=/g' 00_inopinata_transposon_protein-coding_genes > 01_inopinata_transposon_protein-coding_genes_ID

#remove ".t1"

sed -i -e 's/\.t1//g' 01_inopinata_transposon_protein-coding_genes_ID

#remove these genes from the gff

LC_ALL=C fgrep -w -v -f 01_inopinata_transposon_protein-coding_genes_ID $wkdir/39_gene_density/01_prep_gff/inopinata.gff > 02_inopinata_no_transposons.gff &
	#19120 lines
	#for troubleshooting
ln -s $wkdir/39_gene_density/01_prep_gff/inopinata.gff $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/inopinata.gff
	# 21609 lines ; 21609-19120 = 2489 ; this works

rm inopinata.gff

#ok, proceed with the pipeline

#make bedfile for bedtools, add one for summing for gene density

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,"1"}'  02_inopinata_no_transposons.gff > 03_inopinata_no_transposons.bed

#sort bed

bedtools sort -i 03_inopinata_no_transposons.bed > 04_inopinata_no_transposons_sort.bed



#bedtools map ; need to rename chromosomes first

sed -i -e 's/CSP34.Sp34_Chr1/I/g' 04_inopinata_no_transposons_sort.bed
sed -i -e 's/CSP34.Sp34_Chr2/II/g' 04_inopinata_no_transposons_sort.bed
sed -i -e 's/CSP34.Sp34_Chr3/III/g' 04_inopinata_no_transposons_sort.bed
sed -i -e 's/CSP34.Sp34_Chr4/IV/g' 04_inopinata_no_transposons_sort.bed
sed -i -e 's/CSP34.Sp34_Chr5/V/g' 04_inopinata_no_transposons_sort.bed
sed -i -e 's/CSP34.Sp34_ChrX/X/g' 04_inopinata_no_transposons_sort.bed


bedtools map -b 04_inopinata_no_transposons_sort.bed -a $wkdir/windows/inopinata.10kb.windows -c 4 -o sum > 05_inopinata_no_transposons_bedtools_map_10_kb_windows

#replace "." with 0

sed -i -e 's/\./0/g' 05_inopinata_no_transposons_bedtools_map_10_kb_windows

#100kb windows
bedtools map -b 04_inopinata_no_transposons_sort.bed -a $wkdir/windows/inopinata.100kb.windows -c 4 -o sum > 06_inopinata_no_transposons_bedtools_map_100_kb_windows

sed -i -e 's/\./0/g' 06_inopinata_no_transposons_bedtools_map_100_kb_windows

#want to merge that 100_kb_windows set to the rest...

#$wkdir/39_gene_density/06_prep_merge/100_kb_windows/inopinata_perc_n

awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,"inopinata_no_transposon_cds"}' 06_inopinata_no_transposons_bedtools_map_100_kb_windows > 07_inopinata_no_transposons_bedtools_map_100_kb_windows_awk

#prep merge by replacing first tab with ~

sed -e 's/\t/~/' 07_inopinata_no_transposons_bedtools_map_100_kb_windows_awk > 08_inopinata_no_transposons_bedtools_map_100_kb_windows_prep_merge

#merge with repeat count data
perl $wkdir/scripts/merge.pl -k 08_inopinata_no_transposons_bedtools_map_100_kb_windows_prep_merge $wkdir/39_gene_density/06_prep_merge/100_kb_windows/inopinata_perc_n > 09_inopinata_gene_density_perc_n_100kb_windows_merge 2>$wkdir/39_gene_density/07_merge/100_kb_windows/inopinata_no_transposon_error_3-28-19 &

#get data right (percent repetitive, remove spurious "inopinata" column)
awk 'BEGIN {OFS="\t"} {print $1,$2,(($4/100000)*100),$3}' 09_inopinata_gene_density_perc_n_100kb_windows_merge | sed -e 's/~/\t/' > 10_inopinata_gene_density_perc_n_100kb_windows


#combine with previous results
cat $wkdir/39_gene_density/08_cat/all_gene_num_repeat_perc_N_100_kb_windows 10_inopinata_gene_density_perc_n_100kb_windows > 11_all_gene_num_repeat_perc_N_100_kb_windows
	
#now, get positions normalized by chromosome center


mkdir $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species

awk '{print > $5}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/11_all_gene_num_repeat_perc_N_100_kb_windows

rm species

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds

mkdir 13_awk_norm_chr_pos
cd 13_awk_norm_chr_pos
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir inopinata_no_transposon_cds
mkdir nigoni
mkdir remanei


mkdir 13_awk_norm_chr_pos
cd 13_awk_norm_chr_pos
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/briggsae

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/briggsae


cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/elegans

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/elegans

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/inopinata

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/inopinata

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/nigoni

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/nigoni

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/remanei

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/remanei


cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/inopinata_no_transposon_cds

awk '{print > $1}' $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/12_awk_species/inopinata_no_transposon_cds


mkdir $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/briggsae

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7727990-$2)/7727990)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8313577-$2)/8313577)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7289426-$2)/7289426)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8742720-$2)/8742720)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9747579-$2)/9747579)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10770285-$2)/10770285)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/briggsae


cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/elegans

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7536217-$2)/7536217)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7639711-$2)/7639711)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(6891901-$2)/6891901)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8746915-$2)/8746915)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10462090-$2)/10462090)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8859471-$2)/8859471)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/elegans

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/inopinata

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/inopinata



cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/nigoni

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8370129-$2)/8370129)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9614758-$2)/9614758)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(7767739-$2)/7767739)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10195166-$2)/10195166)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11143691-$2)/11143691)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11824229-$2)/11824229)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/nigoni


cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/remanei

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8623773-$2)/8623773)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9967862-$2)/9967862)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(8938925-$2)/8938925)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(12895499-$2)/12895499)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11251229-$2)/11251229)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10750950-$2)/10750950)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/remanei

cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/13_awk_norm_chr_pos/inopinata_no_transposon_cds

awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10297276-$2)/10297276)/2}' I > I.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10058498-$2)/10058498)/2}' II > II.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9718237-$2)/9718237)/2}' III > III.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(10508822-$2)/10508822)/2}' IV > IV.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(11819078-$2)/11819078)/2}' V > V.tmp
awk 'function abs(v) {return v < 0 ? -v : v} {print $0, (abs(9095254-$2)/9095254)/2}' X > X.tmp

cat *.tmp > $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/inopinata_no_transposon_cds



cd $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/14_cat/

cat * > gene_density.tsv
	#formerly "all_gene_dens_norm_chr_pos.tsv"

echo -e "Chr\tBP\tgene_num\tperc_N\tspecies\tnorm_dist_center" | cat - gene_density.tsv > gene_density.tsv.tmp && mv gene_density.tsv.tmp gene_density.tsv

sed -i -e 's/ /\t/g' gene_density.tsv
	#THIS is the data for Figures 6-7; Supplemental Figures 16-17


#############
#############
#############
#############
#Gene length
#############
#############
#############
#############


#ok, im curious about the genomic landscape of gene lengths.......

mkdir 56_gene_length_landscapes




mkdir $wkdir/56_gene_length_landscapes/00_links

#copy gene bed files over
cd $wkdir/39_gene_density/03_bedtools_sort

for i in *; do cp $i $wkdir/56_gene_length_landscapes/00_links_genes_bed/$i; done

#copy the inopinata no-transposon gene bed file over
cp $wkdir/39_gene_density/11_inopinata_remove_transposon_cds/04_inopinata_no_transposons_sort.bed $wkdir/56_gene_length_landscapes/00_links_genes_bed/inopinata_no_transposons.bed

#get gene lengths!
mkdir $wkdir/56_gene_length_landscapes/01_awk

cd $wkdir/56_gene_length_landscapes/00_links

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2}' $i > $wkdir/56_gene_length_landscapes/01_awk/$i; done

#get the genomic landscape of gene lengths with bedtools

mkdir $wkdir/56_gene_length_landscapes/02_bedtools_map/

bedtools map -o mean -c 4 -a $wkdir/windows/briggsae.100kb.windows -b briggsae > $wkdir/56_gene_length_landscapes/02_bedtools_map/briggsae &
bedtools map -o mean -c 4 -a $wkdir/windows/elegans.100kb.windows -b elegans > $wkdir/56_gene_length_landscapes/02_bedtools_map/elegans &
bedtools map -o mean -c 4 -a $wkdir/windows/nigoni.100kb.windows -b nigoni > $wkdir/56_gene_length_landscapes/02_bedtools_map/nigoni &
bedtools map -o mean -c 4 -a $wkdir/windows/inopinata.100kb.windows -b inopinata > $wkdir/56_gene_length_landscapes/02_bedtools_map/inopinata &
bedtools map -o mean -c 4 -a $wkdir/windows/inopinata.100kb.windows -b inopinata_no_transposons > $wkdir/56_gene_length_landscapes/02_bedtools_map/inopinata_no_transposons &
bedtools map -o mean -c 4 -a $wkdir/windows/remanei.100kb.windows -b remanei > $wkdir/56_gene_length_landscapes/02_bedtools_map/remanei &

#add species name

mkdir $wkdir/56_gene_length_landscapes/03_awk

cd $wkdir/56_gene_length_landscapes/02_bedtools_map/

for i in *; do awk 'BEGIN {OFS="\t"} {print $1,$2+1,$4,FILENAME}' $i > $wkdir/56_gene_length_landscapes/03_awk/$i; done

cd $wkdir/56_gene_length_landscapes/03_awk/

cat * > all

echo -e "Chr\tBP\tgene_length\tspecies" | cat - all > all_gene_lengths_8-13-19.tsv

# added the normalized chromosome positions on local machine with previously-existing files
	#this is the data for Supplemental Figure 21-22

#############
#############
#############
#############
#Numbers of insertions per cluster and their size/length (bp) distributions
#############
#############
#############
#############

#count and size distribution of each insertion for each cluster in each assembly




mkdir $wkdir/54_cluster_distributions

#link to gff data
cd $wkdir/54_cluster_distributions
mkdir 00_links

ln -s $wkdir//29_prep_plots/repeatmasker_ii/07_paste/briggsae $wkdir/54_cluster_distributions/00_links/briggsae
ln -s $wkdir//29_prep_plots/repeatmasker_ii/07_paste/elegans $wkdir/54_cluster_distributions/00_links/elegans
ln -s $wkdir//29_prep_plots/repeatmasker_ii/07_paste/inopinata $wkdir/54_cluster_distributions/00_links/inopinata
ln -s $wkdir//29_prep_plots/repeatmasker_ii/07_paste/nigoni $wkdir/54_cluster_distributions/00_links/nigoni
ln -s $wkdir//29_prep_plots/repeatmasker_ii/07_paste/remanei $wkdir/54_cluster_distributions/00_links/remanei

mkdir $wkdir/54_cluster_distributions/01_awk
cd $wkdir/54_cluster_distributions/00_links

#get data right, add 1 for counting
for i in *; do awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $1,$4,$5,$10,$11,($5-$4),FILENAME,"1"}' $i > $wkdir/54_cluster_distributions/01_awk/$i; done &

#get count distributions...
# just get the list of clusters
cd ..
mkdir $wkdir/54_cluster_distributions/02_awk
cd $wkdir/54_cluster_distributions/01_awk

for i in *; do awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $4}' $i > $wkdir/54_cluster_distributions/02_awk/$i; done &

cd $wkdir/54_cluster_distributions/02_awk

for i in *; do echo -e "CLUSTER" | cat - $i > $i.tmp && mv $i.tmp $i; done &


#use count_general.py to count the number of insertions per cluster...

mkdir $wkdir/54_cluster_distributions/03_count_general_clusters

cd $wkdir/54_cluster_distributions/02_awk

for i in *; do python $wkdir/scripts/count_general.py $i $wkdir/03_count_general_clusters/$i; done &

#replace header

cd $wkdir/54_cluster_distributions/03_count_general_clusters

for i in *; do sed -i '1d' $i; done &

for i in *; do echo -e "cluster_id\tcluster_count" | cat - $i > $i.tmp && mv $i.tmp $i; done &

	#cool, the files in $wkdir/54_cluster_distributions/03_count_general_clusters have the count distributions for every cluster in all genomes

#now, what about size/length distributions of all insertions in all genomes?
	#this data is bascally already in $wkdir/54_cluster_distributions/01_awk

mkdir $wkdir/54_cluster_distributions/04_size_distributions

cd $wkdir/54_cluster_distributions/04_size_distributions

mkdir 00_grep_cluster

cd $wkdir/54_cluster_distributions/01_awk

for i in *; do grep "Cluster" $i > $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/$i; done &


mkdir $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster

cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#get positions of each specific cluster
cd briggsae
awk '{print > $4}' $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/briggsae

cd ..
cd elegans
awk '{print > $4}' $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/elegans

cd ..
cd inopinata
awk '{print > $4}' $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/inopinata


cd ..
cd nigoni
awk '{print > $4}' $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/nigoni


cd ..
cd remanei
awk '{print > $4}' $wkdir/54_cluster_distributions/04_size_distributions/00_grep_cluster/remanei

#

mkdir 02_make_histogram
cd 02_make_histogram
mkdir briggsae
mkdir elegans
mkdir inopinata
mkdir nigoni
mkdir remanei

#ok, make the plots

cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster/briggsae

for i in *; do Rscript $wkdir/scripts/make_histogram.R $i; done &
	#make sure to change output directory in R script before running!
	

cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster/elegans

for i in *; do Rscript $wkdir/scripts/make_histogram.R $i; done &
	#make sure to change output directory in R script before running!


cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster/inopinata

for i in *; do Rscript $wkdir/scripts/make_histogram.R $i; done &
	#make sure to change output directory in R script before running!



cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster/nigoni

for i in *; do Rscript $wkdir/scripts/make_histogram.R $i; done &
	#make sure to change output directory in R script before running!



cd $wkdir/54_cluster_distributions/04_size_distributions/01_awk_cluster/remanei

for i in *; do Rscript $wkdir/scripts/make_histogram.R $i; done &
	#make sure to change output directory in R script before running!

#############
#############
#############
#############
#C. inopinata insertion ages
#############
#############
#############
#############


#make the folder

mkdir $wkdir/55_inopinata_insertion_ages

cd  $wkdir/55_inopinata_insertion_ages


#get cluster id # repeatclassifier info

mkdir $wkdir/55_inopinata_insertion_ages/00_grep_cluster_ids

cd $wkdir/55_inopinata_insertion_ages/00_grep_cluster_ids

# this is where the repeat libraries used for the RepeatMasker_II run...

# $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa

grep ">" $wkdir/24_fasta_filter/inopinata_combo_repeats_classified_no_prot.fa | sed -e 's/>//g' > $wkdir/55_inopinata_insertion_ages/00_grep_cluster_ids/inopinata_cluster_ids_with_classification
	#get the ids with classifications

sed -e 's/#/\t/g' inopinata_cluster_ids_with_classification | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' > inopinata_cluster_ids
	#just the ids


#from gff with cluster id info, awk by cluster id to get cluster-id specific gff's

#NOTE, I only want to look at repeats with >100 insertions and not sattelite/simple repeat/low complexity

mkdir $wkdir/55_inopinata_insertion_ages/01_grep_mask_gff

cd $wkdir/55_inopinata_insertion_ages/01_grep_mask_gff

	#ack, forgot strandedness, need to do over!
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,$7,$9}' $wkdir/27_RepeatMasker_II/hard/inopinata_genome.fa.out.gff | grep -v "#" > inopinata_insertions.bed

#clean up bed file

sed -i -e 's/Target //g' inopinata_insertions.bed
sed -i -e 's/ .*//g' inopinata_insertions.bed
sed -i -e 's/"Motif://g' inopinata_insertions.bed
sed -i -e 's/"//g' inopinata_insertions.bed

#exclude simple or low complexity

grep 'Cluster' inopinata_insertions.bed > inopinata_insertions_no_simple_or_low_complexity.bed


#ok, now find clusters with <100 insertions, this data I have. in folder $wkdir/54_cluster_distributions/03_count_general_clusters

grep 'Cluster' $wkdir/54_cluster_distributions/03_count_general_clusters/inopinata | awk ' $2 >= 100 ' > clusters_grtr_or_equal_100_insertions

#looks good

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' clusters_grtr_or_equal_100_insertions > clusters_grtr_or_equal_100_insertions_only_ids

#get just the genomic regions from these clusters


LC_ALL=C fgrep -w -f clusters_grtr_or_equal_100_insertions_only_ids inopinata_insertions_no_simple_or_low_complexity.bed > clusters_grtr_or_equal_100_insertions.bed


#ok, gut check that this worked

wc -l clusters_grtr_or_equal_100_insertions
	#265 clusters_grtr_or_equal_100_insertions

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $5}' clusters_grtr_or_equal_100_insertions.bed | sort | uniq | wc -l
	# 265
	#woo!

#ok, now awk by cluster!!

cd ..
mkdir $wkdir/55_inopinata_insertion_ages/02_awk_bed_by_cluster
cd $wkdir/55_inopinata_insertion_ages/02_awk_bed_by_cluster

awk '{print > $5}' $wkdir/55_inopinata_insertion_ages/01_grep_mask_gff/clusters_grtr_or_equal_100_insertions.bed

#cool, add unique insertion ids...

cd ..
mkdir $wkdir/55_inopinata_insertion_ages/03_nl_inerstion_ids
cd $wkdir/55_inopinata_insertion_ages/02_awk_bed_by_cluster
for i in *; do nl -n ln $i > $wkdir/55_inopinata_insertion_ages/03_nl_inerstion_ids/$i; done

#ok, rearrange

mkdir $wkdir/55_inopinata_insertion_ages/04_awk_move_id_numbers
cd $wkdir/55_inopinata_insertion_ages/03_nl_inerstion_ids

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$1}' $i > $wkdir/55_inopinata_insertion_ages/04_awk_move_id_numbers/$i; done

#replace last tab with underscore

cd $wkdir/55_inopinata_insertion_ages/04_awk_move_id_numbers

for i in *; do sed -i -e 's/\t\([^\t]*\)$/_\1/' $i; done

#bed format

cd ..
mkdir $wkdir/55_inopinata_insertion_ages/05_awk_bed_format
cd $wkdir/55_inopinata_insertion_ages/04_awk_move_id_numbers

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,$5,"score",$4}' $i > $wkdir/55_inopinata_insertion_ages/05_awk_bed_format/$i; done


#ok, now just try to get the insertion sequences sequences with bedfasta

mkdir $wkdir/55_inopinata_insertion_ages/06_bedtools_getfasta_insertion_sequences
cd $wkdir/55_inopinata_insertion_ages/05_awk_bed_format


for i in *; do bedtools getfasta -s -name -fo $wkdir/55_inopinata_insertion_ages/06_bedtools_getfasta_insertion_sequences/$i -fi $wkdir/01_fasta_filter/inopinata_genome.fa -bed $i; done &


#that worked, cool.... now align

cd ..
mkdir $wkdir/55_inopinata_insertion_ages/07_mafft_alignments
cd $wkdir/55_inopinata_insertion_ages/06_bedtools_getfasta_insertion_sequences

#we used mafft version 7.313 for alignments
module load mafft/7.313

for i in *; do mafft --auto --adjustdirection $i > $wkdir/55_inopinata_insertion_ages/07_mafft_alignments/$i; done &

#make trees with untrimmed alignments because automated trimming excludes a bunch of clusters; we used FastTree version 2.1.10

mkdir $wkdir/55_inopinata_insertion_ages/10_fasttree_no_trim/

cd $wkdir/55_inopinata_insertion_ages/07_mafft_alignments/

for i in *; do FastTree -nt -gtr -quiet < $i > /projects/phillipslab/gavincw/repeats_12-18-18/55_inopinata_insertion_ages/10_fasttree_no_trim/$i.nwk &


#trim alignments with trimal option -automated1 on everything ; we used trimal version 1.2
mkdir $wkdir/55_inopinata_insertion_ages/11_trimal_automated1/
cd $wkdir/55_inopinata_insertion_ages/07_mafft_alignments/

for i in *; do trimal -in $i -out $wkdir/55_inopinata_insertion_ages/11_trimal_automated1/$i -automated1; done


#get branch lengths for untrimmed alignments

mkdir $wkdir/55_inopinata_insertion_ages/12_branch_lengths_no_trimal
mkdir $wkdir/55_inopinata_insertion_ages/12_branch_lengths_no_trimal_summary

cd $wkdir/55_inopinata_insertion_ages/10_fasttree_no_trim/

for i in *; do mv $i ${i%.nwk}; done

cd /Users/gavin/genome/repeats_continued_7-25-19/10_fasttree_no_trim/

for i in *; do Rscript $wkdir/scripts/branch_lengths.R $i; done &
	#make sure to change output directory in R script!

#so, just add repeat taxonomy to every thing
	#first, just to the summaries

mkdir $wkdir/55_inopinata_insertion_ages/17_no_trim_branch_length_summary_awk

cd $wkdir/55_inopinata_insertion_ages/12_branch_lengths_no_trimal_summary

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/55_inopinata_insertion_ages/17_no_trim_branch_length_summary_awk/$i; done

cd ..

cd $wkdir/55_inopinata_insertion_ages/17_no_trim_branch_length_summary_awk

sed -i -e 's/Cluster1006/Cluster1006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1006
sed -i -e 's/Cluster1007/Cluster1007\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1007
sed -i -e 's/Cluster100/Cluster100\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster100
sed -i -e 's/Cluster1010/Cluster1010\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1010
sed -i -e 's/Cluster1012/Cluster1012\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1012
sed -i -e 's/Cluster1013/Cluster1013\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1013
sed -i -e 's/Cluster1018/Cluster1018\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster1018
sed -i -e 's/Cluster101/Cluster101\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster101
sed -i -e 's/Cluster1026/Cluster1026\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1026
sed -i -e 's/Cluster107/Cluster107\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster107
sed -i -e 's/Cluster113/Cluster113\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster113
sed -i -e 's/Cluster115/Cluster115\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster115
sed -i -e 's/Cluster1173/Cluster1173\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1173
sed -i -e 's/Cluster1175/Cluster1175\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1175
sed -i -e 's/Cluster1179/Cluster1179\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster1179
sed -i -e 's/Cluster1180/Cluster1180\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1180
sed -i -e 's/Cluster118/Cluster118\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster118
sed -i -e 's/Cluster1210/Cluster1210\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1210
sed -i -e 's/Cluster1211/Cluster1211\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster1211
sed -i -e 's/Cluster1212/Cluster1212\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1212
sed -i -e 's/Cluster1213/Cluster1213\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' Cluster1213
sed -i -e 's/Cluster1226/Cluster1226\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' Cluster1226
sed -i -e 's/Cluster122/Cluster122\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster122
sed -i -e 's/Cluster1230/Cluster1230\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster1230
sed -i -e 's/Cluster1234/Cluster1234\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1234
sed -i -e 's/Cluster1240/Cluster1240\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1240
sed -i -e 's/Cluster1242/Cluster1242\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1242
sed -i -e 's/Cluster1244/Cluster1244\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/g' Cluster1244
sed -i -e 's/Cluster1249/Cluster1249\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1249
sed -i -e 's/Cluster1250/Cluster1250\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1250
sed -i -e 's/Cluster1254/Cluster1254\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1254
sed -i -e 's/Cluster1255/Cluster1255\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster1255
sed -i -e 's/Cluster1258/Cluster1258\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1258
sed -i -e 's/Cluster1378/Cluster1378\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1378
sed -i -e 's/Cluster1381/Cluster1381\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster1381
sed -i -e 's/Cluster1387/Cluster1387\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster1387
sed -i -e 's/Cluster1388/Cluster1388\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster1388
sed -i -e 's/Cluster1390/Cluster1390\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1390
sed -i -e 's/Cluster1414/Cluster1414\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1414
sed -i -e 's/Cluster1417/Cluster1417\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1417
sed -i -e 's/Cluster142/Cluster142\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster142
sed -i -e 's/Cluster1438/Cluster1438\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1438
sed -i -e 's/Cluster1439/Cluster1439\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1439
sed -i -e 's/Cluster1443/Cluster1443\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1443
sed -i -e 's/Cluster1447/Cluster1447\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1447
sed -i -e 's/Cluster1448/Cluster1448\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster1448
sed -i -e 's/Cluster1449/Cluster1449\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster1449
sed -i -e 's/Cluster1451/Cluster1451\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster1451
sed -i -e 's/Cluster1452/Cluster1452\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1452
sed -i -e 's/Cluster1455/Cluster1455\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1455
sed -i -e 's/Cluster1458/Cluster1458\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1458
sed -i -e 's/Cluster1460/Cluster1460\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1460
sed -i -e 's/Cluster1462/Cluster1462\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1462
sed -i -e 's/Cluster1468/Cluster1468\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1468
sed -i -e 's/Cluster1469/Cluster1469\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1469
sed -i -e 's/Cluster1480/Cluster1480\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1480
sed -i -e 's/Cluster150/Cluster150\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster150
sed -i -e 's/Cluster1541/Cluster1541\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster1541
sed -i -e 's/Cluster1591/Cluster1591\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1591
sed -i -e 's/Cluster1593/Cluster1593\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1593
sed -i -e 's/Cluster1594/Cluster1594\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster1594
sed -i -e 's/Cluster1620/Cluster1620\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1620
sed -i -e 's/Cluster1633/Cluster1633\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster1633
sed -i -e 's/Cluster1635/Cluster1635\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1635
sed -i -e 's/Cluster1642/Cluster1642\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1642
sed -i -e 's/Cluster1644/Cluster1644\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1644
sed -i -e 's/Cluster1647/Cluster1647\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1647
sed -i -e 's/Cluster1650/Cluster1650\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1650
sed -i -e 's/Cluster1683/Cluster1683\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1683
sed -i -e 's/Cluster1722/Cluster1722\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster1722
sed -i -e 's/Cluster1744/Cluster1744\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1744
sed -i -e 's/Cluster1750/Cluster1750\tDNA-Merlin\tII\t1\tTIR\tMerlin\tNA/g' Cluster1750
sed -i -e 's/Cluster1752/Cluster1752\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster1752
sed -i -e 's/Cluster1755/Cluster1755\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster1755
sed -i -e 's/Cluster1774/Cluster1774\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1774
sed -i -e 's/Cluster1776/Cluster1776\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster1776
sed -i -e 's/Cluster1800/Cluster1800\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster1800
sed -i -e 's/Cluster1802/Cluster1802\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1802
sed -i -e 's/Cluster1803/Cluster1803\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1803
sed -i -e 's/Cluster1809/Cluster1809\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1809
sed -i -e 's/Cluster1810/Cluster1810\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1810
sed -i -e 's/Cluster1814/Cluster1814\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster1814
sed -i -e 's/Cluster1816/Cluster1816\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1816
sed -i -e 's/Cluster1820/Cluster1820\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1820
sed -i -e 's/Cluster1846/Cluster1846\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1846
sed -i -e 's/Cluster1889/Cluster1889\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1889
sed -i -e 's/Cluster1940/Cluster1940\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1940
sed -i -e 's/Cluster1942/Cluster1942\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1942
sed -i -e 's/Cluster1944/Cluster1944\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1944
sed -i -e 's/Cluster1946/Cluster1946\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1946
sed -i -e 's/Cluster1958/Cluster1958\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1958
sed -i -e 's/Cluster1969/Cluster1969\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster1969
sed -i -e 's/Cluster1979/Cluster1979\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster1979
sed -i -e 's/Cluster1989/Cluster1989\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1989
sed -i -e 's/Cluster1991/Cluster1991\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster1991
sed -i -e 's/Cluster1994/Cluster1994\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster1994
sed -i -e 's/Cluster1997/Cluster1997\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster1997
sed -i -e 's/Cluster2067/Cluster2067\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2067
sed -i -e 's/Cluster2109/Cluster2109\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/g' Cluster2109
sed -i -e 's/Cluster210/Cluster210\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster210
sed -i -e 's/Cluster2121/Cluster2121\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2121
sed -i -e 's/Cluster2124/Cluster2124\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2124
sed -i -e 's/Cluster2125/Cluster2125\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2125
sed -i -e 's/Cluster2129/Cluster2129\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster2129
sed -i -e 's/Cluster2130/Cluster2130\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2130
sed -i -e 's/Cluster2142/Cluster2142\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2142
sed -i -e 's/Cluster2143/Cluster2143\tDNA-PiggyBac\tII\t1\tTIR\tPiggyBac\tNA/g' Cluster2143
sed -i -e 's/Cluster2149/Cluster2149\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2149
sed -i -e 's/Cluster2160/Cluster2160\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2160
sed -i -e 's/Cluster2161/Cluster2161\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2161
sed -i -e 's/Cluster2169/Cluster2169\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2169
sed -i -e 's/Cluster2176/Cluster2176\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2176
sed -i -e 's/Cluster2180/Cluster2180\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2180
sed -i -e 's/Cluster2182/Cluster2182\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster2182
sed -i -e 's/Cluster2187/Cluster2187\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2187
sed -i -e 's/Cluster2191/Cluster2191\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2191
sed -i -e 's/Cluster2192/Cluster2192\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2192
sed -i -e 's/Cluster2260/Cluster2260\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2260
sed -i -e 's/Cluster2267/Cluster2267\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2267
sed -i -e 's/Cluster2298/Cluster2298\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' Cluster2298
sed -i -e 's/Cluster2299/Cluster2299\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2299
sed -i -e 's/Cluster2310/Cluster2310\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2310
sed -i -e 's/Cluster2313/Cluster2313\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2313
sed -i -e 's/Cluster2314/Cluster2314\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2314
sed -i -e 's/Cluster2317/Cluster2317\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster2317
sed -i -e 's/Cluster2342/Cluster2342\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2342
sed -i -e 's/Cluster2344/Cluster2344\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster2344
sed -i -e 's/Cluster2347/Cluster2347\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster2347
sed -i -e 's/Cluster2351/Cluster2351\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2351
sed -i -e 's/Cluster2356/Cluster2356\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2356
sed -i -e 's/Cluster2363/Cluster2363\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster2363
sed -i -e 's/Cluster2389/Cluster2389\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2389
sed -i -e 's/Cluster2476/Cluster2476\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster2476
sed -i -e 's/Cluster2477/Cluster2477\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2477
sed -i -e 's/Cluster2478/Cluster2478\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2478
sed -i -e 's/Cluster2481/Cluster2481\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2481
sed -i -e 's/Cluster2497/Cluster2497\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' Cluster2497
sed -i -e 's/Cluster2498/Cluster2498\tSatellite\tSatellite\tNA\tNA\tNA\tNA\t/g' Cluster2498
sed -i -e 's/Cluster2505/Cluster2505\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster2505
sed -i -e 's/Cluster2518/Cluster2518\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster2518
sed -i -e 's/Cluster2520/Cluster2520\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster2520
sed -i -e 's/Cluster2524/Cluster2524\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster2524
sed -i -e 's/Cluster2527/Cluster2527\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/g' Cluster2527
sed -i -e 's/Cluster2529/Cluster2529\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster2529
sed -i -e 's/Cluster2559/Cluster2559\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2559
sed -i -e 's/Cluster2620/Cluster2620\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2620
sed -i -e 's/Cluster2626/Cluster2626\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2626
sed -i -e 's/Cluster2627/Cluster2627\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/g' Cluster2627
sed -i -e 's/Cluster2629/Cluster2629\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2629
sed -i -e 's/Cluster2643/Cluster2643\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2643
sed -i -e 's/Cluster2656/Cluster2656\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2656
sed -i -e 's/Cluster2665/Cluster2665\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2665
sed -i -e 's/Cluster2667/Cluster2667\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2667
sed -i -e 's/Cluster2671/Cluster2671\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster2671
sed -i -e 's/Cluster2672/Cluster2672\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2672
sed -i -e 's/Cluster2676/Cluster2676\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster2676
sed -i -e 's/Cluster2679/Cluster2679\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2679
sed -i -e 's/Cluster2680/Cluster2680\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2680
sed -i -e 's/Cluster2704/Cluster2704\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2704
sed -i -e 's/Cluster2722/Cluster2722\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster2722
sed -i -e 's/Cluster2786/Cluster2786\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster2786
sed -i -e 's/Cluster2808/Cluster2808\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2808
sed -i -e 's/Cluster2809/Cluster2809\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster2809
sed -i -e 's/Cluster2818/Cluster2818\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/g' Cluster2818
sed -i -e 's/Cluster2829/Cluster2829\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2829
sed -i -e 's/Cluster2830/Cluster2830\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster2830
sed -i -e 's/Cluster2834/Cluster2834\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2834
sed -i -e 's/Cluster2835/Cluster2835\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2835
sed -i -e 's/Cluster2837/Cluster2837\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2837
sed -i -e 's/Cluster2838/Cluster2838\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2838
sed -i -e 's/Cluster2840/Cluster2840\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster2840
sed -i -e 's/Cluster2841/Cluster2841\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster2841
sed -i -e 's/Cluster2843/Cluster2843\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2843
sed -i -e 's/Cluster2846/Cluster2846\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster2846
sed -i -e 's/Cluster2847/Cluster2847\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2847
sed -i -e 's/Cluster2851/Cluster2851\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2851
sed -i -e 's/Cluster2852/Cluster2852\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster2852
sed -i -e 's/Cluster2857/Cluster2857\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2857
sed -i -e 's/Cluster2866/Cluster2866\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2866
sed -i -e 's/Cluster2933/Cluster2933\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster2933
sed -i -e 's/Cluster2962/Cluster2962\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/g' Cluster2962
sed -i -e 's/Cluster2970/Cluster2970\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster2970
sed -i -e 's/Cluster2982/Cluster2982\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2982
sed -i -e 's/Cluster2985/Cluster2985\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2985
sed -i -e 's/Cluster2986/Cluster2986\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2986
sed -i -e 's/Cluster2995/Cluster2995\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster2995
sed -i -e 's/Cluster2996/Cluster2996\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster2996
sed -i -e 's/Cluster2997/Cluster2997\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2997
sed -i -e 's/Cluster2999/Cluster2999\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster2999
sed -i -e 's/Cluster3002/Cluster3002\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster3002
sed -i -e 's/Cluster3003/Cluster3003\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster3003
sed -i -e 's/Cluster3005/Cluster3005\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/g' Cluster3005
sed -i -e 's/Cluster3006/Cluster3006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster3006
sed -i -e 's/Cluster3009/Cluster3009\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster3009
sed -i -e 's/Cluster3029/Cluster3029\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster3029
sed -i -e 's/Cluster3113/Cluster3113\tDNA-MuLE-MuDR\tII\t1\tTIR\tMutator\tMuDR/g' Cluster3113
sed -i -e 's/Cluster3114/Cluster3114\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster3114
sed -i -e 's/Cluster3115/Cluster3115\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster3115
sed -i -e 's/Cluster3116/Cluster3116\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster3116
sed -i -e 's/Cluster3132/Cluster3132\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster3132
sed -i -e 's/Cluster3133/Cluster3133\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster3133
sed -i -e 's/Cluster3134/Cluster3134\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster3134
sed -i -e 's/Cluster3141/Cluster3141\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster3141
sed -i -e 's/Cluster3149/Cluster3149\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/g' Cluster3149
sed -i -e 's/Cluster3164/Cluster3164\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster3164
sed -i -e 's/Cluster3165/Cluster3165\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster3165
sed -i -e 's/Cluster317/Cluster317\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster317
sed -i -e 's/Cluster3217/Cluster3217\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster3217
sed -i -e 's/Cluster37/Cluster37\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster37
sed -i -e 's/Cluster479/Cluster479\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster479
sed -i -e 's/Cluster480/Cluster480\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster480
sed -i -e 's/Cluster481/Cluster481\tDNA-Sola-3\tII\t1\tTIR\tSola\tSola3/g' Cluster481
sed -i -e 's/Cluster483/Cluster483\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster483
sed -i -e 's/Cluster485/Cluster485\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster485
sed -i -e 's/Cluster488/Cluster488\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster488
sed -i -e 's/Cluster489/Cluster489\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/g' Cluster489
sed -i -e 's/Cluster491/Cluster491\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster491
sed -i -e 's/Cluster511/Cluster511\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster511
sed -i -e 's/Cluster514/Cluster514\tDNA-PIF-Harbinger\tII\t1\tTIR\tPIF-Harbinger\tNA/g' Cluster514
sed -i -e 's/Cluster524/Cluster524\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster524
sed -i -e 's/Cluster525/Cluster525\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster525
sed -i -e 's/Cluster527/Cluster527\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster527
sed -i -e 's/Cluster528/Cluster528\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster528
sed -i -e 's/Cluster540/Cluster540\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster540
sed -i -e 's/Cluster542/Cluster542\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster542
sed -i -e 's/Cluster543/Cluster543\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster543
sed -i -e 's/Cluster54/Cluster54\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster54
sed -i -e 's/Cluster600/Cluster600\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster600
sed -i -e 's/Cluster601/Cluster601\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster601
sed -i -e 's/Cluster69/Cluster69\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster69
sed -i -e 's/Cluster717/Cluster717\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster717
sed -i -e 's/Cluster720/Cluster720\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/g' Cluster720
sed -i -e 's/Cluster722/Cluster722\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/g' Cluster722
sed -i -e 's/Cluster724/Cluster724\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster724
sed -i -e 's/Cluster725/Cluster725\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster725
sed -i -e 's/Cluster726/Cluster726\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/g' Cluster726
sed -i -e 's/Cluster727/Cluster727\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/g' Cluster727
sed -i -e 's/Cluster729/Cluster729\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster729
sed -i -e 's/Cluster731/Cluster731\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster731
sed -i -e 's/Cluster738/Cluster738\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster738
sed -i -e 's/Cluster757/Cluster757\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster757
sed -i -e 's/Cluster758/Cluster758\tDNA\tII\tNA\tNA\tNA\tNA\t/g' Cluster758
sed -i -e 's/Cluster774/Cluster774\tUnknown\tNA\tNA\tNA\tNA\tNA/g' Cluster774
sed -i -e 's/Cluster789/Cluster789\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/g' Cluster789
sed -i -e 's/Cluster793/Cluster793\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster793
sed -i -e 's/Cluster794/Cluster794\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster794
sed -i -e 's/Cluster795/Cluster795\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster795
sed -i -e 's/Cluster798/Cluster798\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster798
sed -i -e 's/Cluster800/Cluster800\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster800
sed -i -e 's/Cluster805/Cluster805\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster805
sed -i -e 's/Cluster806/Cluster806\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster806
sed -i -e 's/Cluster80/Cluster80\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster80
sed -i -e 's/Cluster811/Cluster811\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster811
sed -i -e 's/Cluster812/Cluster812\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster812
sed -i -e 's/Cluster813/Cluster813\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster813
sed -i -e 's/Cluster816/Cluster816\tSatellite\tSatellite\tNA\tNA\tNA\tNA\t/g' Cluster816
sed -i -e 's/Cluster856/Cluster856\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster856
sed -i -e 's/Cluster929/Cluster929\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster929
sed -i -e 's/Cluster938/Cluster938\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster938
sed -i -e 's/Cluster961/Cluster961\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/g' Cluster961
sed -i -e 's/Cluster967/Cluster967\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/g' Cluster967
sed -i -e 's/Cluster975/Cluster975\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/g' Cluster975
sed -i -e 's/Cluster992/Cluster992\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster992
sed -i -e 's/Cluster99/Cluster99\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/g' Cluster99
sed -i -e 's/Cluster9/Cluster9\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/g' Cluster9


#put the above in sed.sh and ran it

rm sed.sh

mkdir $wkdir/55_inopinata_insertion_ages/18_no_trim_branch_length_summary_cat
#combine everything
cat * > $wkdir/55_inopinata_insertion_ages/18_no_trim_branch_length_summary_cat/branch_length_no_trim_summary
#add header

cd $wkdir/55_inopinata_insertion_ages/18_no_trim_branch_length_summary_cat/

echo -e "cluster_id\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_superfamily\tmean\tmedian\tmin\tmax\tsd\tiqr\tcount" | cat - branch_length_no_trim_summary > branch_length_no_trim_summary_tmp && mv branch_length_no_trim_summary_tmp branch_lengths_summaries
	#THIS is the data for Supplemental Figures 12-13!!


#ok... now, insertion age by chr position.
#first, as mafft changed insertion ids, fix that.

cd $wkdir/55_inopinata_insertion_ages/12_branch_lengths_no_trimal


for i in *; do sed -i -e 's/_R_//g' $i; done &


#alright, cool.

#ok, first cat the bedfiles for individual clusters/insertions 

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos
mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/00_cat_bedfiles/

cd $wkdir/55_inopinata_insertion_ages/05_awk_bed_format

cat * > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/00_cat_bedfiles/all_insertions_bed

#make sure the appropriate key (cluster id) is first.

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/00_cat_bedfiles/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4,$0}' all_insertions_bed | sort > id_key_all_insertions_bed

#ok, do the same with the branch lengths

cd $wkdir/55_inopinata_insertion_ages/12_branch_lengths_no_trimal

# remove header

for i in *; do sed -i '1d' $i; done &

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/01_cat_branch_lengths/

cat * | sort > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/01_cat_branch_lengths/all_branch_lengths &


#alright, moment of truth. can we bring these files together....

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/02_merge


$wkdir/merge.pl


perl $wkdir/scripts/merge.pl -k -e "no_key" $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/00_cat_bedfiles/id_key_all_insertions_bed $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/01_cat_branch_lengths/all_branch_lengths 2> $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/02_merge/inopinata_insertion_locations_ages.error > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/02_merge/inopinata_insertion_locations_ages &

#wooo that shit worked!

#ok, now let's run the rest of that pipeline.

#get into bed format.

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/03_awk_bed

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/02_merge

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$8}' inopinata_insertion_locations_ages > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/03_awk_bed/inopinata_insertion_ages_bed

#alright, let's sort

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/04_bedtools_sort

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/03_awk_bed

bedtools sort -i inopinata_insertion_ages_bed > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/04_bedtools_sort/inopinata_insertion_ages_bed


#ok, bedtools...

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/05_bedtools_map

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/04_bedtools_sort


bedtools map -o mean -c 4 -a $wkdir/windows/inopinata.10kb.windows -b inopinata_insertion_ages_bed > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/05_bedtools_map/inopinata_insertion_ages_mean_10kb_windows_bed

#remove missing windows

mkdir $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/06_awk

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/05_bedtools_map

awk 'BEGIN {OFS="\t"} $4 != "."' inopinata_insertion_ages_mean_10kb_windows_bed > $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/06_awk/inopinata_insertion_ages_mean_10kb_windows_bed

#remove column three add header...

cd $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/06_awk/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4}' inopinata_insertion_ages_mean_10kb_windows_bed > inopinata_insertion_ages_mean_10kb_windows_bed.tmp

echo -e "Chr\tBP\tinsertion_age" | cat - inopinata_insertion_ages_mean_10kb_windows_bed.tmp > inopinata_insertion_ages_mean_10kb_windows_bed

#get landscape of all insertions using untrimmed alignments



mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/

mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/00_links

cd  $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/00_links

#copy over insertion locations, ages, and ids
cp $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/02_merge/inopinata_insertion_locations_ages inopinata_insertion_locations_ages

mkdir 01_sed

cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/00_links/
#get rid of first underscore to get just cluster ids for joining

sed -e 's/_/\t/1' inopinata_insertion_locations_ages > $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/01_sed/inopinata_insertion_locations_ages

#awk by cluster for adding repeat taxonomy

mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk/

cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk/

awk '{print > $1}' $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/01_sed/inopinata_insertion_locations_ages

#add repeat taxonomy


cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk/

sed -i -e 's/Cluster1006/Cluster1006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1006
sed -i -e 's/Cluster1007/Cluster1007\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1007
sed -i -e 's/Cluster100/Cluster100\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster100
sed -i -e 's/Cluster1010/Cluster1010\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1010
sed -i -e 's/Cluster1012/Cluster1012\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1012
sed -i -e 's/Cluster1013/Cluster1013\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1013
sed -i -e 's/Cluster1018/Cluster1018\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1018
sed -i -e 's/Cluster101/Cluster101\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster101
sed -i -e 's/Cluster1026/Cluster1026\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1026
sed -i -e 's/Cluster107/Cluster107\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster107
sed -i -e 's/Cluster113/Cluster113\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster113
sed -i -e 's/Cluster115/Cluster115\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster115
sed -i -e 's/Cluster1173/Cluster1173\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1173
sed -i -e 's/Cluster1175/Cluster1175\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1175
sed -i -e 's/Cluster1179/Cluster1179\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1179
sed -i -e 's/Cluster1180/Cluster1180\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1180
sed -i -e 's/Cluster118/Cluster118\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster118
sed -i -e 's/Cluster1210/Cluster1210\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1210
sed -i -e 's/Cluster1211/Cluster1211\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1211
sed -i -e 's/Cluster1212/Cluster1212\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1212
sed -i -e 's/Cluster1213/Cluster1213\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster1213
sed -i -e 's/Cluster1226/Cluster1226\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster1226
sed -i -e 's/Cluster122/Cluster122\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster122
sed -i -e 's/Cluster1230/Cluster1230\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1230
sed -i -e 's/Cluster1234/Cluster1234\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1234
sed -i -e 's/Cluster1240/Cluster1240\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1240
sed -i -e 's/Cluster1242/Cluster1242\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1242
sed -i -e 's/Cluster1244/Cluster1244\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/1' Cluster1244
sed -i -e 's/Cluster1249/Cluster1249\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1249
sed -i -e 's/Cluster1250/Cluster1250\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1250
sed -i -e 's/Cluster1254/Cluster1254\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1254
sed -i -e 's/Cluster1255/Cluster1255\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster1255
sed -i -e 's/Cluster1258/Cluster1258\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1258
sed -i -e 's/Cluster1378/Cluster1378\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1378
sed -i -e 's/Cluster1381/Cluster1381\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1381
sed -i -e 's/Cluster1387/Cluster1387\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1387
sed -i -e 's/Cluster1388/Cluster1388\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1388
sed -i -e 's/Cluster1390/Cluster1390\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1390
sed -i -e 's/Cluster1414/Cluster1414\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1414
sed -i -e 's/Cluster1417/Cluster1417\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1417
sed -i -e 's/Cluster142/Cluster142\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster142
sed -i -e 's/Cluster1438/Cluster1438\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1438
sed -i -e 's/Cluster1439/Cluster1439\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1439
sed -i -e 's/Cluster1443/Cluster1443\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1443
sed -i -e 's/Cluster1447/Cluster1447\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1447
sed -i -e 's/Cluster1448/Cluster1448\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster1448
sed -i -e 's/Cluster1449/Cluster1449\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster1449
sed -i -e 's/Cluster1451/Cluster1451\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1451
sed -i -e 's/Cluster1452/Cluster1452\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1452
sed -i -e 's/Cluster1455/Cluster1455\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1455
sed -i -e 's/Cluster1458/Cluster1458\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1458
sed -i -e 's/Cluster1460/Cluster1460\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1460
sed -i -e 's/Cluster1462/Cluster1462\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1462
sed -i -e 's/Cluster1468/Cluster1468\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1468
sed -i -e 's/Cluster1469/Cluster1469\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1469
sed -i -e 's/Cluster1480/Cluster1480\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1480
sed -i -e 's/Cluster150/Cluster150\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster150
sed -i -e 's/Cluster1541/Cluster1541\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster1541
sed -i -e 's/Cluster1591/Cluster1591\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1591
sed -i -e 's/Cluster1593/Cluster1593\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1593
sed -i -e 's/Cluster1594/Cluster1594\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster1594
sed -i -e 's/Cluster1620/Cluster1620\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1620
sed -i -e 's/Cluster1633/Cluster1633\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster1633
sed -i -e 's/Cluster1635/Cluster1635\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1635
sed -i -e 's/Cluster1642/Cluster1642\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1642
sed -i -e 's/Cluster1644/Cluster1644\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1644
sed -i -e 's/Cluster1647/Cluster1647\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1647
sed -i -e 's/Cluster1650/Cluster1650\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1650
sed -i -e 's/Cluster1683/Cluster1683\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1683
sed -i -e 's/Cluster1722/Cluster1722\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1722
sed -i -e 's/Cluster1744/Cluster1744\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1744
sed -i -e 's/Cluster1750/Cluster1750\tDNA-Merlin\tII\t1\tTIR\tMerlin\tNA/1' Cluster1750
sed -i -e 's/Cluster1752/Cluster1752\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1752
sed -i -e 's/Cluster1755/Cluster1755\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster1755
sed -i -e 's/Cluster1774/Cluster1774\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1774
sed -i -e 's/Cluster1776/Cluster1776\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster1776
sed -i -e 's/Cluster1800/Cluster1800\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster1800
sed -i -e 's/Cluster1802/Cluster1802\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1802
sed -i -e 's/Cluster1803/Cluster1803\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1803
sed -i -e 's/Cluster1809/Cluster1809\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1809
sed -i -e 's/Cluster1810/Cluster1810\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1810
sed -i -e 's/Cluster1814/Cluster1814\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1814
sed -i -e 's/Cluster1816/Cluster1816\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1816
sed -i -e 's/Cluster1820/Cluster1820\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1820
sed -i -e 's/Cluster1846/Cluster1846\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1846
sed -i -e 's/Cluster1889/Cluster1889\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1889
sed -i -e 's/Cluster1940/Cluster1940\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1940
sed -i -e 's/Cluster1942/Cluster1942\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1942
sed -i -e 's/Cluster1944/Cluster1944\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1944
sed -i -e 's/Cluster1946/Cluster1946\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1946
sed -i -e 's/Cluster1958/Cluster1958\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1958
sed -i -e 's/Cluster1969/Cluster1969\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1969
sed -i -e 's/Cluster1979/Cluster1979\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1979
sed -i -e 's/Cluster1989/Cluster1989\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1989
sed -i -e 's/Cluster1991/Cluster1991\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1991
sed -i -e 's/Cluster1994/Cluster1994\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1994
sed -i -e 's/Cluster1997/Cluster1997\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster1997
sed -i -e 's/Cluster2067/Cluster2067\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2067
sed -i -e 's/Cluster2109/Cluster2109\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/1' Cluster2109
sed -i -e 's/Cluster210/Cluster210\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster210
sed -i -e 's/Cluster2121/Cluster2121\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2121
sed -i -e 's/Cluster2124/Cluster2124\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2124
sed -i -e 's/Cluster2125/Cluster2125\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2125
sed -i -e 's/Cluster2129/Cluster2129\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2129
sed -i -e 's/Cluster2130/Cluster2130\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2130
sed -i -e 's/Cluster2142/Cluster2142\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2142
sed -i -e 's/Cluster2143/Cluster2143\tDNA-PiggyBac\tII\t1\tTIR\tPiggyBac\tNA/1' Cluster2143
sed -i -e 's/Cluster2149/Cluster2149\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2149
sed -i -e 's/Cluster2160/Cluster2160\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2160
sed -i -e 's/Cluster2161/Cluster2161\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2161
sed -i -e 's/Cluster2169/Cluster2169\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2169
sed -i -e 's/Cluster2176/Cluster2176\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2176
sed -i -e 's/Cluster2180/Cluster2180\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2180
sed -i -e 's/Cluster2182/Cluster2182\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster2182
sed -i -e 's/Cluster2187/Cluster2187\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2187
sed -i -e 's/Cluster2191/Cluster2191\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2191
sed -i -e 's/Cluster2192/Cluster2192\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2192
sed -i -e 's/Cluster2260/Cluster2260\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2260
sed -i -e 's/Cluster2267/Cluster2267\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2267
sed -i -e 's/Cluster2298/Cluster2298\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster2298
sed -i -e 's/Cluster2299/Cluster2299\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2299
sed -i -e 's/Cluster2310/Cluster2310\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2310
sed -i -e 's/Cluster2313/Cluster2313\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2313
sed -i -e 's/Cluster2314/Cluster2314\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2314
sed -i -e 's/Cluster2317/Cluster2317\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster2317
sed -i -e 's/Cluster2342/Cluster2342\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2342
sed -i -e 's/Cluster2344/Cluster2344\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster2344
sed -i -e 's/Cluster2347/Cluster2347\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2347
sed -i -e 's/Cluster2351/Cluster2351\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2351
sed -i -e 's/Cluster2356/Cluster2356\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2356
sed -i -e 's/Cluster2363/Cluster2363\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2363
sed -i -e 's/Cluster2389/Cluster2389\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2389
sed -i -e 's/Cluster2476/Cluster2476\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2476
sed -i -e 's/Cluster2477/Cluster2477\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2477
sed -i -e 's/Cluster2478/Cluster2478\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2478
sed -i -e 's/Cluster2481/Cluster2481\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2481
sed -i -e 's/Cluster2497/Cluster2497\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster2497
sed -i -e 's/Cluster2498/Cluster2498\tSatellite\tSatellite\tNA\tNA\tNA\tNA\t/1' Cluster2498
sed -i -e 's/Cluster2505/Cluster2505\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster2505
sed -i -e 's/Cluster2518/Cluster2518\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2518
sed -i -e 's/Cluster2520/Cluster2520\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster2520
sed -i -e 's/Cluster2524/Cluster2524\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2524
sed -i -e 's/Cluster2527/Cluster2527\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/1' Cluster2527
sed -i -e 's/Cluster2529/Cluster2529\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster2529
sed -i -e 's/Cluster2559/Cluster2559\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2559
sed -i -e 's/Cluster2620/Cluster2620\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2620
sed -i -e 's/Cluster2626/Cluster2626\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2626
sed -i -e 's/Cluster2627/Cluster2627\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/1' Cluster2627
sed -i -e 's/Cluster2629/Cluster2629\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2629
sed -i -e 's/Cluster2643/Cluster2643\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2643
sed -i -e 's/Cluster2656/Cluster2656\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2656
sed -i -e 's/Cluster2665/Cluster2665\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2665
sed -i -e 's/Cluster2667/Cluster2667\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2667
sed -i -e 's/Cluster2671/Cluster2671\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2671
sed -i -e 's/Cluster2672/Cluster2672\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2672
sed -i -e 's/Cluster2676/Cluster2676\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2676
sed -i -e 's/Cluster2679/Cluster2679\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2679
sed -i -e 's/Cluster2680/Cluster2680\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2680
sed -i -e 's/Cluster2704/Cluster2704\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2704
sed -i -e 's/Cluster2722/Cluster2722\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2722
sed -i -e 's/Cluster2786/Cluster2786\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2786
sed -i -e 's/Cluster2808/Cluster2808\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2808
sed -i -e 's/Cluster2809/Cluster2809\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster2809
sed -i -e 's/Cluster2818/Cluster2818\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/1' Cluster2818
sed -i -e 's/Cluster2829/Cluster2829\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2829
sed -i -e 's/Cluster2830/Cluster2830\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster2830
sed -i -e 's/Cluster2834/Cluster2834\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2834
sed -i -e 's/Cluster2835/Cluster2835\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2835
sed -i -e 's/Cluster2837/Cluster2837\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2837
sed -i -e 's/Cluster2838/Cluster2838\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2838
sed -i -e 's/Cluster2840/Cluster2840\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2840
sed -i -e 's/Cluster2841/Cluster2841\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster2841
sed -i -e 's/Cluster2843/Cluster2843\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2843
sed -i -e 's/Cluster2846/Cluster2846\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2846
sed -i -e 's/Cluster2847/Cluster2847\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2847
sed -i -e 's/Cluster2851/Cluster2851\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2851
sed -i -e 's/Cluster2852/Cluster2852\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2852
sed -i -e 's/Cluster2857/Cluster2857\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2857
sed -i -e 's/Cluster2866/Cluster2866\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2866
sed -i -e 's/Cluster2933/Cluster2933\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster2933
sed -i -e 's/Cluster2962/Cluster2962\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2962
sed -i -e 's/Cluster2970/Cluster2970\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2970
sed -i -e 's/Cluster2982/Cluster2982\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2982
sed -i -e 's/Cluster2985/Cluster2985\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2985
sed -i -e 's/Cluster2986/Cluster2986\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2986
sed -i -e 's/Cluster2995/Cluster2995\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2995
sed -i -e 's/Cluster2996/Cluster2996\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster2996
sed -i -e 's/Cluster2997/Cluster2997\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2997
sed -i -e 's/Cluster2999/Cluster2999\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2999
sed -i -e 's/Cluster3002/Cluster3002\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3002
sed -i -e 's/Cluster3003/Cluster3003\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3003
sed -i -e 's/Cluster3005/Cluster3005\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/1' Cluster3005
sed -i -e 's/Cluster3006/Cluster3006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3006
sed -i -e 's/Cluster3009/Cluster3009\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3009
sed -i -e 's/Cluster3029/Cluster3029\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3029
sed -i -e 's/Cluster3113/Cluster3113\tDNA-MuLE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster3113
sed -i -e 's/Cluster3114/Cluster3114\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3114
sed -i -e 's/Cluster3115/Cluster3115\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster3115
sed -i -e 's/Cluster3116/Cluster3116\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3116
sed -i -e 's/Cluster3132/Cluster3132\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3132
sed -i -e 's/Cluster3133/Cluster3133\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3133
sed -i -e 's/Cluster3134/Cluster3134\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3134
sed -i -e 's/Cluster3141/Cluster3141\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3141
sed -i -e 's/Cluster3149/Cluster3149\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA\t/1' Cluster3149
sed -i -e 's/Cluster3164/Cluster3164\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster3164
sed -i -e 's/Cluster3165/Cluster3165\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3165
sed -i -e 's/Cluster317/Cluster317\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster317
sed -i -e 's/Cluster3217/Cluster3217\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3217
sed -i -e 's/Cluster37/Cluster37\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster37
sed -i -e 's/Cluster479/Cluster479\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster479
sed -i -e 's/Cluster480/Cluster480\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster480
sed -i -e 's/Cluster481/Cluster481\tDNA-Sola-3\tII\t1\tTIR\tSola\tSola3/1' Cluster481
sed -i -e 's/Cluster483/Cluster483\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster483
sed -i -e 's/Cluster485/Cluster485\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster485
sed -i -e 's/Cluster488/Cluster488\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster488
sed -i -e 's/Cluster489/Cluster489\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster489
sed -i -e 's/Cluster491/Cluster491\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster491
sed -i -e 's/Cluster511/Cluster511\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster511
sed -i -e 's/Cluster514/Cluster514\tDNA-PIF-Harbinger\tII\t1\tTIR\tPIF-Harbinger\tNA/1' Cluster514
sed -i -e 's/Cluster524/Cluster524\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster524
sed -i -e 's/Cluster525/Cluster525\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster525
sed -i -e 's/Cluster527/Cluster527\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster527
sed -i -e 's/Cluster528/Cluster528\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster528
sed -i -e 's/Cluster540/Cluster540\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster540
sed -i -e 's/Cluster542/Cluster542\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster542
sed -i -e 's/Cluster543/Cluster543\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster543
sed -i -e 's/Cluster54/Cluster54\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster54
sed -i -e 's/Cluster600/Cluster600\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster600
sed -i -e 's/Cluster601/Cluster601\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster601
sed -i -e 's/Cluster69/Cluster69\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster69
sed -i -e 's/Cluster717/Cluster717\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster717
sed -i -e 's/Cluster720/Cluster720\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster720
sed -i -e 's/Cluster722/Cluster722\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster722
sed -i -e 's/Cluster724/Cluster724\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster724
sed -i -e 's/Cluster725/Cluster725\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster725
sed -i -e 's/Cluster726/Cluster726\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster726
sed -i -e 's/Cluster727/Cluster727\tDNA-hAT\tII\t1\tTIR\thAT\tNA\t/1' Cluster727
sed -i -e 's/Cluster729/Cluster729\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster729
sed -i -e 's/Cluster731/Cluster731\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster731
sed -i -e 's/Cluster738/Cluster738\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster738
sed -i -e 's/Cluster757/Cluster757\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster757
sed -i -e 's/Cluster758/Cluster758\tDNA\tII\tNA\tNA\tNA\tNA\t/1' Cluster758
sed -i -e 's/Cluster774/Cluster774\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster774
sed -i -e 's/Cluster789/Cluster789\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster789
sed -i -e 's/Cluster793/Cluster793\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster793
sed -i -e 's/Cluster794/Cluster794\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster794
sed -i -e 's/Cluster795/Cluster795\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster795
sed -i -e 's/Cluster798/Cluster798\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster798
sed -i -e 's/Cluster800/Cluster800\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster800
sed -i -e 's/Cluster805/Cluster805\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster805
sed -i -e 's/Cluster806/Cluster806\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster806
sed -i -e 's/Cluster80/Cluster80\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster80
sed -i -e 's/Cluster811/Cluster811\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster811
sed -i -e 's/Cluster812/Cluster812\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster812
sed -i -e 's/Cluster813/Cluster813\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster813
sed -i -e 's/Cluster816/Cluster816\tSatellite\tSatellite\tNA\tNA\tNA\tNA\t/1' Cluster816
sed -i -e 's/Cluster856/Cluster856\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster856
sed -i -e 's/Cluster929/Cluster929\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster929
sed -i -e 's/Cluster938/Cluster938\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster938
sed -i -e 's/Cluster961/Cluster961\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster961
sed -i -e 's/Cluster967/Cluster967\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster967
sed -i -e 's/Cluster975/Cluster975\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster975
sed -i -e 's/Cluster992/Cluster992\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster992
sed -i -e 's/Cluster99/Cluster99\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster99
sed -i -e 's/Cluster9/Cluster9\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster9

#put the above in sed.sh and ran it

rm sed.sh


#ok, get those things in the right order

mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/03_awk/

cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $9,$10,$11,$1,$12,$3,$4,$5,$6,$7,$14,$15}' $i > $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/03_awk/$i; done

#add header

for i in *; do echo -e "Chr\tBP_start\tBP_end\tCluster_id\tInsertion_id\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tstrand\tage" | cat - $i > $i.tmp; done

for i in *.tmp; do mv $i ${i%.tmp}; done
	#Data for insertion age chromosome position pdf's are ready!




mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/03_R_plot_pdf

cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk
#make chromosomal distribution of insertion ages for each cluster
for i in *; do Rscript $wkdir/scripts/insertion_age_chr_plots.R $i; done &
	#make sure to change output directory in R scripts
cd 03_R_plot_pdf

pdfunite *.pdf inopinata_insertion_ages.pdf
	#THIS is the huge .pdf of insertion age by chromosomal position for each repeat cluster (untrimmed alignments)


#now also get distributions of insertion ages per cluster in histogram form

mkdir $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/04_age_histograms

cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/02_awk


for i in *; do Rscript $wkdir/scripts/age_histograms.R $i; done &
	#remember to change output directory in R scripts


cd $wkdir/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/04_age_histograms

pdfunite *.pdf inopinata_insertion_age_histograms.pdf
	#THIS is the huge .pdf of insertion age histograms for each repeat cluster (untrimmed alignments)




#next, get the insertion ages using trimmed alignments

#however, the number of clusters included must be pruned 


#continuing now, this time just using trimal-trimmed alignments. using the "automated1" option this time. alignments in folder 11_trimal_automated1

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1

#get lengths of all alignments for filtering (really short alignments seem pretty bad for phylogenetics and estimating insertion ages)


mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/00_awk

#this will get the alignment length in bp for each alignment

cd $wkdir/55_inopinata_insertion_ages/11_trimal_automated1

for i in *; do awk '{/>/&&++a||b+=length()}END{print b/a}' $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/00_awk/$i; done &

#ok, get the data in usable form

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/01_awk

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/00_awk/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"}  {print FILENAME,$0}' $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/01_awk/$i; done &


mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/02_cat

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/01_awk/

cat * > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/02_cat/alignment_lengths

#265 alignments
#get the alignments > 49 bp

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/02_cat

mv alignment_lengths 00_alignment_lengths

awk '$2 > 49' 00_alignment_lengths > 01_alignments_grtr_than_49bp

#165 clusters

#get the alignments > 39 bp

awk '$2 > 39' 00_alignment_lengths > 02_alignments_grtr_than_39bp

#175

#get the alignments > 29 bp

awk '$2 > 29' 00_alignment_lengths > 03_alignments_grtr_than_29bp

#189

#ok, do not gain that much. a bunch of alignments are not great apparently?

awk 'BEGIN {FS="\t"} {OFS="\t"}  {print $1}' 01_alignments_grtr_than_49bp > 04_alignments_grtr_than_49bp_cluster_ids

#copy only alignments >49bp to new folder

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/03_cp/

cd $wkdir/55_inopinata_insertion_ages/11_trimal_automated1

xargs -a $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/02_cat/04_alignments_grtr_than_49bp_cluster_ids cp -t $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths/03_cp/

#that worked.

#ok, now exclude sequences that are all gaps


cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/00_alignment_lengths

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/00_linearize_seqs

#linearize sequences

cd 03_cp
    #thanks https://www.biostars.org/p/17680/

for i in *; do sed -e 's/\(^>.*$\)/#\1#/' $i | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/00_linearize_seqs/$i; done &


#ok, now remove sequences that are all gaps

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/01_perl_remove_all_gaps/

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/00_linearize_seqs/


for i in *; do perl -p -e 's/^-*\n//' $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/01_perl_remove_all_gaps/$i; done &

#cool that worked

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/02_awk_remove_empty_seq/

#remove "empty" sequences

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/01_perl_remove_all_gaps/


    #thanks https://www.biostars.org/p/78786/
for i in *; do awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/02_awk_remove_empty_seq/$i; done &

#ok, this shit is now ready for fasttree!

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/02_FastTree/

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/01_remove_sequences_all_gaps/02_awk_remove_empty_seq/

for i in *; do FastTree -nt -gtr -quiet < $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/02_FastTree/$i &


#get branch lengths

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/03_branch_lengths_trimal_automated1

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/02_FastTree/

for i in *; do Rscript $wkdir/scripts/branch_lengths.R $i; done &
	#remember to change output directory in R scripts


cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/03_branch_lengths_trimal_automated1

#clean up data, remove header

for i in *; do sed -i -e 's/_R_//g' $i; done &

for i in *; do sed -i '1d' $i; done &

#combine files
mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/05_cat_insertion_ages

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/03_branch_lengths_trimal_automated1

cat * | sort | uniq > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/05_cat_insertion_ages/all_insertion_ages

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/05_cat_insertion_ages/


awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' all_insertion_ages > insertion_ids

#

#first, all insertion ages by chr position.........

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position

cd 06_global_ages_genomic_position

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/00_links

#just get the inserion ids that I used here



LC_ALL=C fgrep -w -f $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/05_cat_insertion_ages/insertion_ids $wkdir/55_inopinata_insertion_ages/14_branch_lengths_no_trimal_by_chr_pos/00_cat_bedfiles/id_key_all_insertions_bed >  $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/00_links/id_key_trimal_automated1_insertions_bed &

#join positions with ages

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/01_merge

perl $wkdir/scripts/merge.pl -k -e "no_key" $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/00_links/id_key_trimal_automated1_insertions_bed $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/05_cat_insertion_ages/all_insertion_ages 2> $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/inopinata_insertion_locations_ages.error > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/01_merge/inopinata_insertion_locations_ages &

#ok, get ready for bedtools......


mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/02_awk_bed/

cd  $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/01_merge/
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$8}' inopinata_insertion_locations_ages > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/02_awk_bed/inopinata_insertion_ages_bed


#sort 

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position//03_bedtools_sort

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/02_awk_bed/

bedtools sort -i inopinata_insertion_ages_bed > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position//03_bedtools_sort/inopinata_insertion_ages_bed


#ok, bedtools...

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/04_bedtools_map

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/03_bedtools_sort/

bedtools map -o mean -c 4 -a $wkdir/windows/inopinata.10kb.windows -b inopinata_insertion_ages_bed > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/04_bedtools_map/inopinata_insertion_ages_mean_10kb_windows_bed

#remove missing windows

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/05_awk

cd  $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/04_bedtools_map

awk 'BEGIN {OFS="\t"} $4 != "."' inopinata_insertion_ages_mean_10kb_windows_bed > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/05_awk/inopinata_insertion_ages_mean_10kb_windows_bed

#remove column three add header...

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/05_awk

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2+1,$4}' inopinata_insertion_ages_mean_10kb_windows_bed > inopinata_insertion_ages_mean_10kb_windows_bed.tmp

echo -e "Chr\tBP\tinsertion_age" | cat - inopinata_insertion_ages_mean_10kb_windows_bed.tmp > inopinata_insertion_ages_mean_10kb_windows_bed

#ok, that one is done.


#now, more!

#link with taxonomy

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/00_links

cp $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/06_global_ages_genomic_position/01_merge/inopinata_insertion_locations_ages $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/00_links/00_inopinata_insertion_locations_ages

#get rid of first underscore to get just cluster ids for joining

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/00_links

sed -i -e 's/_/\t/1' 00_inopinata_insertion_locations_ages

#get unique cluster ids

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' 00_inopinata_insertion_locations_ages | sort | uniq > 01_unique_cluster_ids

#awk by cluster id

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/01_awk

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/01_awk

awk '{print > $1}' $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/00_links/00_inopinata_insertion_locations_ages

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/01_awk

#put in the repeat taxonomy

sed -i -e 's/Cluster1006/Cluster1006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1006
sed -i -e 's/Cluster1007/Cluster1007\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1007
sed -i -e 's/Cluster100/Cluster100\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster100
sed -i -e 's/Cluster1010/Cluster1010\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1010
sed -i -e 's/Cluster1012/Cluster1012\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1012
sed -i -e 's/Cluster1013/Cluster1013\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1013
sed -i -e 's/Cluster1018/Cluster1018\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1018
sed -i -e 's/Cluster101/Cluster101\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster101
sed -i -e 's/Cluster1026/Cluster1026\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1026
sed -i -e 's/Cluster107/Cluster107\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster107
sed -i -e 's/Cluster113/Cluster113\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster113
sed -i -e 's/Cluster115/Cluster115\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster115
sed -i -e 's/Cluster1173/Cluster1173\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1173
sed -i -e 's/Cluster1175/Cluster1175\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1175
sed -i -e 's/Cluster1179/Cluster1179\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1179
sed -i -e 's/Cluster1180/Cluster1180\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1180
sed -i -e 's/Cluster118/Cluster118\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster118
sed -i -e 's/Cluster1210/Cluster1210\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1210
sed -i -e 's/Cluster1211/Cluster1211\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1211
sed -i -e 's/Cluster1212/Cluster1212\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1212
sed -i -e 's/Cluster1213/Cluster1213\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster1213
sed -i -e 's/Cluster1226/Cluster1226\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster1226
sed -i -e 's/Cluster122/Cluster122\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster122
sed -i -e 's/Cluster1230/Cluster1230\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1230
sed -i -e 's/Cluster1234/Cluster1234\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1234
sed -i -e 's/Cluster1240/Cluster1240\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1240
sed -i -e 's/Cluster1242/Cluster1242\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1242
sed -i -e 's/Cluster1244/Cluster1244\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/1' Cluster1244
sed -i -e 's/Cluster1249/Cluster1249\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1249
sed -i -e 's/Cluster1250/Cluster1250\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1250
sed -i -e 's/Cluster1254/Cluster1254\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1254
sed -i -e 's/Cluster1255/Cluster1255\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster1255
sed -i -e 's/Cluster1258/Cluster1258\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1258
sed -i -e 's/Cluster1378/Cluster1378\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1378
sed -i -e 's/Cluster1381/Cluster1381\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1381
sed -i -e 's/Cluster1387/Cluster1387\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1387
sed -i -e 's/Cluster1388/Cluster1388\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster1388
sed -i -e 's/Cluster1390/Cluster1390\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1390
sed -i -e 's/Cluster1414/Cluster1414\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1414
sed -i -e 's/Cluster1417/Cluster1417\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1417
sed -i -e 's/Cluster142/Cluster142\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster142
sed -i -e 's/Cluster1438/Cluster1438\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1438
sed -i -e 's/Cluster1439/Cluster1439\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1439
sed -i -e 's/Cluster1443/Cluster1443\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1443
sed -i -e 's/Cluster1447/Cluster1447\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1447
sed -i -e 's/Cluster1448/Cluster1448\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster1448
sed -i -e 's/Cluster1449/Cluster1449\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster1449
sed -i -e 's/Cluster1451/Cluster1451\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1451
sed -i -e 's/Cluster1452/Cluster1452\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1452
sed -i -e 's/Cluster1455/Cluster1455\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1455
sed -i -e 's/Cluster1458/Cluster1458\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1458
sed -i -e 's/Cluster1460/Cluster1460\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1460
sed -i -e 's/Cluster1462/Cluster1462\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1462
sed -i -e 's/Cluster1468/Cluster1468\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1468
sed -i -e 's/Cluster1469/Cluster1469\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1469
sed -i -e 's/Cluster1480/Cluster1480\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1480
sed -i -e 's/Cluster150/Cluster150\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster150
sed -i -e 's/Cluster1541/Cluster1541\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster1541
sed -i -e 's/Cluster1591/Cluster1591\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1591
sed -i -e 's/Cluster1593/Cluster1593\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1593
sed -i -e 's/Cluster1594/Cluster1594\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster1594
sed -i -e 's/Cluster1620/Cluster1620\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1620
sed -i -e 's/Cluster1633/Cluster1633\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster1633
sed -i -e 's/Cluster1635/Cluster1635\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1635
sed -i -e 's/Cluster1642/Cluster1642\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1642
sed -i -e 's/Cluster1644/Cluster1644\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1644
sed -i -e 's/Cluster1647/Cluster1647\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1647
sed -i -e 's/Cluster1650/Cluster1650\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1650
sed -i -e 's/Cluster1683/Cluster1683\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1683
sed -i -e 's/Cluster1722/Cluster1722\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1722
sed -i -e 's/Cluster1744/Cluster1744\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1744
sed -i -e 's/Cluster1750/Cluster1750\tDNA-Merlin\tII\t1\tTIR\tMerlin\tNA/1' Cluster1750
sed -i -e 's/Cluster1752/Cluster1752\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1752
sed -i -e 's/Cluster1755/Cluster1755\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster1755
sed -i -e 's/Cluster1774/Cluster1774\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1774
sed -i -e 's/Cluster1776/Cluster1776\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster1776
sed -i -e 's/Cluster1800/Cluster1800\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster1800
sed -i -e 's/Cluster1802/Cluster1802\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1802
sed -i -e 's/Cluster1803/Cluster1803\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1803
sed -i -e 's/Cluster1809/Cluster1809\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1809
sed -i -e 's/Cluster1810/Cluster1810\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1810
sed -i -e 's/Cluster1814/Cluster1814\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster1814
sed -i -e 's/Cluster1816/Cluster1816\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1816
sed -i -e 's/Cluster1820/Cluster1820\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1820
sed -i -e 's/Cluster1846/Cluster1846\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1846
sed -i -e 's/Cluster1889/Cluster1889\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1889
sed -i -e 's/Cluster1940/Cluster1940\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1940
sed -i -e 's/Cluster1942/Cluster1942\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1942
sed -i -e 's/Cluster1944/Cluster1944\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1944
sed -i -e 's/Cluster1946/Cluster1946\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1946
sed -i -e 's/Cluster1958/Cluster1958\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1958
sed -i -e 's/Cluster1969/Cluster1969\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster1969
sed -i -e 's/Cluster1979/Cluster1979\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster1979
sed -i -e 's/Cluster1989/Cluster1989\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1989
sed -i -e 's/Cluster1991/Cluster1991\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster1991
sed -i -e 's/Cluster1994/Cluster1994\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster1994
sed -i -e 's/Cluster1997/Cluster1997\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster1997
sed -i -e 's/Cluster2067/Cluster2067\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2067
sed -i -e 's/Cluster2109/Cluster2109\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/1' Cluster2109
sed -i -e 's/Cluster210/Cluster210\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster210
sed -i -e 's/Cluster2121/Cluster2121\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2121
sed -i -e 's/Cluster2124/Cluster2124\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2124
sed -i -e 's/Cluster2125/Cluster2125\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2125
sed -i -e 's/Cluster2129/Cluster2129\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2129
sed -i -e 's/Cluster2130/Cluster2130\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2130
sed -i -e 's/Cluster2142/Cluster2142\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2142
sed -i -e 's/Cluster2143/Cluster2143\tDNA-PiggyBac\tII\t1\tTIR\tPiggyBac\tNA/1' Cluster2143
sed -i -e 's/Cluster2149/Cluster2149\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2149
sed -i -e 's/Cluster2160/Cluster2160\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2160
sed -i -e 's/Cluster2161/Cluster2161\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2161
sed -i -e 's/Cluster2169/Cluster2169\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2169
sed -i -e 's/Cluster2176/Cluster2176\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2176
sed -i -e 's/Cluster2180/Cluster2180\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2180
sed -i -e 's/Cluster2182/Cluster2182\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster2182
sed -i -e 's/Cluster2187/Cluster2187\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2187
sed -i -e 's/Cluster2191/Cluster2191\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2191
sed -i -e 's/Cluster2192/Cluster2192\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2192
sed -i -e 's/Cluster2260/Cluster2260\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2260
sed -i -e 's/Cluster2267/Cluster2267\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2267
sed -i -e 's/Cluster2298/Cluster2298\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster2298
sed -i -e 's/Cluster2299/Cluster2299\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2299
sed -i -e 's/Cluster2310/Cluster2310\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2310
sed -i -e 's/Cluster2313/Cluster2313\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2313
sed -i -e 's/Cluster2314/Cluster2314\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2314
sed -i -e 's/Cluster2317/Cluster2317\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster2317
sed -i -e 's/Cluster2342/Cluster2342\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2342
sed -i -e 's/Cluster2344/Cluster2344\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster2344
sed -i -e 's/Cluster2347/Cluster2347\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2347
sed -i -e 's/Cluster2351/Cluster2351\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2351
sed -i -e 's/Cluster2356/Cluster2356\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2356
sed -i -e 's/Cluster2363/Cluster2363\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2363
sed -i -e 's/Cluster2389/Cluster2389\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2389
sed -i -e 's/Cluster2476/Cluster2476\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2476
sed -i -e 's/Cluster2477/Cluster2477\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2477
sed -i -e 's/Cluster2478/Cluster2478\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2478
sed -i -e 's/Cluster2481/Cluster2481\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2481
sed -i -e 's/Cluster2497/Cluster2497\tDNA-MULE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster2497
sed -i -e 's/Cluster2498/Cluster2498\tSatellite\tSatellite\tNA\tNA\tNA\tNA/1' Cluster2498
sed -i -e 's/Cluster2505/Cluster2505\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster2505
sed -i -e 's/Cluster2518/Cluster2518\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2518
sed -i -e 's/Cluster2520/Cluster2520\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster2520
sed -i -e 's/Cluster2524/Cluster2524\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2524
sed -i -e 's/Cluster2527/Cluster2527\tDNA-Maverick\tII\t2\tMaverick\tMaverick\tNA/1' Cluster2527
sed -i -e 's/Cluster2529/Cluster2529\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster2529
sed -i -e 's/Cluster2559/Cluster2559\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2559
sed -i -e 's/Cluster2620/Cluster2620\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2620
sed -i -e 's/Cluster2626/Cluster2626\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2626
sed -i -e 's/Cluster2627/Cluster2627\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/1' Cluster2627
sed -i -e 's/Cluster2629/Cluster2629\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2629
sed -i -e 's/Cluster2643/Cluster2643\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2643
sed -i -e 's/Cluster2656/Cluster2656\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2656
sed -i -e 's/Cluster2665/Cluster2665\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2665
sed -i -e 's/Cluster2667/Cluster2667\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2667
sed -i -e 's/Cluster2671/Cluster2671\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2671
sed -i -e 's/Cluster2672/Cluster2672\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2672
sed -i -e 's/Cluster2676/Cluster2676\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2676
sed -i -e 's/Cluster2679/Cluster2679\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2679
sed -i -e 's/Cluster2680/Cluster2680\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2680
sed -i -e 's/Cluster2704/Cluster2704\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2704
sed -i -e 's/Cluster2722/Cluster2722\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2722
sed -i -e 's/Cluster2786/Cluster2786\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2786
sed -i -e 's/Cluster2808/Cluster2808\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2808
sed -i -e 's/Cluster2809/Cluster2809\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster2809
sed -i -e 's/Cluster2818/Cluster2818\tDNA-hAT-Ac\tII\t1\tTIR\thAT\tAc/1' Cluster2818
sed -i -e 's/Cluster2829/Cluster2829\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2829
sed -i -e 's/Cluster2830/Cluster2830\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster2830
sed -i -e 's/Cluster2834/Cluster2834\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2834
sed -i -e 's/Cluster2835/Cluster2835\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2835
sed -i -e 's/Cluster2837/Cluster2837\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2837
sed -i -e 's/Cluster2838/Cluster2838\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2838
sed -i -e 's/Cluster2840/Cluster2840\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster2840
sed -i -e 's/Cluster2841/Cluster2841\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster2841
sed -i -e 's/Cluster2843/Cluster2843\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2843
sed -i -e 's/Cluster2846/Cluster2846\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2846
sed -i -e 's/Cluster2847/Cluster2847\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2847
sed -i -e 's/Cluster2851/Cluster2851\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2851
sed -i -e 's/Cluster2852/Cluster2852\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster2852
sed -i -e 's/Cluster2857/Cluster2857\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2857
sed -i -e 's/Cluster2866/Cluster2866\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2866
sed -i -e 's/Cluster2933/Cluster2933\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster2933
sed -i -e 's/Cluster2962/Cluster2962\tRC-Helitron\tII\t2\tHelitron\tHelitron\tNA/1' Cluster2962
sed -i -e 's/Cluster2970/Cluster2970\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster2970
sed -i -e 's/Cluster2982/Cluster2982\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2982
sed -i -e 's/Cluster2985/Cluster2985\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2985
sed -i -e 's/Cluster2986/Cluster2986\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2986
sed -i -e 's/Cluster2995/Cluster2995\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster2995
sed -i -e 's/Cluster2996/Cluster2996\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster2996
sed -i -e 's/Cluster2997/Cluster2997\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2997
sed -i -e 's/Cluster2999/Cluster2999\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster2999
sed -i -e 's/Cluster3002/Cluster3002\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3002
sed -i -e 's/Cluster3003/Cluster3003\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3003
sed -i -e 's/Cluster3005/Cluster3005\tLTR-DIRS\tI\tNA\tDIRS\tDIRS\tNA/1' Cluster3005
sed -i -e 's/Cluster3006/Cluster3006\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3006
sed -i -e 's/Cluster3009/Cluster3009\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3009
sed -i -e 's/Cluster3029/Cluster3029\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3029
sed -i -e 's/Cluster3113/Cluster3113\tDNA-MuLE-MuDR\tII\t1\tTIR\tMutator\tMuDR/1' Cluster3113
sed -i -e 's/Cluster3114/Cluster3114\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3114
sed -i -e 's/Cluster3115/Cluster3115\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster3115
sed -i -e 's/Cluster3116/Cluster3116\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3116
sed -i -e 's/Cluster3132/Cluster3132\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster3132
sed -i -e 's/Cluster3133/Cluster3133\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3133
sed -i -e 's/Cluster3134/Cluster3134\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3134
sed -i -e 's/Cluster3141/Cluster3141\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster3141
sed -i -e 's/Cluster3149/Cluster3149\tLTR-Gypsy\tI\tNA\tLTR\tGypsy\tNA/1' Cluster3149
sed -i -e 's/Cluster3164/Cluster3164\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster3164
sed -i -e 's/Cluster3165/Cluster3165\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3165
sed -i -e 's/Cluster317/Cluster317\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster317
sed -i -e 's/Cluster3217/Cluster3217\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster3217
sed -i -e 's/Cluster37/Cluster37\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster37
sed -i -e 's/Cluster479/Cluster479\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster479
sed -i -e 's/Cluster480/Cluster480\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster480
sed -i -e 's/Cluster481/Cluster481\tDNA-Sola-3\tII\t1\tTIR\tSola\tSola3/1' Cluster481
sed -i -e 's/Cluster483/Cluster483\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster483
sed -i -e 's/Cluster485/Cluster485\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster485
sed -i -e 's/Cluster488/Cluster488\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster488
sed -i -e 's/Cluster489/Cluster489\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster489
sed -i -e 's/Cluster491/Cluster491\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster491
sed -i -e 's/Cluster511/Cluster511\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster511
sed -i -e 's/Cluster514/Cluster514\tDNA-PIF-Harbinger\tII\t1\tTIR\tPIF-Harbinger\tNA/1' Cluster514
sed -i -e 's/Cluster524/Cluster524\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster524
sed -i -e 's/Cluster525/Cluster525\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster525
sed -i -e 's/Cluster527/Cluster527\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster527
sed -i -e 's/Cluster528/Cluster528\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster528
sed -i -e 's/Cluster540/Cluster540\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster540
sed -i -e 's/Cluster542/Cluster542\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster542
sed -i -e 's/Cluster543/Cluster543\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster543
sed -i -e 's/Cluster54/Cluster54\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster54
sed -i -e 's/Cluster600/Cluster600\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster600
sed -i -e 's/Cluster601/Cluster601\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster601
sed -i -e 's/Cluster69/Cluster69\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster69
sed -i -e 's/Cluster717/Cluster717\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster717
sed -i -e 's/Cluster720/Cluster720\tSimple_repeat\tSimple_repeat\tNA\tNA\tNA\tNA/1' Cluster720
sed -i -e 's/Cluster722/Cluster722\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster722
sed -i -e 's/Cluster724/Cluster724\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster724
sed -i -e 's/Cluster725/Cluster725\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster725
sed -i -e 's/Cluster726/Cluster726\tDNA-TcMar-Tc2\tII\t1\tTIR\tTc1-Mariner\tTc2/1' Cluster726
sed -i -e 's/Cluster727/Cluster727\tDNA-hAT\tII\t1\tTIR\thAT\tNA/1' Cluster727
sed -i -e 's/Cluster729/Cluster729\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster729
sed -i -e 's/Cluster731/Cluster731\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster731
sed -i -e 's/Cluster738/Cluster738\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster738
sed -i -e 's/Cluster757/Cluster757\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster757
sed -i -e 's/Cluster758/Cluster758\tDNA\tII\tNA\tNA\tNA\tNA/1' Cluster758
sed -i -e 's/Cluster774/Cluster774\tUnknown\tNA\tNA\tNA\tNA\tNA/1' Cluster774
sed -i -e 's/Cluster789/Cluster789\tLTR-Pao\tI\tNA\tLTR\tBel-Pao\tPao/1' Cluster789
sed -i -e 's/Cluster793/Cluster793\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster793
sed -i -e 's/Cluster794/Cluster794\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster794
sed -i -e 's/Cluster795/Cluster795\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster795
sed -i -e 's/Cluster798/Cluster798\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster798
sed -i -e 's/Cluster800/Cluster800\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster800
sed -i -e 's/Cluster805/Cluster805\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster805
sed -i -e 's/Cluster806/Cluster806\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster806
sed -i -e 's/Cluster80/Cluster80\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster80
sed -i -e 's/Cluster811/Cluster811\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster811
sed -i -e 's/Cluster812/Cluster812\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster812
sed -i -e 's/Cluster813/Cluster813\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster813
sed -i -e 's/Cluster816/Cluster816\tSatellite\tSatellite\tNA\tNA\tNA\tNA/1' Cluster816
sed -i -e 's/Cluster856/Cluster856\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster856
sed -i -e 's/Cluster929/Cluster929\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster929
sed -i -e 's/Cluster938/Cluster938\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster938
sed -i -e 's/Cluster961/Cluster961\tDNA-TcMar-Tc1\tII\t1\tTIR\tTc1-Mariner\tTc1/1' Cluster961
sed -i -e 's/Cluster967/Cluster967\tDNA-TcMar-Tc4\tII\t1\tTIR\tTc1-Mariner\tTc4/1' Cluster967
sed -i -e 's/Cluster975/Cluster975\tLINE-CR1\tI\tNA\tLINE\tJockey\tCR1/1' Cluster975
sed -i -e 's/Cluster992/Cluster992\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster992
sed -i -e 's/Cluster99/Cluster99\tLINE-RTE-RTE\tI\tNA\tLINE\tRTE\tNA/1' Cluster99
sed -i -e 's/Cluster9/Cluster9\tDNA-TcMar-Mariner\tII\t1\tTIR\tTc1-Mariner\tMariner/1' Cluster9


#put this in sed.sh and ran it

#sweet

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/01_awk

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/02_awk


#get only the columns i want in an order that makes sense

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $9,$10,$11,$1,$12,$3,$4,$5,$6,$7,$14,$15}' $i > $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/02_awk/$i; done

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon/02_awk

#get some headers on there
for i in *; do echo -e "Chr\tBP_start\tBP_end\tCluster_id\tInsertion_id\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tstrand\tage" | cat - $i > $i.tmp; done

for i in *.tmp; do mv $i ${i%.tmp}; done

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1/07_merge_age_insertion_ch_position_taxon


mkdir 03_R_plot_pdf

mkdir 04_age_histograms


#ok, make some plots ; this plots insertion ages by chr position and insertion age histograms FOR EVERY CLUSTER in C. inopinata

for i in *; do Rscript $wkdir/scripts/insertion_age_chr_plots.R $i; done &
	#remember to change output directory in R scripts

for i in *; do Rscript $wkdir/scripts/age_histograms.R $i; done &
	#remember to change output directory in R scripts


cd 03_R_plot_pdf

pdfunite *.pdf inopinata_insertion_ages.pdf
	#THIS is the huge .pdf of insertion age by chromosomal position for each repeat cluster (trimmed alignments)


cd 04_age_histograms

pdfunite *.pdf inopinata_insertion_age_histograms.pdf
	#THIS is the huge .pdf of insertion age histograms for each repeat cluster (untrimmed alignments)

#combine age summary statistics for plotting 

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1//08_cat_summary

mkdir $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1//08_cat_summary/00_sed

cd $wkdir/55_inopinata_insertion_ages/16_branch_lengths_trimal_automated1//04_branch_lengths_trimal_automated1_summary

for i in *; do sed '1d' $i > $wkdir/08_cat_summary/00_sed/$i; done

mkdir $wkdir/08_cat_summary/01_awk

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $i > $wkdir/08_cat_summary/01_awk/$i; done

#used sed.sh to add repeat taxa

cat * > $wkdir/08_cat_summary/branch_lengths_summaries

cd ..

echo -e "cluster_id\trepeat_class\trepeat_subclass\trepeat_order\trepeat_superfamily\trepeat_family\tmean\tmedian\tmin\tmax\tsd\tiqr\tcount" | cat - branch_lengths_summaries > branch_lengths_summaries_tmp && mv branch_lengths_summaries_tmp branch_lengths_summaries
	#THIS is what is used for Supplemental Figures 14-15!
