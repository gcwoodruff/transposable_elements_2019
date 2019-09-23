#this script is by Michelle Stitzer https://github.com/mcstitzer/w22_te_annotation/blob/master/combine_all_TEs.R
#it will generate a disjoined TE gff

library(rtracklayer)
library(data.table)
library(stringr)
library(dplyr)

#command line arguments

args <- commandArgs(TRUE)
rep_file <- args[1]
rep_file_chr <- as.character(args[1]) 


#import the gff
  # te=import.gff3('/home/gavin/genome/genome/repeats_12-18-18/new_gff_idea_8-14-19/B73.structuralTEv2.fulllength.2018-09-19.gff3')

te=import.gff3(rep_file)


### get ready to disjoin
te.dj=disjoin(te, ignore.strand=T)    ### disjoin overlapping TEs so each bp in genome is only covered once.
## relate these disjoined bits to the TE they belong to. 
te.o=findOverlaps(te.dj, te, ignore.strand=T)
mcols(te.dj)=splitAsList(mcols(te)$ID[subjectHits(te.o)],queryHits(te.o))
mcols(te.dj)$ID=sapply(mcols(te.dj)$X, function(x) paste(as.character(x), collapse=','))
mcols(te.dj)$ID=sapply(te.dj$ID, function(x) unlist(str_split(x, ','))[str_count(x, ',')+1])

## more useful is any that are intact
tei=findOverlaps(te.dj, te, type='equal')
te.dj$intact=F
te.dj$intact[queryHits(tei)]=T

## export the gff
  # export.gff3(te.dj, '/home/gavin/genome/genome/repeats_12-18-18/new_gff_idea_8-14-19/B73.allTE.disjoined_gcw_8-15-19.gff3')                 

outfile_title <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/60_disjoined_gff/13_disjoin_gff_args_R/",rep_file_chr, sep = "", collapse = NULL)

export.gff3(te.dj, outfile_title)
