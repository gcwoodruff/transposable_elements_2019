#this extracts the terminal branch lengths of a newick phylogenetic tree file 
#it also prints summary statistics of those branch lengths
library(phytools)

args <- commandArgs(TRUE)

tree_file <- args[1]
tree_file_chr <- as.character(args[1])

tree <- read.tree(file = tree_file)

print(args)

print(tree_file)

print(tree_file_chr)

terminal_branch_lengths <-setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)
#thanks http://blog.phytools.org/2013/10/finding-edge-lengths-of-all-terminal.html

terminal_branch_length_df <- data.frame(inersion_id = names(terminal_branch_lengths), branch_length = terminal_branch_lengths, row.names = NULL)

summary_df <- data.frame(mean = mean(terminal_branch_length_df$branch_length), median = median(terminal_branch_length_df$branch_length), min = min(terminal_branch_length_df$branch_length), max= max(terminal_branch_length_df$branch_length), sd = sd(terminal_branch_length_df$branch_length), iqr = IQR(terminal_branch_length_df$branch_length), count = nrow(terminal_branch_length_df))

file_1 <- paste("/Users/gavin/genome/repeats_continued_7-25-19/12_branch_lengths_no_trimal/",tree_file_chr, sep = "", collapse = NULL)

file_2 <- paste("/Users/gavin/genome/repeats_continued_7-25-19/13_branch_lengths_no_trimal_summary/",tree_file_chr, sep = "", collapse = NULL)

print(file_1)

print(file_2)

write.table(terminal_branch_length_df, file = file_1, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(summary_df, file = file_2, row.names = FALSE, quote = FALSE, sep = "\t")
