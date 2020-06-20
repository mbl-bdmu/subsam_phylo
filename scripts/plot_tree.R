#!/usr/bin/env R
library("ggtree")
args <- commandArgs(trailingOnly = T)
tname <- args[1]
tre <- read.tree(tname)
outpng <- paste(tname, ".png", sep = "")
outpdf <- paste(tname, ".pdf", sep = "")

png(file = outpng)
ggtree(tre, right = F) + geom_treescale(x = 0.05, y = 1) + geom_tiplab(size=1.5, color="purple")
dev.off()

pdf(file = outpdf)
ggtree(tre, right = F) + geom_treescale(x = 0.05, y = 1) + geom_tiplab(size=1.5, color="purple")
dev.off()

# TO DO: add viz of treetime's rerooted phylo with branch lengths 
# scaled to time using timetree's estimated substitution rate under 
# strict clock:branch length = rate of change * time, therefore t=bl/u
