#Get rpkms for our genes
message('loading libraries')
library(magrittr)
library(data.table)
suppressMessages(library(stringr))
library(rtracklayer)
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(edgeR))
suppressMessages(library(tidyverse))

message('...done')

#load arguments
args <- c(
	allfeatcounts = 'feature_counts/all_feature_counts',
	genelengthfile = 'feature_counts/data/OD5P_05_uM_DAC_1/OD5P_05_uM_DAC_1.feature_counts',
	outfile = 'feature_counts/all_rpkms.txt'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

genelengthdf <- genelengthfile%>%fread
allcountdf <- allfeatcounts%>%fread

genelengths <- genelengthdf$Length[match(allcountdf$feature_id,genelengthdf$Geneid)]

allanno <- import('my_gencode.v24lift37.annotation.gtf')

y <- DGEList(counts=allcountdf[,-1],genes=data.frame(gene_id = allcountdf$feature_id,Length=genelengths))
y <- calcNormFactors(y)
RPKM <- rpkm(y)

RPKM%>%as.data.frame%>%
	mutate(gene_id=allcountdf$feature_id)%>%
	mutate(gene_id=allcountdf$feature_id)%>%
	dplyr::select(-matches('_L5|_L7'))%>%
	dplyr::select(gene_id,everything())%>%
	write_tsv(outfile)