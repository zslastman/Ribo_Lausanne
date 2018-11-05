library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)

filter<-dplyr::filter
slice<-dplyr::slice
filter<-dplyr::filter
filter<-dplyr::filter

source('../src/functions.R')

extract_id <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = sampleids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}



#initiaal form when testing
#let's look at the scores for our
cell_lines<-c('OD5P','OMM','ONVC') 
satannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)


#now filter out the weird GRanges columns, for now, and aggregate into one data table 
all_orfs <- satannorfs%>%map(.%>%.$ORFs_tx)

for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character
orfs_dt <- all_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')


# orfs_dt%<>%group_by(sample)%>%mutate(fdr =  p.adjust(pval, method='fdr'))
# orfs_dt %<>% dplyr::filter(fdr<0.05)%>%ungroup

orfs_dt%>%nrow

#TODO - eveyrthing in linc

n_genes_translated <- orfs_dt %>% group_by(sample)%>%dplyr::summarise(n_genes=n_distinct(gene_id))

readthroughs <- satannorfs%>%map_df(.%>%.$ORFs_readthroughs%>%length)%>%stack%>%set_colnames(c('n_readthroughs','sample'))
readthroughs$sample%<>%as.character
readthroughs$n_readthroughs%<>%as.numeric

genetypehits <- satannorfs%>%lapply(
	.%>%.$ORFs_tx%>%.[T,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame%>%distinct%>%group_by(gene_biotype)%>%tally%>%
	mutate(gene_biotype = ifelse(gene_biotype%>%str_detect('pseudogene'),'pseudogene',gene_biotype))%>%
	group_by(gene_biotype)%>%summarise(n=sum(n))%>%
	dplyr::filter(gene_biotype %in% c('antisense','lincRNA','protein_coding','pseudogene')))%>%
	bind_rows(.id='sample')%>%
	spread(gene_biotype,n)

satann_summary_table <- genetypehits%>%left_join(readthroughs)

satann_summary_table%>%write_tsv('satann_summary.tsv'%T>%{message(normalizePath(.))})

i

satannorfs[[1]]$ORFs_tx[,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame