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


#pie charts of orf types
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


ccdsorftypenames <- c('CCDS_ORFS','uORFs','non CCDS_ORFs','NMD','dORFs')
nccdsorftypenames <- c('Processed _transcript','other','lincRNA','Antisense','Retained Intron','Pseudogene')

sample_i = 'OD5P'

orfs_dt%<>%group_by(sample)%>%
	mutate(ccdscat = case_when(
		compatible_biotype %in%c('nonsense_mediated_decay') ~ 'NMD',
		ORF_category_Tx_compatible %in%c('dORF','overl_dORF') ~ 'dORFs',
		ORF_category_Tx_compatible %in%c('uORF') ~ 'uORF',
		ORF_category_Tx_compatible %in%c('novel') ~ 'non CCDS_ORFs',
		TRUE  ~ 'CCDS_ORFs',
))

orfs_dt%<>%group_by(sample)%>%
	mutate(nccdscat = case_when(
		ccdscat != 'non CCDS_ORFs' ~ 'NA',
		transcript_biotype %in% c('lincRNA') ~ 'lincRNA',
		transcript_biotype %in%c('antisense') ~ 'Antisense',
		transcript_biotype %in%c('retained intron') ~ 'Retained Intron',
		transcript_biotype %in%c('processed_transcript') ~ 'Processed Transcript',
		transcript_biotype %>%str_detect('pseudogene') ~ 'Pseudogene',
		TRUE  ~ 'NA'
))

sampdt <- orfs_dt%>%filter(sample==sample_i)
pdf('../plots/merged_summary2/sample_i')
pie%>%args
pie(
	sampdt$ccdscat%>%table,
	sampdt$ccdscat%>%table%>%names,
	col=c('red','darkyellow','green','blue')
	)

dev.off()

#
orfs_dt[1,]%>%t
orfs_dt$ORF_category_Gen%>%table
orfs_dt$compatible_biotype%>%table
orfs_dt$ORF_category_Tx%>%table
orfs_dt$ORF_category_Tx_compatible%>%table
orfs_dt$transcript_biotype%>%table
orfs_dt$gene_biotype%>%table
orfs_dt$gene_biotype%>%table


