library(Biostrings)
library(rtracklayer)
library(tidyverse)
library(magrittr)
library(data.table)
library(GenomicFeatures)
source('/fast/groups/ag_ohler/work/dharnet_m/Ribo_Lausanne/src/functions.R')
select<-dplyr::select

# datafiles <- allbigwigs

# make_track<-function(file, trackname, sizefactor, vtr,isneg=FALSE,...){

#   selection =   GRanges(vtr$transcript_id,IRanges(1,sum(width(vtrexons))))

#   if(str_detect(trackname,'total')) isneg=TRUE
#   if(isneg){ 
#     importfunction = import.bw.neg
#     # ylims = c(-30,0)
#   }else{
#     importfunction= function(file,selection) import.bw(file,sel=selection)
#     # ylims = c(0,10)
#   }
#   DataTrack(file,name = trackname,chromosome=vtr$transcript_id,stream=TRUE,importFunction = importfunction)

# }

# ribotracks <- datafiles%>%
#   data_frame(file=.)%>%
#   mutate(trackname =extract_id(file,sampleids),
#   	strand=ifelse(str_detect(file,'neg'),'-','+')
#   )%>%
#   split(.,.$strand)%>%
#   map(.%>%{map2(.$file,.$trackname,.f = make_track)}%>%setNames(.,.$trackname))


options(ucscChromosomeNames=FALSE)



#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=2)%>%
	set_colnames(c("MS-identified peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid
#

lncRNA_table$gene_id%>%n_distinct

#get exon info
annofile<-'./my_gencode.v24lift37.annotation.gtf'
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)
transcripts<-read_compressed_gfile(annofile,'transcript')

#We will first search all of our ORFs for peptides of interest.
cell_lines<-c('OD5P','OMM','ONVC') 
satannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
#load all our satann data
if(!exists('satannorfs')) { 
	satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)
}
#sample orf df
sample_orf_dt<-satannorfs%>%imap(~ data_frame(sample=.y,ORF_id_tr=.x$ORFs_tx$ORF_id_tr,gene_id=.x$ORFs_tx$gene_id))%>%bind_rows
llprotseqs <- satannorfs%>%
	# .[names(.)%>%str_detect('OD5P')]%>%
	map(.%>%.$ORFs_tx%>%{setNames(.$Protein,.$ORF_id_tr)})%>%Reduce(f=c)%>%unique

#
cryptic_peps <- '/fast/groups/ag_ohler/work/dharnet_m/Ribo_Lausanne/ext_data/Cryptic_Peptides_Type_I.txt'%>%fread(header=FALSE)
cryptic_peps%<>%set_colnames(c('MS-identified peptide','ORF_id_tr'))

peptides <- c(lncRNA_table$"MS-identified peptide",cryptic_peps[[1]])



#now search all our ORFs for each peptide
#get a df with peptide/start/end/width/ORF_id_tr
allprotseqs <- satannorfs%>%
	.[names(.)%>%str_detect('OD5P')]%>%
	map(.%>%.$ORFs_tx%>%{setNames(.$Protein,.$ORF_id_tr)})%>%Reduce(f=c)%>%unique


peptide_hits_df <- map(peptides,function(peptide) vmatchPattern(peptide,allprotseqs)%>%unlist%>%as.data.frame)%>%
	setNames(peptides)%>%
	bind_rows(.id='peptide')
peptide_hits_df%<>%rename('names'='ORF_id_tr')
#offset the start and end to transcript coordinates
peptide_hits_df%<>%mutate(orfstart=as.numeric(str_extract(ORF_id_tr,regex('(?<=_)\\d+'))))
peptide_hits_df%<>%mutate(orfend=as.numeric(str_extract(ORF_id_tr,regex('\\d+$'))))
peptide_hits_df%<>%mutate(start=((start-1)*3) +1 + orfstart -1 )
peptide_hits_df%<>%mutate(end=((end)*3) + orfstart -1 )
peptide_hits_df%<>%dplyr::select(-matches(c('width')))

stop()

peptide_hits_df <- map(peptides,function(peptide) vmatchPattern(peptide,allprotseqs)%>%unlist%>%as.data.frame)%>%
	setNames(peptides)%>%
	bind_rows(.id='peptide')
peptide_hits_df%<>%rename('names'='ORF_id_tr')
#offset the start and end to transcript coordinates
peptide_hits_df%<>%mutate(orfstart=as.numeric(str_extract(ORF_id_tr,regex('(?<=_)\\d+'))))
peptide_hits_df%<>%mutate(orfend=as.numeric(str_extract(ORF_id_tr,regex('\\d+$'))))
peptide_hits_df%<>%mutate(start=((start-1)*3) +1 + orfstart -1 )
peptide_hits_df%<>%mutate(end=((end)*3) + orfstart -1 )
peptide_hits_df%<>%dplyr::select(-matches(c('width')))

stop()













#now get which samples our ORF is in.55
peptide_hits_df%<>%left_join(sample_orf_dt)

#ensure each peptide matches the listed gene
g2t_df <- data_frame(gene_id=transcripts$gene_id,transcript_id=transcripts$transcript_id)%>%distinct
peptide_hits_df <- peptide_hits_df%>%mutate(transcript_id=str_extract(ORF_id_tr,'^[^_]+'))%>%left_join(g2t_df)


genome='mine'
genomefile='pipeline/my_hg38.fa'


#######First get list of wigs
# allbigwigs<-Sys.glob('pipeline/bigwigs/*OD5P*/*/*chr.bw')
allbigwigs<- Sys.glob('../pipeline/riboqc/data/*/_P_sites_*.bw')%>%grep(v=TRUE,inv=TRUE,pat='uniq')
#now as a pairlist
bigwigpairlist <- allbigwigs%>%
	data_frame(file=.)%>%
	# mutate(base=file%>%basename%>%str_replace('pos|neg',''))%>%
	mutate(base=file%>%dirname%>%basename)%>%
	# mutate(strand=file%>%basename%>%str_extract('pos|neg'))%>%
	mutate(strand=file%>%basename%>%str_extract('plus|minus'))%>%
	mutate(strand = case_when(
		strand=='plus'~'+',
		strand=='minus'~'-'
	)) %>% 
	arrange(desc(strand))%>%
	{split(.,.$base)}%>%{map(.,~split(.$file,.$strand))%>%map(rev)}

stopifnot(bigwigpairlist%>%names%>%n_distinct%>%`>`(3))
stopifnot(bigwigpairlist%>%.[[1]]%>%names%>%`==`(c('+','-')))
cellnames <- allbigwigs%>%dirname%>%basename%>%str_extract('[^_]+')%>%unique
#this is ugly but the samples are badly named, no time to fix
cellnames%<>%grep(v=TRUE,inv=TRUE,patt='OMM475')


satannorfs[[1]]$ORFs_tx$ORF_category_Tx_compatible%>%table
ORF_stats_df<-satannorfs%>%.[TRUE]%>%map(.%>%.$ORFs_tx%>%mcols%>%.[c('TrP_pNpM','ORF_id_tr','pval','pval_uniq','ORF_category_Tx_compatible')]%>%as.data.frame%>%
	set_colnames(c('normalized_periodic_expr','ORF_id_tr','pval','pval_uniq','ORF_category_Tx_compatible')))%>%bind_rows(.id='sample')
ORF_stats_df%>%head(2)
##Also create an extra info sheet including a) the expression level and b) the 
peptide_info_df<-peptide_hits_df%>%
	left_join(ORF_stats_df,by=c('sample','ORF_id_tr'))%>%
	select(-orfstart,-orfend)

# peptide_info_df%>%group_by(peptide,ORF_id_tr)%>%nest%>%filter(n()>1)%>%slice(1:2)%>%unnest
# peptide_info_df[1,]%>%t


peptideinfofile<-'../tables/peptide_info_df.tsv'
peptide_info_df%>%write_tsv(peptideinfofile)
normalizePath(peptideinfofile)


peptide_info_df%>%.$gene_id%>%n_distinct
peptide_info_df%>%filter(peptide %in% lncRNA_table[[1]])%>%.$gene_id%>%n_distinct
peptide_info_df%>%filter(peptide %in% lncRNA_table[[1]])%>%.$peptide%>%n_distinct

11/55



lncRNA_table$peptide = lncRNA_table[[1]]
lncRNA_table%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%.$found%>%mean

lncRNA_table%>%filter(`Gene Type`=='lncRNA')%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%group_by(gene_id)%>%summarise(found=any(found))%>%.$found%>%mean

