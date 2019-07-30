library(Biostrings)
library(rtracklayer)
library(tidyverse)
library(magrittr)
library(data.table)
library(here)
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


#get exon info
annofile<-here('pipeline/my_gencode.v24lift37.annotation.gtf')

if(!is(exons,'GRangesList')){
	exonsgr<-mymemoise(read_compressed_gfile)(annofile,'exon')
	names(exonsgr)<-exonsgr$transcript_id
	exonsdb<-makeTxDbFromGRanges(exonsgr)
	exons<-exonsBy(exonsdb)
}
transcripts<-mymemoise(read_compressed_gfile)(annofile,'transcript')


sampleparams <- fread(here('pipeline/sample_parameter.csv'))
allsatannfiles <- Sys.glob(here('pipeline/SaTAnn/*/*final_SaTAnn_res*'))%T>%{stopifnot(file.exists(.))}


#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel(here('ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx'),sheet=2)%>%
	set_colnames(c("peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid
#

#Also load the cryptic peptides form OD5P
cryptic_peps <- here('ext_data/Cryptic_Peptides_Type_I.txt')%>%fread(header=FALSE)
cryptic_peps%<>%set_colnames(c('peptide','ORF_id_tr'))
cryptic_peps%>%head

#Merge both tables
od5ppeptidetbl <- bind_rows(lncRNA_table%>%mutate(type=`Gene Type`),cryptic_peps%>%mutate(type='cryptic'))
od5ppeptidetbl$cell_line<-'OD5P'
od5ppeptidetbl%<>%select(peptide=peptide,gene_id=gene_id,gene_name=`Gene Name`,everything())




###Also read in the peptides from other cell lines
nonod5p_lncfiles <- Sys.glob(here('ext_data/*_lncRNAs.xlsx'))
nonod5p_lncfiles %<>% setNames(nonod5p_lncfiles%>%basename%>%str_extract('[^_]+'))
nonod5p_lnc_peptides <- nonod5p_lncfiles%>%map(~readxl::read_excel(.))%>%bind_rows(.id='cell_line')
nonod5p_lnc_peptides%<>%select(peptide=Peptide_Sequence,gene_id=Gene_ID,gene_name=Gene_name,everything())


nonod5p_lnc_peptides$cell_line%<>%str_replace('0','O')
nonod5p_lnc_peptides$type='lncRNA'
nonod5p_lnc_peptides$cell_line%<>%recode('OMM745'='OMM')


allpeptides <- bind_rows(od5ppeptidetbl,nonod5p_lnc_peptides)

library(readxl)

peptides <- c(
	lncRNA_table$"MS-identified peptide",
	cryptic_peps[[1]]
	)

#We will first search all of our ORFs for peptides of interest.
cell_lines<- allpeptides$cell_line%>%unique


satannfiles <- lapply(cell_lines,function(cell_line){
	Sys.glob(here(str_interp('pipeline/SaTAnn/*${cell_line}*/*final_SaTAnn_res*')))%T>%{stopifnot(file.exists(.))}%>%setNames(.,.)
})%>%flatten_chr






# #load all our satann data
# satannorfs <- 
# 	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
# 	satannfiles%>%
# 	setNames(.,basename(dirname(.)))%>%
# 	mclapply(load_objs)

# satannorfs%<>%map('SaTAnn_results')

#sample orf df

# sample_orf_dt<-satannorfs%>%imap(~ data_frame(sample=.y,ORF_id_tr=.x$ORFs_tx$ORF_id_tr,gene_id=.x$ORFs_tx$gene_id))%>%bind_rows

satannfastas <- lapply(cell_lines,function(cell_line){
	Sys.glob(here(str_interp('pipeline/SaTAnn/*${cell_line}*/*fasta')))%T>%{stopifnot(file.exists(.))}%>%setNames(.,.)
})%>%setNames(cell_lines)%>%enframe('cell_line','file')%>%unnest

#Now get the cell line names in line (bad files names necessitate this annoying kludge)
satannfastas$cell_line%<>%recode(OMM='OMM475')
allpeptides$cell_line%<>%recode(OMM='OMM475')
stopifnot(satannfastas$cell_line %in% allpeptides$cell_line)

#get info on orfs from the sample files
sample_orf_dt <- satannfastas$file%>%split(satannfastas$cell_line)%>%map(function(files){
	files<-files%>%setNames(files%>%dirname)
	mclapply(files,function(file){
		fread(str_interp('grep -e ">" ${file}'),header=F)%>%
		set_colnames(c('ORF_id_tr','gene_type','gene_id','ORF_category_Gen','ORF_category_Tx_compatible'))%>%
		mutate(ORF_id_tr=str_replace(ORF_id_tr,'>',''))%>%mutate()
	})%>%bind_rows(.id='sample')
})%>%bind_rows(.id='cell_line')
sample_orf_dt%<>%as_tibble

cell_line_prots <- satannfastas$file%>%split(satannfastas$cell_line)%>%lapply(.%>%readAAStringSet%>%unique)


# for(cell_line in cell_lines){

		# .[names(.)%>%str_detect('OD5P')]%>%
	

#now search all our ORFs for each peptide
#get a df with peptide/start/end/width/ORF_id_tr


cell_line_hits <- mclapply(unique(satannfastas$cell_line),mymemoise(function(cell_line){
	cell_line_protseqs <- cell_line_prots[[cell_line]]


	peptides <- allpeptides%>%filter(cell_line==cell_line)%>%.$peptide

	peptide_hits_df <- map(peptides,function(peptide) vmatchPattern(peptide,cell_line_protseqs)%>%unlist%>%as.data.frame)%>%
		setNames(peptides)%>%
		bind_rows(.id='peptide')
	peptide_hits_df$names%<>%str_extract('[^|]+')
	peptide_hits_df%<>%rename('names'='ORF_id_tr')
	#offset the start and end to transcript coordinates
	peptide_hits_df%<>%mutate(orfstart=as.numeric(str_extract(ORF_id_tr,regex('(?<=_)\\d+'))))
	peptide_hits_df%<>%mutate(orfend=as.numeric(str_extract(ORF_id_tr,regex('\\d+$'))))
	peptide_hits_df%<>%mutate(start=((start-1)*3) +1 + orfstart -1 )
	peptide_hits_df%<>%mutate(end=((end)*3) + orfstart -1 )
	peptide_hits_df%<>%dplyr::select(-matches(c('width')))
	peptide_hits_df$cell_line=cell_line

	message('.')
	peptide_hits_df
}))


cell_line_hits%<>%setNames(unique(satannfastas$cell_line))

cell_line_hits%<>%bind_rows(.id='cell_line')

#cell_line_hits%>%write_tsv(here('tables/peptide_hits_alllines.tsv'))







#now get which samples our ORF is in.55
cell_line_hits%<>%.[c('cell_line','peptide','start','end','ORF_id_tr','orfstart','orfend')]
cell_line_hits%<>%as_tibble

cell_line_hits%<>%left_join(sample_orf_dt%>%distinct,by=c('cell_line','ORF_id_tr'))
cell_line_hits$sample%<>%basename
#ensure each peptide matches the listed gene
g2t_df <- data_frame(gene_id=transcripts$gene_id,transcript_id=transcripts$transcript_id)%>%distinct
cell_line_hits <- cell_line_hits%>%mutate(transcript_id=str_extract(ORF_id_tr,'^[^_]+'))%>%left_join(g2t_df)

cell_line_hits

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
ORF_stats_df%<>%as_tibble
ORF_stats_df%>%head(2)
##Also create an extra info sheet including a) the expression level and b) the 

cell_line_hits<-cell_line_hits%>%
	safe_left_join(ORF_stats_df%>%select(sample,ORF_id_tr,normalized_periodic_expr,pval,pval_uniq),by=c('sample','ORF_id_tr'))%>%
	select(-orfstart,-orfend)

# peptide_info_df%>%group_by(peptide,ORF_id_tr)%>%nest%>%filter(n()>1)%>%slice(1:2)%>%unnest
# peptide_info_df[1,]%>%t


peptideinfofile<-here('tables/peptide_info_df.tsv')
cell_line_hits%>%write_tsv(peptideinfofile)
normalizePath(peptideinfofile)

cell_line_hits%>%filter(cell_line!='OD5P')%>%write_tsv(here('tables/peptide_info_df_newlines.tsv'))



# peptide_info_df%>%.$gene_id%>%n_distinct
# peptide_info_df%>%filter(peptide %in% lncRNA_table[[1]])%>%.$gene_id%>%n_distinct
# peptide_info_df%>%filter(peptide %in% lncRNA_table[[1]])%>%.$peptide%>%n_distinct

# 11/55



# lncRNA_table$peptide = lncRNA_table[[1]]
# lncRNA_table%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%.$found%>%mean

# lncRNA_table%>%filter(`Gene Type`=='lncRNA')%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%group_by(gene_id)%>%summarise(found=any(found))%>%.$found%>%mean

