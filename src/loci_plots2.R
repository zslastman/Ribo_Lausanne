
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


options(ucscChromosomeNames=F)



#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=2)%>%
	set_colnames(c("MS-identified peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid
#
cellnames <- allbigwigs%>%dirname%>%basename%>%str_extract('[^_]+')%>%unique
#this is ugly but the samples are badly named, no time to fix
cellnames%<>%grep(v=TRUE,inv=TRUE,patt='OMM475')


#get exon info
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)
transcripts<-read_compressed_gfile(annofile,'transcript')

#We will first search all of our ORFs for peptides of interest.
cell_lines<-c('OD5P','OMM','ONVC') 
satannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
#load all our satann data
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

#sample orf df
sample_orf_dt<-satannorfs%>%imap(~ data_frame(sample=.y,ORF_id_tr=.x$ORFs_tx$ORF_id_tr,gene_id=gene_id))%>%bind_rows
[satannorfs[[1]][['ORFs_tx']]
allprotseqs <- satannorfs%>%
	# .[names(.)%>%str_detect('OD5P')]%>%
	map(.%>%.$ORFs_tx%>%{setNames(.$Protein,.$ORF_id_tr)})%>%Reduce(f=c)%>%unique

#

peptides <- lncRNA_table$"MS-identified peptide"

#now search all our ORFs for each peptide
#get a df with peptide/start/end/width/ORF_id_tr
peptide_hits_df <- map(peptides,function(peptide) vmatchPattern(peptide,allprotseqs)%>%unlist%>%as.data.frame)%>%
	setNames(peptides)%>%
	bind_rows(.id='peptide')
peptide_hits_df%<>%rename('names'='ORF_id_tr')
#offset the start and end to transcript coordinates
peptide_hits_df%<>%mutate(orfstart=as.numeric(str_extract(ORF_id_tr,regex('(?<=_)\\d+'))))
peptide_hits_df%<>%mutate(orfend=as.numeric(str_extract(ORF_id_tr,regex('\\d+$'))))
peptide_hits_df%<>%mutate(start=((start-1)*3) +1 + orfstart -1 )
peptide_hits_df%<>%mutate(end=((end)*3) + orfstart -1 )
peptide_hits_df%<>%select(-matches(c('width')))

#now get which samples our ORF is in.
peptide_hits_df%<>%left_join(sample_orf_dt)

peptide_hits_df%>%group_by(peptide,start,end,ORF_id_tr)%>%nest(.key='samples')
peptide_hits_df%>%group_by(peptide,start,end,ORF_id_tr)%>%nest(.key='samples')

#ensure each peptide matches the listed gene
g2t_df <- data_frame(gene_id=transcripts$gene_id,transcript_id=transcripts$transcript_id)%>%distinct
peptide_hits_df <- peptide_hits_df%>%mutate(transcript_id=str_extract(ORF_id_tr,'^[^_]+'))%>%left_join(g2t_df)


genome='mine'
genomefile='pipeline/my_hg38.fa'


#######First get list of wigs
# allbigwigs<-Sys.glob('pipeline/bigwigs/*OD5P*/*/*chr.bw')
allbigwigs<- Sys.glob('../pipeline/riboqc/data/*/_P_sites_*.bw')%>%grep(v=T,inv=T,pat='uniq')
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


