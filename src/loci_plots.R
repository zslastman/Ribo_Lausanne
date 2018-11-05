library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
source('../src/functions.R')


#load some utility functions

#load excel data, parse
stopifnot(c("ere_Peptides", "lncRNA_peptides") == readxl::excel_sheets('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx'))




#read in the data on the repeat associated peptides
ere_table <- readxl::read_excel('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=1)
colnames(ere_table) <- c("MS_identified_peptide", "ID", "seqnames", "start", "end", "repName", 
"strand", "class.fam", "unMerged", "FPKM RNA seq OD5PCTRL", "FPKM RNA seq OD5PDAC")


#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=2)%>%
	set_colnames(c("MS-identified peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid



#import genome as fasta file object
genome <- ('../../genomes/hg38.fa'%T>%{stopifnot(file.exists(.)) })%>%FaFile
annofile <- '../pipeline/my_gencode.v24lift37.annotation.gtf'%T>%{stopifnot(file.exists(.)) }

# transcriptsother<-read_compressed_gfile('annotation/gencode.v21.annotation.gtf','transcript')

ere_gr <- DT2GR(ere_table,seqinf=seqinfo(genome))



genes<-read_compressed_gfile(annofile,'gene')
start_codons<-read_compressed_gfile(annofile,'start_codon')
transcripts<-read_compressed_gfile(annofile,'transcript')
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)


#now get the relevant genes for our lincRNAs
stopifnot(all(lncRNA_table$gene_id %in% genes$gene_id))
lincgenes <- subset(genes, gene_id %in% lncRNA_table$gene_id)





# allbigwigs<-Sys.glob('pipeline/bigwigs/*OD5P*/*/*chr.bw')
allbigwigs<- Sys.glob('../pipeline/riboqc/data/*/_P_sites_*.bw')%>%grep(v=T,inv=T,pat='uniq')


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


#############Now individual Plots

#TODO this logic can be split up a bit to work on bam files say
# stopmut_transcripts <- fread('stopmut_transcripts.txt')[[1]]
stopmut_transcripts <- NULL

linctranscripts <- transcripts %>% subset(gene_id %in% lincgenes$gene_id)

ranges2plot <- c(
	ere_gr%>%setNames(.,.$ID)%>%.[,NULL]%>%split(.,names(.)),
	exons[linctranscripts$transcript_id]%>%unlist%>%.[,NULL]%>%split(.,names(.)),
	exons[stopmut_transcripts]
)

ribofilenames <- names(bigwigpairlist)


extract_oneof <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = ids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}
sampledf <- fread('../pipeline/sample_parameter.csv')
sampleids <- sampledf %>%filter(assay=='ribo')%>%.[[1]]
extract_id <- purrr::partial(extract_oneof,ids=sampleids)
stopifnot(all(!is.na(extract_id(ribofilenames))))

sampleids %>%setdiff( ribofilenames)

# rangename='ENST00000374922'
# ribofilename=ribofilenames[[1]]


# for (ribofilename in ribofilenames[1]){

# 	riboprofdata <- get_riboprofdata(ranges2plot['ENST00000374922.7'],bigwigpairlist[[ribofilename]])
# 	ribofileid <- extract_id(ribofilename,sampleids)


# 	riboprofdata%>%group_by(tid)%>%summarise(sum(score))

# 	# for (rangename in names(ranges2plot)){
# 	for (rangename in 'ENST00000374922.7'){



# 		range = ranges2plot[[rangename]]
# 		rangename%<>%str_replace_all('[\\[\\]]','')

# 		plottitle <- str_interp("Riboseq Read Profile for:\n${rangename} = ${range}\nfile:${ribofileid}")		

# 		rangeriboprofdata <- riboprofdata %>%dplyr::filter(tid == rangename)
# 		if(nrow(rangeriboprofdata)<4){
# 			message(str_interp('skipping ${plottitle}'))
# 			next
# 		}
# 		transprofplot <-
# 			rangeriboprofdata %>%
# 			ggplot(aes(fill=frame,color=frame,x=pos,y=score))+
# 			geom_bar(stat='identity')+
# 			coord_cartesian(xlim=c(0,sum(width(ranges2plot[rangename]))))+
# 			ggtitle(plottitle)+
# 			theme_bw()
        
# 		rangeseq <- range%>%getSeq(x=FaFile('../pipeline/my_hg38.fa'))%>%Reduce(f=c)
# 		startpos <- rangeseq%>%str_locate('ATG')%>%.[,1]%>%{data_frame(pos=.)}%>%mutate(type='start_codon',frame = pos %%3)
# 		stoppos <- rangeseq%>%str_locate('TAG|TAA|TGA')%>%.[,1]%>%{data_frame(pos=.)}%>%mutate(type='start_codon',frame = pos %%3)
		


# 		#name the plot file
# 		transprofplotfile <- str_interp('plots/loci_riboprofiles/${rangename}/${ribofileid}/${rangename}_${ribofileid}_prof.svg')
# #		transprofplotfile%<>%str_replace_all('[^a-zA-z0-9_\\./]+','')
# 		transprofplotfile%>%dirname%>%dir.create(rec=TRUE,showWarnings=FALSE)
# 		#print the plot
		
# 		svglite(transprofplotfile);multiplot(transprofplot,layout=mylayout);dev.off()

# 		message(transprofplotfile%>%normalizePath)
# 	}
# }


#initiaal form when testing
#let's look at the scores for our
cell_lines<-c('OD5P','OMM','ONVC') 
satannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

# all_linc_orfs <- satannorfs%>%map(.%>%.$ORFs_tx%>%subset(clipid(seqnames) %in% clipid(linctranscripts$transcript_id)))
all_linc_orfs <- satannorfs%>%map(.%>%.$ORFs_tx%>%subset((seqnames) %in% (linctranscripts$transcript_id)))

#all_linc_orfs

#TODO_ what are these Granges columsn???
#now filter out the weird GRanges columns, for now, and aggregate into one data table 
for (i in seq_along(all_linc_orfs))all_linc_orfs[[i]]$Protein%<>%as.character
linc_orfs_dt <- all_linc_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')


# linc_orfs_dt%<>%group_by(sample)%>%mutate(fdr =  p.adjust(pval, method='fdr'))
linc_orfs_dt %<>% dplyr::filter(pval<0.05)%>%ungroup

linc_orfs_dt%>%nrow

#TODO - eveyrthing in linc

linc_orfs_dt_with_linctable <- linc_orfs_dt%>%inner_join(lncRNA_table,by='gene_id')
linc_orfs_dt_with_linctable%<>%mutate(contains_peptide = str_detect(Protein,`MS-identified peptide`))

best_linc_orfs<-linc_orfs_dt_with_linctable%>%
	group_by(sample,gene_id)%>%
	arrange(desc(contains_peptide),pval)%>%
	dplyr::slice(1)


OD5Pbest_linc_orfs<-best_linc_orfs
OD5Pbest_linc_orfs<-best_linc_orfs%>%subset(sample%>%str_detect('OD5P'))

ranges_with_goodorf <- ranges2plot %>% .[
	match(OD5Pbest_linc_orfs$transcript_id%>%unique,
		ranges2plot%>%names)
]


# orfs2plot <- best_linc_orfs[match(names(ranges_with_goodorf),best_linc_orfs$transcript_id),]
#for our grange in transcript space
# transcript_seqinfo <- ranges_with_goodorf%>%width%>%sum%>%{Seqinfo(seqnames=names(.),seqlength=.)}
# orftrack<-OD5Pbest_linc_orfs%>%DT2GR(seqinf=transcript_seqinfo)%>%.[,NULL]%>% 
# 	# c(GRanges(c('ENST00000433514:1-20:+','ENST00000433514:2-21:+','ENST00000433514:3-22:+')))%>%
# 	unique%>%
# 	{AnnotationTrack(name='ORFs',.,feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue')}

transcript_seqinfo <- exons[unique(linc_orfs_dt$transcript_id)]%>%width%>%sum%>%{Seqinfo(seqnames=names(.),seqlength=.)}

# orftrack<-linc_orfs_dt_with_linctable%>%DT2GR(seqinf=transcript_seqinfo)%>%.[,NULL]%>% 
# 	unique%>%
# 	{Gviz::AnnotationTrack(
# 		name='ORFs',.,
# 		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
# 		id=rep('foo',length(.))}

#fraction of our lincs that have orfs in them
linctranscripts$gene_id%>%clipid%>%unique%>%is_in(best_linc_orfs$gene_id)%>%{list(sum(.),mean(.))}

#let's check up on if the annotaiton diff matters
linctranscripts$gene_id%>%clipid%>%unique%>%is_in(best_linc_orfs%>%dplyr::filter(contains_peptide)%>%.$gene_id)%>%{list(sum(.),mean(.))}

#fraction of our lincs that have orfs in them
linctranscripts$gene_id%>%clipid%>%unique%>%is_in(OD5Pbest_linc_orfs$gene_id)%>%{list(sum(.),mean(.))}

#let's check up on if the annotaiton diff matters
linctranscripts$gene_id%>%clipid%>%unique%>%is_in(OD5Pbest_linc_orfs%>%dplyr::filter(contains_peptide)%>%.$gene_id)%>%{list(sum(.),mean(.))}


#This code will insert single nucleotide variants into our 


# testexon <- mappedback%>%resize(6)%>%.[,NULL]%>%setNames('test')
# testlocation <- GRanges('test:1-3:+')
# mapFromTranscripts(testlocation,testexon)

# mappedback <- OD5Pbest_linc_orfs%>%DT2GR(seqinf=transcript_seqinfo)%>%invertStrand%>%mapFromTranscripts(ranges_with_goodorf)
# mappedback%>%resize(3,'start')%>%{Rsamtools::getSeq(x=FaFile('../pipeline/my_hg38.fa'),.)}%>%Reduce(f='c')
# mappedback%>%{Rsamtools::getSeq(x=FaFile('../pipeline/my_hg38.fa'),.)}%>%Reduce(f='c')
# mappedback%>%mapToTranscripts(ranges_with_goodorf)

# GRanges('ENST00000433514:1-10')%>%mapFromTranscripts(ranges_with_goodorf)

# #

# test<-OD5Pbest_linc_orfs[1,]%>%t
# testtr <- test['seqnames',]
# linc_orfs_dt_with_linctable%>%subset(seqnames==testtr)%>%.[c('start','end')]




#now modify the peptides with our vcf data


# allsatannorfs <- 
# 	Sys.glob('SaTAnn/*Final_ORFs*')%>%
# 	setNames(.,extract_id(.))%>%
# 	mclapply(load_objs)

# #now ust the tx granges
# all_orfs <- allsatannorfs%>% map(. %>%.$ORFs_tx)

# #TODO_ what are these Granges columsn???
# #now filter out the weird GRanges columns, for now, and aggregate into one data table 
# for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character
# orfs_dt <- all_orfs%>%map(.%>%{
# 	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
# 	.[,issimplelist]%>%GR2DT
# })%>%bind_rows(.id='sample')

# orfs_dt%<>%group_by(sample)%>%mutate(fdr =  p.adjust(pval, method='fdr'))

# #orfs_dt %<>% dplyr::filter(fdr<0.05)
# orfs_dt %<>% dplyr::filter(pval<0.05)

# trstrands<-allsatannorfs%>%lapply(.%>%.$ORFs_gen%>%{data_frame(ORF_id_tr=names(.),transcript_strand=as.character(strand(.)))})%>%
# 	bind_rows(.id='sample')%>%
# 	distinct(ORF_id_tr,transcript_strand)

# stopifnot(all(trstrands$ORF_id_tr %in% orfs_dt$ORF_id_tr))
# stopifnot(all(orfs_dt$ORF_id_tr %in% trstrands$ORF_id_tr))

# orfs_dt%<>%left_join(trstrands)

# orfs_dt%<>%mutate(strand=transcript_strand)

# orfs_gr <- orfs_dt %>%
# 	mutate(seqnames = clipid(seqnames))%>%
# 	DT2GR(seqinf=transcript_seqinfo[unique(.$seqnames)])

