library(magrittr)
library(tidyverse)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)
library(purrr)

filter<-dplyr::filter
slice<-dplyr::slice
filter<-dplyr::filter
filter<-dplyr::filter

source("/fast/groups/ag_ohler/dharnet_m/Ribo_Lausanne/src/functions.R")

#load arguments
args <- c(
	vcffile = '../ext_data/vcfs_october_2018/vcfs/0D5P_M11_240517/M11_mq_240517.vcf',
	orfgenomifile = 'SaTAnn/OD5P_05_uM_DAC_1/SaTAnn_Final_ORFs_files',
	genome = 'my_hg19.fa',
	outputfile = 'modified_fastas/tmp'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

outputfolder<-dirname(outputfile)
outputfolder <- paste0(outputfolder,'/')
outputfolder%>%dir.create(showWarn=F,rec=TRUE)

genome%T>%{stopifnot(file.exists(.))}
genome = Rsamtools::FaFile(genome)

#now we need to apply this function to each of the ORFs found in relevant cell lines

#############Get our vcf data
vcfgr<-
	vcffile%>%
	fread%>%
	{colnames(.)%<>%str_replace('#','');colnames(.)[1:2]<-c('seqnames','start');.}%>%
	mutate(width=nchar(ALT),seqnames=paste0('chr',seqnames))%>%
	DT2GR(seqinf=genome%>%seqinfo)

# vcfgrs %<>% map(. %>% .[.$INFO=='SNP'])
# vcfgrs %<>% map(. %>% .[nchar(.$REF)==1])
# vcfgrs %<>% map(. %>% .[nchar(.$ALT)==1])
vcfgr %<>% setNames(.,.$ID)

vcfpullseqs <- vcfgr%>%getSeq(x=genome)%>%as.character

vcfrefannoseq <- vcfgr$REF
vcfaltannoseq <- vcfgr$ALT

#AHA - sometimes the variant is just the new genome's 
table(vcfpullseqs==vcfrefannoseq)
stopifnot(all((vcfpullseqs==vcfrefannoseq)||(vcfpullseqs==vcfaltannoseq)))

#map the vcfs onto the relevant transcripts, then work out the overlaps, then modify teh strings

# get_transcriptvcfgr<-function(x,exons){
# 	out = mapToTranscripts(x,exons)
# 	out$ALT = x$ALT[out$xHits]
# 	out
# }
# #finally map to the relevant exons
# vcfgr_tx <- vcfgrs%>%map(get_transcriptvcfgr,exons)

# vcfgr_tx%<>%map(.%>%{
# 	isneg <- as.vector(strand(.)=='-')
# 	.$ALT[isneg] <- as.character(reverseComplement(DNAStringSet(.$ALT[isneg])))
# 	.$REF[isneg] <- as.character(reverseComplement(DNAStringSet(.$REF[isneg])))
# 	.
# })


##########This function produces modified DNA sequences from gr of exons, a vcf, and a genome file


# satann_trgr<- read_compressed_gfile('groupedsatan/OD5P.gtf','sequence_feature')
# satann_trgr%>%mcols%>%.[1,]
# vcfname <-basename(vcffile)

# satandata <- load_objs(orfobjectfile)

# satandata%>%.$ORFs_gen%>%export('orfs_genomic.gtf')


# ENSP00000334393.3      rs75062661_0:T141A

genomicorffiles <- 'groupedsatan/*.genomic.gtf'%>%Sys.glob
fastafiles <- genomicorffiles%>%str_replace('genomic.gtf','fasta')%>%setNames(.,genomicorffiles)


# vcfgr<-vcfgr['rs75062661']
# satandata$ORFs_tx[]


for(genomicorffile in genomicorffiles){

	orfsgenome <- import(genomicorffile)%>%setNames(.,.$ID)

	seqs<-import(fastafiles[[genomicorffile]])%>%head(1000)
	seqorfids<-names(seqs)%>%str_extract(regex('.*?(?=\\|)'))

	orfsgenome <- orfsgenome%>%split(.,names(.))%>%.[seqorfids]	
	#get the modified fasta headers
	modnames <- injectSNPsHeader(orfsgenome,vcfgr,genome)
	stop()
	names(seqs)%<>%paste0('|',str_split_fixed(modnames,'\\|',3)%>%.[,2])

	# dnafile <- file.path(outputfolder,paste0('cell_line','_','celllinefile.fa'))%T>%message
	# writeXStringSet(seqs,dnafile)
	protfile <- file.path(outputfolder,paste0(basename(genomicorffile),'_',basename(vcffile),'modified.fasta'))%T>%message
	writeXStringSet(seqs,protfile)

	import(protfile)%>%names
}







