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

source("/fast_new/work/groups/ag_ohler/dharnet_m/Ribo_Lausanne/src/functions.R")

#load arguments
args <- c(
	genomicorffile = 'groupedsatan/OD5P.genomic.gtf',
	genome = 'my_hg19.fa',
	vcffiles = '',
	outputfile = 'modified_fastas/tmp'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

# save.image()
# stop(getwd())

outputfolder<-dirname(outputfile)
outputfolder <- paste0(outputfolder,'/')
outputfolder%>%dir.create(showWarn=F,rec=TRUE)

genome = Rsamtools::FaFile(genome)

fastafile <- genomicorffile%>%str_replace('genomic.gtf','fasta')


#now we need to apply this function to each of the ORFs found in relevant cell lines

#############Get our vcf data
vcffiles<-str_split(vcffiles,',')%>%unlist
vcfgr<-
	vcffiles%>%
	lapply(.%>%
	fread%>%
	{colnames(.)%<>%str_replace('#','');colnames(.)[1:2]<-c('seqnames','start');.}%>%
	mutate(width=nchar(ALT),seqnames=paste0('chr',seqnames))
	)%>%
	bind_rows%>%
	distinct(seqnames,start,width,ALT,ID,.keep_all=TRUE)%>%
	DT2GR(seqinf=genome%>%seqinfo)

vcfgr %<>% setNames(.,.$ID)

names(vcfgr) <- names(vcfgr)%>%data_frame(a=.)%>%group_by(a)%>%transmute(b=paste0(a,'_',seq_along(a)))%>%.$b
# vcfgrs %<>% map(. %>% .[.$INFO=='SNP'])
# vcfgrs %<>% map(. %>% .[nchar(.$REF)==1])
# vcfgrs %<>% map(. %>% .[nchar(.$ALT)==1])


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

# genomicorffiles <- 'groupedsatan/*.genomic.gtf'%>%Sys.glob
# fastafile <- genomicorffiles%>%str_replace('genomic.gtf','fasta')%>%setNames(.,genomicorffiles)
# stopifnot(file.exi)

# vcfgr<-vcfgr['rs75062661']
# satandata$ORFs_tx[]

vcfgr$indel <- !(width(vcfgr)==1 & (nchar(vcfgr$ALT)==1) )

# for(genomicorffile in genomicorffiles){
orfsgenome <- rtracklayer::import(genomicorffile)%>%setNames(.,.$ID)

#read in our protein sequences as AAStringset objects
seqs<-	Biostrings::readAAStringSet(fastafile)

#get the names of these sequences
seqorfids<-names(seqs)%>%str_extract(regex('.*?(?=\\|)'))

#orfgrs in same order as our sequences
stopifnot(setequal(seqorfids,names(orfsgenome)))
orfsgenome <- orfsgenome%>%split(.,names(.))%>%.[seqorfids]	

#get the modified fasta headers
message('getting snp headers')
snpheaders<-injectSNPsHeader(orfsgenome,vcfgr%>%subset(!indel),genome)
#getting indel headers
message('getting indel headers')
indelheaders<-injectIndelsHeader(orfsgenome,vcfgr%>%subset(indel),genome,exons)

#now merge these headers
modnames <- full_join(
  snpheaders%>%str_split_fixed('\\|',3)%>%.[,-3,drop=F]%>%as.data.frame,
  indelheaders%>%str_split_fixed('\\|',3)%>%.[,-3,drop=F]%>%as.data.frame,
  by='V1'
)%>%
{.[match(seqorfids,.$V1),]}%>%
{paste0(.$V1,'|',.$V2.x,',',.$V2.y)}%>%str_replace('\\|,$','|')


#and name the sequence object
names(seqs) %<>%paste0('|',str_split_fixed(modnames,'\\|',3)%>%.[,2])
names(seqs) %<>%str_replace(',$','')
#finally, export our new protein fasta with it's modified
message('exporting')

vcfnames <- vcffiles%>%map(basename)%>%map(tools::file_path_sans_ext)%>%paste0(collapse='_')

# protfile <- file.path(outputfolder,paste0(basename(genomicorffile),'_',vcfnames,'modified.fasta'))%T>%message
protfile <- outputfile%T>%message

writeXStringSet(seqs,protfile)
# import(protfile)%>%names
# }