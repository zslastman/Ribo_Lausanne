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

source("/fast/groups/ag_ohler/dharnet_m/Ribo_Lausanne/functions.R")

#load arguments
args <- c(
	vcffile = '../ext_data/vcfs_october_2018/vcfs/0D5P_M11_240517/M11_mq_240517.vcf',
	orfobjectfile = 'SaTAnn/OD5P_05_uM_DAC_1/SaTAnn_Final_ORFs_files',
	genome = 'my_hg19.fa',
	outputfile = ''
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
orfsgenome<-satandata[1:3000]
injectSNPsIndels <- function(orfsgenome,vcfgr,genome,snpstringsep=','){
	
	rs2187121

	stopifnot(is(orfsgenome,'GRanges'))
	stopifnot(is(vcfgr,'GRanges'))
	vcfgr_snps <- vcfgr[nchar(vcfgr$ALT)==1]
	vcfgr$ALT<-DNAStringSet(vcfgr$ALT)
	orfsgenome <- orfsgenome%>%split(.,names(.))
	allorfnames <- names(orfsgenome)
	mutsonorfs <- mapToTranscripts(vcfgr,orfsgenome)
	orfsmutted <- orfsgenome[mutsonorfs$transcriptsHits]
	# posorfdnaseq <- orfsgenome[unique(seqnames(mutsonorfs))]%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
	posorfdnaseq <- orfsmutted%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
	mutsonorfs$ALT = vcfgr$ALT[mutsonorfs$xHits]
	msnames <- as.character(seqnames(mutsonorfs))
	#now modify the dna strings
	stopifnot(all(msnames %in% names(posorfdnaseq)))
#	nonindelinds<-which(!mutsonorfs_isindel)
	mutlocs <- start(mutsonorfs)
	startaainds <- ceiling(mutlocs/3)
	startbpinds <- ((startaainds-1)*3)+1
	codseqs<-subseq(DNAStringSet(posorfdnaseq[msnames]),startbpinds,width=3)
	refaas <- translate(codseqs)
	subseq(codseqs,(((mutlocs)-1)%%3)+1,width=1) <- mutsonorfs$ALT
	altaas <- translate(codseqs)
	orfmodstrings <- data_frame(
		orfname=names(orfsgenome)[mutsonorfs$transcriptsHits],
		snpstring=paste0(refaas,startaainds,altaas))
	#collapse
	orfmodstrings <- orfmodstrings%>%group_by(orfname)%>%summarise(snpstring=paste0(snpstring,collapse=snpstringsep))

	data_frame(orfname=allorfnames)%>%left_join(orfmodstrings)%>%replace_na(replace=list(snpstring=''))%>%head%>%do.call(paste0,.)


	for(i in which(!mutsonorfs_isindel)){
		mutloc<-start(mutsonorfs[i])
		startaaind <- ceiling(mutloc/3)
		startbpind <- ((startaaind-1)*3)+1
		codseq <- dnaseq[startbpind:(startbpind+2)]
		refaa <- translate(codseq)
		codseq[(((mutloc)-1)%%3)+1] <- mutsonorfs$ALT[[i]]
		altaa <- translate(codseq)
		paste0(refaa,startaaind,altaa)

		names(posorfdnaseq)[msnames[[i]]]

		aaref <- translate()

		dnaseq<-posorfdnaseq[[msnames[i]]]
		aamut <- translate(dnaseq[startbpind:endbpind])


		protseq
		#codonends ceiling(1:7/3)*3
		ceiling(start(mutsonorfs)[i]/3)*3+1
		floor(start(mutsonorfs)[i]/3)+1
		mutcodons <- dnaseq[codstart:codend]
		mutinds <- start(mutsonorfs[i]):end(mutsonorfs[i])
		
		posorfdnaseq[[msnames[i]]][]
		posorfdnaseq[[msnames[i]]][start(mutsonorfs[i]):end(mutsonorfs[i])]<-mutsonorfs$ALT[[i]]
	}

	for(i in which(mutsonorfs_isindel)){
		seq <- posorfdnaseq[[msnames[i]]]
		posorfdnaseq[[msnames[i]]] <- c(
			seq[1:(start(mutsonorfs[i])-1)],
			mutsonorfs[i]$ALT[[1]],
			seq[(end(mutsonorfs[i])+1):length(seq)]
		)		
	}

	negorfs <- orfsgenome%>%subset(strand=='-')%>%split(.,names(.))
	mutsonorfs <- mapToTranscripts(vcfgr,negorfs)

	# negorfdnaseq <- negorfs[unique(seqnames(mutsonorfs))]%>%lapply(.%>%getSeq(x=genome)%>%rev%>%Reduce(f=c))
	negorfdnaseq <- negorfs%>%lapply(.%>%getSeq(x=genome)%>%rev%>%Reduce(f=c))
	mutsonorfs$ALT = reverseComplement(DNAStringSet(vcfgr$ALT[mutsonorfs$xHits]))

	msnames <- as.character(seqnames(mutsonorfs))
	stopifnot(all(msnames %in% names(negorfdnaseq)))
	#now modify the dna strings
	mutsonorfs_isindel = ! ((nchar(mutsonorfs$ALT)==1) & (width(mutsonorfs)==1))
	for(i in which(!mutsonorfs_isindel)) negorfdnaseq[[msnames[i]]][start(mutsonorfs[i]):end(mutsonorfs[i])]<-mutsonorfs$ALT[[i]]
	for(i in which(mutsonorfs_isindel)){
		seq <- negorfdnaseq[[msnames[i]]]
		negorfdnaseq[[msnames[i]]] <- c(
			seq[1:(start(mutsonorfs[i])-1)],
			mutsonorfs[i]$ALT[[1]],
			seq[(end(mutsonorfs[i])+1):length(seq)]
		)		
	}
	return(list(
			seqs = DNAStringSet(c(posorfdnaseq,negorfdnaseq)[unique(names(orfsgenome))])
		)
	)
}

vcfname <-basename(vcffile)

for(celllinefile in celllinefiles){
	satandata <- load_objs(orfobjectfile)$ORFs_gen
	snpinject <- injectSNPsIndels(satandata,vcfgr,mygenome)
	
	satandata%>%getSeq(x=genome)%>%

	stop()
	seqs <- snpinject$seqs

	dnafile <- file.path(outputfolder,paste0('cell_line','_','celllinefile.fa'))%T>%message
	writeXStringSet(seqs,dnafile)
	protfile <- file.path(outputfolder,paste0('cell_line','_','celllinefile.prot.fa'))%T>%message
	writeXStringSet(translate(seqs),protfile)

}





