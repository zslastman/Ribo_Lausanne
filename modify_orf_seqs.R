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

fafile <- 'pipeline/my_hg38.fa'%T>%{stopifnot(file.exists(.))}
mygenome = FaFile(fafile)

args <- commandArgs(trailing=T)%T>%message
satanfiles <- Sys.glob(paste0(args[1],'/*/*Final_ORFs*'))
stopifnot(length(satanfiles)>0)
# orfsgenome[1]
# vcfgr[1]
#now we need to apply this function to each of the ORFs found in relevant cell lines

#############Get our vcf data
vcfs <- Sys.glob('ext_data/*.vcf')
#may god forgive me
vcfgrs<-vcfs%>%map(.%>%fread%>%{colnames(.)%<>%str_replace('#','');colnames(.)[1:2]<-c('seqnames','start');.}%>%mutate(width=nchar(ALT),seqnames=paste0('chr',seqnames))%>%DT2GR(FaFile('../genomes/hg19.fa')%>%seqinfo))
chainurl = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
chainfile = basename(chainurl)
# system(str_interp('wget ${chainurl} -O ${chainfile}'))
# system(str_interp('gunzip ${chainfile}'))
chainfile = str_replace(chainfile,'.gz','')
chain <- import.chain(chainfile)

# vcfgrs %<>% map(. %>% .[.$INFO=='SNP'])
# vcfgrs %<>% map(. %>% .[nchar(.$REF)==1])
# vcfgrs %<>% map(. %>% .[nchar(.$ALT)==1])
vcfgrs %<>% map(.%>%setNames(.,.$ID))
vcfgrs %<>% map(liftOver,chain)
vcfgrs %<>% map(unlist)
vcfgrs	%<>%map(~.[names(.)%>%table%>%keep(~.==1)%>%names])

vcfpullseqs <- vcfgrs%>%map(unlist)%>%.[[1]]%>%getSeq(x=FaFile('../genomes/hg38.fa'))%>%as.character

vcfrefannoseq <- vcfgrs%>%map(unlist)%>%.[[1]]%>%.$REF
vcfaltannoseq <- vcfgrs%>%map(unlist)%>%.[[1]]%>%.$ALT

#AHA - sometimes the variant is just the new genome's 
# table(vcfpullseqs==vcfrefannoseq)
stopifnot(all((vcfpullseqs==vcfrefannoseq)||(vcfpullseqs==vcfaltannoseq)))

cell_lines <- vcfs%>%str_extract('(?<=/).*?(?=_)')

names(vcfgrs) <- cell_lines

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

injectSNPsIndels <- function(orfsgenome,vcfgr,genome){
	stopifnot(is(orfsgenome,'GRanges'))
	stopifnot(is(vcfgr,'GRanges'))
	
	vcfgr$ALT<-DNAStringSet(vcfgr$ALT)
	posorfs <- orfsgenome%>%subset(strand=='+')%>%split(.,names(.))
	mutsonorfs <- mapToTranscripts(vcfgr,posorfs)
	# posorfdnaseq <- posorfs[unique(seqnames(mutsonorfs))]%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
	posorfdnaseq <- posorfs%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
	mutsonorfs$ALT = vcfgr$ALT[mutsonorfs$xHits]
	msnames <- as.character(seqnames(mutsonorfs))
	#now modify the dna strings
	stopifnot(all(msnames %in% names(posorfdnaseq)))
	mutsonorfs_isindel = ! ((nchar(mutsonorfs$ALT)==1) & (width(mutsonorfs)==1))

	for(i in which(!mutsonorfs_isindel)){
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



# testmut<-vcfgr[1]%>%{GRanges('chrY',IRanges(end(orfsgenome[137073])-1,end(orfsgenome[137073])-1), ALT='C' )}
# testorf<-orfsgenome[137073]
# injectSNPsIndels(testorf,testmut,genome)

cell_lines_short <- cell_lines%>%
	map(str_extract,c('OD5P','OMM','ONVC')) %>%
	map_chr(keep,Negate(is.na))

names(vcfgrs) <- cell_lines_short

satancell_lines <-satanfiles%>%
	map(str_extract,c('OD5P','OMM','ONVC'))%>%
	map_chr(keep,Negate(is.na))

modified_seq_dir <- 'Modified_seqs'%T>%dir.create(showWarn=T)
for(cell_line in cell_lines_short){
	vcfgr <- vcfgrs[[cell_line]]
	celllinefiles <- satanfiles[satancell_lines==cell_line]
	for(celllinefile in celllinefiles){
		satandata <- load_objs(celllinefile)$ORFs_gen
		snpinject <- injectSNPsIndels(satandata,vcfgr,mygenome)
		seqs <- snpinject$seqs

		dnafile <- file.path(modified_seq_dir,paste0('cell_line','_','celllinefile.fa'))%T>%message
		writeXStringSet(seqs,dnafile)
		protfile <- file.path(modified_seq_dir,paste0('cell_line','_','celllinefile.prot.fa'))%T>%message
		writeXStringSet(translate(seqs),protfile)

	}
}
