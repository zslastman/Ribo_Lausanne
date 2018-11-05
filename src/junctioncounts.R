

#This script will calculate a set of regions, say all utrs, along with some data about those, and an annotation of which are positive and negate. It will then use the data as well as some caclulated statistics like length and gc content to try and construct a realistic background set for the positive sequences
#can probably jsut generate the sequences directly actually.
message('loading libraries')
library(magrittr)
library(data.table)
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
library(rtracklayer)
suppressMessages(library(tidyverse))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))
library(parallel)
library(BiocParallel)
register(MulticoreParam(8))

CHUNKSIZE=1e5

message('...done')

#load arguments
args <- c(
	uniquejuncrds = 'junctions/uniqueintrons.rds',
	bf = 'star/data/OD5P_ctrl_3/OD5P_ctrl_3.bam',
	outputfile='junctioncounts/OD5P_ctrl_3/OD5P_ctrl_3.junctioncounts.tsv'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

outputfolder<-dirname(outputfile)
outputfolder <- paste0(outputfolder,'/')
outputfolder%>%dir.create(showWarn=F,rec=TRUE)


uniqueintrons<-readRDS(uniquejuncrds)

# # junctionbams <- Sys.glob('star.bak/data/*/*bam')%>%
# # 	discard(~str_detect(.,'transcript'))%>%
# # 	# .[[3]]
# # 	identity
# # junctionbams[[3]]
# bamnames <- junctionbams%>%basename%>%tools::file_path_sans_ext()%>%setNames(junctionbams)
# junctionbams%<>%setNames(bamnames)

# juncgr<-uniqueintrons[,NULL]%>%resize(1)

# bamcounts <- mclapply(mc.cores=1,
# 	junctionbams[1],juncgr,FUN=function(junctionbam,juncgr){
# # bamcounts <- mclapply(junctionbams,juncgr,FUN=function(junctionbam,uniqueintrons){
# 	param <- ScanBamParam(
# 		flag=scanBamFlag(isSecondaryAlignment=FALSE),what=c("mapq"),
# 		which=juncgr,
# 		tag = "MD"
# 	)
# 	reads <- readGAlignments(BamFile(file=junctionbam, yieldSize=NA),param = param)
# 	#summary gr of the junctions and their scores
# 	readjuncs <- summarizeJunctions(reads)
# 	#match predicted juctnions to those in reads
# 	juncmatch <- match(juncgr%>%{strand(.)='*';.},readjuncs)
# 	#assign count
# 	ifelse(
# 		as.vector(strand(juncgr)=='+'),
# 		readjuncs$plus_score[juncmatch],
# 		readjuncs$minus_score[juncmatch]
# 	)
# })
# library(GenomicAlignments)


juncgr <- uniqueintrons

param <- ScanBamParam(
	flag=scanBamFlag(isSecondaryAlignment=FALSE)
)

# pbamcounts<-mclapply(junctionbams,function(bf){
bf <- BamFile(file=bf, yieldSize=CHUNKSIZE)
YIELD <- function(x, ...) readGAlignments(x, param=param)
MAP <- function(reads){
	cat('.')
	#summary gr of the junctions and their scores
	readjuncs <- summarizeJunctions(reads)
	#match predicted juctnions to those in reads
	juncmatch <- match(juncgr%>%{strand(.)='*';.},readjuncs)
	#assign count
	out <- ifelse(
		as.vector(strand(juncgr)=='+'),
		readjuncs$plus_score[juncmatch],
		readjuncs$minus_score[juncmatch]
	)
	out[is.na(out)] <- 0
	out
}
REDUCE <- `+`
DONE <- function(value) length(value) == 0L

bamcounts <- GenomicFiles::reduceByYield(bf, YIELD, MAP, REDUCE, DONE, parallel=T)

# uids <- fread('junctions/id_unique.txt')

data.frame(bamcounts) %>%
	rownames_to_column %>% 
	set_colnames(c('id_in_unique',basename(bf$path))) %>%
	# inner_join(uids, by='id_in_unique') %>%
	write_tsv(outputfile)

# })

# stopifnot(bamcounts%>%map_lgl(is.null)%>%any%>%not)
# tbamcounts%>%map_lgl(is.null)%>%which

# bamcounts->bamcountsbak

# bamcounts %<>% simplify2array
# bamcounts%<>%apply(2,function(x){x[is.na(x)]<-0;x})
# bamcounts%<>%as_data_frame
# bamcounts%<>%mutate(id_in_unique=juncgr$id_in_unique)


# stop()


# fintrons <- mergedjunctionintrons%>%subset(ID%>%str_detect('junction_diff_hits'))

# jmatch <- match(hits$id_in_unique,juncgr$id_in_unique)
# hitcounts <- bamcounts[jmatch,]

# ###First deal with our hits
# hittable<-juncfiles%>%str_subset('hits')%>%fread
# hittable$ID%<>%as.character
# hits$ID<-hits$ID%>%str_replace('.*_','')


# hittable%<>%
# 	left_join(distinct(data_frame(ID=hits$ID,id_in_unique=hits$id_in_unique)),by='ID')%>%
# 	left_join(distinct(bamcounts),by='id_in_unique')%>%
# 	select(-id_in_unique)%>%
# 	identity

# hittable$anyOD5P <- hittable%>%select(one_of(bamnames%>%str_subset('OD5P')))%>%as.matrix%>%array_branch(1)%>%map_lgl(~any(.>0))
# hittable$any <- hittable%>%select(one_of(bamnames))%>%as.matrix%>%array_branch(1)%>%map_lgl(~any(.>0))

# hittable%>%write_tsv('spliced_peptide_hits_with_od5pcounts.tsv'%T>%message)
# mergedjunctionexons%>%
# 	subset(str_detect(ID,'junction_diff_hits'))%>%
# 	export('../junctionfiles/junction_diff_hits.gtf')


# hittable<-read_tsv('spliced_peptide_hits_with_od5pcounts.tsv')


#now output counts in each 

# maxjunc <- readjuncs%>%.[which.max(.$score)]

# juncs%>%subsetByOverlaps(resize(maxjunc,width(maxjunc)+300,'center'))


# end(fpexons%>%subset(strand=='+')) - start(tpexons%>%subset(strand=='+'))
# end(fpexons%>%subset(strand=='-')) - start(tpexons%>%subset(strand=='-'))
# ifelse(juncdat$gene_strand=='-',end(tpexons),start(tpexons))- ifelse(juncdat$gene_strand=='-',start(fpexons),end(fpexons))

# 

# library(Rsamtools)

# fpexons <- ifelse(juncdat[['gene_strand']]=='-',juncex2,juncex1)
# tpexons <- ifelse(juncdat[['gene_strand']]=='-',juncex1,juncex2)

# juncexseq1<-getSeq(x=FaFile('../../genomes/hg19.fa'),juncex1)%>%setNames(NULL)
# juncexseq2<-getSeq(x=FaFile('../../genomes/hg19.fa'),juncex2)%>%setNames(NULL)


# # junctiongr <- GRanges(juncdat[[1]],IRanges(nchar(fpexons),w=rep(2,length(fpexons))))



# str(
#  mapply(F=paste0,
#  	fpexons,tpexons
#  	)
#  )




# translate(juncexseq1)