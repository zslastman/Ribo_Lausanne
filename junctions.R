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

message('...done')

#load arguments
getwd()
args <- c(
	juncfiledir = '../ext_data/junctionfiles/',
  	GFF = 'my_gencode.v24lift37.annotation.gff3',
   	v19anno='../annotation/gencode.v19.annotation.gtf',
	v22anno='../annotation/gencode.v22.annotation.gtf',
	v24anno='../annotation/gencode.v24lift37.annotation.gtf',
  	outputfolder = 'junctions'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])


dir.create(outputfolder,showWarn=T)

read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
	f=tempfile();
	stopifnot(file.exists(annofile))
	catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
	system(str_interp('${catbin} ${annofile} |   grep -e "\t${annotype}\t" > ${f}'));
	out = import(f,format=fformat) 
	file.remove(f)
	out
}


library(rtracklayer)
dplyr::select->select
dplyr::filter->filter

allanno <- GenomicFeatures::makeTxDbFromGFF(GFF)
introns<-intronsByTranscript(allanno,use=T)

gffbase <- (tools::file_path_sans_ext(basename(GFF)))

introns%>%
	unlist%>%
	{data_frame(chromosome=as.character(seqnames(.)),start=start(.),end=end(.))}%>%
	write_tsv(col_names=FALSE,paste0('junction_exons/',gffbase,'.tsv')) 


####Now the junction files
genes <- read_compressed_gfile(GTF,'gene')

# juncfile <- 'ref_metadata.filt.tsv'
juncfiles <- list.files(rec=TRUE,patt='.*.tsv$',juncfiledir,full=TRUE)
juncids <- juncfiles%>%basename%>%tools::file_path_sans_ext()
juncfile<-juncfiles%>%str_subset('hits')
juncfile<-juncfiles%>%sample(1)

genechrs <- 
	c(v19anno,v22anno,v24anno)%>%
	map(.%>%{
		annofile = .
		str_interp('cat ${annofile} | perl -lne \'/(.*?)\\t.*gene_id\\W+([\\w\\.]+)/;print $1,"\\t",$2\'')%>%
		fread(header=F)%>%
		set_colnames(c('seqnames','gene_id'))
	})%>%
	bind_rows%>%
	distinct(seqnames,gene_id)


testj<-mclapply(juncfiles%>%str_subset('hits'),function(juncfile){
	message(juncfile)
	juncid <- juncfile%>%basename%>%tools::file_path_sans_ext()
	juncdat <- fread(juncfile)
	if(!'output_id' %in% colnames(juncdat)) juncdat[['output_id']] <- juncdat[['ID']]

	if(!'gene_id' %in% colnames(juncdat)) juncdat[['gene_id']] <- juncdat[['gene_name']]
	juncdat%>%colnames
	juncdat %<>% dplyr::filter(!str_detect(exons_coor,'\\.;\\.'))


	# juncdat$gene_id %<>%str_replace('\\.\\d+$','')
	juncdat%<>%left_join(genechrs,by=c('gene_id'))
	stopifnot(!any(is.na(juncdat$seqnames)))
	

	coordmat <- juncdat[['exons_coor']]%>%str_split_fixed(patt=';',n=4)
	coordmat %<>% apply(2,str_replace,'\\.','NA')
	coordmat %<>% apply(2,as.numeric)
	coordmatgoodrange <- (coordmat[,1]<coordmat[,2]) & (coordmat[,3]<coordmat[,4])
	coordmat %<>% .[coordmatgoodrange,]
	juncdat %<>% .[coordmatgoodrange,]
	

	message('make grs')

	fpexons <- GRanges(juncdat$seqnames,IRanges(coordmat[,1],coordmat[,2]),strand=juncdat$gene_strand)
	start(fpexons) <- start(fpexons)+1
	tpexons <- GRanges(juncdat$seqnames,IRanges(coordmat[,3],coordmat[,4]),strand=juncdat$gene_strand)
	start(tpexons) <- start(tpexons)+1
	junctranscript <- GRanges(juncdat$seqnames,
		IRanges(pmin(coordmat[,1],coordmat[,3]),pmax(coordmat[,2],coordmat[,4])),strand=juncdat$gene_strand)
	start(junctranscript) <- start(junctranscript)+1



	juncintrons <- pgap(fpexons,tpexons)

	tsvfile = paste0(outputfolder,juncid,'.tsv')
	juncintrons%>%{data_frame(chromosome=as.character(seqnames(.)),start=start(.),end=end(.))}%>%
		write_tsv(col_names=FALSE,tsvfile) 

	juncintrons$type = 'intron'
	juncintrons$gene_id <- juncdat$gene_id
	juncintrons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	juncintrons$ID <- paste0(juncid,'_',juncdat$output_id)

	fpexons$type <- 'exon'
	tpexons$type <- 'exon'
	fpexons$exon_number = 1
	tpexons$exon_number = 2
	fpexons$gene_id <- juncdat$gene_id
	tpexons$gene_id <- juncdat$gene_id
	fpexons$Parent <- paste0(juncid,'_',juncdat$output_id)
	tpexons$Parent <- paste0(juncid,'_',juncdat$output_id)
	fpexons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	tpexons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	fpexons$ID <-  paste0(juncid,'_',juncdat$output_id,'_1')
	tpexons$ID <- paste0(juncid,'_',juncdat$output_id,'_2')
	junctranscript$type <- 'transcript'
	junctranscript$exon_number <- NA
	junctranscript$gene_id <- juncdat$gene_id
	junctranscript$Parent <- juncdat$gene_id
	junctranscript$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	junctranscript$ID <- paste0(juncid,'_',juncdat$output_id)
	junctranscript%>%mcols%>%colnames
	fpexons%>%mcols%>%colnames

	allgtf = c(fpexons,tpexons,junctranscript)

	export(allgtf,paste0(outputfolder,juncid,'.gtf')) 
	
	list(
		allgtf,
		juncintrons,
		tsvfile
	)
})


mergedjunctionexons <- junctiongrlist%>%map(1)%>%Reduce(f=c)%>%sort
mergedjunctionintrons <- junctiongrlist%>%map(2)%>%Reduce(f=c)%>%sort
tsvfiles <- junctiongrlist%>%map(3)

dir.create(showWarnings = F,rec=TRUE,'junction_exons')
export(mergedjunctionexons,paste0('junction_exons/','all','.gtf')) 
export(mergedjunctionintrons,paste0('junction_exons/','allintrons','.gtf')) 


mcols(mergedjunctionexons)[,c('gene_id','transcript_id')]%>%as.data.frame%>%write_tsv('junction_gene_transcript_map.tsv')


# system("find junction_exons/ | grep -v 'all.gtf' | xargs cat > junction_exons/all.gtf")
alltsvstr <- paste0(tsvfiles,collapse=' ')
system(str_interp("cat ${alltsvstr}| sort | uniq | grep -ve start > junction_exons/all.tsv")

getwd()

source('../functions.R')

alljsdt <- mergedjunctionexons%>%GR2DT

juncs<-fread('junction_exons/all.tsv')
juncs <- juncs%>%subset(V1=='chr1')%>%arrange(V2)%>%as_data_frame%>%distinct
juncs <- GRanges(juncs[[1]],IRanges(juncs[[2]],juncs[[3]]))


junctionbams <- Sys.glob('star/data/*/*bam')%>%
	discard(~str_detect(.,'transcript'))%>%
	# .[[3]]
	identity
bamnames <- junctionbams%>%basename%>%tools::file_path_sans_ext()%>%setNames(junctionbams)
junctionbams%<>%setNames(bamnames)

testjunctables<-Sys.glob('junction_exons/junction*.tsv')%>%setNames(.,basename(.))


hits<-mergedjunctionintrons%>%subset(ID%>%str_detect('junction_diff_hits'))



uniqueintrons <- mergedjunctionintrons%>%unique
uniqueintrons$id_in_unique <- seq_along(uniqueintrons)
mergedjunctionintrons$id_in_unique <- match(mergedjunctionintrons,uniqueintrons)
juncgr<-uniqueintrons



# juncgr<-test
# junctionbam <- junctionbams[33]
# juncgr <- testj[[1]][[2]]
#ignore the new samples,, for now....
junctionbams<-junctionbams[junctionbams%>%str_detect('ctrl\\dL')%>%not]

bamcounts <- mclapply(mc.cores=20,
	junctionbams,juncgr,FUN=function(junctionbam,juncgr){
# bamcounts <- mclapply(junctionbams,juncgr,FUN=function(junctionbam,uniqueintrons){
	param <- ScanBamParam(
		flag=scanBamFlag(isSecondaryAlignment=FALSE),what=c("mapq"),
		which=juncgr%>%resize(1),
		tag = "MD"
	)
	reads <- readGAlignments(BamFile(file=junctionbam, yieldSize=NA),param = param)
	#summary gr of the junctions and their scores
	readjuncs <- summarizeJunctions(reads)
	#match predicted juctnions to those in reads
	juncmatch <- match(juncgr%>%{strand(.)='*';.},readjuncs)
	#assign count
	ifelse(
		as.vector(strand(juncgr)=='+'),
		readjuncs$plus_score[juncmatch],
		readjuncs$minus_score[juncmatch]
	)
})
stopifnot(bamcounts%>%map_lgl(is.null)%>%any%>%not)
tbamcounts%>%map_lgl(is.null)%>%which

bamcounts->bamcountsbak

bamcounts %<>% simplify2array
bamcounts%<>%apply(2,function(x){x[is.na(x)]<-0;x})
bamcounts%<>%as_data_frame
bamcounts%<>%mutate(id_in_unique=juncgr$id_in_unique)


stop()


fintrons <- mergedjunctionintrons%>%subset(ID%>%str_detect('junction_diff_hits'))

jmatch <- match(hits$id_in_unique,juncgr$id_in_unique)
hitcounts <- bamcounts[jmatch,]




###First deal with our hits
hittable<-juncfiles%>%str_subset('hits')%>%fread
hittable$ID%<>%as.character
hits$ID<-hits$ID%>%str_replace('.*_','')



hittable%<>%
	left_join(distinct(data_frame(ID=hits$ID,id_in_unique=hits$id_in_unique)),by='ID')%>%
	left_join(distinct(bamcounts),by='id_in_unique')%>%
	select(-id_in_unique)%>%
	identity

hittable$anyOD5P <- hittable%>%select(one_of(bamnames%>%str_subset('OD5P')))%>%as.matrix%>%array_branch(1)%>%map_lgl(~any(.>0))
hittable$any <- hittable%>%select(one_of(bamnames))%>%as.matrix%>%array_branch(1)%>%map_lgl(~any(.>0))

hittable%>%write_tsv('spliced_peptide_hits_with_od5pcounts.tsv'%T>%message)
mergedjunctionexons%>%
	subset(str_detect(ID,'junction_diff_hits'))%>%
	export('../junctionfiles/junction_diff_hits.gtf')


hittable<-read_tsv('spliced_peptide_hits_with_od5pcounts.tsv')


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