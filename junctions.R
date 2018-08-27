# sshbih
# filter<-dplyr::filter
# screen -x ril
read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
	f=tempfile();
	stopifnot(file.exists(annofile))
	catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
	system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
	out = import(f,format=fformat) 
	file.remove(f)
	out
}
library(rtracklayer)
dplyr::select->select
dplyr::filter->filter

# genes <- read_compressed_gfile('../annotation/gencode.v19.annotation.gtf','gene')
# juncfile <- 'ref_metadata.filt.tsv'
juncfiles <- list.files(rec=TRUE,patt='.*.tsv$','../junctionfiles/',full=TRUE)
juncids <- juncfiles%>%basename%>%tools::file_path_sans_ext()

for(juncfile in rev(juncfiles)){
	message(juncfile)
	juncid <- juncfile%>%basename%>%tools::file_path_sans_ext()
	juncdat <- fread(juncfile)
	if(!'output_id' %in% colnames(juncdat)) juncdat[['output_id']] <- juncdat[['ID']]
	if(!'gene_id' %in% colnames(juncdat)) juncdat[['gene_id']] <- juncdat[['gene_name']]
	juncdat%>%colnames
	juncdat %<>% dplyr::filter(!str_detect(exons_coor,'\\.;\\.'))


	seqexonmap<-data_frame(seqnames=seqnames(genes)%>%as.vector,gene_id=genes$gene_id)%>%distinct
	juncdat%<>%left_join(seqexonmap,by=c('gene_id'))

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
	# junctranscript <- GRanges(juncdat$seqnames,
	# 	IRanges(pmin(coordmat[,1],coordmat[,3]),pmax(coordmat[,2],coordmat[,4])),strand=juncdat$gene_strand)
	# start(junctranscript) <- start(junctranscript)+1



	juncintrons <- pgap(fpexons,tpexons)

	# juncintrons%>%{data_frame(chromosome=as.character(seqnames(.)),start=start(.),end=end(.))}%>%
	# 	write_tsv(col_names=FALSE,paste0('junction_exons/',juncid,'.tsv')) 

	juncintrons$type = 'intron'
	juncintrons$gene_id <- juncdat$gene_id
	juncintrons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	juncintrons$ID <- paste0(juncid,'_',juncdat$output_id)

	# fpexons$type <- 'exon'
	# tpexons$type <- 'exon'
	# fpexons$exon_number = 1
	# tpexons$exon_number = 2
	# fpexons$gene_id <- juncdat$gene_id
	# tpexons$gene_id <- juncdat$gene_id
	# fpexons$Parent <- paste0(juncid,'_',juncdat$output_id)
	# tpexons$Parent <- paste0(juncid,'_',juncdat$output_id)
	# fpexons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	# tpexons$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	# fpexons$ID <-  paste0(juncid,'_',juncdat$output_id,'_1')
	# tpexons$ID <- paste0(juncid,'_',juncdat$output_id,'_2')
	# junctranscript$type <- 'transcript'
	# junctranscript$exon_number <- NA
	# junctranscript$gene_id <- juncdat$gene_id
	# junctranscript$Parent <- juncdat$gene_id
	# junctranscript$transcript_id <- paste0(juncid,'_',juncdat$output_id)
	# junctranscript$ID <- paste0(juncid,'_',juncdat$output_id)
	# junctranscript%>%mcols%>%colnames
	# fpexons%>%mcols%>%colnames
	# if(!exists('mergedjunctionexons')) 	mergedjunctionexons <- sort(c(fpexons,tpexons,junctranscript))
	# mergedjunctionexons <- sort(c(mergedjunctionexons,fpexons,tpexons,junctranscript))

	if(!exists('mergedjunctionintrons')) 	mergedjunctionintrons <- c(juncintrons)
	mergedjunctionintrons <- c(mergedjunctionintrons,juncintrons)
}

dir.create(showWarnings = F,rec=TRUE,'junction_exons')
# export(mergedjunctionexons,paste0('junction_exons/',juncid,'.gtf')) 

mcols(mergedjunctionexons)[,c('gene_id','transcript_id')]%>%as.data.frame%>%write_tsv('junction_gene_transcript_map.tsv')


system("find junction_exons/ | grep -v all.gtf | xargs cat > junction_exons/all.gtf")
system("find junction_exons/ | grep -v all.tsv | xargs cat | sort | uniq > junction_exons/all.tsv")

getwd()

source('../functions.R')

alljsdt <- mergedjunctionexons%>%GR2DT


juncs<-fread('junction_exons/all.tsv')
juncs <- juncs%>%subset(V1=='chr1')%>%arrange(V2)%>%as_data_frame%>%distinct
juncs <- GRanges(juncs[[1]],IRanges(juncs[[2]],juncs[[3]]))


junctionbams <- Sys.glob('junction_star/data/*/*bam')%>%
	discard(~str_detect(.,'transcript'))%>%
	# .[[3]]
	identity

testjunctables<-Sys.glob('junction_exons/junction*.tsv')%>%setNames(.,basename(.))

mergedjunctionintrons%>%length%>%divide_by(1e6)
mergedjunctionexons%>%unique%>%length%>%divide_by(1e6)

#we'll also shift these one to left and right to allow for that
uniqueintrons <- mergedjunctionintrons%>%unique
uniqueintrons$id_in_unique <- seq_along(uniqueintrons)
mergedjunctionintrons$id_in_unique <- match(mergedjunctionintrons,uniqueintrons)


#this turns out to p. much unnecessary
uniqueintrons_shift <- c(
	uniqueintrons%>%{.$shift=0;.},
	uniqueintrons%>%shift(1)%>%{.$shift=1;.},
	uniqueintrons%>%shift(-1)%>%{.$shift=-1;.}
)


for(junctionbam in (junctionbams%>%str_subset('OD5P')%>%.[1:2])){
	
	bamnames<-bamnames[[junctionbam]]

	#get the reads as a gapped alignment object


	param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")
	reads <- readGAlignments(BamFile(file=junctionbam, yieldSize=NA),param = param)


	#summary gr of the junctions and their scores
	readjuncs <- summarizeJunctions(reads)

	#match predicted juctnions to those in reads
	# juncmatch <- match(uniqueintrons_shift%>%(strand(.)='*';.),readjuncs)
	# #assign count
	# mcols(uniqueintrons_shift)[[bamname]]<-readjuncs$score[juncmatch]

	#match predicted juctnions to those in reads
	juncmatch <- match(uniqueintrons%>%{strand(.)='*';.},readjuncs)
	#assign count
	mcols(uniqueintrons)[[bamname]] <- ifelse(
		as.vector(strand(uniqueintrons)=='+'),
		readjuncs$plus_score[juncmatch],
		readjuncs$minus_score[juncmatch]
	)

}
hits<-mergedjunctionintrons%>%subset(ID%>%str_detect('junction_diff_hits'))
test < uniqueintrons%>%subset(id_in_unique%in%(test$id_in_unique))

bamnames <- junctionbams%>%basename%>%tools::file_path_sans_ext()%>%setNames(junctionbams)
junctionbams%<>%setNames(bamnames)

juncgr<-uniqueintrons
juncgr<-test

bamcounts <- mclapply(
	junctionbams,
	bamnames,juncgr,FUN=function(junctionbam,bamnames,juncgr){
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


bamcounts %<>% simplify2array
colnames(bamcounts)
bamcounts%<>%apply(2,function(x){x[is.na(x)]<-0;x})
bamcounts%<>%as_data_frame
bamcounts%<>%mutate(id_in_unique=juncgr$id_in_unique)



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