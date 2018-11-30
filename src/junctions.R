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
outputfolder <- paste0(outputfolder,'/')

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
GFF%>%file.exists
allanno <- GenomicFeatures::makeTxDbFromGFF(GFF)
introns<-intronsByTranscript(allanno,use=T)

gffbase <- (tools::file_path_sans_ext(basename(GFF)))

introns%>%
	unlist%>%
	{data_frame(chromosome=as.character(seqnames(.)),start=start(.),end=end(.))}%>%
	write_tsv(col_names=FALSE,paste0(outputfolder,gffbase,'.tsv')) 



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


junctiongrlist<-mclapply(juncfiles,function(juncfile){
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

dir.create(showWarnings = F,rec=TRUE,outputfolder)
export(mergedjunctionexons,paste0(outputfolder,'all','.gtf')) 
export(mergedjunctionintrons,paste0(outputfolder,'allintrons','.gtf')) 



mcols(mergedjunctionexons)[,c('gene_id','transcript_id')]%>%as.data.frame%>%write_tsv('junction_gene_transcript_map.tsv')


# system("find junction_exons/ | grep -v 'all.gtf' | xargs cat > junction_exons/all.gtf")
alltsvstr <- paste0(tsvfiles,collapse=' ')
system(str_interp("cat ${alltsvstr}| sort | uniq | grep -ve start > junction_exons/all.tsv"))


source('../functions.R')





testjuncs <- import('junctions/junction_diff_hits.gtf')%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%lapply(gaps)%>%lapply('[',2)%>%GRangesList%>%unlist

uniqueintrons <- mergedjunctionintrons%>%unique
uniqueintrons$id_in_unique <- seq_along(uniqueintrons)
mergedjunctionintrons$id_in_unique <- match(mergedjunctionintrons,uniqueintrons)

#save these 
mcols(mergedjunctionintrons)[,c('ID','id_in_unique')]%>%as.data.frame%>%write_tsv('junctions/id_unique.txt')

uniqueintrons%>%saveRDS('junctions/uniqueintrons.rds')

jcfiles <- Sys.glob("junctioncounts/*/*.junctioncounts.tsv")%>%grep(v=T,inv=T,patt='L5|L7')
jcfile <- 'junctioncounts/OD5P_05_uM_DAC_1/OD5P_05_uM_DAC_1.junctioncounts.tsv'

hitcounts<-mclapply(jcfiles,function(jcfile) fread(str_interp("grep -e hit junctions/id_unique.txt | cut -f 2 | awk '{ print(\"^\"$1\"\\\\s\")}' | grep -f - ${jcfile}")))

hitcounts%>%Reduce(f=partial(left_join,by='V1'))%>%dplyr::select(-V1)%>%rowSums%>%`!=`(0)%>%table

hitcounts[jcfiles%>%str_detect('OD5P')]%>%Reduce(f=partial(left_join,by='V1'))%>%dplyr::select(-V1)%>%rowSums%>%`!=`(0)%>%table
