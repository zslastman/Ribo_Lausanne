library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(Gviz)
library(rtracklayer)
genomestr<-'hg19'

options(ucscChromosomeNames=F)
source("/fast_new/work/groups/ag_ohler/dharnet_m/Ribo_Lausanne/src/functions.R")

sampledf <- fread('../pipeline/sample_parameter.csv')
sampleids <- sampledf %>%filter(assay=='ribo')%>%.[[1]]
extract_id <- purrr::partial(extract_oneof,ids=sampleids)

plotfolder <- 'plots/orf_riboprofiles/' %T>% dir.create

#get exon data
annofile <- '../pipeline/my_gencode.v24lift37.annotation.gtf'%T>%{stopifnot(file.exists(.)) }
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)

#get bigwigpairlists

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
ribofilenames<-names(bigwigpairlist)
stopifnot(all(!is.na(extract_id(ribofilenames))))




#get our table of links - from which we need the orf names, and samples, and peptide sequences
lncs2_table <- readxl::read_excel('../ext_data/OD5P_ctrl_case_sequence_search.xlsx')
lncs2_table%<>%rename('HLA-peptide','"MS-identified peptide')
lncs2_table$gene_id <- str_extract(lncs2_table$header,'ENSG[^|]+')
lncs2_table$orfname <-lncs2_table$header%>%str_split('(>|\\|)')%>%map_chr(2)
lncs2_table$tr <- lncs2_table$orfname%>%str_split('_')%>%map_chr(1)

#
orfnames <- lncs2_table$orfname
orfsamples <- lncs2_table$header%>%str_split('(>|\\|)')%>%map_chr(7)
orftranscripts <- lncs2_table$tr

#make a track of all orfs
allsatannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
allsatannorfs <- 
	allsatannfiles%>%
	setNames(.,extract_id(.))%>%
	mclapply(load_objs)
#get orfs near ours
all_orfs_tr_orfs <- allsatannorfs%>%map(.%>%.$ORFs_tx%>%subset((seqnames) %in% (orftranscripts)))
#as a data.table
for (i in seq_along(all_orfs_tr_orfs))all_orfs_tr_orfs[[i]]$Protein%<>%as.character
all_orfs_tr_orfs <- all_orfs_tr_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')
#width of the transcripts
transcript_seqinfo <- exons[unique(orftranscripts)]%>%width%>%sum%>%{Seqinfo(seqnames=names(.),seqlength=.)}
#
orftrack<-all_orfs_tr_orfs%>%
	DT2GR(seqinf=transcript_seqinfo)%>%
	unique%>%
	{Gviz::AnnotationTrack(
		name='ORFs',.,
		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
		id='.',showFeatureId=TRUE,genome=genomestr)}

displayPars(orftrack) <- list(height=2)

#connect orfs to their peptide info
peptide_loc<-all_orfs_tr_orfs%>%inner_join(
	lncs2_table%>%dplyr::select(peptide=`HLA-peptide`,orfname),
	by=c('ORF_id_tr'='orfname')
)
#use this to change their start and end
peptide_loc[c('start','end')]<-peptide_loc%$%{
	map2(Protein, peptide,str_locate)}%>%
	map(as.data.frame)%>%
	bind_rows%>%
	mutate(start=((start-1)*3) + 1,end=end*3)%>%
	set_colnames(c('start','end'))
peptidetrack<-peptide_loc%>%
	DT2GR(seqinf=transcript_seqinfo)%>%
	unique%>%
	{Gviz::AnnotationTrack(
		name='peptide',.,
		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
		id='.',showFeatureId=TRUE,genome=genomestr)}



get_riboproftrack<- function(exons_tr,bigwigpair,rangename){
	stopifnot(c('+','-') %in% names(bigwigpair))

	# for(i in ls()) assign(i,get(i),envir=.GlobalEnv)
	profilegrange <- 
		suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
		# suppressWarnings(lapply(bigwigpair,import,which = unlist(testwindow)))
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		{suppressWarnings(Reduce(f=c,.))}%>%
		# subsetByOverlaps(testwindow)%>%
		subset(score>0)
	# profilegrange<<-profilegrange
	#now map our prifle data to the exons space
	message(sum(profilegrange$score))
	gr <- suppressWarnings(mapToTranscripts(profilegrange,exons_tr))
	# (mapToTranscripts(testwindow,exons_tr))
	if(length(gr)==1) gr <- c(gr,gr)
	gr$score<-profilegrange$score[gr$xHits];
	# gr%>%subset(score>30)
	# profilegrange%>%subset(score>30)

	datname = bigwigpair%>%extract_id%>%.[[1]]
	
	scoremat <- rep(0,length(gr)*3)%>%matrix(ncol=3)
	scoremat[matrix( c(1:length(gr),(start(gr)%%3)+1 ) ,ncol=2 ) ] <- gr$score
	mcols(gr) = scoremat
	rgbvect <- c('red','green','blue')%>%sort

	if(length(gr)==0) {groupvect <- NULL}else{
		# groupvect <- groupvect[(start(gr)%%3)+1]
		groupvect <- paste0('frame ',1:3)
	}

	DataTrack(gr,name=	datname,chr=rangename,groups=groupvect,col = rgbvect, fill=rgbvect, 
		cex.title=0.3,legend=FALSE,genome=genomestr,chromosome=rangename)
}


getylims<-function(ribotracks,plotgenomewindow){
	inwindranges <- ribotracks%>%map(~ .@range %>% overlapsAny(plotgenomewindow)%>%which)
	ymax <- map2(map(ribotracks,~.@data%>%colSums),inwindranges, ~ .x[.y])%>%unlist%>%sort%>%nth(-2)%>%{.*2}
	c(0,ymax)
}

plot_orf_riboprofile <- function(orf_trspace,orfexons_gen,orf_bigwigpairs,peptidetrack){

	#
	orftr <- names(orf_trspace)%>%str_split('_')%>%map_chr(1)

	#get the orfs inside our viewing windows
	orfexons <- orfexons_gen[[orftr]]
	exontrack<-
		IRanges(c(1,cumsum(width(orfexons))[-length(orfexons)]+1),cumsum(width(orfexons)))%>%
		GRanges(orftr,.)%>%
		{AnnotationTrack(.,name='',col='black',fill='yellow',strand='*',
			genome=genomestr,id=' '	,showFeatureId=TRUE,
			chr=orftr)
		}
	#define a viewing wndow around our orf
	plotgenomewindow <- orf_trspace %>% {
		orf_trspace <- .
		stopifnot(length(orf_trspace)>0)
		plotstart <- start(orf_trspace)%>%min
		plotend <- end(orf_trspace)%>%max
		orfswidth <- plotend-plotstart
		plotstart <- max(0,plotstart - (orfswidth))
		plotend <- min(sum(width(orfexons)),plotend + (orfswidth))

		GRanges(orftr,IRanges(plotstart,plotend))%>%split(.,seqnames(.))
	}

	#get the tracks of data in tr space
	ribotracks <- orf_bigwigpairs%>%
		map(~get_riboproftrack(GRangesList(orfexons)%>%setNames(orftr),.,orftr))

	stopifnot(0 != ribotracks%>%map(~.@range)%>%GRangesList%>%unlist%>%length)
	stopifnot(0 != exontrack@range%>%subsetByOverlaps(plotgenomewindow)%>%length)
	stopifnot(0 != orftrack@range%>%subsetByOverlaps(plotgenomewindow)%>%length)
	# browser()
	#now plot it all
	genomelocstring <- 
		GRanges(
			orfexons%>%seqnames%>%.[1]%>%as.character,
			IRanges(orfexons%>%start%>%.[1],orfexons%>%end%>%.[1])
		)	
	plottitle <- paste0('Riboseq Profile:\n',names(orf_trspace)[1],'\n',orftr,': ',genomelocstring)


	plotTracks(main=plottitle,cex.main=0.5,
		from=start(plotgenomewindow)[[1]],to=end(plotgenomewindow)[[1]],#zoomed in on the orf in question
		c(
			ribotracks, # plot the riboseq signal in tr space
			exontrack, # plot exons in tr space
			orftrack, #plot the orf (and nearby orfs) in tr space
			peptidetrack,
			GenomeAxisTrack(range=GRanges(orftr,IRanges(0,sum(width(orfexons)))),id=orftr)
		),
		
		ylim=getylims(ribotracks,plotgenomewindow),
		# ylim=c(0,3),
		# type='hist',
		chr=orftr,
		genome=genomestr
	)
}


#for each orf
# for(i in seq_along(orfnames)){
for(i in orfnames%>%str_detect('ENST00000457649_85_144')%>%which){

	#get our orf name and sample
	orfname <- orfnames[[i]]
	orfsample <- orfsamples[[i]]

	orfsamples2plot<-sampledf%>%
		group_by(cell_line)%>%
		filter(any(sample_id==orfsample))%>%
		filter(assay=='ribo')%>%.$sample_id

	orf_bigwigpairs <- bigwigpairlist[orfsamples2plot]
	#get the grange object for our orf
	orf_trspace <- allsatannorfs[[orfsample]]$ORFs_tx[orfname]
	
	#get the exons for our transcript
	orftr <- orfname%>%str_split('_')%>%map_chr(1)
	orfexons_gen <- exons[orftr]

	#name te file
	transprofplotfile <- str_interp('${plotfolder}${orfname}_riboprof.pdf')

	#run the plotting code
	pdf(transprofplotfile)
	plot_orf_riboprofile(orf_trspace,orfexons_gen,orf_bigwigpairs,peptidetrack)
	dev.off()
	message(normalizePath(transprofplotfile))

}
