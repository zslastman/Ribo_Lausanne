


datafiles <- allbigwigs

ribotracks <- datafiles%>%
  data_frame(file=.)%>%
  mutate(trackname =extract_id(file,sampleids),
  	strand=ifelse(str_detect(file,'neg'),'-','+')
  )%>%
  split(.,.$strand)%>%
  map(.%>%{map2(.$file,.$trackname,.f = make_track)}%>%setNames(.,.$trackname))


options(ucscChromosomeNames=F)


get_riboproftrack<- function(exons_tr,bigwigpair){
	
	message( names(bigwigpair))
	stopifnot(c('+','-') %in% names(bigwigpair))

	for(i in ls()) assign(i,get(i),envir=.GlobalEnv)

	profilegrange <- 
		suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		{suppressWarnings(Reduce(f=c,.))}%>%
		subset(score>0)
	#now map our prifle data to the exons space
	gr <- suppressWarnings(mapToTranscripts(profilegrange,exons_tr))
	gr$score<-profilegrange$score[gr$xHits];
	DataTrack(gr[,'score'],genome=genome,col=c('red','green','blue')[(start(gr)%%2)+1])
}



genome='mine'
genomefile='pipeline/my_hg38.fa'
make_startstoptrack<-function(range,fafile='pipeline/my_hg38.fa'){
	rangeseq <- range%>%getSeq(x=FaFile(fafile))%>%Reduce(f=c)

	Mp<-rangeseq%>%
		str_locate('ATG')%>%
		{GRanges(seqnames(range)[1],IRanges(.[,1],.[,2]))}%>%
		{AnnotationTrack(.,col=c('red','green','blue')[(start(.)%%2)+1],genome=genome,strand='+')}

	Mn<-rangeseq%>%
		str_locate(as.character(reverseComplement(DNAString('ATG'))))%>%
		{GRanges(seqnames(range)[1],IRanges(.[,1],.[,2]))}%>%
		{AnnotationTrack(.,col=c('red','green','blue')[(start(.)%%2)+1],genome=genome,strand='-')}

	stpp<-rangeseq%>%
		str_locate('TAA|TAG|TGA')%>%
		{GRanges(seqnames(range)[1],IRanges(.[,1],.[,2]))}%>%
		{AnnotationTrack(.,col=c('red','green','blue')[(start(.)%%2)+1],genome=genome,strand='+')}

	stpn<-rangeseq%>%
		str_locate('TTA|CTA|TCA')%>%
		{GRanges(seqnames(range)[1],IRanges(.[,1],.[,2]))}%>%
		{AnnotationTrack(.,col=c('red','green','blue')[(start(.)%%2)+1],genome=genome,strand='-')}

	return(list(Mp,Mn,stpp,stpn))

}


for (ribofilename in ribofilenames[1]){
	ribofileid <- extract_id(ribofilename,sampleids)
	# for (rangename in names(ranges2plot)){
	for (rangename in 'ENST00000374922.7'){

		range = ranges2plot[[rangename]]
		rangelength = sum(width(range))
		rangename%<>%str_replace_all('[\\[\\]]','')

		plottitle <- str_interp("Riboseq Read Profile for:\n${rangename} = ${range}\nfile:${ribofileid}")		


		codtrack<-make_startstoptrack(range,genomefile)
		exontrack<-
			GRanges(rangename,IRanges(c(1,cumsum(width(range))[-length(range)]+1),cumsum(width(range))))%>%
			AnnotationTrack(col='yellow',strand='*',genome=genome)



		#name the plot file
		transprofplotfile <- str_interp('plots/loci_riboprofiles/${rangename}/${ribofileid}/${rangename}_${ribofileid}_prof.svg')
#		transprofplotfile%<>%str_replace_all('[^a-zA-z0-9_\\./]+','')
		transprofplotfile%>%dirname%>%dir.create(rec=TRUE,showWarnings=FALSE)
		#print the plot

		svglite(transprofplotfile)
		plotTracks(
			c(
				map(bigwigpairlist[1:2],~get_riboproftrack(range,.)),
				codtrack,
				exontrack
			),type='hist',
			chr=seqnames(range),from=0,to=rangelength
		)
		dev.off()
		message(transprofplotfile%>%normalizePath)
	}
}

load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}
clipid <- . %>% str_replace('\\.\\d+$','')

library(assertthat)
#let's look at the scores for our 
satannorfs <- 
	Sys.glob('SaTAnn/*_Mar5newst_Final_ORFs_files')%>%
	setNames(.,extract_id(.))%>%
	mclapply(load_objs)

all_linc_orfs <- satannorfs%>%map(.%>%.$ORFs_tx%>%subset(clipid(seqnames) %in% clipid(linctranscripts$transcript_id)))
#TODO_ what are these Granges columsn???
#now filter out the weird GRanges columns, for now, and aggregate into one data table 
linc_orfs_dt <- all_linc_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')

#TODO - eveyrthing in linc
lncRNA_table$gene_id%<>%clipid
lncRNA_table%>%inner_join(lncRNA_table,by='gene_id')
lncRNA_table

#this guy's gene id isn't in the talbe even though it matches a transcript who's gene id is
#versioninf issue.
	missing<-linc_orfs_dt$gene_id%>%str_replace('\\.\\d+$','')%>%setdiff(lncRNA_table$gene_id%>%str_replace('\\.\\d+$',''))




linc_orfs_dt%>%subset(gene_id==missing)%>%.$transcript_id%>%is_in(linctranscripts$transcript_id)

linctranscripts%>%subset(gene_id==missing)

best_linc_orfs<-linc_orfs_dt%>%
	group_by(sample,gene_id)%>%
	dplyr::slice(which.min(pval))

#fraction of our lincs that have orfs in them
linctranscripts$gene_id%>%unique%>%is_in(best_linc_orfs$gene_id)%>%{list(sum(.),mean(.))}




