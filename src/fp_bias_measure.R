library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)
library(stringr)
library(data.table)
source('../src/functions.R')
select<-dplyr::select
summarise<-dplyr::summarise
filter<-dplyr::filter
slice<-dplyr::select



NGENES<- 1e3
annotation <- import('my_gencode.v24lift37.annotation.gtf')
sample_ids <- fread('sample_parameter.csv')%>%head(30)%>%.$sample_id
bamfiles <- Sys.glob('star/data/*/*.bam')%>%grep(v=T,inv=T,patt='transcript')
cdslengths <- annotation%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%width%>%sum
shortgenes <- names(which(cdslengths<1e3))

fpbiaslist <- list()
for (sample_id in sample_ids[1]){
	try({
		bamfile <- str_subset(bamfiles,sample_id)
	message(bamfile)

	topgenes <- fread('feature_counts/all_feature_counts')%>%
		select(feature_id,count=!!sample_id)%>%
		filter(!feature_id  %in% shortgenes)%>%
		arrange(desc(count))%>%
		head(NGENES)%>%.$feature_id

	topcds <- annotation%>%subset(type=='CDS')%>%subset(gene_id %in% topgenes)%>%.[,c('transcript_id','gene_id')]
	bamseqs<-seqinfo(BamFile(bamfile))@seqnames
	topcds<-keepSeqlevels(topcds,intersect(bamseqs,seqnames(topcds)),pruning='coarse')

	readsontopgene<-readGappedReads(bamfile,param=ScanBamParam(which=topcds))%>%as('GRanges')%>%resize(1)

	topcdscounts<-countOverlaps(topcds,readsontopgene)

	besttrs<-data_frame(count=topcdscounts,transcript_id=topcds$transcript_id,gene_id=topcds$gene_id)%>%
		group_by(gene_id,transcript_id)%>%
		summarise(count=sum(count))%>%
		group_by(gene_id)%>%
		dplyr::slice(which.max(count))%>%
		.$transcript_id

	topcds%<>%split(.,.$transcript_id)
	topcds <- topcds[besttrs]



	#map the reads to the best transcirpts for the metaplots
	readsontrs<-mapToTranscripts(readsontopgene,topcds[besttrs])
	#count numbers in the first 250 and next 250
	startcounts<-GRanges(besttrs,IRanges(50,250))%>%countOverlaps(readsontrs)
	starpluscounts <- GRanges(besttrs,IRanges(250,550))%>%countOverlaps(readsontrs)

	(startcounts / starpluscounts)%>%keep(is.finite)%>%keep(Negate(is.nan))%>%mean

	readcov <- readsontrs%>%coverage
	cdslengths <- readcov%>%lapply(length)%>%unlist

	readcov<-readcov[cdslengths>1e3]
	cdslengths <- readcov%>%lapply(length)%>%unlist

	#make sets of genes based on gene length

	botquart<-quantile(cdslengths,0.25)
	topquart<-quantile(cdslengths,0.75)
	
	absvshort<-readcov[cdslengths<botquart]%>%sapply(function(.)as.vector(.)%>%c(.,rep(NA,1e3))%>%head(1e3))
	absvlong<-readcov[cdslengths>topquart]%>%sapply(function(.)as.vector(.)%>%c(.,rep(NA,1e3))%>%head(1e3))


	posdensdf<-bind_rows(list(
		`Bottom 25% Length` = data_frame(pos = 1:nrow(absvshort), meanNormalizedRiboDensity=absvshort%>%sweep(.,2,STATS=colSums(.,na.rm=T),FUN='/')%>%rowMeans(na.rm=T)),
		`Top 25% Length` = data_frame(pos = 1:nrow(absvlong), meanNormalizedRiboDensity=absvlong%>%sweep(.,2,STATS=colSums(.,na.rm=T),FUN='/')%>%rowMeans(na.rm=T))
	),.id='gene_set')

	fpbiaslist %<>% append(
		posdensdf%>%mutate(seg=floor((pos)/250))%>%filter(seg%in%0:1)%>%group_by(seg)%>%summarise(segmean = sum(meanNormalizedRiboDensity)) %>% {.$segmean[1]/.$segmean[2]}
		)

	# posdensdf%>%filter(gene_set==unique(gene_set)[1])%>%.[[3]]%>%rollmax(6)%>%diff%>%plot
	stop()

	plots<-ggpubr::ggarrange(nrow=2,
		posdensdf%>%
		ggplot(aes(y=meanNormalizedRiboDensity,x=pos,color=gene_set))+geom_smooth()+
		scale_x_continuous(name='Distance From AUG')+theme_minimal()+ggtitle(str_interp('Riboseq Coverage Over first 1kb ${sample_id}')),
		qplot(x=identity(cdslengths),fill=I('darkblue'),geom='histogram')+
		geom_vline(aes(xintercept=identity(botquart),size=I(2)))+geom_vline(aes(xintercept=identity(topquart),size=I(2)))+
		ggtitle(str_interp('CDS Length Quartiles '))+
		theme_minimal()+scale_x_log10(name='CDS length (bp)')
	)
	plots
	dev.off()
	ggsave(plots,file=str_interp('../plots/Coverage_Absdist_${sample_id}.pdf'))
	ggsave(plots,file=str_interp('../plots/Coverage_Absdist_${sample_id}.svg'))
	})
}

fpbiaslist%>%unlist%>%setNames(sample_ids)%>%enframe('sample','fp_bias')%>%write_tsv('fp_bias.txt')







