library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)

#load some utility functions
source("/fast/groups/ag_ohler/dharnet_m/Ribo_Lausanne/functions.R")
#load excel data, parse
stopifnot(c("ere_Peptides", "lncRNA_peptides") == readxl::excel_sheets('ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx'))

#read in the data on the repeat associated peptides
ere_table <- readxl::read_excel('ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=1)
colnames(ere_table) <- c("MS_identified_peptide", "ID", "seqnames", "start", "end", "repName", 
"strand", "class.fam", "unMerged", "FPKM RNA seq OD5PCTRL", "FPKM RNA seq OD5PDAC")


#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel('ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=2)%>%
	set_colnames(c("MS-identified peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))


#import genome as fasta file object
genome <- ('../genomes/hg38.fa'%T>%{stopifnot(file.exists(.)) })%>%FaFile
annofile <- 'annotation/gencode.v22.annotation.gtf'%T>%{stopifnot(file.exists(.)) }

ere_gr <- DT2GR(ere_table,seqinf=seqinfo(genome))

read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
	f=tempfile();
	stopifnot(file.exists(annofile))
	catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
	system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
	out = import(f,format=fformat) 
	file.remove(f)
	out
}

genes<-read_compressed_gfile(annofile,'gene')
start_codons<-read_compressed_gfile(annofile,'start_codon')
transcripts<-read_compressed_gfile(annofile,'transcript')
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)


#now get the relevant genes for our lincRNAs
stopifnot(all(lncRNA_table$gene_id %in% genes$gene_id))
lincgenes <- subset(genes, gene_id %in% lncRNA_table$gene_id)


#let's play with maptToTrnascripts

# testexons <- exons%>%head(10e3)%>%subset(transcript_id!=last(transcript_id))


# # rpls <- genes%>%subset(gene_name%>%str_detect('RPL'))
# rpls19 <- genes%>%subset(gene_name%>%str_detect('RPL'))

# import('bigwigs/star_OD5P_05_uM_DAC_1_Aligned.sortedByCoord.out.bam_Ribo_coverage_plus.bw',which = rpls)%>%{strand(.)<-'+';.} %>%
# 	subsetByOverlaps(rpls)%>%.$score%>%sum

# import('bigwigs/star_OD5P_05_uM_DAC_1_Aligned.sortedByCoord.out.bam_Ribo_coverage_plus.bw',which = rpls)%>%{strand(.)<-'+';.} %>%
# 	subsetByOverlaps(rpls19)%>%.$score%>%sum




##########Producing metaprofiles to 


allbigwigs<-Sys.glob('bigwigs/*_Ribo_coverage_*.bw')


bigwigpairlist <- allbigwigs%>%
	data_frame(file=.)%>%
	mutate(base=file%>%basename%>%str_replace('plus|minus',''))%>%
	mutate(strand=file%>%basename%>%str_extract('plus|minus'))%>%
	mutate(strand = case_when(
		strand=='plus'~'+',
		strand=='minus'~'-'
	)) %>% 
	arrange(desc(strand))%>%
	{split(.,.$base)}%>%{map(.,~split(.$file,.$strand))%>%map(rev)}

stopifnot(bigwigpairlist%>%names%>%n_distinct%>%`>`(3))
stopifnot(bigwigpairlist%>%.[[1]]%>%names%>%`==`(c('+','-')))

startcod_windows<-start_codons%>%
	# subset(gene_name%>%str_detect('RPL'))%>%
	resize(101,'center')

#TODO this logic can be split up a bit to work on bam files say
get_5p_profiles<- function(startcod_windows,bigwigpair){

	
	stopifnot(c('+','-') %in% names(bigwigpair))
	
	# for(i in ls()) assign(i,get(i),envir=.GlobalEnv)
	
	fp_profile_data <- 
		lapply(bigwigpair,FUN=import,which = startcod_windows)%>%
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		Reduce(f=c,.)%>%
		mergeByOverlaps(startcod_windows[,'transcript_id'])%>%
		mergeGR2DT%>%
		transmute(tid=transcript_id,pos = start-start.1+1,score = score)%>%
		group_by(tid)%>%
		# dplyr::filter(sum(score)>20)%>%
		mutate(score = score /sum(score))%>%
		mutate(pos=pos-50)%>%
		group_by(pos)%>%
		summarize(score=mean(score))%>%
		mutate(frame=as.factor(pos %% 3))
	
	return(fp_profile_data)
}

pairnum = 1
fp_profile_data <- get_5p_profiles(startcod_windows[,'transcript_id'],bigwigpairlist[[pairnum]])

rangenames <- 'Start Codon Regions'
ribofilenames <- names(bigwigpairlist)[pairnum]
fp_profplot <- 
	ggplot(fp_profile_data,aes(fill=frame,x=pos,y=score))+
	geom_bar(stat='identity')+
	coord_cartesian(xlim=-50:50)+
	ggtitle(str_interp("5' Metaprofile for ranges:\n${rangenames}\nfile:${ribofilenames}"))
	theme_bw()

fpproffile <- 'plots/fp_profplot_lorenzRibo_rpls.svg'%T>%{normalizePath(.)%>%message}

svglite(fpproffile);print(fp_profplot);dev.off()





# sitesovergenes <- 
# 	lapply(bigwigpairlist[[1]],import,which = unlist(genes)%>%head(1000))%>%
# 	{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
# 	Reduce(f=c,.)%>%
# 	subset(score>0)


summarize_ov_score <- function(genes,sitesovergenes,myfun=sum){
	data_frame(gene=queryHits(ov),score=sitesovergenes$score[subjectHits(ov)])%>%
	left_join(data_frame(gene=seq_along(genes)),.)%>%
	group_by(gene)%>%
	dplyr::summarise(score=myfun(score))%>%
	arrange()
	.$score
}

# genes$score <- summarize_ov_score(genes,sitesovergenes)

# genes$score%>%table





#############Now individual Plots

#TODO this logic can be split up a bit to work on bam files say
get_riboprofdata<- function(exons_tr,bigwigpair){
	
	message( names(bigwigpair))
	stopifnot(c('+','-') %in% names(bigwigpair))

	for(i in ls()) assign(i,get(i),envir=.GlobalEnv)




	seqlevs <-list(profilegrange,exons_tr)%>%unique%>%as.character

	shared_seqinfo <- suppressWarnings(intersect(seqinfo(BigWigFile(bigwigpair[[1]])),seqinfo(exons_tr)))
	
	trseqinfo <- Seqinfo(seqnames=names(exons_tr),seqlengths=as.vector(sum(width(exons_tr))))

	profilegrange <- 
		suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		Reduce(f=c,.)%>%
		subset(score>0)


	#now map our prifle data to the exons space
	rle <- suppressWarnings(mapToTranscripts(profilegrange,exons_tr))
	rle$score<-profilegrange$score[rle$xHits];
	seqinfo(rle)<-trseqinfo[seqinfo(rle)@seqnames]
	rle <- coverage(rle, weight='score')

	rle %>%
		#selecting the transcirpt we want
		lapply(FUN=.%>%
		#turn that into a dataframe
			{	
				pos = which(.!=0)
				data_frame(pos=pos, score = as.vector(.[pos]))
			}
		)%>%
		bind_rows(.id='tid')%>%
		mutate(frame=as.factor(pos %% 3))
}


#pick a set of files and a transcript
mytrans <- transcripts%>%subset(gene_name%>%str_detect('RPL'))%>%sample(1)
pairnum = 1
#get the profile data
riboprofdata <- get_riboprofdata(exons[mytrans$transcript_id],bigwigpairlist[[pairnum]])
#check there's data
stopifnot(riboprofdata$score%>%sum%>%`>`(0))
#now print the plot
rangenames <- mytrans$transcript_name
ribofilenames <- names(bigwigpairlist)[pairnum]
#point to focus on so we can see the frame
spoint <- riboprofdata%>%{.$pos[.$score>median(.$score)][1]}
#
transprofplot <-
	riboprofdata %>%
	# mutate(score= log10(score+(0.1*(min(score))))) %>%
	ggplot(aes(fill=frame,color=frame,x=pos,y=score))+
	geom_bar(stat='identity')+
	coord_cartesian(xlim=c(spoint-50,spoint+50))+
	ggtitle(str_interp("Riboseq Read Profile for:\n${mytrans}\nfile:${ribofilenames}"))+
	theme_bw()

lociplotfile <- 'plots/loci_riboprofiles/test.svg'%T>%{normalizePath(.)%>%message}

svglite(lociplotfile);print(transprofplot);dev.off()
















# readxl::read_excel('ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx')

# vcf <- readLines('ext_data/OD5P_mq_050718.vcf',10)


# df$ensembl_gene_id <- gsub('\\..+$', '', df$ensembl_gene_id)

# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- df$ensembl_gene_id

# symbol <- getBM(filters = "ensembl_gene_id",
#                 attributes = c("ensembl_gene_id","hgnc_symbol"),
#                 values = genes, 
#                 mart = mart)

# df <- merge(x = symbol, 
#               y = df, 
#               by.x="ensembl_gene_id",
#               by.y="ensembl_gene_id")



# ##______
# #this code will lift our vcfs over to the 38 annotation
# GRanges(
# 	c('chr1:1-10:+',
# 	'chr1:100-110:+',
# 	'chr1:201-210:+',
# 	'chr1:301-310:+')
# )%>%{.$transcript_name='a';.}%>%
# {.$exon_='a';.}

# chainurl = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
# chainfile = basename(chainulr)
# system(str_interp('wget ${chainurl} -O ${chainfile}'))
# chain <- import.chain(chainfile)
# vcfgrs %<>% map(liftOver,chain)