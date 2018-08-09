library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#very similiar to above but reutrns the spec value as welll as the ftest
take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
     if(length(x)<25){
          remain<-50-length(x)
          x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
     }
     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     #this calculates a frequency object from a numeric vector using the 
     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)

     #this selects only the frequency range surrounding 0.33 - the one we are concerned with - but isn't used
     resSpec2<-dropFreqs(resSpec1,0.29,0.39)
     
     freq_max_3nt<-resSpec1$freq[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
     
     Fmax_3nt<-resSpec1$mtm$Ftest[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]

     spect_3nt<-resSpec1$spec[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]

     return(c(Fmax_3nt,spect_3nt))
     
}



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
	f=tempfile();
	stopifnot(file.exists(annofile))
	catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
	system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
	out = import(f,format=fformat) 
	file.remove(f)
	out
}


get_riboprofdata<- function(exons_tr,bigwigpair){
	
	message( names(bigwigpair))
	stopifnot(c('+','-') %in% names(bigwigpair))

	for(i in ls()) assign(i,get(i),envir=.GlobalEnv)



	profilegrange <- 
		suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		Reduce(f=c,.)%>%
		subset(score>0)

	seqlevs <-list(profilegrange,exons_tr)%>%unique%>%as.character

	shared_seqinfo <- suppressWarnings(intersect(seqinfo(BigWigFile(bigwigpair[[1]])),seqinfo(exons_tr)))
	
	trseqinfo <- Seqinfo(seqnames=names(exons_tr),seqlengths=as.vector(sum(width(exons_tr))))

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



genes<-read_compressed_gfile(annofile,'gene')
start_codons<-read_compressed_gfile(annofile,'start_codon')
transcripts<-read_compressed_gfile(annofile,'transcript')
exons<-read_compressed_gfile(annofile,'exon')
names(exons)<-exons$transcript_id
exons<-makeTxDbFromGRanges(exons)%>%exonsBy(use.names=TRUE)


#now get the relevant genes for our lincRNAs
stopifnot(all(lncRNA_table$gene_id %in% genes$gene_id))
lincgenes <- subset(genes, gene_id %in% lncRNA_table$gene_id)




allbigwigs<-Sys.glob('pipeline/bigwigs/*/*/*chr.bw')


bigwigpairlist <- allbigwigs%>%
	data_frame(file=.)%>%
	mutate(base=file%>%basename%>%str_replace('pos|neg',''))%>%
	mutate(strand=file%>%basename%>%str_extract('pos|neg'))%>%
	mutate(strand = case_when(
		strand=='pos'~'+',
		strand=='neg'~'-'
	)) %>% 
	arrange(desc(strand))%>%
	{split(.,.$base)}%>%{map(.,~split(.$file,.$strand))%>%map(rev)}

stopifnot(bigwigpairlist%>%names%>%n_distinct%>%`>`(3))
stopifnot(bigwigpairlist%>%.[[1]]%>%names%>%`==`(c('+','-')))

#############Now individual Plots

#TODO this logic can be split up a bit to work on bam files say


linctranscripts <- transcripts %>% subset(gene_id %in% lincgenes$gene_id)

ranges2plot <- c(
	ere_gr%>%setNames(.,.$ID)%>%.[,NULL]%>%split(.,names(.)),
	exons[linctranscripts$transcript_id]%>%unlist%>%.[,NULL]%>%split(.,names(.))
)

ribofilenames <- names(bigwigpairlist)


extract_id <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = sampleids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}
sampleids <- fread('pipeline/sample_parameter.csv')[[1]]
extract_id(ribofilenames,sampleids)

ribofilenames%>%str_extract(regex('(?<=star_).*(?=DAC|ctrl.*)',ignore.case=T))


rangename='ENST00000374922.7'
ribofilename=ribofilenames[[1]]


for (ribofilename in ribofilenames[1]){

	riboprofdata <- get_riboprofdata(ranges2plot['ENST00000374922.7'],bigwigpairlist[[ribofilename]])
	ribofileid <- extract_id(ribofilename,sampleids)


	riboprofdata%>%group_by(tid)%>%summarise(sum(score))

	# for (rangename in names(ranges2plot)){
	for (rangename in 'ENST00000374922.7'){



		range = ranges2plot[[rangename]]
		rangename%<>%str_replace_all('[\\[\\]]','')

		plottitle <- str_interp("Riboseq Read Profile for:\n${rangename} = ${range}\nfile:${ribofileid}")		

		rangeriboprofdata <- riboprofdata %>%dplyr::filter(tid == rangename)
		if(nrow(rangeriboprofdata)<4){
			message(str_interp('skipping ${plottitle}'))
			next
		}
		transprofplot <-
			rangeriboprofdata %>%
			ggplot(aes(fill=frame,color=frame,x=pos,y=score))+
			geom_bar(stat='identity')+
			coord_cartesian(xlim=c(0,sum(width(ranges2plot[rangename]))))+
			ggtitle(plottitle)+
			theme_bw()
        
		rangeseq <- range%>%getSeq(x=FaFile('pipeline/my_hg38.fa'))%>%Reduce(f=c)
		startpos <- rangeseq%>%str_locate('ATG')%>%.[,1]%>%{data_frame(pos=.)}%>%mutate(type='start_codon',frame = pos %%3)
		stoppos <- rangeseq%>%str_locate('TAG|TAA|TGA')%>%.[,1]%>%{data_frame(pos=.)}%>%mutate(type='start_codon',frame = pos %%3)
		


		#name the plot file
		transprofplotfile <- str_interp('plots/loci_riboprofiles/${rangename}/${ribofileid}/${rangename}_${ribofileid}_prof.svg')
#		transprofplotfile%<>%str_replace_all('[^a-zA-z0-9_\\./]+','')
		transprofplotfile%>%dirname%>%dir.create(rec=TRUE,showWarnings=FALSE)
		#print the plot
		
		svglite(transprofplotfile);multiplot(transprofplot,layout=mylayout);dev.off()

		message(transprofplotfile%>%normalizePath)
	}
}


