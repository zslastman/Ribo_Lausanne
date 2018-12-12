
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
library(multitaper)
register(MulticoreParam(8))

#for a set of exonsspectestdf
library(rtracklayer)
library(GenomicFeatures)
library(multitaper)
library(ORFik)

message('...done')

#load arguments
args <- c(
	juncgtf = 'junctions/ref_metadata.filt_REF.C3N-02289.filt_L1.gtf',
	pluspsites = 'riboqc/data/OD5P_05_uM_DAC_1/_P_sites_plus.bw',
	negpsites = 'riboqc/data/OD5P_05_uM_DAC_1/_P_sites_plus.bw',
	outputfile='junctionspec/ref_metadata.filt_REF.C3N-02289.filt_L1.tsv'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

outputfolder<-dirname(outputfile)
outputfolder <- paste0(outputfolder,'/')
outputfolder%>%dir.create(showWarn=F,rec=TRUE)





n <- 300
x <- rnorm(n)


#this function takes in a numeric vector x (psites) and parameters for the multitaper FFT
#and precalculated slepian values, and returns the maximum F value (for pvalue calcs)
#and the magnitudes of the FFT at the closest frequency to 0.33
take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
     if(length(x)<25){
          remain<-50-length(x)
          x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
     }
     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
          
     closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
     
     freq_max_3nt<-resSpec1$freq[closestfreqind]
     Fmax_3nt<-resSpec1$mtm$Ftest[closestfreqind]
     spect_3nt<-resSpec1$spec[closestfreqind]
     return(c(Fmax_3nt,spect_3nt))
     
}
#this takes in a numeric vector of psite values (i.e. the number of psites at a given bp),
#and possibly the number of slepians and bandwidth for the multitaper test, and returnst the
#pvalue of the FFT.
ftestvect<-function(psit,k=24,bw=12){
	psit <- as.vector(psit)
	if(sum(psit!=0) < 3) return(NA)
	if(sum(psit) < 10) return(NA)
		
	slepians_values<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
	vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = slepians_values)
	pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
	return(c(pval))
}





#get our exons

exons <- import(juncgtf)
exons<-exons%>%subset(type=='exon')
exons<-exons[exons$ID%>%duplicated%>%not,]
exons<-exons[order(exons$ID)]
#get fp and tp exons#verify we have two exons in our sets
stopifnot(exons$ID%>%table%>%is_in(1)%>%all)
stopifnot(exons$transcript_id%>%table%>%is_in(2)%>%all)

isfpexon <- exons$ID%>%str_detect('_1$')
stopifnot(mean(isfpexon)==0.5)

fpexons <- exons[isfpexon]
tpexons <- exons[!isfpexon]



# 'grep start_codon my_gencode.v24lift37.annotation.gtf > startcods.gtf'%>%system
# 'grep stop_codon my_gencode.v24lift37.annotation.gtf > stopcods.gtf'%>%system
startcods<-rtracklayer::import('startcods.gtf')
stopcods<-rtracklayer::import('stopcods.gtf')

# startcods<-fread(('grep start_codon my_gencode.v24lift37.annotation.gtf'))


#function to clip our exons at the fp or tp end with a second set of grs
clipGR <- function(gr1,gr2,threeprime=FALSE){
	
	gr1$.__matchid = seq_along(gr1)
	ovdf <-mergeByOverlaps(gr1,gr2)
	ovdf = ovdf[str_replace(ovdf$gr1$gene_id,'\\.\\d+$','') == ovdf$gr2$gene_id,]

	gr1
	stopifnot(all(strand(ovdf$gr1)==strand(ovdf$gr2)))

	ispos = as.vector(strand(ovdf$gr1)%in% c('+','*'))
	if(threeprime) ispos = ! ispos#if we're doing 3' then just reverse this


	start(ovdf$gr1) <- as.vector(ifelse(ispos,start(ovdf$gr2),start(ovdf$gr1)))
	end(ovdf$gr1) <- as.vector(ifelse(!ispos,end(ovdf$gr2),end(ovdf$gr1)))

	gr1[ovdf$gr1$.__matchid] <- ovdf$gr1
	
	gr1$.__matchid <- NULL

	gr1
}

clippedfpexons = clipGR (fpexons , startcods,threeprime=FALSE)
clippedtpexons = clipGR (fpexons , startcods,threeprime=TRUE)
clippedexons <- c(clippedfpexons,clippedtpexons) 
#use any overlapping exons to hit them


#

#load up our psites as bigwig objects

poscovgr <- import(pluspsites,which=subset(exons,strand=='+'))
negcovgr <- import(negpsites,which=subset(exons,strand=='-'))

GRViews<-function(rle,gr){
  stopifnot(is(rle,'RleList'))
  stopifnot(is(gr,'GenomicRanges'))

  gr = sort(gr) #sort
  
  v=Views(rle[seqlevels(gr)],as( gr[,NULL],'IntegerRangesList') )[seqlevels(gr)] #get views
  
  v=v[sapply(v,function(x){length(x)!=0})] #get a vector, leaving out empty chromosomes
  
  return(v)
}

ftestvect<-function(psit,k=24,bw=12){
	psit <- as.vector(psit)
	if(sum(psit!=0) < 3) return(NA)
	if(sum(psit) < 10) return(NA)
		
	sl<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
	vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = sl)
	pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
	return(c(pval))
}

names(exons) <- exons$ID
stopifnot(!anyDuplicated(exons$ID))

pexons <- subset(exons,strand=='+')
nexons <- subset(exons,strand=='-')

pviews <- GRViews(coverage(poscovgr,weight='score'),pexons)
nviews <- GRViews(coverage(negcovgr,weight='score'),nexons)


pspec <- pviews %>%
	lapply( function(x) x[sum(x)>10])%>%
	lapply(head,10)%>%map(function(x)viewApply(x,ftestvect)) 
nspec <- nviews %>%
	lapply(function(x) x[sum(x)>10])%>%
	lapply(head,10)%>%map(function(x)viewApply(x,ftestvect)) 

pval_df <-
	rbind(
		stack(unlist2(pspec))%>%dplyr::select(ID=ind,p_spec=values),
		stack(unlist2(nspec))%>%dplyr::select(ID=ind,p_spec=values)
	)%>%mutate(fp=str_detect(ID,'_1$')) %>%
	mutate(ID=str_replace(ID,'_\\d$',''))%>%
	spread(fp,p_spec)%>%
	dplyr::select(ID,fp_spec=`TRUE`,tp_spec=`FALSE`)

pval_df %>%write_tsv(outputfile)






# # #create a granges list ocntaining the paired exons as transcripts, and the individual exons as transcripts
# # exonsbytrid <- exons%>%subset(type=='exon')%>%split(.,.$transcript_id)
# # exonslist<-c(exonsbytrid,exons%>%subset(type=='exon')%>%split(.,.$ID))

# # segmentseqs<-xscat(
# # 	getSeq(clippedfpexons,x=FaFile('my_hg19.fa')),
# # 	getSeq(clippedtpexons,x=FaFile('my_hg19.fa'))
# # )

# # exonslist%>%length
# # segmentseqs%>%length

# # ORFik::findORFs(segmentseqs,longestORF=TRUE)


# # addedstarts<-'ATGCATGCATG'
# # addedstops<-'TAACTAACTAA'

# # inseq <- 'TAAGGGGGG'
# # nchar(inseq)
# # nchar(addedstarts)+nchar(inseq)


# #if a reading frame starts between 1 and 11, it begins before the exon somewhere
# #9 is the 1st reading frame,1 is the second, 5 is the third
# #thus nchar(start)+nchar(seq)+3 is the first open reading frame, +4 is the 2nd, + 8 is the third

# #get the longest open reading frame, clip the start and stop back, use resulting range to clip the vector

# test <- ORFik::findORFs(paste0(addedstarts,addedstops))
# 14-11


# seqlens <- nchar(segmentseqs)
# addstartlen <- nchar(addedstarts)
# longestorfs <- ORFik::findORFs(xscat(addedstarts,segmentseqs,addedstops),longest=TRUE)%>%unlist
# start(longestorfs)	= pmax(start(longestorfs), addstartlen+1) - addstartlen
# end(longestorfs)	= pmin(end(longestorfs), addstartlen+seqlens) - addstartlen

# # ORFik::findORFs(paste0(addedstarts,addedstops))
# # ORFik::findORFs('ATGCATGCATGCCTAACTAACTAA')
# # ORFik::findORFs('ATGCATGCATGCCCTAACTAACTAA')

# #load our psites into the exon space, mapping over teh scores at each
# trseqinfo <- Seqinfo(seqnames=names(exonslist),seqlengths=as.vector(sum(width(exonslist))))
# mappedcov <- mapToTranscripts(covgr, exonslist)
# seqinfo(mappedcov)<-trseqinfo[seqinfo(mappedcov)@seqnames]
# mappedcov$score <- covgr$score[mappedcov$xHits]
# #get our vectors of psite coverage, and then test those
# spectests<-mappedcov%>%
# 	coverage(weight='score')%>%
# 	lapply(as.vector)%>%
# 	map(ftestvect)
# #format results into a data frame
# spectestdf<-spectests%>%
# 	simplify2array%>%t%>%as.data.frame%>%
# 	rownames_to_column%>%
# 	set_colnames(c('ID','F','spec','pval'))%>%
# 	arrange(ID)
# #write the results into a file
# write_tsv(spectestdf,outfile)














# # (mapToTranscripts(testwindow,exons_tr))
# gr$score<-profilegrange$score[gr$xHits];
# # gr%>%subset(score>30)
# # profilegrange%>%subset(score>30)

# datname = bigwigpair%>%extract_id%>%.[[1]]

# scoremat <- rep(0,length(gr)*3)%>%matrix(ncol=3)
# scoremat[matrix( c(1:length(gr),(start(gr)%%3)+1 ) ,ncol=2 ) ] <- gr$score
# mcols(gr) = scoremat
# rgbvect <- c('red','green','blue')%>%sort

# if(length(gr)==0) {groupvect <- NULL}else{
# 	# groupvect <- groupvect[(start(gr)%%3)+1]
# 	groupvect <- paste0('frame ',1:3)
# } 





# k=24
# bw=12

# sl<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)

# take_Fvals_spect(nuccounts,k,bw,sl)

# take_Fvals_spect(c(rep(0,3),nuccounts,rep(0,3)),k,bw,sl)

# take_Fvals_spect



# ## Read some sample ecg data
# ecg <- read.table('http://www.indyrad.iupui.edu/public/mmiller3/sample-ecg-1kHz.txt')
# names(ecg) <- c('t','ecg')

# ecg$t <- ecg$t/1000  # convert from ms to s

# par(mfrow=c(2,2))

# ## Plot the ecg:
# plot(ecg ~ t, data=ecg, type='l', main='ECG data sampled at 1 kHz', xlab='Time [s]')

# ## Calculate fft(ecg):
# ecg$fft <- fft(ecg$ecg)

# ## Plot fft(ecg):
# #plot(ecg$fft, type='l')

# ## Plot Mod(fft(ecg)):
# plot(Mod(ecg$fft), type='l', log='y', main='FFT of ecg vs index')

# ## Find the sample period:
# delta <- ecg$t[2] - ecg$t[1]

# ## Calculate the Nyquist frequency:
# f.Nyquist <- 1 / 2 / delta

# ## Calculate the frequencies.  (Since ecg$t is in seconds, delta
# ## is in seconds, f.Nyquist is in Hz and ecg$freq is in Hz)
# ## (Note: I may be off by 1 in indexing here ????)
# ecg$freq <- f.Nyquist*c(seq(nrow(ecg)/2), -rev(seq(nrow(ecg)/2)))/(nrow(ecg)/2)

# ## Plot fft vs frequency
# plot(Mod(fft) ~ freq, data=ecg, type='l', log='y', main='FFT of ECG vs frequency', xlab='Frequency [Hz]')

# ## Now let's look at some artificial data:
# x <- seq(100000)/1000  # pretend we're sampling at 1 kHz

# ## We'll put in two frequency components, plus a dc offset
# f1 <- 5  # Hz
# f2 <- 2  # Hz
# y <- 0.1*sin(2*pi*f1*x) + sin(2*pi*f2*x) + 50
# fft.y <- fft(y)

# delta <- x[2] - x[1]
# f.Nyquist <- 1 / 2 / delta
# f <- f.Nyquist*c(seq(length(x)/2), -rev(seq(length(x)/2)))/(length(x)/2)

# par(mfrow=c(2,2))
# plot(x,y, type='l', xlim=c(0,20))
# plot(f, Mod(fft.y), type='l', log='y')

# ## Now let's zoom in and mark the points were I expect to see peaks:
# plot(f, Mod(fft.y), type='l', log='y', xlim=c(-10,10))
# rug(c(-f1, -f2, 0, f1, f2), col='red', side=3)

# dev.off()





# #this is similiar to the above, but smooshes teh functions together (not really sure
# #why they were seperate), and returns the actual mtm object.
# #this contains the power spectrum ( .$spec ) as well as Fvalues and other stuff
# #see ?multitaper::spec.mtm
# get_mmt_spec<-function(psit,n_tapers=24,time_bw=12,slepians_values){
# 	psit <- as.vector(psit)
# 	if(sum(psit!=0) < 3) return(NA)
# 	if(sum(psit) < 10) return(NA)
		
# 	slepians_values<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=n_tapers,nw=time_bw)
	
#      if(length(psit)<25){
#           remain<-50-length(psit)
#           psit<-c(rep(0,as.integer(remain/2)),psit,rep(0,remain%%2+as.integer(remain/2)))
#      }
#      if(length(psit)<1024/2){padding<-1024}
#      if(length(psit)>=1024/2){padding<-"default"}

#      resSpec1 <- spec.mtm(as.ts(psit), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
     
#      resSpec1
     
# }

# # txtplot(get_mmt_spec(rep(c(1000,0,0),20))$spec)

# txtplot(get_mmt_spec(rep(c(5,0,0),1))$spec)



# fft_specdense <-function(psit){
# 	len <- length(psit)
# 	ft <- fft(psit)
# 	ft*(Conj(ft))
# }

# resSpec1$spec


# tvect <- rep(c(3,0,0),12)
# tvects1 <- rep(c(3,0,0),12)[-1]
# tvects2 <- rep(c(3,0,0),12)[-c(1:2)]
# normvect <- rnorm(36)
# normtvect <- tvect + normvect


# spec.mtm()

