
suppressMessages(library("purrr"))
suppressMessages(library("data.table"))
suppressMessages(library("assertthat"))
suppressMessages(library("magrittr"))
suppressMessages(library("Rsamtools"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("multitaper"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("foreach"))
suppressMessages(library("doMC"))

source('../src/functions.R')


#This script is intended to get TrP values for all the ORFs in a given

library(GenomicFeatures)

GTF_annotation <- load_objs('riboqc/my_gencode.v24lift37.annotation.annoout')[[1]]

###Objective here is to seperate some of the mess above - not sure best way to get that psite vector... probably using 
fasta = '../export/ORFs/Fastas/OD5P/OD5P.fasta'
sentorfs <- system(str_interp('grep -Poe "ENST[0-9_]+" ${fasta}'),intern=T)
allmyorfob <- sentorfs%>%str_replace('_',':')%>%str_replace('_','-')%>%GRanges
names(allmyorfob)<-sentorfs


allresfiles <- Sys.glob('riboqc/data/*/*to_SaTAnn*')

allresobs<-allresfiles%>%mclapply(load_objs)

names(allresobs)<-basename(dirname(allresfiles))

#Get TrP values for all 
for(libname in names(allresobs)){try({
     #get the psite object
     to_SaTAnn <- allresobs[[libname]]$to_SaTAnn

     psites <- to_SaTAnn$P_sites_all

     junc_psites <- to_SaTAnn$junctions
     txs <- as.character(unique(allmyorfob@seqnames))
     #First map to transcripts
     mappedpsites <- mapToTranscripts(psites,GTF_annotation$exons_txs[txs])
     mappedpsites$score <- psites$score[mappedpsites$xHits]
     mcols(mappedpsites)%<>%.[,'score',drop=F]
     #now map to ORFs.
     orfpsites <- mapToTranscripts(mappedpsites,split(allmyorfob,names(allmyorfob)))
     orfpsites$score <- mappedpsites$score[orfpsites$xHits]
     orfpsitescov <- orfpsites%>%coverage(weight='score')
     #vectors with ORF's psites
     orfpsitescov%<>%lapply(as.vector)

     #Get the 
     bw=12
     nw=12
     tapers=24

     #Now get the TrP value for all the ORFs
     # orftrpres <- orfpsitescov%>%head(10)%>%lapply(tapers=24,bw=12,FUN=safely(function(ps,tapers,bw){
     message('get TrP')
     orftrpres <- orfpsitescov%>%mclapply(mc.cores=detectCores(),tapers=24,bw=12,FUN=safely(function(ps,tapers,bw){
               # browser()
               # ps <- as.vector(trcov[[tr]][start:end])
               testres <-  take_Fvals_spect(x = ps , slepians_values = dpss(n=max(50,length(ps)),k=tapers,nw=bw),n_tapers=tapers,time_bw=bw)
               pval<-pf(q=testres[1],df1=2,df2=(2*24)-2,lower.tail=F)
               list(TrP=testres[2],pval=pval,P_sites=sum(ps))
          }))
     message('saving')
     orftrpres%>%saveRDS(paste0('SaTAnn/alltrps',libname,basename(fasta),'.rds'))

     message(libname)

})}

#Now get the TrPs for all our detected ORFs in OD5P
