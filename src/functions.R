
library(memoise)
library(stringr)

###memoise
mycache=memoise::cache_filesystem(here::here('R_cache'))
mymemoise <- function(f) if(!is.memoised(f)) memoise(f,cache=mycache) else f


#' safe_left_join
#' @description left join that fails if a row in x is either duplicated or
#'   unmatched.
#' @param x table to join
#' @param y table to join
#' @param by a character vector of column names to join by.
#' @param verbose Default is TRUE.
#' @export


safe_left_join = function (x, y, by = NULL, verbose = TRUE) {
  rows_start = nrow(x)

  if (is.null(by)) {
    by = intersect(names(x), names(y))
  } else {
    by = as.character(by)
  }

  y[["..1.."]] = 1
  x = left_join(x, y, by)

  if (nrow(x) > rows_start) {
    stop("Rows have been duplicated in 'safe' left join")
  }

  if (any(ind <- is.na(x[["..1.."]]))) {
    sample = sample(which(ind), min(10, sum(ind)))
    examples = distinct(x[sample, by, drop = FALSE])
    if (verbose) print(examples)
    stop(sprintf("Failed to match %d rows in x.", sum(ind)))
  }

  x[["..1.."]] = NULL

  x

}

counts_to_tpm <- function(counts, featureLength, meanFragmentLength, featureName){
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
   featureLength - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  ###  featureNameAdjusted <- featureName[idx]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
   denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the column names from the original matrix add feature names to rownames (or extra column?)
  colnames(tpm) <- colnames(counts)
  ###rownames(tpm) <- featureNameAdjusted
  return(tpm)
}

strandis_grl <- function(grl,selstr){
  strandchars <- c('+','-','*')
  stopifnot(is(grl,'GRangesList'))
  stopifnot(selstr %in% c('+','-'))
  strtab<-table(strand(grl))
  strtab <- strtab[,strandchars,drop=FALSE]!=0
  #ensure our list has only one strand each
  stopifnot(all(rowSums(strtab[,1:3,drop=F])==1))
  #get the logical vector
  grl_isstr <- strtab[,selstr]
  #return logical vector
  grl_isstr
}

get_orf_utrs <- function(orfsgenome,exons){
  trlens <- sum(width(exons))
  
  orf_utrs <- 
    names(orfsgenome)%>%
    # paste(sep='_',names(exons)[2],'100','200')  %>%
    unique%>%
    str_split_fixed('_',3)%>%
    as.data.frame%>%
    set_colnames(c('seqnames','start','end'))%>%
    as('GRanges')%>%
    resize(1,'end')%>%
    strandshift(1)%>%
    #get them as transcript ranges
    {end(.) <- trlens[as.character(.@seqnames)];.}%>%
    setNames(.,str_replace_all(as.character(.),'-|:','_'))

  orf_utrs_genome<-orf_utrs%>%pmapFromTranscripts(exons[seqnames(orf_utrs)])%>%.[width(.)>0]
  names(orf_utrs_genome) <- names(orf_utrs)
  orf_utrs_genome
}


injectSNPsHeader <- function(orfsgenome,vcfgr,genome,snpstringsep=','){
  #this stops e.g. CTG getting translated to M
  mygeneticcode<- GENETIC_CODE
  attr(mygeneticcode,'alt_init_codons')<-c('ATG')[0]

  stopifnot(is(orfsgenome,'GRangesList'))
  allorfnames <- names(orfsgenome)
  stopifnot(is(vcfgr,'GRanges'))

  stopifnot(all(nchar(vcfgr$ALT)==1 & (width(vcfgr)==1)))
  vcfgr$ALT<-DNAStringSet(vcfgr$ALT)
  
  mutsonorfs <- mapToTranscripts(vcfgr,orfsgenome)
  orfsmutted <- orfsgenome[mutsonorfs$transcriptsHits]

  #make sure our grlist is properly sorted and get the sequences
  orfsmutted%<>%sort
  orfsmutted[strandis_grl(  orfsmutted,'-')]%<>%lapply(rev)%>%GRangesList
  orfdnaseq <- orfsmutted%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))

  mutsonorfs$ALT = ifelse(strand(mutsonorfs)=='+',vcfgr$ALT[mutsonorfs$xHits],reverseComplement(vcfgr$ALT[mutsonorfs$xHits]))
  msnames <- as.character(seqnames(mutsonorfs))
  
  #now modify the dna strings
  stopifnot(all(msnames %in% names(orfdnaseq)))
# nonindelinds<-which(!mutsonorfs_isindel)
  mutlocs <- start(mutsonorfs)
  startaainds <- ceiling(mutlocs/3)
  startbpinds <- ((startaainds-1)*3)+1
  codseqs<-subseq(DNAStringSet(orfdnaseq[msnames]),startbpinds,width=3)
  refaas <- translate(codseqs,mygeneticcode)


  subseq(codseqs,(((mutlocs)-1)%%3)+1,width=1) <- mutsonorfs$ALT
  altaas <- translate(codseqs,mygeneticcode)
  synon= refaas==altaas
  mean(synon)



  orfmodstrings <- data_frame(
    orfname=names(orfsgenome)[mutsonorfs$transcriptsHits],
    snpname=names(mutsonorfs),
    snpstring=paste0(refaas,startaainds,altaas))

  orfmodstrings <- orfmodstrings[!synon,]
  #collapse
  names(mutsonorfs) %<>% paste0(.,'_',start(mutsonorfs)%%3)
  orfmodstrings <- orfmodstrings%>%group_by(orfname)%>%summarise(snpstring=paste0(paste0(snpname,':',snpstring),collapse=snpstringsep))

  modfastaheaders <- data_frame(orfname=allorfnames)%>%
    left_join(orfmodstrings)%>%
    replace_na(replace=list(snpstring=''))%>%
    do.call(partial(paste,sep="|"),.)%>%
    paste0('|')

  stopifnot(identical(modfastaheaders%>%str_split_fixed('\\|',n=3)%>%.[,1],allorfnames))

  modfastaheaders
}

# # vcfgrbak<-vcfgr
# vcfgr<-
# vcfgr$ALT<-'G'

# injectSNPsHeader(orfsgenome,orfsgenome[[1]][1]%>%resize(4)%>%resize(1,'end')%>%{.$ALT<-'G';names(.)<-'foo';.},genome)

# negtestorf<-negorfsgenome[head(which(elementNROWS(negorfsgenome)==2),1)]


getSeq_grlist <- function(grlist,x,...){
  #make sure our grlist is properly sorted and get the sequences
  nms <- names(grlist)
  grlist%<>%sort
  grlist[strandis_grl(  grlist,'-')]%<>%lapply(rev)%>%GRangesList
  grlistdnaseq <- grlist%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
  DNAStringSet(grlistdnaseq[nms])
}



injectIndelsHeader <- function(orfsgenome,vcfgr,genome,exons,snpstringsep=','){
 #this stops e.g. CTG getting translated to M
  mygeneticcode<- GENETIC_CODE
  attr(mygeneticcode,'alt_init_codons')<-c('ATG')[0]

  stopifnot(is(orfsgenome,'GRangesList'))
  allorfnames <- names(orfsgenome)
  stopifnot(is(vcfgr,'GRanges'))

  stopifnot(!any(nchar(vcfgr$ALT)==1 & (width(vcfgr)==1)))
  vcfgr$ALT<-DNAStringSet(vcfgr$ALT)
  
  mutsonorfs <- mapToTranscripts(vcfgr,orfsgenome)
  if(length(mutsonorfs)==0){
    nullout <-
      data_frame(orfname=allorfnames)%>%
        mutate(snpstring='')%>%
        do.call(partial(paste,sep="|"),.)%>%
        paste0('|')
    return(nullout)

  }
  orfsmutted <- orfsgenome[mutsonorfs$transcriptsHits]

  orfdnaseq <- getSeq_grlist(orfsmutted)


  mutsonorfs$ALT = ifelse(strand(mutsonorfs)=='+',vcfgr$ALT[mutsonorfs$xHits],reverseComplement(vcfgr$ALT[mutsonorfs$xHits]))
  msnames <- as.character(seqnames(mutsonorfs))
  
  #now modify the dna strings
  stopifnot(all(msnames %in% names(orfdnaseq)))
  
  mutlocs <- start(mutsonorfs)
  mutends <- end(mutsonorfs)
  
  isstart <- mutlocs < 4
  prevaainds <- ceiling(mutlocs/3) - 1+ isstart
  lastaainds <- ceiling(mutends/3)
  startbpinds <- ((prevaainds-1)*3)+1
  lastbpinds <- ((lastaainds-1)*3)+3
  #browser()

  dnaseqs<-subseq(DNAStringSet(orfdnaseq[msnames]),startbpinds,lastbpinds)
  refaas <- translate(dnaseqs,mygeneticcode)
  #insert the mutation sequence
  subseq(dnaseqs,((mutlocs%%3))+3-(3*isstart),width=width(mutsonorfs)) <- mutsonorfs$ALT
  frameshift<- nchar(dnaseqs)%%3 !=0

  if(any(frameshift)){
    frameshiftutrs <- get_orf_utrs(orfsgenome[msnames[frameshift]],exons)
    frameshiftutrseqs <- DNAStringSet(getSeq_grlist(frameshiftutrs,genome))
  }else{
    frameshiftutrseqs <- letters[0]
  }
  dnaseqs[frameshift] <- paste0(
    dnaseqs[frameshift],
    subseq(DNAStringSet(orfdnaseq[msnames][frameshift]),lastbpinds[frameshift]+1,-1),
    frameshiftutrseqs
  )

  #translate after trimming down to mod 3
  dnaseqs <- DNAStringSet(substr(dnaseqs,1,(nchar(dnaseqs)%>%{floor(./3)*3})))
  altaas <- translate(dnaseqs,mygeneticcode)

  #clip after first stop codon
  aastringends <- altaas%>%str_locate('\\*')%>%.[,1]
  aastringends[is.na(aastringends)] <- nchar(altaas)[is.na(aastringends)]
  altaas <- substr(altaas,1,aastringends)

  orfmodstrings <- data_frame(
    orfname=names(orfsgenome)[mutsonorfs$transcriptsHits],
    snpname=names(mutsonorfs),
    snpstring=paste0(refaas,prevaainds,altaas))

  orfmodstrings <- orfmodstrings[refaas!=altaas,]
  
  #collapse
  names(mutsonorfs) %<>% paste0(.,'_',start(mutsonorfs)%%3)
  orfmodstrings <- orfmodstrings%>%group_by(orfname)%>%summarise(snpstring=paste0(paste0(snpname,':',snpstring),collapse=snpstringsep))

  modfastaheaders <- data_frame(orfname=allorfnames)%>%
    left_join(orfmodstrings)%>%
    replace_na(replace=list(snpstring=''))%>%
    do.call(partial(paste,sep="|"),.)%>%
    paste0('|')

  stopifnot(identical(modfastaheaders%>%str_split_fixed('\\|',n=3)%>%.[,1],allorfnames))

  modfastaheaders
}


# snpheaders<-injectIndelsHeader(negtestorf,
#   c(
#     negtestorf[[1]][1]%>%resize(12,'end')%>%resize(2,'start')%>%{.$ALT<-'C';names(.)<-'foo';.},
#     negtestorf[[1]][2]%>%resize(8,'start')%>%resize(2,'end')%>%{.$ALT<-'AAATGATGATG';names(.)<-'foo';.},
#     NULL  
#   ),

#   genome,exons)

# indelheaders<-injectSNPsHeader(negtestorf,negtestorf[[1]][2]%>%resize(4,'start')%>%resize(1,'end')%>%{.$ALT<-'C';names(.)<-'foo';.},genome)

# full_join(
#   snpheaders%>%str_split_fixed('\\|',3)%>%.[,-3,drop=F]%>%as.data.frame,
# indelheaders%>%str_split_fixed('\\|',3)%>%.[,-3,drop=F]%>%as.data.frame
# )


injectSNPsIndels <- function(orfsgenome,vcfgr,genome,snpstringsep=','){
  

  stopifnot(is(orfsgenome,'GRanges'))
  stopifnot(is(vcfgr,'GRanges'))
  vcfgr_snps <- vcfgr[nchar(vcfgr$ALT)==1]
  vcfgr$ALT<-DNAStringSet(vcfgr$ALT)
  orfsgenome <- orfsgenome%>%split(.,names(.))
  allorfnames <- names(orfsgenome)
  mutsonorfs <- mapToTranscripts(vcfgr,orfsgenome)
  orfsmutted <- orfsgenome[mutsonorfs$transcriptsHits]
  # posorfdnaseq <- orfsgenome[unique(seqnames(mutsonorfs))]%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
  posorfdnaseq <- orfsmutted%>%lapply(.%>%getSeq(x=genome)%>%Reduce(f=c))
  mutsonorfs$ALT = vcfgr$ALT[mutsonorfs$xHits]
  msnames <- as.character(seqnames(mutsonorfs))
  #now modify the dna strings
  stopifnot(all(msnames %in% names(posorfdnaseq)))
# nonindelinds<-which(!mutsonorfs_isindel)
  mutlocs <- start(mutsonorfs)
  startaainds <- ceiling(mutlocs/3)
  startbpinds <- ((startaainds-1)*3)+1
  codseqs<-subseq(DNAStringSet(posorfdnaseq[msnames]),startbpinds,width=3)
  refaas <- translate(codseqs)
  subseq(codseqs,(((mutlocs)-1)%%3)+1,width=1) <- mutsonorfs$ALT
  altaas <- translate(codseqs)
  orfmodstrings <- data_frame(
    orfname=names(orfsgenome)[mutsonorfs$transcriptsHits],
    snpstring=paste0(refaas,startaainds,altaas))
  #collapse
  orfmodstrings <- orfmodstrings%>%group_by(orfname)%>%summarise(snpstring=paste0(snpstring,collapse=snpstringsep))

  data_frame(orfname=allorfnames)%>%left_join(orfmodstrings)%>%replace_na(replace=list(snpstring=''))%>%head%>%do.call(paste0,.)


  for(i in which(!mutsonorfs_isindel)){
    mutloc<-start(mutsonorfs[i])
    startaaind <- ceiling(mutloc/3)
    startbpind <- ((startaaind-1)*3)+1
    codseq <- dnaseq[startbpind:(startbpind+2)]
    refaa <- translate(codseq)
    codseq[(((mutloc)-1)%%3)+1] <- mutsonorfs$ALT[[i]]
    altaa <- translate(codseq)
    paste0(refaa,startaaind,altaa)

    names(posorfdnaseq)[msnames[[i]]]

    aaref <- translate()

    dnaseq<-posorfdnaseq[[msnames[i]]]
    aamut <- translate(dnaseq[startbpind:endbpind])


    protseq
    #codonends ceiling(1:7/3)*3
    ceiling(start(mutsonorfs)[i]/3)*3+1
    floor(start(mutsonorfs)[i]/3)+1
    mutcodons <- dnaseq[codstart:codend]
    mutinds <- start(mutsonorfs[i]):end(mutsonorfs[i])
    
    posorfdnaseq[[msnames[i]]][]
    posorfdnaseq[[msnames[i]]][start(mutsonorfs[i]):end(mutsonorfs[i])]<-mutsonorfs$ALT[[i]]
  }

  for(i in which(mutsonorfs_isindel)){
    seq <- posorfdnaseq[[msnames[i]]]
    posorfdnaseq[[msnames[i]]] <- c(
      seq[1:(start(mutsonorfs[i])-1)],
      mutsonorfs[i]$ALT[[1]],
      seq[(end(mutsonorfs[i])+1):length(seq)]
    )   
  }

  negorfs <- orfsgenome%>%subset(strand=='-')%>%split(.,names(.))
  mutsonorfs <- mapToTranscripts(vcfgr,negorfs)

  # negorfdnaseq <- negorfs[unique(seqnames(mutsonorfs))]%>%lapply(.%>%getSeq(x=genome)%>%rev%>%Reduce(f=c))
  negorfdnaseq <- negorfs%>%lapply(.%>%getSeq(x=genome)%>%rev%>%Reduce(f=c))
  mutsonorfs$ALT = reverseComplement(DNAStringSet(vcfgr$ALT[mutsonorfs$xHits]))

  msnames <- as.character(seqnames(mutsonorfs))
  stopifnot(all(msnames %in% names(negorfdnaseq)))
  #now modify the dna strings
  mutsonorfs_isindel = ! ((nchar(mutsonorfs$ALT)==1) & (width(mutsonorfs)==1))
  for(i in which(!mutsonorfs_isindel)) negorfdnaseq[[msnames[i]]][start(mutsonorfs[i]):end(mutsonorfs[i])]<-mutsonorfs$ALT[[i]]
  for(i in which(mutsonorfs_isindel)){
    seq <- negorfdnaseq[[msnames[i]]]
    negorfdnaseq[[msnames[i]]] <- c(
      seq[1:(start(mutsonorfs[i])-1)],
      mutsonorfs[i]$ALT[[1]],
      seq[(end(mutsonorfs[i])+1):length(seq)]
    )   
  }
  return(list(
      seqs = DNAStringSet(c(posorfdnaseq,negorfdnaseq)[unique(names(orfsgenome))])
    )
  )
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
load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


extract_oneof <- function(strings,ids){
  
  matchlist <- map(strings,~str_extract(pattern = ids,string = .))
  
  matchnum <- matchlist%>%map(~sum(!is.na(.)))
  
  stopifnot(all(matchnum < 2 ))
  stopifnot(all(matchnum > 0))

  matches <- matchlist%>%map_chr(keep,Negate(is.na))

  matches
}

DT2GR = function(dt,seqinf=si,checksi=TRUE){

  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(c('seqnames','start')%in%colnames(dt))
  stopifnot(c('width')%in%colnames(dt)|c('end')%in%colnames(dt))
  if(checksi){stopifnot(dt[['seqnames']] %in% seqlevels(seqinf))
  }else{seqinf=NULL}
  
  hasend=FALSE
  haswidth=FALSE

  if('end' %in% colnames(dt) ){
    stopifnot (dt[['end']] %>% `>`(0) %>%all)
    hasend=TRUE
  }
  if('width' %in% colnames(dt) ){
    stopifnot (dt[['width']] %>% `>`(0) %>%all)
    haswidth=TRUE
  }
  
  stopifnot(dt[['start']] %>% is.numeric)
  stopifnot(hasend|haswidth )
  
  if(haswidth & ! hasend ){
    dt[['end']]  = dt[['start']]+dt[['width']]-1 
  } 
  if(hasend ){

  } 

  #strand
  if(! 'strand' %in% colnames(dt)){
    dt[['strand']] ='*'
  }

  stopifnot(dt[['strand']] %in% c('+','-','*'))
  



  mdatcols = colnames(dt) %>% setdiff(c('seqnames','start','width','strand','end')) 
  #create basic granges object
  if(checksi){
    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']],seqinfo=seqinf)
  }else{    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])}

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

    stopifnot(all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')))

  gr
}


mergeGR2DT = function(mergegr){
  grcols = mergegr%>%vapply(is,TRUE,'GRanges')%>%which
  mergedt=cbind(
    mergegr[[grcols[[1]]]]%>%GR2DT,
    mergegr[[grcols[[2]]]]%>%GR2DT
  )
  cnames = colnames(mergedt)
  cnames[duplicated(cnames)]%<>%paste0('.1')
  setnames(mergedt,cnames)
  mergedt
}

GR2DT = function(gr){
  if(is.data.table(gr)) {
    warning('already a data table')
    return(gr)
  }
  #all columns must be vectors
  for(i in seq_len(ncol(mcols(gr)))){
    mcols(gr)[[i]] = as.vector(mcols(gr)[[i]])
  }

  dt = as.data.frame(gr,row.names=NULL,optional=FALSE)%>%as.data.table
  dt$strand= as.character(dt$strand)
  # setkey(dt,seqnames,strand,start,end)

  stopifnot( all(colnames( mcols(gr)) %in% colnames(dt) )  )
  dt
}



load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}

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

strandshift<-function(gr,dist){ 
  sign=ifelse(strand(gr)=='-',-1,1)
  GenomicRanges::shift(gr,dist*sign)
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

