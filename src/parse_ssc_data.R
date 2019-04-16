##

#Read in the scrnaseq data
scrnaseq <- fread('../ext_data/lognormalized_sscrnaseq.csv')

stopifnot(scrnaseq%>%colnames%>%head(3)%>%`==`(c("ensemblID", "geneName", "in_list")))
scrnaseqmat <- scrnaseq[,-c(1:3)]%>%as.matrix%>%set_colnames(NULL)
scrnaseqinf <- scrnaseq[,1:3]
#check if the number I get matches the manuscript


fracdetected<-scrnaseqmat%>%`>`(0)%>%apply(1,mean)%>%setNames(scrnaseqinf$ensemblID)%>%enframe('ensemblID','detected')
#looks like fraction is correct
stopifnot(fracdetected%>%filter(geneName=='LINC00520')%>%pluck('detected')%>%between(0.72,0.79))


#Now, save the object
dir.create('../data')
fracdetected%>%saveRDS('../data/fracdetected.rds')

