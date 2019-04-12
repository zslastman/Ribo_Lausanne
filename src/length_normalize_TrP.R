satanngtfs <- Sys.glob('SaTAnn/*/SaTAnn_Detected_ORFs.gtf')
satannfastas <- Sys.glob('groupedsatan/*.fasta')

satanfastaorfs <- map(satannfastas,.%>%readLines%>%str_subset('>')%>%str_extract(regex('(?<=>).*?(?=\\|)')))
names(satanfastaorfs) <- satannfastas%>%basename%>%str_replace('.fasta','')

source('../src/functions.R')
library(data.table)
library(GenomicRanges)
library(tidyverse)
slice<-dplyr::slice

for(satanngtffile in satanngtfs){
	
	message(satanngtffile)
	satanngtf <- rtracklayer::import(satanngtffile)


	satanngtf%>%
		subset(!is.na(Iso_pct_P_sites))%>%
		GR2DT%>%
		group_by(gene_id)%>%
		slice(which.max(Iso_pct_TrP)) %>%
		select(ORF_id,TrP_pNpM,P_sites)%>%
		write_tsv(str_replace(satanngtffile,'.gtf$','.TrPvals.tsv'))

}

'ENST00000443026_1427_1546'

Sys.glob('SaTAnn/*/*.TrPvals.tsv')%>%
	setNames(.,.)%>%map_df(fread,.id='file')%>%
	mutate(sample=basename(dirname(file)))%>%
	select(sample,gene_id,ORF_id,TrP_pNpM)%>%
	spread(sample,TrP_pNpM)%>%
	write_tsv('SaTAnn/all_TrPs.tsv')

normalizePath('SaTAnn/all_TrPs.tsv')

allorfexprdf<-fread('SaTAnn/all_TrPs.tsv')
 orfids <- allorfexprdf%>%.$ORF_id 
 orftrs <- orfids %>%str_extract(regex('.*(?=_\\d+_\\d+)'))

library(rtracklayer)
exons <- read_compressed_gfile('my_gencode.v24lift37.annotation.gtf','exon')

orfexons <- exons%>%split(.,.$transcript_id)%>%.[orftrs]

bamfiles <- Sys.glob('star/data/*/*.bam')
bamfile <- bamfiles[1]
library(GenomicAlignments)


orftrovs <- mclapply(mc.cores = 10,bamfiles,function(bamfile){
	message('.')
	summarizeOverlaps(orfexons,bamfile)
})
assay(orftrovs[[1]])%>%rownames%>%setdiff(names(orfexons),.)

orftrcounts<-orftrovs[T]%>%setNames(basename(dirname(bamfiles[T])))%>%lapply(.%>%assay%>%as.data.frame%>%set_colnames('count')%>%rownames_to_column('tr'))%>%bind_rows(.id='sample')%>%
	mutate(tr = str_replace(tr,'\\.\\d+',''))%>%
	group_by(sample,tr)%>%
	summarise(count=sum(count,na.rm=T))%>%
	spread(sample,count)





satannobfiles <- Sys.glob('SaTAnn/*/SaTAnn_Final_ORFs_files')
#get raw psites from the satann objects
rawpsitedf <- satannobfiles %>% setNames(.,.)%>%map_df(.id='file',~ load_objs(.)$ORFs_tx%>%{mcols(.[intersect(orfids,names(.))])[,c('P_sites_raw','ORF_id_tr')]}%>%as.data.frame)
rawpsitedf%<>%mutate(sample = basename(dirname(file)))
rawpsitedf$file<-NULL

rawpsitedf%<>%spread(sample,P_sites_raw)

for(cell_line in c('OD5P','ONVC','OMM')){
	stopifnot(cell_line %in% colnames(allorfexprdf))
	rawpsitedf[[cell_line]] <- apply(rawpsitedf%>%select(matches(cell_line)),1,sum)
	
}