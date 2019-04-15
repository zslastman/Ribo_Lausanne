# library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)
library(parallel)

filter<-dplyr::filter
slice<-dplyr::slice
filter<-dplyr::filter
filter<-dplyr::filter

source('../src/functions.R')

extract_id <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = sampleids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}





#initiaal form when testing
#let's look at the scores for our
cell_lines<-c('OD5P','OMM','ONVC') 
satannfiles <- Sys.glob('SaTAnn/*/*Final_ORFs*')%T>%{stopifnot(length(.)>0)}
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

#first create gtf for all the 

for(i in seq_along(satannorfs)){
	rtracklayer::export(satannorfs[[i]]$ORFs_gen,paste0(dirname(satannfiles[i]),'/orfs_genomic.gtf'))
}

#now filter out the weird GRanges columns, for now, and aggregate into one data table 
all_orfs <- satannorfs%>%map(.%>%.$ORFs_tx)

for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character
orfs_dt <- all_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')
#TODO - eveyrthing in linc

n_genes_translated <- orfs_dt %>% group_by(sample)%>%dplyr::summarise(n_genes=n_distinct(gene_id))

n_genes_translated%>% filter(sample%>%str_detect('OD5P'))

readthroughs <- satannorfs%>%map_df(.%>%.$ORFs_readthroughs%>%length)%>%stack%>%set_colnames(c('n_readthroughs','sample'))
readthroughs$sample%<>%as.character
readthroughs$n_readthroughs%<>%as.numeric

#number of gene hits by gene category
genetypehits <- satannorfs%>%lapply(
	.%>%.$ORFs_tx%>%.[T,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame%>%distinct%>%group_by(gene_biotype)%>%tally%>%
	mutate(gene_biotype = ifelse(gene_biotype%>%str_detect('pseudogene'),'pseudogene',gene_biotype))%>%
	group_by(gene_biotype)%>%summarise(n=sum(n))%>%
	dplyr::filter(gene_biotype %in% c('antisense','lincRNA','protein_coding','pseudogene')))%>%
	bind_rows(.id='sample')%>%
	spread(gene_biotype,n)

#get number of genes with N_truncations
ntruncgenefrac <- satannorfs%>%map(~.$ORFs_tx%>%mcols%>%
	.[,c('gene_id','ORF_category_Tx_compatible')]%>%
	as.data.frame%>%
	group_by(gene_id)%>%
	summarise(n_trunc=any(ORF_category_Tx_compatible%in%'N_truncation') & (!any(ORF_category_Tx_compatible%in%c('ORF_annotated'))) )%>%
	.$n_trunc%>%mean)
ntruncgenefrac%<>%enframe('sample','frac_genes_ntrunc')
ntruncgenefrac$frac_genes_ntrunc%<>%unlist
#number of short ORFs detected
shorfrac <- satannorfs%>%map_dbl(~mean(width(.$ORFs_tx) < 250))%>%enframe('sample','shortORF_frac')

fp_bias <- fread('fp_bias.txt')

satann_summary_table <- genetypehits%>%left_join(readthroughs)%>%left_join(ntruncgenefrac)%>%left_join(shorfrac)%>%left_join(fp_bias)


cortest<-cor.test(satann_summary_table$fp_bias,satann_summary_table$frac_genes_ntrunc)
corpval<-cortest$p.value%>%round(5)
corest<-cortest$estimate%>%round(5)

#Plot number of 
pdf(h=4,w=6,'../plots/n_trunc_genes_vs_fp_bias.svg')
qplot(satann_summary_table$fp_bias^(-1),satann_summary_table$frac_genes_ntrunc)+
	scale_x_continuous(name="Read 5' ends AUG-250/250-500")+
	scale_y_continuous(name='Genes with only N-truncations')+
	geom_smooth(method='lm')+
	theme_minimal()+
	ggtitle(str_interp("N-truncated Genes as a function of 5' bias\nRho = ${corest}, p = ${corpval}"))
dev.off()


cortest<-cor.test(satann_summary_table$shortORF_frac,satann_summary_table$frac_genes_ntrunc)
corpval<-cortest$p.value%>%formatC(format = "e", digits = 2)
corest<-cortest$estimate%>%round(2)

pdf(h=4,w=6,'../plots/shortORF_frac.svg'%T>%{normalizePath(.)%>%message})
qplot(satann_summary_table$shortORF_frac,satann_summary_table$frac_genes_ntrunc)+
	scale_x_continuous(name="Read 5' ends AUG-250/250-500")+
	scale_y_continuous(name='Fraction of ORFs that are less than 250bp long')+
	geom_smooth(method='lm')+
	theme_minimal()+
	ggtitle(str_interp("Short ORF Detection as a function of 5' bias\nRho = ${corest}, p = ${corpval}"))
dev.off()


satann_summary_table$sample

mergesamples<-satann_summary_table$sample%>%grep(v=T,inv=T,patt='_')

pdf(h=4,w=6,'/fast_new/work/groups/ag_ohler/dharnet_m/Ribo_Lausanne/plots/fp_bias_batches.svg')
satann_summary_table%>%
	filter(!sample %in% mergesamples)%>%
	mutate(new = sample%>%str_detect('ctrl\\d+(B|L)'))%>%
	mutate(sample_prep_location=ifelse(str_detect(sample,'B$'),'Berlin','Lausanne'))%>%
	mutate(Batch=ifelse(!new,'Batch1','Batch2'))%>%
	arrange(Batch,sample_prep_location)%>%
	dplyr::select(sample,sample_prep_location,fp_bias,Batch)%>%{
			ggplot(.,aes(x=sample,fill=sample_prep_location,y=fp_bias,group=Batch))+
			scale_x_discrete(limits=rev(factor(.$sample,unique(.$sample))))+
			scale_y_continuous(name = "Read 5' ends AUG-250/250-500")+
			stat_identity(geom='bar')+theme_minimal()+theme(axis.text.x=element_text(angle=45,vjust=0.5,size=8))+coord_flip()
	}
dev.off()

satann_summary_table%>%write_tsv('satann_summary.tsv'%T>%{message(normalizePath(.))})

satann_summary_table%>%as.data.frame

satannorfs[[1]]$ORFs_tx[,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame
satannorfs[[1]]$ORFs_tx$ORF_category_Tx_compatible%>%unique




#for each group of samples,
#we want to merge ORFs that are wholly contained within the samples.

#get the nocall ranges for that group of samples
#use find

duplicated(c(1,1,2))

allgorfs<-map(satannorfs,'ORFs_gen')%>%GRangesList%>%unlist(use.names = F)
allgorfs <- allgorfs[!duplicated(paste0(names(allgorfs),as.character(allgorfs)))]
allgorfs$ORF_id_tr <- names(allgorfs)
groupsamples <- names(all_orfs)%>%str_subset('OD5P')

cell_line_i<-'OD5P'
for(cell_line_i in c('OD5P','ONVC','OMM')){

	groupsamples <- c(groupsamples%>%str_subset('^[^_]+$'),groupsamples)%>%unique

	grouporfs <- all_orfs[groupsamples]%>% map(~.[,NULL])%>% GRangesList %>% unlist %>% unique
		
	ovs <- grouporfs %>% GenomicRanges::findOverlaps(type='within') %>% as.data.frame %>% filter(!queryHits==subjectHits)

	grouporfs <- grouporfs[! 1:length(grouporfs) %in% ovs$queryHits]

	mergedorfdt <- orfs_dt%>%inner_join(names(grouporfs)%>%str_split_fixed('\\.',2)%>%as.data.frame%>%set_colnames(c('sample','ORF_id_tr')))

	allsi<-all_orfs%>%map(~.[,NULL])%>%GRangesList%>%unlist%>%seqinfo

	# outfasta <- 

	protseqs <- mergedorfdt%>%{AAStringSet(setNames(.$Protein,.$ORF_id_tr))}
	names(protseqs)<-paste(mergedorfdt$ORF_id_tr,mergedorfdt$gene_biotype,mergedorfdt$gene_id,mergedorfdt$ORF_category_Gen,mergedorfdt$ORF_category_Tx_compatible,sep="|")

	protseqs%>%writeXStringSet(str_interp('groupedsatan/${cell_line_i}.fasta'))

	mergedorfdt %>% DT2GR(allsi)%>%export(str_interp('groupedsatan/${cell_line_i}.gtf'))
	#also export record of genomic locations
	stopifnot(all(mergedorfdt$ORF_id_tr %in% allgorfs$ORF_id_tr))
	#now export genomic coordinates
	allgorfs%>%subset(ORF_id_tr %in% mergedorfdt$ORF_id_tr)%>%export(str_interp('groupedsatan/${cell_line_i}.genomic.gtf'))

}


# diffgorfs <- 	setdiff(mergedorfdt$ORF_id_tr, allgorfs$ORF_id_tr)%>%.[666]

# satannorfs[[mergedorfdt%>%subset(ORF_id_tr==diffgorfs)%>%.$sample]]$ORFs_gen[diffgorfs]


diffgorfs%in%allgorfs$ORF_id_tr



