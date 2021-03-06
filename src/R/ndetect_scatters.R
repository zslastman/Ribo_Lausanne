#Get rpkms for our genes
message('loading libraries')
library(magrittr)
library(assertthat)
library(data.table)
library(ggExtra)
suppressMessages(library(stringr))
library(rtracklayer)
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(edgeR))
suppressMessages(library(tidyverse))



source('/fast/groups/ag_ohler/work/dharnet_m/Ribo_Lausanne/src/functions.R')

annofile <- 'riboqc/my_gencode.v24lift37.annotation.annoout'
lncRNA_peptide_file <- '../ext_data/20190415_lncRNA_RNAseqmax_FINAL.txt'
fasta <-  '../export/ORFs/Fastas/OD5P/OD5P.fasta'

message('...done')

if(!exists("GTF_annotation")) GTF_annotation<-load_obj(annofile)

#get the orfs in the fasta file I sent
sentorfs <- system(str_interp('grep -Poe "ENST[0-9_]+" ${fasta}'),intern=T)
sent_gids <- system(str_interp('grep -Poe "ENSG[0-9_]+" ${fasta}'),intern=T)
sent_pcoding <- sent_gids %in% system(str_interp('grep -e "protein_coding" ${fasta} | grep -Poe "ENSG[0-9_]+" '),intern=T)

sentgtf <- '../export/ORFs/Fastas/OD5P/OD5P.gtf'%>%import
sentseqs <- import(fasta)
all(sent_gids %in% sentgtf$gene_id)
all(sentorfs %in% sentgtf$ORF_id_tr)

#The real question here is, are there ORFs with substantial numbers of riboseq reads that we miss in the analysis?

#read in the data on the lincRNA and ribonovel peptides
lncRNA_table <- fread(lncRNA_peptide_file)
lncRNA_table%<>%set_colnames(c("peptide","transcript_id",  "Gene Name","gene_id"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid

ribonovelhits <-'../ext_data/riboseq_peptides.txt'%>%fread
ribonovelhits%<>%set_colnames(c('id','peptide','gene_id','gene_name','rnaseq_detected'))
ribonovelhits$gene_id%<>%clipid

#many of the genes for this aren't in our feature_counts_all....
#allcountdf[feature_id %in% lncRNA_table$gene_id,]
#but they are in the annotation used for SaTAnnn
peptide_gids <- c(lncRNA_table$gene_id,ribonovelhits$gene_id) 
stopifnot(lncRNA_table$gene_id %in% names(GTF_annotation$genes))
stopifnot(peptide_gids %in% names(GTF_annotation$genes))

#So where are they lost..
peptidegenetranscripts<-GTF_annotation$txs_gene[peptide_gids]
g2tx <- peptidegenetranscripts%>%unlist%>%{data_frame(gid=names(.),trid=.$tx_name)}%>%distinct(gid,trid)
peptidegeneexons <- GTF_annotation$exons_txs[GTF_annotation$txs_gene[peptide_gids]%>%unlist%>%.$tx_name]
peptidegeneexons<-peptidegeneexons%>%unlist%>%{split(.,g2tx$gid[match(names(.),g2tx$trid)])}

# genegr<-lincgenetranscripts
genegr<- peptidegeneexons


#Get the matrix of raw counts for these
allbams <- Sys.glob('star/data/*/*.bam')
allbams %<>% str_subset('OD5P')
#Now get the matrix of raw counts

ovdata<-summarizeOverlaps(features=genegr,reads = allbams,inter.feature=FALSE)
countmat <- assay(ovdata)
colnames(countmat)%<>%str_replace('.bam$','')

#And psite counts
allresfiles <- Sys.glob('riboqc/data/*/*to_SaTAnn*')
if(!exists('allpsiteobs')) allpsiteobs<-allresfiles%>%mclapply(load_objs)
names(allpsiteobs)<-allresfiles%>%dirname%>%basename


getCounts<-function(features,reads){
	counts <- rep(0,length(features))
	ov<-findOverlaps(features,reads,ignore.strand=FALSE)

	summreads <- reads$score[subjectHits(ov)]%>%split(queryHits(ov))%>%map_dbl(sum)

	counts[as.numeric(names(summreads))] <- summreads
	counts

}

ribosamps <- colnames(countmat)%>%str_subset('L[57]',neg=T)

psitecountmat <- lapply(allpsiteobs%>%map('to_SaTAnn')%>%map('P_sites_all')%>%.[ribosamps],getCounts,features=genegr)
psitecountmat%<>%simplify2array
psitecountmat%<>%set_rownames(names(genegr))



gids2plot <- lncRNA_table$gene_id%>%unique

stopifnot(gids2plot%in%rownames(psitecountmat))
stopifnot(gids2plot%in%rownames(countmat))

##Now we need to test how many genes have nested stuff in them....
complex_locus <- genegr[gids2plot]%>%countOverlaps(GTF_annotation$genes)%>%`>`(1)




###Now also read in the scRNAseq data
scdetectfile<-'../data/fracdetected.rds'
if(!file.exists(scdetectfile)) source('./src/parse_ssc_data.R')
#load the data
fracdetected<-readRDS(scdetectfile)
table(lncRNA_table$gene_id %in% fracdetected$ensemblID)
stopifnot(sum(unique(lncRNA_table$gene_id) %in% fracdetected$ensemblID)==35)
#create named vector with cell fractions
scRNAseq_fractionall <-peptide_gids%>%
	match(fracdetected$ensemblID)%>%
	fracdetected$detected[.]%>%
	replace_na(0)%>%
	setNames(.,peptide_gids)
	
scRNAseq_fraction<-scRNAseq_fractionall[gids2plot]

`Gene has ORF` = gids2plot %in% sent_gids

sentseqs <- readAAStringSet(fasta)
lncRNA_table$peptide_hit <- lncRNA_table$peptide%>%lapply(vmatchPattern,subject=sentseqs)%>%map_dbl(.%>%unlist%>%length)
table(lncRNA_table$peptide_hit!=0)

# mycanonseq <- readA5AStringSet('my_gencode.v24lift37.annotation.protein.fa')
# peptide_hits_canon <- lncRNA_table$peptide%>%lapply(vmatchPattern,subject=mycanonseq)%>%map_dbl(.%>%unlist%>%length)
# table(peptide_hits!=0)



# #checking with fuzzy match
# peptide_hits_fuzzy <- lncRNA_table$peptide[`Gene has%>%lapply(vmatchPattern,max.mismatch=1,subject=sentseqs)%>%map_dbl(.%>%unlist%>%length)
# table(peptide_hits_fuzzy==0)



# peptideinfofile<-'../tables/peptide_info_df.tsv'
# peptide_info_df<-read_tsv(peptideinfofile)
# lncRNA_table$peptide = lncRNA_table[[1]]
#  <- lncRNA_table%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%{setNames(.$found,.$gene_id)}
sentseqs<- readAAStringSet(fasta)
lncRNA_table$`peptidehit` <- lncRNA_table$peptide%>%mclapply(vmatchPattern,subject=sentseqs)%>%map_lgl(.%>%unlist%>%length%>%`>`(0))
`peptidehit` <- lncRNA_table%>%group_by(gene_id)%>%summarise(`peptidehit`=any(`peptidehit`))%>%{setNames(.$`peptidehit`,.$gene_id)}%>%.[gids2plot]

stopifnot(lncRNA_table$`Peptide in ORF`%>%table == (c(61,16)))
stopifnot((scRNAseq_fraction[gids2plot])%>%`==`(0)%>%not%>%table == (c(36,35)))

#label outliers
labels2plot <- gids2plot%>%setNames(.,.)
labels2plot[`Gene has ORF`]<- ''
labels2plotscRNAseq <-labels2plot
labels2plotscRNAseq[scRNAseq_fraction < 0.1 ]<-''
labels2plotscRNAseq[psitecountmat[gids2plot,ribosamps]%>%rowSums%>%`<`(130)]<-''
labels2plotspsites <- labels2plot
labels2plotspsites[scRNAseq_fraction < 0.1  ]<-''
labels2plotspsites[psitecountmat[gids2plot,ribosamps]%>%rowSums%>%`<`(130) ]<-''


plotlist<-list
library(ggrepel)
psiterawreadscatterfile<-here('./plots/psite_rawread_hitscatter.pdf')
pdf(psiterawreadscatterfile)
plot1<-qplot(
	x=countmat[gids2plot,ribosamps]%>%rowSums%>%add(1),
	xlab='Raw read Counts OD5P all Ribo Samples',
	y=psitecountmat[gids2plot,ribosamps]%>%rowSums%>%add(1),ylab='Raw Psite Counts OD5P all Ribo samples',
	log='xy',
	shape=`Peptide in ORF`,
	color = `Gene has ORF`,
	label=labels2plotspsites,
	fill = `Gene has ORF`,guide=FALSE)+
	theme_bw()+
	geom_text_repel(alpha=I(0.5),show.legend=FALSE)+
	geom_abline(slope=1,linetype=2)+
	ggtitle('Psites vs. Reads for all\nlincRNA hit peptide attributed Genes\nexons only')
ggMarginal(plot1,groupColour=TRUE,groupFill=TRUE)
dev.off();message(normalizePath(psiterawreadscatterfile))
 

psite_scrnaseq<-here('./plots/psite_scrnaseq_hitscatter.pdf')
pdf(psite_scrnaseq)
plot2<-qplot(
	x=scRNAseq_fraction[gids2plot],
	xlab='Fraction of Cells with scRNAseq Reads',
	y=psitecountmat[gids2plot,ribosamps]%>%rowSums%>%add(1),ylab='Raw Psite Counts OD5P all Ribo samples',
	log='y',
	shape=`Peptide in ORF`,
	color = `Gene has ORF`,
	label=labels2plotscRNAseq,
	fill = `Gene has ORF`)+
	geom_text_repel(alpha=I(0.5),show.legend=FALSE)+
	theme_bw()+
	ggtitle('Psites vs. scRNAseq for all\nlincRNA hit peptide attributed Genes')
# ggMarginal(plot)
print(plot2)
dev.off();message(normalizePath(psite_scrnaseq))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcolor2 <- gg_color_hue(2)[2]



novel_nchlapdf<-ribonovelhits%>%filter(gene_id %in% sent_gids[!sent_pcoding])

##Now also do this for the /ibo-only genes
psite_riboscrnaseq<-here('plots/nc_riboonly_psite_scrnaseq_hitscatter.pdf')
pdf(psite_riboscrnaseq)
plot3<-qplot(
	x=scRNAseq_fractionall[novel_nchlapdf$gene_id],
	xlab='Fraction of Cells with scRNAseq Reads',
	y=psitecountmat[novel_nchlapdf$gene_id,ribosamps]%>%rowSums%>%add(1),ylab='Raw Psite Counts OD5P all Ribo samples',
	log='y',
	# label=labels2plotscRNAseq,
	shape=I('triangle'),
	color=I(ggcolor2)
	)+
	# geom_text_repel(alpha=I(0.5),show.legend=FALSE)+
	theme_bw()+
	ggtitle('Psites vs. scRNAseq for all Peptides found using Riboseq-translatome\n')
# ggMarginal(plot)
print(plot3)
# print(plot(1))
dev.off()

##Now also do this for the /ibo-only genes
psite_riboscrnaseq<-here('plots/riboonly_psite_scrnaseq_hitscatter.pdf')
pdf(psite_riboscrnaseq)
plot3<-qplot(
	x=scRNAseq_fractionall[ribonovelhits$gene_id],
	xlab='Fraction of Cells with scRNAseq Reads',
	y=psitecountmat[ribonovelhits$gene_id,ribosamps]%>%rowSums%>%add(1),ylab='Raw Psite Counts OD5P all Ribo samples',
	log='y',
	# label=labels2plotscRNAseq,
	shape=I('triangle'),
	color=ribonovelhits$gene_id %in% sent_gids[!sent_pcoding]
	)+
	scale_color_manual(name='Non-coding',values=c('black',ggcolor2))+
	# geom_text_repel(alpha=I(0.5),show.legend=FALSE)+
	theme_bw()+
	ggtitle('Psites vs. scRNAseq for all Peptides found using Riboseq-translatome\n')
# ggMarginal(plot)
print(plot3)
# print(plot(1))
dev.off()


stop()









################################################################################
########Addressing gencode discrepency (v22 vs 24 )
################################################################################
	

gencodeold<-import('../ext_data/annotation/gencode.v22.annotation.gtf')
gencodeold$gene_id%<>%str_replace('\\.\\d+','')
stopifnot(mean(gencodeold$gene_id %in% (GTF_annotation$txs_gene%>%names))>.98)
newgenes <- GTF_annotation$txs_gene%>%names%>%setdiff(gencodeold$gene_id)
lostgenes <- gencodeold$gene_id%>%setdiff(GTF_annotation$txs_gene%>%names)

ribonovelhits <-'../ext_data/riboseq_peptides.txt'%>%fread
ribonovelhits%>%head
ribonovelhits$Gene_ID%>%str_replace('\\.\\d+','') %>% is_in( c(newgenes,lostgenes))
ribonovelhits$Gene_ID%>%str_replace('\\.\\d+','') %>% is_in( GTF_annotation$txs_gene%>%names)

lncRNA_table%>%colnames
lncRNA_table%>%filter(!`Gene has ORF`)%>%.$gene_id %>%is_in( lostgenes)



#why isn't this guy getting detected??
testgene<-psitecountmat[,ribosamps]%>%rowSums%>%`>`(1000)%>%which%>%names
stopifnot(testgene=='ENSG00000261645')
testgene<-'ENSG00000261645'
# lncRNA_table$peptide %5n% sents
sentseqs <- import(fasta)
#scatterplot showing the P-site counts, and 


psitecountmat[,ribosamps]%>%rowSums%>%add(1)


# allcounts <- mclapply(allbams,summarizeOverlaps,features = (orfgr%>%head))


#we could also look at the TEs of these genes...


# #load arguments%>%length
# args <- c(
# 	allfeatcounts = 'feature_counts/all_feature_counts',
# 	genelengthfile = 'feature_counts/data/OD5P_05_uM_DAC_1/OD5P_05_uM_DAC_1.feature_counts',
# = 'feature_counts/all_rpkms.txt'
# )
# = 'feature_counts/all_rpkms.txt'

# for(i in names(args)) assign(i,args[i])
# genelengthdf <- genelengthfile%>%fread
# allcountdf <! allfeatcounts%>%fread
# genelengths <- genelengthdf$Length[match(allcountdf$feature_id,genelengthdf$Geneid)]
# # 5
# #allanno <- import('my_gencode.v24lift37.annotation.gtf')
# debug(counts_to_t!m)


fraglengths<-colnames(allcountdf[,-1])%>%{ifelse(str_detect(.,'_L[57]'),101,29)}
lncRNA_table$peptide %5n% sents
sentseqs <- import(fasta)
 counts_to_tpm(as.matrix(allcountdf[,-1]),genelengths,fraglengths,allcountdf[[1]])














# y <- DGEList(counts=allcountdf[,-1],genes=data.frame(gene_id = allcountdf$feature_id,Length=genelengths))
# y <- calcNormFactors(y)
# RPKM <- rpkm(y)

# RPKM%>%length%>%as.data.frame%>%
# 	mutate(gene_id=allcountdf$feature_id)%>%
# 	mutate(gene_id=allcountdf$feature_id)%>%
# 	dplyr::select(-matches('_L5|_L7'))%>%
# 	dplyr::select(gene_id,everything())%>%
# 	write_tsv(outfile)



# # so for the highly transcribed-but-not-detected gene, it looks like it's not detected because the P-sites are attributed to a different gene....
# #' gene detected by satann: 'ENSG00000214391' 

# testpeptide <- lncRNA_table%>%filter(gene_id=='ENSG00000261645')%>%.[[1]]
# #the psites their TRANSCRIPTS overlap, is identical
# GTF_annotation$txs_gene['ENSG00000214391']%>%mergeByOverlaps(to_SaTAnn$P_sites_all)%>%.$score
# GTF_annotation$txs_gene['ENSG00000261645']%>%mergeByOverlaps(to_SaTAnn$P_sites_all)%>%.$score


# stestpeptideprote5ins

# str_detect(prote!ns,testpeptide)

# GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%export('nested.gtf')
# GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%export('nestee.gtf')

# nestedcounts<-GTF_annotation$exons_txs[GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%.$tx_name]%>%{summarizeOverlaps(features=.,reads = allbams)}
# nesteecounts<-GTF_annotation$exons_txs[GTF_annotation$txs_gene['ENSG00000261645']%>%unlist%>%.$tx_name]%>%{summarizeOverlaps(features=.,reads = allbams)}
# assay(nesteecounts)
# 0

# stopifnot(rownames(countmat)==rownames(psitecountmat))

# normalizePath('nested.gtf')
# options(error=NULL)
# #



# peptide_info_df%>%filter(gene_id=='ENSG00000261645')




# #Next weird case....
# highpsitenames <- countmat[,ribosamps]%>%rowSums%>%keep(~ log10(.) > 3.5)%>%sort%>%names

# plotdf<-enframe%>%length(countmat[,ribosamps]%>%rowSums,'gene','reads')%>%cbind(psites=psitecountmat[,ribosamps]%>%rowSums)%>%
# 	mutate(called=`3bp Periodicity`,complex_locus=complex_locus)%>%
# 	# filter(called==FALSE)%>%
# 	arrange(desc(reads))

# testgene<-plotdf%>%filter(called==FALSE)%>%.$gene%>%.[2]

# testgene=='ENSG00000220157')
# #This one is just losing a shitload of it's Psites to RiboseQC...
# testgene=='ENSG00000220157')

# psitecountmat[testgene,]
# readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=testgenegr))

# import!'riboqc/data/OD5P_05_uM_DAC_1/_P_sites_plus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))
# #There are VERY few psites for this thing....
# import('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))$score%>%sum
# #to_satann object agrees with above... did I somehow subset it or something????
# psites <- allres!bs[[ribosamps[1]]]$to_SaTAnn$P_sites_all%>%subsetByOverlaps(testgenegr)
# psites%>%.$score%>%sum
# #NOPE - the new riboqc pipeline agrees with it....
# import('new_riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))
# #maybe if I tate oly the exons?
# testgeneexons <- GTF_annotation$exons_txs[unlist(GTF_annotation$txs_gene[testgene])$tx_name]
# readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=unlist(testgeneexons)))
# readsingene%>%subset(strand=='+')
# readsingene%>%subset(strand=='-')
# #No.... maybe if I limit it to transcripts???
# testgeneexons <- GTF_annotation$cds_txs[unlist(GTF_annotation$txs_gene[testgene])$tx_name]
# readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=unlist(testgeneexons)))


# readsingene%>%subsetByOverlaps(psites)%>%width%>%table
# #Okay at least all the psites overlap stuff in the gene
# psites%>%subsetByOverlaps(readsingene)

# allbams[1]

# normalizePath(allbams[1])
# normalizePath('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw')

# fread('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_calcs')
# readsingene%>%subset(width%in%c(29,28))
# stopifnot(rownames(countmat)==rownames(psitecountmat))
# countmat[
# highpsitenames,ribosamps]%>%rowSums
# testgene<-psitecountmat[,ribosamps]%>%rowSums%>%keep(~ . >1000)%>%which%>%names
# stopifnot(testgene=='ENSG00000261645')
# testgene<-'ENSG00000261645'



# testgene<-plotdf%>%filter(called==FALSE)%>%.$gene%>%.[3]
# stopifnot(testgene=='ENSG00000215548')
# testgenegr <- GTF_annotation$gene[testgene]

# #So this one just looks like it's actually detected. Whhy not in fasta??
# 'ENSG00000215548'


# ###New run?

# newpoolfasta = lapply(Sys.glob('SaTAnn.bak/*OD5P*/SaTAnn_Protein_sequences.fasta'),readAAStringSet)%>%Reduce(f=c)
# sentidsnew <- newpoolfasta%>%names%>%str_extract('[^_]+')%>%unique
# idtbl<-GTF_annotation$txs_gene%>%unlist%>%{data_frame(gene_id=names(.),transcript_id=.$tx_name)}
# sentgenesnew<-data_frame(transcript_id=sentids)%>%left_join(idtbl)
# table(lncRNA_table$gene_id %in% sentgenesnew$gene_id)

