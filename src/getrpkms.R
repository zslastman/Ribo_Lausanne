#Get rpkms for our genes
message('loading libraries')
library(magrittr)
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

message('...done')

GTF_annotation<-load_obj('riboqc/my_gencode.v24lift37.annotation.annoout')

#get the orfs in the fasta file I sent
fasta = '../export/ORFs/Fastas/OD5P/OD5P.fasta'
sentorfs <- system(str_interp('grep -Poe "ENST[0-9_]+" ${fasta}'),intern=T)
sent_gids <- system(str_interp('grep -Poe "ENSG[0-9_]+" ${fasta}'),intern=T)
sentgtf <- '../export/ORFs/Fastas/OD5P/OD5P.gtf'%>%import
all(sent_gids %in% sentgtf$gene_id)
all(sentorfs %in% sentgtf$ORF_id_tr)

#The real question here is, are there ORFs with substantial numbers of riboseq reads that we miss in the analysis?

#read in the data on the lincRNA peptides
lncRNA_table <- readxl::read_excel('../ext_data/20180719_OD5P_lncRNA_RE_peptides_combiner.xlsx',sheet=2)%>%
	set_colnames(c("MS-identified peptide", "gene_id", "Gene Name", 
"Gene Type", "FPKM RNAseq OD5PCTRL", "FPKM RNAseq DAC"))
clipid <- . %>% str_replace('\\.\\d+$','')
lncRNA_table$gene_id%<>%clipid
#many of the genes for this aren't in our feature_counts_all....
#allcountdf[feature_id %in% lncRNA_table$gene_id,]
#but they are in the annotation used for SaTAnnn
stopifnot(lncRNA_table$gene_id %in% names(GTF_annotation$genes))

#So where are they lost..
lincgenetranscripts<-GTF_annotation$txs_gene[lncRNA_table$gene_id]
lincgeneexons <- GTF_annotation$exons_txs[GTF_annotation$txs_gene[lncRNA_table$gene_id]%>%unlist%>%.$tx_name]


#Get the matrix of raw counts for these
allbams <- Sys.glob('star/data/*/*.bam')
allbams %<>% str_subset('OD5P')
ovdata<-summarizeOverlaps(features=lincgenetranscripts,reads = allbams,inter.feature=FALSE)
countmat <- assay(ovdata)
colnames(countmat)%<>%str_replace('.bam$','')

#And psite counts
allresfiles <- Sys.glob('riboqc/data/*/*to_SaTAnn*')
allresobs<-allresfiles%>%mclapply(load_objs)
names(allresobs)<-allresfiles%>%dirname%>%basename


getCounts<-function(features,reads){
	counts <- rep(0,length(features))
	ov<-findOverlaps(features,reads,ignore.strand=FALSE)

	summreads <- reads$score[subjectHits(ov)]%>%split(queryHits(ov))%>%map_dbl(sum)

	counts[as.numeric(names(summreads))] <- summreads
	counts

}

ribosamps <- colnames(countmat)%>%str_subset('L[57]',neg=T)

psitecountmat <- lapply(allresobs%>%map('to_SaTAnn')%>%map('P_sites_all')%>%.[ribosamps],getCounts,features=lincgenetranscripts)
psitecountmat%<>%simplify2array
psitecountmat%<>%set_rownames(names(lincgenetranscripts))

##Now we need to test how many genes have nested stuff in them....
complex_locus <- lincgenetranscripts%>%countOverlaps(GTF_annotation$genes)%>%`>`(1)



psiterawreadscatterfile<-'../plots/psite_rawread_hitscatter.pdf'
#dirname(psiterawreadscatterfile)%>%dir.create
`Gene has ORF` = rownames(countmat) %in% sent_gids

peptideinfofile<-'../tables/peptide_info_df.tsv'
peptide_info_df<-read_tsv(peptideinfofile)
lncRNA_table$peptide = lncRNA_table[[1]]
`Peptide in ORF` <- lncRNA_table%>%mutate(found = peptide %in% peptide_info_df$peptide)%>%{setNames(.$found,.$gene_id)}

`Peptide in ORF`%>%table

pdf(psiterawreadscatterfile)
plot<-qplot(
	x=countmat[,ribosamps]%>%rowSums%>%add(1),
	xlab='Raw read Counts OD5P all Ribo Samples',
	y=psitecountmat[,ribosamps]%>%rowSums%>%add(1),ylab='Raw Psite Counts OD5P all Ribo samples',
	log='xy',
	shape=`Peptide in ORF`,
	color = `Gene has ORF`,
	fill = `Gene has ORF`)+
	theme_bw()+
	geom_abline(slope=1,linetype=2)+
	ggtitle('Psites vs. Reads for all\nlincRNA hit peptide attributed Genes')
ggMarginal(plot,groupColour=TRUE,groupFill=TRUE)
dev.off();message(normalizePath(psiterawreadscatterfile))



#why isn't this guy getting detected??
testgene<-psitecountmat[,ribosamps]%>%rowSums%>%`>`(1000)%>%which%>%names
stopifnot(testgene=='ENSG00000261645')
testgene<-'ENSG00000261645'



#scatterplot showing the P-site counts, and 


psitecountmat[,ribosamps]%>%rowSums%>%add(1)


# allcounts <- mclapply(allbams,summarizeOverlaps,features = (orfgr%>%head))


#we could also look at the TEs of these genes...


#load arguments
args <- c(
	allfeatcounts = 'feature_counts/all_feature_counts',
	genelengthfile = 'feature_counts/data/OD5P_05_uM_DAC_1/OD5P_05_uM_DAC_1.feature_counts',
	outfile = 'feature_counts/all_rpkms.txt'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

genelengthdf <- genelengthfile%>%fread
allcountdf <- allfeatcounts%>%fread

genelengths <- genelengthdf$Length[match(allcountdf$feature_id,genelengthdf$Geneid)]

#allanno <- import('my_gencode.v24lift37.annotation.gtf')
debug(counts_to_tpm)


fraglengths<-colnames(allcountdf[,-1])%>%{ifelse(str_detect(.,'_L[57]'),101,29)}

TPMdf <- counts_to_tpm(as.matrix(allcountdf[,-1]),genelengths,fraglengths,allcountdf[[1]])



















y <- DGEList(counts=allcountdf[,-1],genes=data.frame(gene_id = allcountdf$feature_id,Length=genelengths))
y <- calcNormFactors(y)
RPKM <- rpkm(y)

RPKM%>%as.data.frame%>%
	mutate(gene_id=allcountdf$feature_id)%>%
	mutate(gene_id=allcountdf$feature_id)%>%
	dplyr::select(-matches('_L5|_L7'))%>%
	dplyr::select(gene_id,everything())%>%
	write_tsv(outfile)


'LHRP3'


#' so for the highly transcribed-but-not-detected gene, it looks like it's not detected because the P-sites are attributed to a different gene....
#' gene detected by satann: 'ENSG00000214391' 
#' gene attributed to peptide 'ENSG00000261645'
testpeptide <- lncRNA_table%>%filter(gene_id=='ENSG00000261645')%>%.[[1]]

#the psites their TRANSCRIPTS overlap, is identical
GTF_annotation$txs_gene['ENSG00000214391']%>%mergeByOverlaps(to_SaTAnn$P_sites_all)%>%.$score
GTF_annotation$txs_gene['ENSG00000261645']%>%mergeByOverlaps(to_SaTAnn$P_sites_all)%>%.$score

stestpeptideproteins

str_detect(proteins,testpeptide)

GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%export('nested.gtf')
GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%export('nestee.gtf')

nestedcounts<-GTF_annotation$exons_txs[GTF_annotation$txs_gene['ENSG00000214391']%>%unlist%>%.$tx_name]%>%{summarizeOverlaps(features=.,reads = allbams)}
nesteecounts<-GTF_annotation$exons_txs[GTF_annotation$txs_gene['ENSG00000261645']%>%unlist%>%.$tx_name]%>%{summarizeOverlaps(features=.,reads = allbams)}
assay(nesteecounts)
0

stopifnot(rownames(countmat)==rownames(psitecountmat))

normalizePath('nested.gtf')
options(error=NULL)
#



peptide_info_df%>%filter(gene_id=='ENSG00000261645')




#Next weird case....
highpsitenames <- countmat[,ribosamps]%>%rowSums%>%keep(~ log10(.) > 3.5)%>%sort%>%names

plotdf<-enframe(countmat[,ribosamps]%>%rowSums,'gene','reads')%>%cbind(psites=psitecountmat[,ribosamps]%>%rowSums)%>%
	mutate(called=`3bp Periodicity`,complex_locus=complex_locus)%>%
	# filter(called==FALSE)%>%
	arrange(desc(reads))

testgene<-plotdf%>%filter(called==FALSE)%>%.$gene%>%.[2]

stopifnot(testgene=='ENSG00000220157')
#This one is just losing a shitload of it's Psites to RiboseQC...
countmat[testgene,]
psitecountmat[testgene,]
testgenegr <- GTF_annotation$gene[testgene]
#plenty of reads on it...
readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=testgenegr))

import('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_plus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))
#There are VERY few psites for this thing....
import('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))$score%>%sum
#to_satann object agrees with above... did I somehow subset it or something????
psites <- allresobs[[ribosamps[1]]]$to_SaTAnn$P_sites_all%>%subsetByOverlaps(testgenegr)
psites%>%.$score%>%sum
#NOPE - the new riboqc pipeline agrees with it....
import('new_riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw',which=testgenegr%>%keepSeqlevels(testgenegr@seqnames))
#maybe if I tate oly the exons?
testgeneexons <- GTF_annotation$exons_txs[unlist(GTF_annotation$txs_gene[testgene])$tx_name]
readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=unlist(testgeneexons)))
readsingene%>%subset(strand=='+')
readsingene%>%subset(strand=='-')
#No.... maybe if I limit it to transcripts???
testgeneexons <- GTF_annotation$cds_txs[unlist(GTF_annotation$txs_gene[testgene])$tx_name]
readsingene<-readGAlignments(allbams[1],param=ScanBamParam(which=unlist(testgeneexons)))


readsingene%>%subsetByOverlaps(psites)%>%width%>%table
#Okay at least all the psites overlap stuff in the gene
psites%>%subsetByOverlaps(readsingene)

allbams[1]

normalizePath(allbams[1])
normalizePath('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_minus.bw')

fread('riboqc/data/OD5P_05_uM_DAC_1/_P_sites_calcs')
readsingene%>%subset(width%in%c(29,28))
stopifnot(rownames(countmat)==rownames(psitecountmat))
countmat[
highpsitenames,ribosamps]%>%rowSums
testgene<-psitecountmat[,ribosamps]%>%rowSums%>%keep(~ . >1000)%>%which%>%names
stopifnot(testgene=='ENSG00000261645')
testgene<-'ENSG00000261645'



testgene<-plotdf%>%filter(called==FALSE)%>%.$gene%>%.[3]
stopifnot(testgene=='ENSG00000215548')
testgenegr <- GTF_annotation$gene[testgene]

#So this one just looks like it's actually detected. Whhy not in fasta??
'ENSG00000215548'



