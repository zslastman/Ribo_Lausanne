suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,DESeq2))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,xtail))

args <- c(
  countfile='feature_counts/all_feature_counts',
  outdir= 'xtail'
)

args[] = commandArgs(trailingOnly=TRUE)
for(i in names(args)) assign(i,args[i])

featurecountsagg <- data.table::fread(countfile)

sampledf = fread('sample_parameter.csv')%>%
  filter(!str_detect(sample_id,'test'))%>%
  group_by(cell_line)%>%
  filter(!sample_id%>%str_detect('125_'))%>%
  filter(!sample_id%>%str_detect('ctrl\\d(L|B)'))%>%
  filter(any(assay=='ribo'),any(assay=='total'))

cell_lines<-sampledf$cell_line%>%unique

featurecountsagg = featurecountsagg%>%select(feature_id,one_of(sampledf$sample_id))

xtailfiles = file.path(paste0(outdir,'/xtail_',cell_lines,'.txt'))%>%setNames(cell_lines)

for(cell_linei in cell_lines){

    mrnasamps <- sampledf%>%filter(assay=='total',cell_line==cell_linei)%>%.$sample_id
    ribosamps <- sampledf%>%filter(assay=='ribo',cell_line==cell_linei)%>%.$sample_id

    #samples to use
    xtailcounts =   
      featurecountsagg%>%
      dplyr::select(matches(paste0('feature_id|')),one_of(mrnasamps),one_of(ribosamps))
    
    #pick the largest libraries, so we have matching pairs
    libsizes <- colSums(xtailcounts[,-1])
    ribosamps <- c(
        libsizes[ribosamps]%>%.[str_detect(names(.),regex(ignore=T,'ctrl'))]%>%sort%>%tail(2)%>%names,
        libsizes[ribosamps]%>%.[str_detect(names(.),regex(ignore=T,'DAC'))]%>%sort%>%tail(2)%>%names
    )
    mrnasamps <- c(
        libsizes[mrnasamps]%>%.[str_detect(names(.),regex(ignore=T,'ctrl'))]%>%sort%>%tail(2)%>%names,
        libsizes[mrnasamps]%>%.[str_detect(names(.),regex(ignore=T,'DAC'))]%>%sort%>%tail(2)%>%names
    )
    stopifnot(length(c(ribosamps,mrnasamps))==8)
    conditionvect <- c('control','control','treat','treat')

    xtailcounts = xtailcounts[!apply(xtailcounts,1,.%>%is.na%>%any),]
    mrnatab <- xtailcounts%>%select(one_of(mrnasamps))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)
    ribotab <- xtailcounts%>%select(one_of(ribosamps))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)

    xtailres <- xtail::xtail(
      mrnatab,
      ribotab,
      conditionvect,
      threads=8
    )

    xtailtable <- xtailres$resultsTable%>%rownames_to_column%>%set_colnames(
      c("feature_id","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", "E145_log2TE", 
      "E13_log2TE", "log2FC_TE_v2", "pvalue_v2", "log2fc", 
      "p_value", "adj_p_value")
    )

    xtailfiles[cell_linei]%>%dirname%>%dir.create(show=T)
    write_tsv(xtailtable,xtailfiles[cell_linei])
  
}