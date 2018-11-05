import sys
import os
import collections
import subprocess
import pandas as pd

RiboDiffscript = ("../../cortexomics/Applications/RiboDiff/scripts/TE.py")
ribodiffconda = '/fast/users/harnettd_c/miniconda3/envs/ribodiff/'
(countfile, 
  dispdiff, 
  outdir) =(list((
  'feature_counts/all_feature_counts',1,'ribodiff'))+sys.argv[1:])[-3:]

sampledata = pd.read_csv('sample_parameter.csv')

contrasts = collections.OrderedDict()

countdata = pd.read_table(countfile,sep=' ')

cell = sampledata['cell_line'].unique()[0]
for cell in sampledata['cell_line'].unique():
  celldata = sampledata[sampledata['cell_line']==cell]

  tgroups = celldata.dosage.unique()

  ctrlgrp = list(filter(lambda x: 'ctrl' in x.lower(),tgroups))
  tgrps = list(filter(lambda x: '05_um_dac' in x.lower(),tgroups))

  if not len(ctrlgrp) == 1 : next

  ctrlgrp = ctrlgrp[0]

  if not len(tgroups) > 0 : next
  
  for tgrp in tgrps:
    contrastdata = (celldata[celldata['dosage'].isin([ctrlgrp,tgrp])]).sort_values(['dosage','assay'],ascending=False)

    contrastdata['Data_Type'] = contrastdata.assay.replace('ribo','Ribo-Seq').replace('total','RNA-Seq')
    contrastdata['Conditions'] = contrastdata.dosage
    contrastdata['Samples'] = contrastdata.sample_id
    
    os.makedirs(os.path.join(outdir,cell+tgrp), exist_ok=True)
    
    sampletablefile =os.path.join(outdir,cell+tgrp,'ribodiff_samples.tsv')
      
    contrastdata.loc[:,['Samples','Data_Type','Conditions']].to_csv(sampletablefile,sep=',',index=False)
    
    countfile =os.path.join(outdir,cell+tgrp,'ribodiff_counts.tsv')
    countcols = [countdata.columns[0]]+list(contrastdata.sample_id)
    countdata.loc[:,countcols].to_csv(countfile,sep='\t',index=False)
  
    ribodiffresfile = os.path.join(outdir,cell+tgrp,'ribodiffres.tsv')

    ribodiffcmd = ("source activate  "+ribodiffconda+" ;" +
      "python "+ RiboDiffscript +
      " -d"+ str(dispdiff)+
      " -e  "+sampletablefile+
      " -c "+countfile+
      " -o "+ribodiffresfile)

    subprocess.call(ribodiffcmd,shell=True)



import glob
ribodiffresfiles = map(os.path.abspath,glob.glob(outdir+'/*/ribodiffres*tsv'))
allribodiffres = pd.concat([pd.read_table(i).assign(contrast=i) for i in ribodiffresfiles])
allribodiffres.contrast = allribodiffres.contrast.str.split('/',expand=True)[8]
allribodiffres.columns = ['feature_id','disper','disperRNA','p_value','adj_p_value','TE1','TE2','log2fc','unnamed','contrast']

allribodiffres.query('adj_p_value < 0.05')

# #read in 
# ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
# ribodiffcols=cols(
#   col_character(),
#   col_double(),
#   col_double(),
#   col_double(),
#   col_double(),
#   col_double(),
#   col_double()
# )
# ribodiffcols$cols%<>%setNames(ribodiffcolsnms)

# ribodiffcontrastobs <- riboseqresfiles%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)

# exprtablessizefactdf<-data_frame(sample_id=sample_list,sizefact=sizeFactors(dds)[sample_list])
# ribobasemeans<-ribocountsfiles%>%map(.%>%{quietly(read_tsv)(.)$result}%>%gather(sample_id,count,-feature_id)%>%
#     left_join(sizefactdf,by='sample_id')%>%
#     group_by(feature_id
#       )%>%summarize(base_mean=sum(count/sizefact))
# )
# ribodiffcontrastobs<- map2(ribodiffcontrastobs,ribobasemeans,left_join)
# for(i in seq_along(ribodiffcontrastobs)){ribodiffcontrastobs[[i]]$log2fc_se<-NA}
# for(i in seq_along(ribodiffcontrastobs)){ribodiffcontrastobs[[i]]$stat<-NA}

# names(ribodiffcontrastobs) = paste0('TE_',names(ribodiffcontrastobs))
# my_contrast_objects %<>% {append(.[setdiff(names(.),names(ribodiffcontrastobs))],ribodiffcontrastobs)}




# ALPHA=0.05
# riboseqresfiles%>%map(.%>%read_tsv%>%filter(padj<ALPHA)%>%.$geneIDs%>%unlist%>%unique)

# inputfolder<-'exprdata'
# outputfolder <- 'exprdata_filt'
# ribodifffolder = 'ribodiff'
# exprtablefiles <- Sys.glob(paste0(inputfolder,'/*'))
# exprtablefiles %<>% str_subset('tsv$|txt$')
# exprtables <- map(exprtablefiles,read_tsv)

# ribodifffiltgenes <- Sys.glob(file.path(ribodifffolder,'riboseqres_*'))%>%
#   map(.%>%read_tsv%>%.[[1]])

# for (exprtablefile in exprtablefiles){
#   exprtable<-read_tsv(exprtablefile)
#   if('gene' in colnames(exprtable))     gcol = 'gene'
#   if('gene_name' in colnames(exprtable))     gcol = 'gene_name'
#   passesfilt <- exprtable[[gcol]]%in%ribodifffiltgenes
#   assert_that(!all(passesfilt==FALSE))
#   message(paste0('Got ',sum(passesfilt),' genes passing filter for ',exprtablefile))
#   exprtablefilt <- exprtable[passesfilt,]
#   write_tsv(file.path(outputf
#     older,basename(exprtablefile)))
# }
