# the qc step at the moment assumes that the data is single ended
# #
# shell.executable("bash")
# shell.prefix("set -e  pipefail;")
# user set parameter
TMPDIR = '../tmp'
SCRIPTDIR = '../git/rna_seq/scripts'

# #reference genome
REF_orig = '../genomes/hg38.fa'

import os


# # used by star
# STARINDEX = "/fast/projects/cubit/0.12.0/static_data/precomputed/STAR/2.4.1d/GENCODE/M12/GRCm38/p5.chr_scaff/50/"
# # used by 
# GTF_orig = '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gtf'
# GTF_cdsfilt = 'static_local/gencode.vM12.annotation_cdsfilt.gtf'
# CDSGTF = 'static_local/gencode.vM12.annotation.cds.gtf'

# # used by infer_experiment
# BED = 'static_local/gencode.vM12.annotation.bed'
# GFF = '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3'

REF = 'my_'+os.path.splitext(os.path.split(REF_orig)[1])[0]+'.fa'


#STAR uses the rsem index
STARINDEX = 'rsemref'

# used by infer_experiment
GFF_orig = 'annotation/gencode.v21.annotation.gtf.gz'

ANNOBASE = 'my_'+os.path.splitext(os.path.split(GFF_orig)[1])[0]

GFF = ANNOBASE+'.gff3'
GFF_mod = ANNOBASE+'.mod.gff3'

BED = ANNOBASE+'.bed'

# used by 
GTF = ANNOBASE+'.gtf'
CDSGTF = ANNOBASE+'.cdsfilt.gtf'
#used to make indices
RNAFASTA = ANNOBASE+'.transcript.fa'
CDSFASTA = ANNOBASE+'.cds.fa'

# need to think on this.... we start with a genome and a gff3 file. We add the extra transcript to our gff3 file, (optionally eliminating all others) 
#then we create the transcript fasta file from our gff3 file, and the gtf file. Then, we use these to make salmon, and rsem/star indices

# htseq parameter
HTSEQ_MODE = 'union'


# used by qc
SeQC_GTF = 'gencode.vM12.chr_scaff.annotation.no_genes.protein_coding.rRNA.gtf'
SeQC_REF = REF


RSEMINDEXFOLD="rsemref"

SAMPLELINES = [line.strip().split(',') for line in open("sample_parameter.csv").readlines()]

# #switch for testmode
# if(config.get('test',0)): 
#   print('\n\n--------testmode on -------------\n\n')
#   SAMPLELINES = SAMPLELINES[0:2]  
#   origSAMPLE = [ entry[SAMPLELINES[0].index('sample_id')] for entry in SAMPLELINES[1:]]
#   SAMPLES = ['test']


SAMPLES = [ entry[SAMPLELINES[0].index('sample_id')] for entry in SAMPLELINES[1:]]



#get our info as dictionaries
def tab_to_dict(SAMPLELINES,valcol):
  valind = SAMPLELINES[0].index(valcol)
  vals   = [ entry[valind] for entry in SAMPLELINES[1:]]
  return dict(zip(SAMPLES,vals))

#sample - info dictionaries
LIBRARY_DICT          = tab_to_dict(SAMPLELINES,'library_layout')
READ_PATTERN_DICT     = tab_to_dict(SAMPLELINES,'read_pattern')
PROTOCOL_DICT         = tab_to_dict(SAMPLELINES, 'protocol')
FRAG_LENGTH_MEAN_DICT = tab_to_dict(SAMPLELINES, 'fragment_length_mean')
FRAG_LENGTH_SD_DICT   = tab_to_dict(SAMPLELINES, 'fragment_length_sd')

print(SAMPLELINES)

ASSAY_DICT            = tab_to_dict(SAMPLELINES, 'assay')


#for f in $(echo input/*); do for q in $( echo ${f}/* ); do echo $f $q; done; done | sed 's/input\///' > pipeline/sample_file.txt
SAMPLEFASTQLINES = [line.strip().split(' ') for line in open("sample_file.txt").readlines()]
FASTQS = [l[1] for l in SAMPLEFASTQLINES]
FASTQSAMPLES = [l[0] for l in SAMPLEFASTQLINES]
FASTQSAMPLEDICT = dict(zip(FASTQS,FASTQSAMPLES))
SAMPLEFASTQDICT = {v:[i for i in FASTQSAMPLEDICT.keys() if FASTQSAMPLEDICT[i] == v ] for k,v in FASTQSAMPLEDICT.items()}
assert FASTQSAMPLES in SAMPLES
