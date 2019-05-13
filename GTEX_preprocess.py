#! python

# Graphing library imports


# logging.info('logging info')
# logging.warning('logging.warning')

import os, sys, glob, logging
import inspect, urllib

import pandas, scipy, numpy
from pprint import pprint

import time, datetime

import seaborn as sns

import plotly
import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.graph_objs as go
from plotly.offline import iplot, init_notebook_mode

# Using plotly + cufflinks in offline mode
import cufflinks
import inspect
__file__ = os.path.abspath(inspect.getfile(lambda: None))
import pathlib
script_dir = pathlib.Path().resolve()


SIMPLE_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(level=logging.INFO, format=SIMPLE_FORMAT, datefmt=DATEFMT)
logging.info("logging")

# Gene transcripts per million data
cache_dir = os.path.join(pathlib.Path().resolve(), 'data_cache')

# Detailed data file:
# This file is quite large ~ 1GB
gtex_tpm_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
gtex_tpm_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz')
# Using unzipped version
gtex_tpm_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct')

gtex_tpm_med_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')

# Download to a cache
SIMPLE_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(level=logging.INFO, format=SIMPLE_FORMAT, datefmt=DATEFMT)
logging.info("logging")

# GTEx data from https://gtexportal.org/home/datasets
# Gene transcripts per million data
cache_dir = os.path.join(pathlib.Path().resolve(), 'data_cache')
if not os.path.exists(cache_dir):
    logging.info('Createing {}'.format(cache_dir))
    os.makedirs(cache_dir, mode=0o777)

# Detailed data file:
# This file is quite large ~ 1GB
gtex_tpm_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
gtex_tpm_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz')
# Using unzipped version
# gtex_tpm_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct')
if not os.path.exists(gtex_tpm_fn):
    print("downloading {}".format(gtex_tpm_link))
    urllib.request.urlretrieve(gtex_tpm_link, gtex_tpm_fn)

# Just looking at median values by tissue (SMTSD) greatly reduces file size
gtex_tpm_med_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
gtex_tpm_med_fn = os.path.join(cache_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz')
if not os.path.exists(gtex_tpm_med_fn):
    print("downloading {}".format(gtex_tpm_med_link))
    urllib.request.urlretrieve(gtex_tpm_med_link, gtex_tpm_med_fn)

# Sample annotation data

# Subject Phenotype data
# Age, sex and Hardy Scale death circumstances
gtex_pheno_link = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'
gtex_pheno_fn = os.path.join(cache_dir, 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt')

# Main sample data of interest is:
# 'SMTS': 'Tissue Type, area from which the tissue sample was taken.  This is a parent value to SMTSD.'
# 'SMTSD': 'SMTS Detailed'
gtex_attr_link = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt'
gtex_attr_fn = os.path.join(cache_dir, 'GTEx_v7_Annotations_SampleAttributesDS.txt')
gtex_attr_desc_link = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx'
gtex_attr_desc_fn = os.path.join(cache_dir, 'GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx')
if not os.path.exists(gtex_attr_fn):
    logging.info('downloading {}'.format(gtex_attr_link))
    urllib.request.urlretrieve(gtex_attr_link, gtex_attr_fn)
if not os.path.exists(gtex_attr_desc_fn):
    logging.info('downloading {}'.format(gtex_attr_desc_link))
    urllib.request.urlretrieve(gtex_attr_desc_link, gtex_attr_desc_fn)

# Load tissue specific detatils
gtex_attr = pandas.read_csv(gtex_attr_fn, sep='\t')
tissue_types = sorted(gtex_attr['SMTS'].unique())
tissue_dict = dict()
for tissue in tissue_types:
    tissue_dict[tissue] = dict()
    tissue_dict[tissue]['samples'] = gtex_attr[gtex_attr['SMTS'] == tissue]['SAMPID'].to_list()
    tissue_dict[tissue]['subtypes'] = dict()
    for subtype in gtex_attr[gtex_attr['SMTS'] == tissue]['SMTSD'].unique():
        # Get all samples by subtype
        tissue_dict[tissue]['subtypes'][subtype] = gtex_attr[gtex_attr['SMTSD'] == subtype]['SAMPID'].to_list()

    sorted(gtex_attr[gtex_attr['SMTS'] == tissue]['SMTSD'].unique())
tissue_subtypes = sorted(gtex_attr['SMTSD'].unique())

# log number of tissue samples
logging.info('%s separate tissue types - total samples %s', len(tissue_dict), len(gtex_attr['SAMPID']))
for tissue in sorted(tissue_dict.keys()):
    logging.info('%s (n=%s)', tissue, len(tissue_dict[tissue]['samples']))
    if len(tissue_dict[tissue]['subtypes']) > 1:
        for subtype in sorted(tissue_dict[tissue]['subtypes'].keys()):
            logging.info('\t%s (n=%s)', subtype, len(tissue_dict[tissue]['subtypes'][subtype]))

logging.info('Loading the full GTEX table...')
gtex = pandas.read_csv(gtex_tpm_fn, sep='\t', skiprows=2)
logging.info('\t finished loading gtex')
logging.info(gtex.info())
# Set up multindex values
logging.info('resetting index')
logging.info(gtex.columns)
logging.info(gtex.index)

sample_cols = [c for c in gtex.columns if c not in ('Description', 'Name')]
gtex['ensid'] = gtex['Name'].apply(lambda x: x.split('.')[0])
gtex['ens_ver'] = gtex['Name'].apply(lambda x: x.split('.')[-1])
gtex.rename(columns={'Description': 'Symbol'}, inplace=True)
gtex.drop('Name', inplace=True, axis=1)
gtex.set_index(['Symbol', 'ensid', 'ens_ver'], inplace=True)
logging.info('creating column lists')
gtex_attr.set_index('SAMPID', inplace=True)

def subtype_rename(subtype):
    if ' - ' in subtype:
        if subtype.startswith('Kidney - '):
            return 'Kidney'
        elif subtype.startswith('Nerve - '):
            return 'Nerve'
        else:
            return subtype.split(' - ')[-1]
    return subtype
tissue_list = [gtex_attr.loc[col]['SMTS'] for col in sample_cols]
gtex_attr['subtype'] = gtex_attr['SMTSD'].apply(subtype_rename)
subtype_list = [gtex_attr.loc[col]['subtype'] for col in sample_cols]

logging.info('adding multiindex columns')
gtex.columns = [tissue_list, subtype_list, sample_cols]
gtex.columns.names = ['tissue', 'subtype', 'sample']
logging.info(gtex.info())
logging.info(gtex.describe())
logging.info('applying safelog2 transform')
# applying log2 transformation
safelog2 = lambda x: numpy.log2(x) if x else numpy.nan
gtex = gtex.applymap(safelog2)
logging.info(gtex.info())
logging.info(gtex.describe())

if not os.path.exists(os.path.join(cache_dir, 'preprocess')):
    os.makedirs(os.path.join(cache_dir, 'preprocess'))

logging.info('Calculating %s', 'gtex_log2tpm.describe.tsv')
out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_all.describe.tsv')
desc = gtex.T.describe().T
desc['skew'] = gtex.apply(lambda row: row.skew(), axis=1)
desc['kurtosis'] = gtex.apply(lambda row: row.kurtosis(), axis=1)
desc.to_csv(out_fn, sep='\t', header=True, index=True)
logging.info('wrote %s', os.path.abspath(out_fn))

for tissue in sorted(set(tissue_list)):
    logging.info("Calculating '%s' statistics", tissue)
    gt = gtex.loc[:, tissue]
    subtypes = sorted(set([i[0] for i in gt.columns.to_flat_index()]))
    for subtype in subtypes:
        out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_{}__{}.describe.tsv'.format(tissue, subtype))
        logging.info('Calculating %s: %s', tissue, subtype)
        gts = gtex.loc[:, gtex.columns.get_level_values('subtype') == subtype]        
        desc = gts.T.describe().T
        desc['skew'] = gts.apply(lambda row: row.skew(), axis=1)
        desc['kurtosis'] = gts.apply(lambda row: row.kurtosis(), axis=1)
        desc.to_csv(out_fn, sep='\t', header=True, index=True)
        logging.info('wrote %s', os.path.abspath(out_fn))
    if len(subtypes) > 1:
        out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_{}__{}.describe.tsv'.format(tissue, 'all'))
        desc = gt.T.describe().T
        desc['skew'] = gt.apply(lambda row: row.skew(), axis=1)
        desc['kurtosis'] = gt.apply(lambda row: row.kurtosis(), axis=1)
        desc.to_csv(out_fn, sep='\t', header=True, index=True)
        logging.info('wrote %s', os.path.abspath(out_fn))


