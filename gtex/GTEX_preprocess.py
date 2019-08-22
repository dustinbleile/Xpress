#! python
import shutil
import os, sys, glob, logging
import gzip
import inspect, urllib


import pandas, scipy, numpy
from pprint import pprint

import time, datetime

import seaborn as sns

import inspect
__file__ = os.path.abspath(inspect.getfile(lambda: None))
import pathlib
script_dir = pathlib.Path().resolve()

SIMPLE_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(level=logging.INFO, format=SIMPLE_FORMAT, datefmt=DATEFMT)
logging.info("logging")


class ExpressionData:
    def __init__(self, **kwargs):
        self.compendium_group = kwargs.pop('compendium_group', '')
        self.compendium_group_ver = kwargs.pop('compendium_group_ver', '')


class GtexData(ExpressionData):
    """
    GTEx data helper
    """
    GTEX_VERSION = 'v7'

    GTEX_TPM_LINK = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
    GTEX_TPM_BASENAME = 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'

    GTEX_TPM_MEDIAN_LINK = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
    GTEX_TPM_MEDIAN_BASENAME = 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct'

    GTEX_PHENO_LINK = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'
    GTEX_ATTR_LINK = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt'
    GTEX_ATTR_DESC_LINK = 'https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx'

    def __init__(self, cache_dir=None, gtex_tpm_fn=''):
        # Gene transcripts per million data
        self.cache_dir = os.path.join(os.path.dirname(__file__), 'data_cache') if not cache_dir else cache_dir
        self.gtex_tpm_fn = gtex_tpm_fn if gtex_tpm_fn else os.path.join(self.cache_dir, self.GTEX_TPM_BASENAME)
        self.gtex_tpm_med_fn = os.path.join(self.cache_dir, self.GTEX_TPM_MEDIAN_BASENAME)
        # Main sample data of interest is:
        # 'SMTS': 'Tissue Type, area from which the tissue sample was taken.  This is a parent value to SMTSD.'
        # 'SMTSD': 'SMTS Detailed'
        self.gtex_attr_fn = os.path.join(self.cache_dir, 'GTEx_v7_Annotations_SampleAttributesDS.txt')
        self.gtex_attr_desc_fn = os.path.join(self.cache_dir, 'GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx')
        # Subject Phenotype data
        # Age, sex and Hardy Scale death circumstances
        self.gtex_pheno_fn = os.path.join(self.cache_dir, 'GTEx_v7_Annotations_SubjectPhenotypesDS.txt')

        self.tissue_dict = None
        self.GTEX7_ROW_TABLE = None
        self.GTEX7_COL_TABLE = None
        self.cache_data()

    def tissues(self):
        """
        All possible samples, that may not exist in the current table
        """
        if self.tissue_dict is not None:
            return self.tissue_dict

        # Load tissue specific detatils
        gtex_attr = pandas.read_csv(self.gtex_attr_fn, sep='\t')
        tissue_types = sorted(gtex_attr['SMTS'].unique())
        tissue_dict = dict()
        for tissue in tissue_types:
            tissue_dict[tissue] = dict()
            tissue_dict[tissue]['samples'] = gtex_attr[gtex_attr['SMTS'] == tissue]['SAMPID'].tolist()
            tissue_dict[tissue]['subtypes'] = dict()
            for subtype in gtex_attr[gtex_attr['SMTS'] == tissue]['SMTSD'].unique():
                # Get all samples by subtype
                tissue_dict[tissue]['subtypes'][subtype] = gtex_attr[gtex_attr['SMTSD'] == subtype]['SAMPID'].tolist()

            sorted(gtex_attr[gtex_attr['SMTS'] == tissue]['SMTSD'].unique())

        # log number of tissue samples
        logging.info('%s separate tissue types - total samples %s', len(tissue_dict), len(gtex_attr['SAMPID']))
        for tissue in sorted(tissue_dict.keys()):
            logging.info('%s (n=%s)', tissue, len(tissue_dict[tissue]['samples']))
            if len(tissue_dict[tissue]['subtypes']) > 1:
                for subtype in sorted(tissue_dict[tissue]['subtypes'].keys()):
                    logging.info('\t%s (n=%s)', subtype, len(tissue_dict[tissue]['subtypes'][subtype]))
        self.tissue_dict = tissue_dict
        return self.tissue_dict

    def gtex_gene_table(self, gene_names, comparators=None):
        """
        gene_names list of gene_name and alts - use first value found.
        """
        sample_list = []
        if comparators:
            comparators = [comparators] if isinstance(comparators, str) else comparators
            tables = []
            for comp in comparators:
                sample_list = []
                if '__' not in comp:
                    sample_list.extend(self.tissues()[comp]['samples'])
                else:
                    tissue, subtype = comp.split('__')
                    sample_list.extend(self.tissues()[tissue]['subtypes'][subtype])
                part_table = self.gtex_tpm_partial_table_load(gene_names, sample_list)
                part_table.rename(columns={0: comp + '_GTEXv7(TPM)'}, inplace=True)
                tables.append(part_table)
            table = pandas.concat(tables, sort=False)
        else:
            table = self.gtex_tpm_partial_table_load(gene_names, sample_list)
        return table

    def gtpm_row_by_genename(self, gene_name):
        if self.GTEX7_ROW_TABLE is None:
            logging.info('Loading reference GTEX7 row table')
            self.GTEX7_ROW_TABLE = pandas.read_csv(self.gtex_tpm_fn, skiprows=2, sep='\t', usecols=range(3))
            self.GTEX7_ROW_TABLE['ensid'] = self.GTEX7_ROW_TABLE['Name'].apply(lambda x: x.split('.')[0])  # stable version ENSG
        sub_table = self.GTEX7_ROW_TABLE[self.GTEX7_ROW_TABLE['Description'] == gene_name]
        if sub_table.empty:
            sub_table = self.GTEX7_ROW_TABLE[self.GTEX7_ROW_TABLE['ensid'] == gene_name]
        if sub_table.empty:
            logging.error("Missing '%s row in GTEX TPM table'", gene_name)
            return None
        return [2 + 1 + i for i in sub_table.index][0]

    def gtpm_keep_rows(self, gene_list=None):
        row_list = [self.gtpm_row_by_genename(gene) for gene in gene_list]
        return [2] + [gene for gene in row_list if gene]

    def gtex_tpm_partial_table_load(self, gene_list, sample_list=None):
        start_time = time.time()
        if isinstance(gene_list, str):
            gene_list = [gene_list]
        if gene_list:
            keep_rows = self.gtpm_keep_rows(gene_list)
            skiprows_func = lambda x: x not in keep_rows
            logging.info('loading %s GTEX TPM rows', len(keep_rows))
        else:
            logging.info('loading all rows')
            skiprows_func = lambda x: x in [0, 1]

        if sample_list:
            if self.GTEX7_COL_TABLE is None:
                logging.info('Loading reference GTEX7 column table')
                self.GTEX7_COL_TABLE = pandas.read_csv(self.gtex_tpm_fn, skiprows=lambda x: x not in [2, 3], sep='\t')
            logging.info('Limiting load to %s potential samples', len(sample_list))
            # Not all samples are in tables
            sample_list = [s for s in sample_list if s in self.GTEX7_COL_TABLE.columns]
            logging.info('Limiting load to %s existing samples', len(sample_list))
            columns = ['Name', 'Description'] + sample_list
        else:
            sample_list = None

        gtpm_partial = pandas.read_csv(
            self.gtex_tpm_fn,
            skiprows=skiprows_func,  # Only selected gene rows
            sep='\t',
            header=0,
            usecols=sample_list,  # Number of samples limit
        )

        end_time = time.time()
        logging.info("Gtex TPM partial table took {} minutes".format((end_time - start_time)/60))
        return gtpm_partial

    def cache_data(self):
        """
        Get all big data files into cache
        """
        # GTEx data from https://gtexportal.org/home/datasets
        # Gene transcripts per million data
        if not os.path.exists(self.cache_dir):
            logging.info('Creating {}'.format(self.cache_dir))
            os.makedirs(self.cache_dir, mode=0o777)

        # Detailed data file:
        # This file is quite large even zipped ~ 1GB
        if not os.path.exists(self.gtex_tpm_fn):
            gzip_temp = self.gtex_tpm_fn + '.gz'
            if not os.path.exists(gzip_temp):
                logging.info('Getting large GTEx TPM file - take some time')
                logging.info('\tDownloading from: %s', self.GTEX_TPM_LINK)
                logging.info('\tDownloading to : %s', gzip_temp)
                urllib.request.urlretrieve(self.GTEX_TPM_LINK, gzip_temp)
                if gzip_temp:
                    logging.info('\tFinished Downloading: %s', self.GTEX_TPM_LINK)
                else:
                    logging.error('\tDownloading Error!  No file %s', gzip_temp)
            if os.path.exists(gzip_temp):
                logging.info('Found zipped GTEX TPM file: %s', gzip_temp)
                logging.info('\tUnzipping file')
                with gzip.open(gzip_temp, "rb") as f_in, open(self.gtex_tpm_fn, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                logging.info('\tFinished unzipping %s', self.gtex_tpm_fn)
                if os.path.exists(self.gtex_tpm_fn):
                    logging.info('\tRemoving temp file: %s', gzip_temp)
                    os.remove(gzip_temp)
        else:
            logging.info('Found: %s', self.gtex_tpm_fn)

        # Just looking at median values by tissue (SMTSD) greatly reduces file size
        if not os.path.exists(self.gtex_tpm_med_fn):
            gzip_temp = self.gtex_tpm_med_fn + '.gz'
            if not os.path.exists(gzip_temp):
                logging.info('Getting GTEx median TPM by tissue File.')
                logging.info('\tDownloading from: %s', self.GTEX_TPM_MEDIAN_LINK)
                logging.info('\tDownloading to : %s', gzip_temp)
                urllib.request.urlretrieve(self.GTEX_TPM_MEDIAN_LINK, gzip_temp)
                if gzip_temp:
                    logging.info('\tFinished Downloading: %s', self.GTEX_TPM_MEDIAN_LINK)
                else:
                    logging.error('\tDownloading Error!  No file %s', gzip_temp)
            if os.path.exists(gzip_temp):
                logging.info('Found zipped GTEX Median TPM file: %s', gzip_temp)
                logging.info('\tUnzipping file')
                with gzip.open(gzip_temp, "rb") as f_in, open(self.gtex_tpm_med_fn, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                logging.info('\tFinished unzipping %s', self.gtex_tpm_med_fn)
                if os.path.exists(self.gtex_tpm_med_fn):
                    logging.info('\tRemoving temp file: %s', gzip_temp)
                    os.remove(gzip_temp)
        else:
            logging.info('Found: %s', self.gtex_tpm_med_fn)

        #  Sample annotation data
        if not os.path.exists(self.gtex_attr_fn):
            logging.info('downloading {}'.format(self.GTEX_ATTR_LINK))
            urllib.request.urlretrieve(self.GTEX_ATTR_LINK, self.gtex_attr_fn)
        else:
            logging.info('Found: %s', self.gtex_attr_fn)
        if not os.path.exists(self.gtex_attr_desc_fn):
            logging.info('downloading {}'.format(self.GTEX_ATTR_DESC_LINK))
            urllib.request.urlretrieve(self.GTEX_ATTR_DESC_LINK, self.gtex_attr_desc_fn)
        else:
            logging.info('Found: %s', self.gtex_attr_desc_fn)
        if not os.path.exists(self.gtex_pheno_fn):
            logging.info('downloading {}'.format(self.GTEX_PHENO_LINK))
            urllib.request.urlretrieve(self.GTEX_PHENO_LINK, self.gtex_pheno_fn)
        else:
            logging.info('Found: %s', self.gtex_pheno_fn)

    def low_median_tissue_expression_table(self, cutoff=0.01):
        """
        Use median expression table to find genes with most values below cutoff.
        """
        gm = pandas.read_csv(self.gtex_tpm_med_fn, sep='\t', skiprows=2)
        gmt = gm.set_index(['Description', 'gene_id']).T
        gmd = gmt.describe()
        low_tissue_expression = gmd.T[gmd.T['75%'] <= 0.01]
        return low_tissue_expression

    def low_expression_genes_table(self, gtex=None, cutoff=0.01, sample_fraction=0.75, sample_list=None):
        """
        Get lists of genes with low expression.

        All genes where at least <sample_fraction> of genes are expressing below <cutoff>
        """
        if gtex is None:
            # load full table of all genes
            gtex = self.gtex_gene_table(None)

        desc_cols = [c for c in gtex.columns if not c.startswith('GTEX')]
        if sample_list:
            gtex = gtex[desc_cols + sample_list]
        else:
            sample_list = [c for c in gtex.columns if c not in desc_cols]
        num_samples = len(sample_list)
        sample_cutoff = min(1, int(round(num_samples * sample_fraction)))
        filter_func = lambda row: sample_cutoff > len([sample for sample in sample_list if row[sample] < cutoff])
        ft = gtex[gtex.apply(filter_func, axis=1)]
        return ft[desc_cols]

    def preprocess(self, cache_dir=None):
        cache_dir = self.cache_dir if not cache_dir else cache_dir
        if not os.path.exists(os.path.join(cache_dir, 'preprocess')):
            os.makedirs(os.path.join(cache_dir, 'preprocess'))

        logging.info('Finding low median expression set')
        low = self.low_median_tissue_expression_table()
        low_fn = os.path.join(cache_dir, 'preprocess', 'low_median_tissue_expression.tsv')
        logging.info('Saving low median expression file: %s', low_fn)
        low.to_csv(low_fn, sep='\t', header=True, index=True)
        low = pandas.read_csv(low_fn, sep='\t')

        logging.info('Filtering out low expression genes before applying log2 transform')
        low_symbols = sorted(low['Description'].tolist())

        logging.info('Loading the full GTEX table...')
        gtex = self.gtex_gene_table(None)
        logging.info('\t finished loading gtex of %s genese', len(gtex))
        logging.info('Removing %s low tissue expression genes', len(low_symbols))
        gtex = gtex[~(gtex['Description'].isin(low_symbols))]

        logging.info(gtex.info())
        # Set up multindex values
        logging.info('resetting index')
        logging.info(gtex.columns)
        logging.info(gtex.index)

        gtex['ensid'] = gtex['Name'].apply(lambda x: x.split('.')[0])
        gtex['ens_ver'] = gtex['Name'].apply(lambda x: x.split('.')[-1])
        gtex.rename(columns={'Description': 'Symbol'}, inplace=True)
        gtex.drop('Name', inplace=True, axis=1)
        gtex.set_index(['Symbol', 'ensid', 'ens_ver'], inplace=True)

        logging.info(gtex.info())
        logging.info(gtex.describe())
        logging.info('applying safelog2 transform')
        # applying log2 transformation
        safelog2 = lambda x: numpy.log2(x) if x else numpy.nan
        gtex = gtex.applymap(safelog2)

        log2fn = os.path.join(cache_dir, 'preprocess', 'GTEX_log2.tsv')
        logging.info('Saving log2 version: %s', log2fn)
        gtex.to_csv(os.path.join(cache_dir, 'preprocess', 'GTEX_log2.tsv'), index=True, header=True, sep='\t')

        logging.info('Calculating %s', 'gtex_log2tpm.describe.tsv')
        out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_all.describe.tsv')
        desc = gtex.T.describe().T
        desc['skew'] = gtex.apply(lambda row: row.skew(), axis=1)
        desc['kurtosis'] = gtex.apply(lambda row: row.kurtosis(), axis=1)
        desc.to_csv(out_fn, sep='\t', header=True, index=True)
        logging.info('wrote %s', os.path.abspath(out_fn))

        td = self.tissues()
        for tissue in sorted(td.keys()):
            logging.info("Calculating '%s' statistics", tissue)
            subtypes = sorted(td[tissue]['subtypes'].keys())
            for subtype in subtypes:
                out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_{}__{}.describe.tsv'.format(tissue, subtype))
                logging.info('Calculating %s: %s', tissue, subtype)
                samples = [s for s in td[tissue]['subtypes'][subtype] if s in gtex.columns]
                gts = gtex[samples]
                desc = gts.T.describe().T
                desc['skew'] = gts.apply(lambda row: row.skew(), axis=1)
                desc['kurtosis'] = gts.apply(lambda row: row.kurtosis(), axis=1)
                desc.to_csv(out_fn, sep='\t', header=True, index=True)
                logging.info('wrote %s', os.path.abspath(out_fn))
            if len(subtypes) > 1:
                out_fn = os.path.join(cache_dir, 'preprocess', 'gtex_log2tpm_{}__{}.describe.tsv'.format(tissue, 'all'))
                samples = [s for s in td[tissue]['samples'] if s in gtex.columns]
                gt = gtex[samples]
                desc = gt.T.describe().T
                desc['skew'] = gt.apply(lambda row: row.skew(), axis=1)
                desc['kurtosis'] = gt.apply(lambda row: row.kurtosis(), axis=1)
                desc.to_csv(out_fn, sep='\t', header=True, index=True)
                logging.info('wrote %s', os.path.abspath(out_fn))


if __name__ == '__main__':
    gtex = GtexData()
    if os.path.exists('/gsc/resources/expression/GTEx/v7'):
        gtex.preprocess(cache_dir='/gsc/resources/expression/GTEx/v7')
    else:
        gtex.preprocess()
