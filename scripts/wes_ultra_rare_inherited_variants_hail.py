import hail as hl

import seaborn as sns
import matplotlib.pyplot as plt

from typing import Union
from tenacity import retry, after_log, before_sleep_log, retry_if_exception_type, stop_after_attempt, wait_exponential
from firecloud import fiss
from firecloud.errors import FireCloudServerError
import firecloud.api as fapi
import re
import sys
import time
import math
import ast

import pandas as pd
import numpy as np
import os

import logging
from logging import INFO, DEBUG
logger = logging.getLogger()
logger.setLevel(INFO)

pd.set_option('display.float_format', '{:.10f}'.format)

import matplotlib.pyplot as plt
import seaborn as sns

import argparse

def parse_args():
    p = argparse.ArgumentParser()

    add = p.add_argument
    add("--filt-mt-uri", required=True)
    add("--ped-uri", required=True)

    add("--gnomad-af-threshold", type=float, default=0.001)
    add("--cohort-ac-threshold", type=int, default=20)
    add("--cohort-af-threshold", type=float, default=0.001)
    add("--mem", type=float)

    return p.parse_args()


args = parse_args()

filt_mt_uri = args.filt_mt_uri
ped_uri = args.ped_uri
gnomad_af_threshold = args.gnomad_af_threshold
cohort_ac_threshold = args.cohort_ac_threshold
cohort_af_threshold = args.cohort_af_threshold
mem = args.mem

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.executor.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'}, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

# Functions

# Parent-aware TDT annotations
#
# This draws heavily from Jack Kosmicki's function, 
# https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test,
# and, when summed, should give the same output values
#
# Requires trio dataset and, like Hail's TDT function, does not cover the Y
#
# Also note:
# 1) Uses sexes from ped and assumes fathers are male and mothers are female
# 2) Does not allow fathers to be het when they should be hemizygous
# 3) To match Jack's function, requires fathers to have a genotype even when considering regions where 
#    proband is hemizygous
def parent_aware_t_u_annotations_v4(td):

    # First decide copy state
    td = td.annotate_entries(autosomal_copy_state = ( hl.case()
            .when(td.locus.in_autosome() | td.locus.in_x_par() | (td.is_female == True), True)
            .when(td.locus.in_x_nonpar() & (td.is_female == False), False)
            .default(hl.missing('bool')) ) )
    # Note: the above uses the "is_female" from the ped and not from the dataset itself
    
    # Now annotate t & u values
    td = td.annotate_entries(
        
        t_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_var()) )) ), 1, 0),
        
        t_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_var()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_var() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),
        # I could consider removing any reference at all to father's genotype in this last line    
        
        u_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_ref()) )) ), 1, 0),
        
        u_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_ref()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_ref() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),   
        # Again, could consider removing any reference at all to father's genotype in this last line
        
        t_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0),

        u_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0)        
    )
        
    return (td)

# Load

filt_mt = hl.read_matrix_table(filt_mt_uri)
filt_mt = filt_mt.filter_rows(filt_mt.gnomad_non_neuro_AF>gnomad_af_threshold, keep=False)

filt_mt = filt_mt.filter_rows((filt_mt.info.cohort_AC<=cohort_ac_threshold) |
                              (filt_mt.info.cohort_AF<=cohort_af_threshold))

ped_df = pd.read_csv(ped_uri, sep='\t')

cropped_ped_uri = f"{ped_uri.split('.ped')[0]}_cropped.ped"
ped_df.iloc[:,:6].to_csv(cropped_ped_uri, index=False, sep='\t')
cropped_ped_uri

pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

complete_trio_samples = [trio.s for trio in pedigree.complete_trios()] + \
    [trio.pat_id for trio in pedigree.complete_trios()] + \
    [trio.mat_id for trio in pedigree.complete_trios()]

pedigree = pedigree.filter_to(complete_trio_samples)

td = hl.trio_matrix(filt_mt, pedigree, complete_trios=True)

td = parent_aware_t_u_annotations_v4(td)

td = td.annotate_entries(total_t_from_parents = td.t_from_dad + td.t_from_mom,
                        total_u_from_parents = td.u_from_dad + td.u_from_mom)

# Original function for site-level
tdt_table_filtered = hl.transmission_disequilibrium_test(filt_mt, pedigree)

td = td.annotate_rows(tdt=tdt_table_filtered[td.row_key])

td_mt_uri = f"{filt_mt_uri.split('.mt')[0]}.tdt.mt"
td = td.checkpoint(td_mt_uri, overwrite=True)

td = hl.read_matrix_table(td_mt_uri)

# Ultra-rare inherited

# Ultra-rare inherited
inh_td = td.filter_entries((td.total_t_from_parents==1) & (td.total_u_from_parents==0))
inh_td = inh_td.filter_rows(hl.agg.count_where((hl.is_defined(inh_td.proband_entry.GT)) |
                  (hl.is_defined(inh_td.mother_entry.GT)) |
                  (hl.is_defined(inh_td.father_entry.GT)))>0)

inh_td_mt_uri = f"{td_mt_uri.split('.mt')[0]}.inherited.mt"
inh_td = inh_td.checkpoint(inh_td_mt_uri, overwrite=True)

inh_output_uri = f"{td_mt_uri.split('.mt')[0]}.inherited.tsv.gz"
inh_td.entries().flatten().export(inh_output_uri)

# inh_df = inh_td.entries().flatten().to_pandas()

# inh_df[['REF','ALT']] = inh_df['alleles'].apply(pd.Series)
# inh_df['ID'] = inh_df[['locus','REF','ALT']].astype(str).apply(':'.join, axis=1)

# inh_output_uri = f"{td_mt_uri.split('.mt')[0]}.inherited.tsv.gz"
# inh_df.to_csv(inh_output_uri, sep='\t', index=False)