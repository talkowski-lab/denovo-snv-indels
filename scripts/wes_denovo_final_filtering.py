import datetime
import pandas as pd
import numpy as np
import pandas as pd
import argparse
import ast
import sys
import os

def str2bool(v):
    return ast.literal_eval(v.capitalize())

parser = argparse.ArgumentParser(description="WES de novo final filtering")
parser.add_argument("--de-novo-merged", required=True)
parser.add_argument("--cohort-prefix", required=True)
parser.add_argument("--vqslod-cutoff-snv", type=int, required=True)
parser.add_argument("--vqslod-cutoff-indel", type=int, required=True)
parser.add_argument("--af-threshold", type=float, required=True)
parser.add_argument("--ad-alt-threshold", type=int, required=True)
parser.add_argument("--cores", required=True)
parser.add_argument("--mem", type=float, required=True, help="Memory in GB")
parser.add_argument("--single-variant", type=str2bool, required=True)

args = parser.parse_args()

de_novo_merged = args.de_novo_merged
cohort_prefix = args.cohort_prefix
vqslod_cutoff_snv = args.vqslod_cutoff_snv
vqslod_cutoff_indel = args.vqslod_cutoff_indel
MAF_thresh = args.af_threshold
AD_alt_threshold = args.ad_alt_threshold
cores = args.cores
mem = int(np.floor(args.mem))
single_variant = args.single_variant

df = pd.read_csv(de_novo_merged, sep='\t')
df['VarKey'] = df[['ID', 'proband.s']].astype(str).agg(':'.join, axis=1)
df['Consequence'] = df.Consequence.replace({np.nan: '[]'}).apply(ast.literal_eval).agg(','.join)

# Filter to only hets
df = df[~df['proband_entry.GT'].isin(['1/1','1|1'])]

# Filter out LOW confidence
df = df[df.confidence!='LOW']

if 'VQSLOD' in df.columns:
    # Filter SNPs on VQSLOD
    df = df[(df.isIndel)|((df.VQSLOD >= vqslod_cutoff_snv) | df.VQSLOD.isna())] 

    # Filter indels on VQSLOD
    df = df[(df.isSNV)|((df.VQSLOD >= vqslod_cutoff_indel) | df.VQSLOD.isna())] 

# Set frequency threshold

# Filter on dataset AF
df = df[df.cohort_AF<= MAF_thresh]

# Filter on gnomAD
df = df[(df.gnomad_non_neuro_AF.isna())|(df.gnomad_non_neuro_AF<=MAF_thresh)]

# Filter on AD_alt
df['proband_entry.AD_alt'] = df['proband_entry.AD'].apply(ast.literal_eval).str[1].astype(int)
df = df[df['proband_entry.AD_alt']>=AD_alt_threshold]

if single_variant:
    # Pick one variant per gene per sample
    df['SAMPLE_GENE'] = df[['proband.s', 'SYMBOL']].astype(str).agg('-'.join, axis=1)
    df = pd.concat([df[df.isCoding].sort_values(['SAMPLE_GENE','csq_score'], ascending=False).groupby('SAMPLE_GENE').head(1).reset_index(drop=True),
                    df[~df.isCoding]])

df.to_csv(cohort_prefix+'_de_novo_filtered_final.tsv', sep='\t', index=False)