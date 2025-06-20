# from wgs_bagging_pu_clean_clean.ipynb
import gzip
import os
import glob
import io
import ast
import warnings
import seaborn as sns
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import sklearn.metrics
import sklearn.model_selection

from baggingPU import BaggingClassifierPU
from sklearn.metrics import RocCurveDisplay
from google.cloud import storage
from sklearn.ensemble import RandomForestClassifier
from io import StringIO

plt.rcParams.update({'font.size': 10})

import sys
import pandas as pd
import numpy as np
import sklearn.metrics
import sklearn.model_selection
import ast
import pickle 

from sklearn.ensemble import RandomForestClassifier

vcf_metrics_tsv = sys.argv[1]
ultra_rare_inherited_tsv = sys.argv[2]
ultra_rare_parents_tsv = sys.argv[3]
cohort_prefix = sys.argv[4]
var_type = sys.argv[5]
variant_features = sys.argv[6].split(',')
sample_features = sys.argv[7].split(',')
vqslod_cutoff = float(sys.argv[8])
rep_regions = sys.argv[9]
metric = sys.argv[10]  # ['roc-auc', 'accuracy', 'f1', 'fp_fn_ratio', 'neg_log_loss', 'recall', 'precision']
n_estimators_rf = int(sys.argv[11])
n_bags = int(sys.argv[12])
filter_pass_before = ast.literal_eval(sys.argv[13].capitalize())
n_jobs = int(sys.argv[14])

# check for empty sample_features
if sample_features==["",""]:
    sample_features = []

# check for empty variant_features
if variant_features==["",""]:
    variant_features = []

def fp_fn_ratio(y, y_pred):
    FP = ((y==0) & (y_pred==1)).sum()
    TP = ((y==1) & (y_pred==1)).sum()
    FN = ((y==1) & (y_pred==0)).sum()
    TN = ((y==0) & (y_pred==0)).sum()
    FN_rate = FN / (FN + TN)
    FP_rate = FP / (FP + TP)
    if np.isnan(FP_rate):
        FP_rate = FN_rate
    return float(FP_rate / FN_rate)

metrics_to_funcs = {'roc-auc': sklearn.metrics.roc_auc_score, 'accuracy': True, 'f1': sklearn.metrics.f1_score,
                    'fp_fn_ratio': fp_fn_ratio, 'neg_log_loss': sklearn.metrics.log_loss, 
                    'recall': sklearn.metrics.recall_score, 'precision': sklearn.metrics.precision_score}

def load_variants(vcf_metrics_tsv, ultra_rare_inherited_tsv, var_type, use_random_gq=False, pl_filter=True): 
    ultra_rare = pd.read_csv(ultra_rare_inherited_tsv, sep='\t').replace({'.': np.nan, None: np.nan})
    ultra_rare['Indel_type'] = ultra_rare.apply(lambda x: 'Insertion' if (len(x.ALT) - len(x.REF)) > 0 else 'Deletion', axis=1)
    ultra_rare.loc[ultra_rare.TYPE=='SNV', 'Indel_type'] = 'SNV'
    ultra_rare = ultra_rare[ultra_rare.TYPE==var_type]
    # ultra_rare = ultra_rare[~ultra_rare.SAMPLE.isin(outlier_samples)]

    final_output = pd.read_csv(vcf_metrics_tsv, sep='\t').replace({'.': np.nan, None: np.nan})
    final_output['Indel_type'] = final_output.apply(lambda x: 'Insertion' if (len(x.ALT) - len(x.REF)) > 0 else 'Deletion', axis=1)
    final_output.loc[final_output.TYPE=='SNV', 'Indel_type'] = 'SNV'
    final_output = final_output[final_output.TYPE==var_type]
    # final_output = final_output[~final_output.SAMPLE.isin(outlier_samples)]

    ultra_rare = ultra_rare[(~ultra_rare.AD_father.isna()) & (~ultra_rare.AD_mother.isna())].reset_index(drop=True)

    ultra_rare['VarKey'] = ultra_rare.ID + ':' + ultra_rare.SAMPLE

    # ultra_rare['AB_mother'] = ultra_rare.AD_mother.apply(ast.literal_eval).apply(min) / ultra_rare.AD_mother.apply(ast.literal_eval).apply(sum)
    # ultra_rare['AB_father'] = ultra_rare.AD_father.apply(ast.literal_eval).apply(min) / ultra_rare.AD_father.apply(ast.literal_eval).apply(sum)
    # ultra_rare['AB_sample'] = ultra_rare.AD_sample.apply(ast.literal_eval).apply(min) / ultra_rare.AD_sample.apply(ast.literal_eval).apply(sum)
    if pl_filter:
        ultra_rare = ultra_rare[~ultra_rare.PL_sample.isna()]
        ultra_rare = ultra_rare[~ultra_rare.PL_sample.str.contains('None')]
        ultra_rare[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = ultra_rare.PL_sample.replace({np.nan: '[0, 0, 0]'}).str.strip('][').str.split(', ', expand=True).astype(int)
    ultra_rare['GQ_mean'] = ultra_rare[['GQ_sample', 'GQ_mother', 'GQ_father']].mean(axis=1)

    # https://gatk.broadinstitute.org/hc/en-us/articles/360040507131-FilterVcf-Picard-
    # Allele balance is calculated for heterozygotes as the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote genotypes at the same locus are measured independently
    
    for samp in ['sample', 'mother', 'father']:
        final_output[f"DPC_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).apply(sum)
        final_output[f"AB_{samp}"] = final_output[f"AD_{samp}"].apply(ast.literal_eval).str[1] / final_output[f"DPC_{samp}"]
    
    # repetitive from step6, but if final_output isn't the output of step06
    final_output = final_output[(final_output.PL_sample!='.')&(~final_output.PL_sample.str.contains('\.'))]
    try:
        final_output[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = final_output.PL_sample.str.split(",", expand=True).astype(int)
    except:  # NEW 12/12/2024
        final_output[['PL_sample_0.0', 'PL_sample_0.1', 'PL_sample_1.1']] = final_output.PL_sample.replace({np.nan: '[0, 0, 0]'}).str.strip('][').str.split(', ', expand=True).astype(int)
    # NEW 12/12/2024 removed VAF
    # parent_metrics = ['GQ', 'AB', 'DPC', 'VAF']
    parent_metrics = ['GQ', 'AB', 'DPC']
    for metric in parent_metrics:
        final_output[f'{metric}_parent'] = final_output[[f'{metric}_mother', f'{metric}_father']].min(axis=1)

    ultra_rare['GQ_hom'] = ultra_rare.apply(lambda row: row.GQ_mother if row.GT_mother in ['0/0', '0|0']
                                            else row.GQ_father, axis=1)
    # ultra_rare['AB_parent'] = ultra_rare.apply(lambda row: row.AB_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.AB_father, axis=1)
    # ultra_rare['DPC_parent'] = ultra_rare.apply(lambda row: row.DPC_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.DPC_father, axis=1)
    # ultra_rare['VAF_parent'] = ultra_rare.apply(lambda row: row.VAF_mother if row.GT_mother in ['0/0', '0|0'] 
    #                                             else row.VAF_father, axis=1)
    if use_random_gq:
        ultra_rare['GQ_parent'] = ultra_rare[['GQ_hom', 'GQ_random']].min(axis=1)
    else:
        ultra_rare['GQ_parent'] = ultra_rare[['GQ_mother', 'GQ_father']].min(axis=1)
    ultra_rare['AB_parent'] = ultra_rare[['AB_mother', 'AB_father']].min(axis=1)
    ultra_rare['DPC_parent'] = ultra_rare[['DPC_mother', 'DPC_father']].min(axis=1)
    ultra_rare['VAF_parent'] = ultra_rare[['VAF_mother', 'VAF_father']].min(axis=1)
    
    final_output['GQ_parent'] = final_output[['GQ_mother', 'GQ_father']].min(axis=1)
    final_output['AB_parent'] = final_output[['AB_mother', 'AB_father']].min(axis=1)
    final_output['DPC_parent'] = final_output[['DPC_mother', 'DPC_father']].min(axis=1)
    # final_output['VAF_parent'] = final_output[['VAF_mother', 'VAF_father']].min(axis=1)  # NEW 12/12/2024 commented out

    final_output_raw = final_output.copy();
    ultra_rare_raw = ultra_rare.copy();
    return final_output, ultra_rare, final_output_raw, ultra_rare_raw

def filter_variants(final_output, ultra_rare, final_output_raw, ultra_rare_raw, filter_pass=True):
    if (~final_output['VQSLOD'].isna()).sum()!=0:
        final_output = final_output[~final_output['VQSLOD'].isna()].reset_index(drop=True)
        final_output['VQSLOD'] = final_output.VQSLOD.astype(float)
    if (final_output.VQSLOD>vqslod_cutoff).sum()!=0:
        final_output = final_output[final_output.VQSLOD>vqslod_cutoff]

    final_output = final_output[final_output.LEN<=50]
    if filter_pass:
        final_output = final_output[final_output.FILTER=='PASS']

    ultra_rare = ultra_rare[ultra_rare.POLYX<=10]
    ultra_rare = ultra_rare[ultra_rare.LEN<=50]
    try:
        ultra_rare = ultra_rare[~ultra_rare.VQSLOD.isna()]
        ultra_rare = ultra_rare[ultra_rare.VQSLOD>vqslod_cutoff]
    except:
        pass
    ultra_rare = ultra_rare[(~ultra_rare.VAF_sample.isna())&
                            (~ultra_rare.VAF_mother.isna())&
                            (~ultra_rare.VAF_father.isna())].reset_index(drop=True)

    ultra_rare['label'] = 1  # positive
    final_output['label'] = 0
    # overlapping variants set to unlabeled
    ultra_rare = ultra_rare[~ultra_rare.VarKey.isin(np.intersect1d(ultra_rare.VarKey, final_output.VarKey))]
  
    merged_output = pd.concat([ultra_rare, final_output]).reset_index(drop=True)
    merged_output['var_type'] = merged_output.label.map({1: 'ultra-rare', 0: 'unlabeled'})

    return final_output, ultra_rare, merged_output

def BaggingPU(X, y, kf, model, n_bags, n_jobs):
    y_pred_bag = np.zeros(y.size)
    output_bag = np.zeros(y.size)
    classifiers = []
    
    for i, (train_index, test_index) in enumerate(kf.split(X)):
        print(f"----------------- SPLIT {i+1}/{kf.get_n_splits()} -----------------")
        X_train = X.iloc[train_index,:]
        y_train = y.iloc[train_index]
        X_test = X.iloc[test_index,:]
        y_test = y.iloc[test_index]  
#         print(sum(y_train==0))
        bc = BaggingClassifierPU(
            model, 
            n_estimators = n_bags, n_jobs = n_jobs, 
#             max_samples = sum(y_train==0)  # all unlabeled
        )
        bc.fit(X_train, y_train)
        y_pred_bag[test_index] = bc.predict(X_test)
        output_bag[test_index] = bc.predict_proba(X_test)[:,1]
        classifiers.append(bc)
    return y_pred_bag, output_bag, classifiers


def runBaggingPU_RF(X, y, model, merged_output, features, suffix, n_bags=10, n_jobs=-1):
    kf = sklearn.model_selection.KFold(n_splits=5)
    y_pred_bag, output_bag, classifiers_bag = BaggingPU(X, y, kf, model, n_bags=n_bags, n_jobs=n_jobs)

    print('---- {} ----'.format('Bagging PU'))
    print(sklearn.metrics.confusion_matrix(y, y_pred_bag))

    results = pd.DataFrame({'label': y, 
                            'VarKey': merged_output.iloc[X.index].VarKey,
                            f'predict_proba_bag{suffix}': output_bag,
                            f'pred_bag{suffix}': y_pred_bag.astype(int)
                            })
    return results, classifiers_bag

def get_importances_oob_scores(X, y, merged_output, features, suffix, n_estimators_rf=100, n_bags=10, oob_score=True, n_jobs=-1):
    model = RandomForestClassifier(warm_start=True, oob_score=oob_score, n_estimators=n_estimators_rf)
    results, estimators = runBaggingPU_RF(X, y, model, merged_output, features, suffix, n_bags=n_bags, n_jobs=n_jobs)

    importances = pd.DataFrame(np.array([[estimators[j].estimators_[i].feature_importances_ 
                for i in range(len(estimators[j].estimators_))] 
                    for j in range(len(estimators))]).reshape((len(estimators)*len(estimators[0]),len(features))),
                        columns=features)
    oob_scores = pd.DataFrame(np.array([[estimators[j].estimators_[i].oob_score_ 
                for i in range(len(estimators[j].estimators_))] 
                    for j in range(len(estimators))]).T)
    
    return results, estimators, importances, oob_scores

def get_optimized_results(features, results, suffix, metric='fp_fn_ratio'):
    results[f"features{suffix}"] = ', '.join(features)
    results['metric'] = metric

    opt_auc_scores = pd.DataFrame()

    for p_thr in np.arange(0, 0.5, 0.05): 
        y_pred = [1 if x>p_thr else 0 for x in results[f"predict_proba_bag{suffix}"]]
        opt_auc_scores.loc[p_thr, 'roc_auc_score'] = sklearn.metrics.roc_auc_score(results.label, y_pred)

    best_p_thr = round(opt_auc_scores.roc_auc_score.idxmax(), 2) 

    results[f'pred_bag_optimized{suffix}'] = [1 if x>best_p_thr else 0 for x in results[f"predict_proba_bag{suffix}"]]
    results['p_thr'] = best_p_thr
    return results

def runBaggingPU_level_features(merged_output_subset, features, n_estimators_rf, n_bags, suffix='', n_jobs=-1):
    merged_output_subset = merged_output_subset[~merged_output_subset[features].isna().any(axis=1)].reset_index(drop=True)
    X = merged_output_subset[features].sample(frac=1)
    y = merged_output_subset['label'].astype(int).loc[X.index]
    results, estimators_optimized, importances_optimized, oob_scores_optimized = get_importances_oob_scores(X, y, merged_output_subset, features, 
                                                                                                          suffix=suffix, n_estimators_rf=n_estimators_rf, n_bags=n_bags, oob_score=metrics_to_funcs[metric],
                                                                                                          n_jobs=n_jobs)  
    results_optimized = get_optimized_results(features, results, suffix)
        
    merged_output_subset.index = merged_output_subset.VarKey
    results_optimized.index = results_optimized.VarKey

    # filter out non-PASS
    # TODO: hard-coded PASS only filter?
    # merged_results = results_optimized[(~results_optimized.VarKey.isin(merged_output_subset[merged_output_subset.FILTER!='PASS'].VarKey))
    #                                   |(results_optimized.label==1)]
    merged_results = results_optimized.copy()
    return merged_results, estimators_optimized, importances_optimized, oob_scores_optimized  # 12/12/2024 NEW: return all

# sample-level training (ultra-rare inherited)
final_output, ultra_rare, final_output_raw, ultra_rare_raw = load_variants(vcf_metrics_tsv, ultra_rare_inherited_tsv, var_type)

final_output, ultra_rare, merged_output = filter_variants(final_output, ultra_rare, 
                                                          final_output_raw, ultra_rare_raw, filter_pass=filter_pass_before)
# variant-level training (ultra-rare only in one parent)
final_output_var, ultra_rare_var, final_output_var_raw, ultra_rare_var_raw = load_variants(vcf_metrics_tsv, ultra_rare_parents_tsv, var_type, pl_filter=False)

final_output_var, ultra_rare_var, merged_output_var = filter_variants(final_output_var, ultra_rare_var, 
                                                          final_output_var_raw, ultra_rare_var_raw, filter_pass=False)

if var_type=='Indel':
    # Big indels: LEN>3
    # sample-level
    print("---------------------- Running Large Indels (LEN>3) sample-level ----------------------")
    merged_output_big_indels = merged_output[merged_output.LEN>3].reset_index(drop=True)
    big_sample_features = [feature for feature in sample_features if merged_output_big_indels[feature].isna().all()==False]
    big_indels, big_indels_estimators_optimized, big_indels_importances_optimized, big_indels_oob_scores_optimized = runBaggingPU_level_features(merged_output_big_indels, big_sample_features, n_estimators_rf, n_bags, 
                                            suffix='_sample_level', n_jobs=n_jobs)
    # variant-level
    print("---------------------- Running Large Indels (LEN>3) variant-level ----------------------")
    passes_sample_level = big_indels[(big_indels['pred_bag_optimized_sample_level']==1)].VarKey

    merged_output_big_indels_var = merged_output_var[(merged_output_var.LEN>3)
                                                    & ((merged_output_var.label==1) | (merged_output_var.VarKey.isin(passes_sample_level)))].reset_index(drop=True)

    big_variant_features = [feature for feature in variant_features if merged_output_big_indels_var[feature].isna().all()==False]
    big_indels_var, big_indels_var_estimators_optimized, big_indels_var_importances_optimized, big_indels_var_oob_scores_optimized = runBaggingPU_level_features(merged_output_big_indels_var, big_variant_features, n_estimators_rf, n_bags, 
                                            suffix='_variant_level', n_jobs=n_jobs)

    # small indels: LEN<=3
    # sample-level
    print("---------------------- Running Small Indels (LEN<=3) sample-level ----------------------")
    merged_output_small_indels = merged_output[merged_output.LEN<=3].reset_index(drop=True)
    small_sample_features = [feature for feature in sample_features if merged_output_small_indels[feature].isna().all()==False]
    small_indels, small_indels_estimators_optimized, small_indels_importances_optimized, small_indels_oob_scores_optimized = runBaggingPU_level_features(merged_output_small_indels, small_sample_features, n_estimators_rf, n_bags, 
                                            suffix='_sample_level', n_jobs=n_jobs)
    # variant-level
    print("---------------------- Running Small Indels (LEN<=3) variant-level ----------------------")
    passes_sample_level = small_indels[(small_indels['pred_bag_optimized_sample_level']==1)].VarKey

    merged_output_small_indels_var = merged_output_var[(merged_output_var.LEN<=3)
                                                    & ((merged_output_var.label==1) | (merged_output_var.VarKey.isin(passes_sample_level)))].reset_index(drop=True)
    small_variant_features = [feature for feature in variant_features if merged_output_small_indels_var[feature].isna().all()==False]
    small_indels_var, small_indels_var_estimators_optimized, small_indels_var_importances_optimized, small_indels_var_oob_scores_optimized = runBaggingPU_level_features(merged_output_small_indels_var, small_variant_features, n_estimators_rf, n_bags, 
                                            suffix='_variant_level', n_jobs=n_jobs)

    # merge with step06 output
    small_passes_variant_level = small_indels_var[(small_indels_var['pred_bag_optimized_variant_level']==1)].VarKey
    big_passes_variant_level = big_indels_var[(big_indels_var['pred_bag_optimized_variant_level']==1)].VarKey
    passes_variant_level = pd.concat([small_passes_variant_level, big_passes_variant_level])

if var_type=='SNV':
    # sample-level
    print("---------------------- Running SNVs sample-level ----------------------")
    sample_features = [feature for feature in sample_features if merged_output[feature].isna().all()==False]
    if len(sample_features) > 0:
        all_snvs, all_snvs_estimators_optimized, all_snvs_importances_optimized, all_snvs_oob_scores_optimized = runBaggingPU_level_features(merged_output, sample_features, n_estimators_rf, n_bags, 
                                                suffix='_sample_level', n_jobs=n_jobs)
    else: 
        all_snvs = merged_output.copy()
        all_snvs['pred_bag_optimized_sample_level'] = 1    
    # variant-level
    print("---------------------- Running SNVs variant-level ----------------------")
    passes_sample_level = all_snvs[(all_snvs['pred_bag_optimized_sample_level']==1)].VarKey
    merged_output_var = merged_output_var[((merged_output_var.label==1) | (merged_output_var.VarKey.isin(passes_sample_level)))].reset_index(drop=True)
    variant_features = [feature for feature in variant_features if merged_output_var[feature].isna().all()==False]
    if len(variant_features) > 0:
        all_snvs_var, all_snvs_var_estimators_optimized, all_snvs_var_importances_optimized, all_snvs_var_oob_scores_optimized = runBaggingPU_level_features(merged_output_var, variant_features, n_estimators_rf, n_bags, 
                                            suffix='_variant_level', n_jobs=n_jobs)
    else:
        all_snvs_var = merged_output_var.copy()
        all_snvs_var['pred_bag_optimized_variant_level'] = 1    
    passes_variant_level = all_snvs_var[(all_snvs_var['pred_bag_optimized_variant_level']==1)].VarKey

# Save
final_output = final_output_raw.copy()

rep_reg = pd.read_csv(rep_regions, sep='\t', header=None)
final_output['repetitive_region'] = final_output.VarKey.isin(rep_reg[3])
# final_output['multiallelic'] = (final_output.DPC_sample!=final_output.DP_sample)\
#                 |(final_output.DPC_mother!=final_output.DP_mother)\
#                 |(final_output.DPC_father!=final_output.DP_father)

final_output = final_output[(final_output.VarKey.isin(passes_variant_level))
                            | (final_output.TYPE!=var_type)]

base_filename = os.path.basename(vcf_metrics_tsv).split('.tsv.gz')[0]
final_output.to_csv(f"{base_filename}_pu_{var_type}.tsv", sep='\t', index=False)

# 12/12/2024 NEW: feature importance plots
if var_type=='Indel': 
    importances = {'small Indels, sample-level': small_indels_importances_optimized, 
        'small Indels, variant-level': small_indels_var_importances_optimized,
        'big Indels, sample-level': big_indels_importances_optimized,
        'big Indels, variant-level': big_indels_var_importances_optimized}
    fig, ax = plt.subplots(2, 2, figsize=(8, 8));
    for i, (title, importance_df) in enumerate(importances.items()):
        row_n, col_n = i//2, i%2
        sns.violinplot(importance_df, ax=ax[row_n][col_n]);
        ax[row_n][col_n].set_title(title);
        ax[row_n][col_n].set_xticks(ticks=ax[row_n][col_n].get_xticks(), labels=ax[row_n][col_n].get_xticklabels(), 
                        rotation=45, horizontalalignment='right');
if var_type=='SNV':
    importances = {'SNVs, sample-level': all_snvs_importances_optimized, 
        'SNVs, variant-level': all_snvs_var_importances_optimized}
    fig, ax = plt.subplots(1, 2, figsize=(8, 4));
    for i, (title, importance_df) in enumerate(importances.items()):
        sns.violinplot(importance_df, ax=ax[i]);
        ax[i].set_title(title);
        ax[i].set_xticks(ticks=ax[i].get_xticks(), labels=ax[i].get_xticklabels(), 
                        rotation=45, horizontalalignment='right');

plt.tight_layout();
plt.savefig(f"{base_filename}_pu_{var_type}_feature_importances.png");
