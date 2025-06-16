
import multiprocessing as mp
import os

import pandas as pd
import polars as pl
import xgboost as xgb

from ppm.constants import (
    CRYPTIC_STRATA,
    TRAIN_FEATURES,
    N_CV_GROUPS,
    SPLICED_FEATURES,
    TRANSCRIPT_FEATURES
)
from ppm.preprocess import gather_positive_samples, add_features
from ppm.preprocess_spliced import gather_positive_samples_spliced

def score_multi_mappers(config):
    for pep_len in range(7, 16):
        if config.model == 'spliced':
            preprocess_spliced_multi_mappers(config, pep_len)
        if config.model == 'cryptic':
            preprocess_cryptic_multi_mappers(config, pep_len)

    if config.model == 'spliced':
        score_spliced_multi_mappers(config)
    elif config.model == 'cryptic':
        score_cryptic_multi_mappers(config)

def score_cryptic_multi_mappers(config):
    strat_dfs = []
    feature_set = TRAIN_FEATURES[config.model] + TRANSCRIPT_FEATURES[config.cell_line] 

    for stratum in CRYPTIC_STRATA:
        strat_dfs.append(pd.concat([
            pd.read_csv(f'{config.output_folder}/mmDatasets/{stratum}_{pep_len}.csv')
            for pep_len in range(7,16) if os.path.exists(
                f'{config.output_folder}/mmDatasets/{stratum}_{pep_len}.csv'
            )
        ]))
    total_pep_df = pd.concat(strat_dfs)

    for model_idx in range(N_CV_GROUPS):
        clf = xgb.XGBClassifier()
        clf.load_model(f'{config.output_folder}/models/clf_combined_{model_idx}.json')
        print(feature_set)
        total_pep_df[f'score_{model_idx}'] = clf.predict_proba(total_pep_df[feature_set])[:,1]

    total_pep_df['meanScore'] = total_pep_df[
        [f'score_{model_idx}' for model_idx in range(N_CV_GROUPS)]
    ].mean(axis=1)
    total_pep_df = total_pep_df.sort_values(by='meanScore', ascending=False)
    total_pep_df['stratum'] = total_pep_df['stratum'].apply(lambda x: CRYPTIC_STRATA[x])
    unique_pep_df = total_pep_df.drop_duplicates(subset=['peptide'])

    unique_pep_df.to_csv(f'{config.output_folder}/mm_scored.csv', index=False)

def score_spliced_multi_mappers(config):
    feature_set = SPLICED_FEATURES + TRANSCRIPT_FEATURES[config.cell_line]

    total_pep_dfs = []
    for pep_len in range(7,16):
        if len(os.listdir(f'{config.output_folder}/mmDatasets/{pep_len}')):
            total_pep_df = pd.concat([
                pd.read_parquet(f'{config.output_folder}/mmDatasets/{pep_len}/{parq_file}')
                for parq_file in os.listdir(f'{config.output_folder}/mmDatasets/{pep_len}')
            ])
            total_pep_dfs.append(total_pep_df)
    total_pep_df = pd.concat(total_pep_dfs)

    for model_idx in range(10):
        clf = xgb.XGBClassifier()
        clf.load_model(f'{config.output_folder}/models/clf_spliced_{model_idx}.json')
        total_pep_df[f'score_{model_idx}'] = clf.predict_proba(
            total_pep_df[feature_set]
        )[:,1]

    total_pep_df['meanScore'] = total_pep_df[
        [f'score_{model_idx}' for model_idx in range(10)]
    ].mean(axis=1)
    total_pep_df = total_pep_df.sort_values(by='meanScore', ascending=False)
    unique_pep_df = total_pep_df.drop_duplicates(subset=['peptide'])
    unique_pep_df.to_csv(f'{config.output_folder}/mm_scored.csv', index=False)

def preprocess_cryptic_multi_mappers(config, pep_len):
    for stratum in CRYPTIC_STRATA:
        pos_pep_df = gather_positive_samples(stratum, True, config, pep_len, True)
        if not pos_pep_df.shape[0]:
            continue

        # Additional features for cryptic and canonical peptides.
        pos_pep_df = add_features(pos_pep_df)
        pos_pep_df = pos_pep_df.with_columns(
            pl.lit(CRYPTIC_STRATA.index(stratum)).alias('stratum')
        )

        pos_pep_df.write_csv(f'{config.output_folder}/mmDatasets/{stratum}_{pep_len}.csv')

def preprocess_spliced_multi_mappers(config, pep_len):
    if not os.path.exists(f'{config.output_folder}/mmDatasets/{pep_len}'):
        os.mkdir(f'{config.output_folder}/mmDatasets/{pep_len}')

    gather_positive_samples_spliced(
        config, pep_len, is_mm=True,
    )
