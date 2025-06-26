
import os

import pandas as pd
import xgboost as xgb

from ppm.constants import (
    PROTEOMICS_FEATURES,
    TRAIN_FEATURES,
    TRANSCRIPT_FEATURES,
)
from ppm.train_utils import add_shap_values, run_hpt
N_ITERATIONS = 10
CV_SIZE = 10



def train_spliced_models(config):
    """ Functions to train models to distinguish peptides identified in IP. vs random.
    """
    feature_set = TRAIN_FEATURES['spliced'] + TRANSCRIPT_FEATURES[config.cell_line] + PROTEOMICS_FEATURES[config.cell_line]
    # total_pep_df = pd.read_parquet(f'background_analsis/trainingDatasets/{cell_line}/spliced.parquet')
    total_pep_dfs = []
    for pep_len in [9, 10, 11, 12]:
        total_pep_df = pd.concat([
            pd.read_parquet(f'{config.output_folder}/trainingDatasets//{pep_len}/{parq_file}')
            for parq_file in os.listdir(f'{config.output_folder}/trainingDatasets/{pep_len}')
        ])
        total_pep_dfs.append(total_pep_df)
    total_pep_df = pd.concat(total_pep_dfs)
    total_pep_df['peptideCount'] = total_pep_df.groupby('peptide')['sr1'].transform('count')
    total_pep_df = create_cv_groups(total_pep_df)

    params = run_hpt(total_pep_df, feature_set)
    scored_pep_df, feature_df_1 = run_cv_training(
        total_pep_df, feature_set, 'prediction_xgb0', 0, params=params, total_pep_df=total_pep_df
    )


    for idx in range(N_ITERATIONS):
        unique_pep_df = scored_pep_df.sort_values(
            by=f'prediction_xgb{idx}', ascending=False
        ).drop_duplicates(subset=['peptide'])
        unique_pep_df = unique_pep_df.reset_index(drop=True)
        params = run_hpt(unique_pep_df, feature_set)

        scored_pep_df, feature_df_2 = run_cv_training(
            unique_pep_df, feature_set, f'prediction_xgb{idx+1}', idx+1, params=params,
            save_folder=config.output_folder, total_pep_df=scored_pep_df, explain=(idx==(N_ITERATIONS-1))
        )
    
    feature_df_2.to_csv(f'{config.output_folder}/models/importances.csv', index=False)
    unique_pep_df = scored_pep_df.sort_values(
        by=f'prediction_xgb{idx+1}', ascending=False
    ).drop_duplicates(subset=['peptide'])

    unique_pep_df.to_csv(
        f'{config.output_folder}/unique_peps_scored.csv', index=False,
    )


def create_cv_groups(total_pep_df):
    """ Function to divide the training data into cross validation groups.
    """
    pos_peps = total_pep_df[total_pep_df['label'] == 1][['peptide']].drop_duplicates()
    neg_peps = total_pep_df[total_pep_df['label'] == 0][['peptide']].drop_duplicates()
    cv_peps = []
    for peptide_df in [pos_peps, neg_peps]:
        peptide_df = peptide_df.sample(frac=1, random_state=42).reset_index(drop=True)
        peptide_df['cvGroup'] = peptide_df.index % CV_SIZE
        cv_peps.append(peptide_df)

    total_pep_df = pd.merge(total_pep_df, pd.concat(cv_peps), how='inner', on='peptide')
    return total_pep_df

def run_cv_training(
        unique_pep_df, feature_set, result_col, idx,
        params={}, save_folder=None, total_pep_df=None, explain=False
    ):
    test_dfs = []
    importances = {'feature': feature_set}
    for i in range(CV_SIZE):
        clf = xgb.XGBClassifier(**params)
        train_df = unique_pep_df[unique_pep_df['cvGroup'] != i]
        if f'prediction_xgb{idx-1}' in unique_pep_df.columns:
            print(f'starting with {train_df.shape[0]}')
            cut_off = train_df[f'prediction_xgb{idx-1}'].quantile(0.1)
            train_df = train_df[
                (train_df['label'] == 0) |
                (train_df[f'prediction_xgb{idx-1}'] > cut_off)
            ]
            print(f'reduced with cut off {cut_off} to {train_df.shape[0]}')
        if total_pep_df is None:
            test_df = unique_pep_df[unique_pep_df['cvGroup'] == i]
        else:
            test_df = total_pep_df[total_pep_df['cvGroup'] == i]

        clf.fit(train_df[feature_set], train_df['label'])

        importances[f'model_{i}'] = clf.feature_importances_
        test_df[result_col] = clf.predict_proba(test_df[feature_set])[:,1]
    
        if explain:
            test_df = add_shap_values(test_df, clf, feature_set)
        test_dfs.append(test_df)

        if save_folder is not None:
            clf.save_model(f'{save_folder}/models/clf_spliced_{i}.json')

    feat_imp_df = pd.DataFrame(importances)

    return pd.concat(test_dfs), feat_imp_df


