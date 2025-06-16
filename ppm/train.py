
import os

import pandas as pd
from sklearn.model_selection import RandomizedSearchCV, LeaveOneGroupOut
import xgboost as xgb

from ppm.constants import (
    TRAIN_FEATURES,
    CRYPTIC_STRATA,
    TRANSCRIPT_FEATURES,
    N_CV_GROUPS,
)
from ppm.train_utils import add_shap_values, run_hpt

def train_models(config):
    """ Functions to train models to distinguish peptides identified in IP. vs random.
    """
    if config.model == 'cryptic':
        feature_set = TRAIN_FEATURES['cryptic'] + TRANSCRIPT_FEATURES[config.cell_line]
        strata = CRYPTIC_STRATA
    else:
        feature_set = TRAIN_FEATURES['canonical'] + TRANSCRIPT_FEATURES[config.cell_line]
        strata = ['canonical']

    strat_dfs = []
    for stratum in strata:
        strat_dfs.append(pd.concat([
            pd.read_csv(f'{config.output_folder}/trainingDatasets/{stratum}_{pep_len}.csv')
            for pep_len in [
                9,
                10,
                11,
                12,
            ] if os.path.exists(
                f'{config.output_folder}/trainingDatasets/{stratum}_{pep_len}.csv'
            )
        ]))
    total_pep_df = pd.concat(strat_dfs)

    total_pep_df = create_cv_groups(total_pep_df)

    params = run_hpt(total_pep_df, feature_set)
    print(params)
    scored_pep_df = run_cv_training(
        total_pep_df, feature_set, 'prediction_xgb1', config, params=params
    )

    unique_pep_df = scored_pep_df.sort_values(
        by='prediction_xgb1', ascending=False
    ).drop_duplicates(subset=['peptide'])
    unique_pep_df = unique_pep_df.reset_index(drop=True)
    params = run_hpt(unique_pep_df, feature_set)
    print(params)

    unique_pep_df, all_scored_pep_df = run_cv_training(
        unique_pep_df, feature_set, 'prediction_xgb2', config, params=params,
        save_key='combined', explain=True, score_df=scored_pep_df,
    )

    unique_pep_df.to_csv(
        f'{config.output_folder}/unique_peps_scored.csv', index=False,
    )
    all_scored_pep_df.to_csv(
        f'{config.output_folder}/all_peps_scored.csv', index=False,
    )



def create_cv_groups(total_pep_df):
    """ Function to divide the training data into cross validation groups.
    """
    pos_peps = total_pep_df[total_pep_df['label'] == 1][['peptide']].drop_duplicates()
    neg_peps = total_pep_df[total_pep_df['label'] == 0][['peptide']].drop_duplicates()
    cv_peps = []
    for peptide_df in [pos_peps, neg_peps]:
        peptide_df = peptide_df.sample(frac=1, random_state=42).reset_index(drop=True)
        peptide_df['cvGroup'] = peptide_df.index % N_CV_GROUPS
        cv_peps.append(peptide_df)

    antigen_df = total_pep_df[['proteinID']].drop_duplicates()

    antigen_df = antigen_df.sample(frac=1, random_state=42).reset_index(drop=True)
    antigen_df['cvGroup'] = antigen_df.index % N_CV_GROUPS

    total_pep_df = pd.merge(total_pep_df, antigen_df, how='inner', on='proteinID')
    return total_pep_df

def run_cv_training(
        all_df, feature_set, result_col, config,
        params={}, save_key=None, explain=False, score_df=None,
    ):
    test_dfs = []
    test_score_dfs = []
    importances = {'feature': feature_set}

    for i in range(N_CV_GROUPS):
        clf = xgb.XGBClassifier(**params)
        train_df = all_df[all_df['cvGroup'] != i]
        test_df = all_df[all_df['cvGroup'] == i]
        if score_df is not None:
            test_score_df = score_df[score_df['cvGroup'] == i]
        clf.fit(train_df[feature_set], train_df['label'])

        test_df[result_col] = clf.predict_proba(test_df[feature_set])[:,1]

        if explain:
            test_df = add_shap_values(test_df, clf, feature_set)

        test_dfs.append(test_df)
        if score_df is not None:
            test_score_dfs.append(test_score_df)

        if save_key is not None:
            importances[f'model_{i}'] = clf.feature_importances_
            clf.save_model(f'{config.output_folder}/models/clf_{save_key}_{i}.json')


    if save_key is not None:
        importance_df = pd.DataFrame(importances)
        print(importance_df)
        importance_df.to_csv(f'{config.output_folder}/models/importances.csv', index=False)
    if score_df is None:
        return pd.concat(test_dfs)

    return pd.concat(test_dfs), pd.concat(test_score_dfs)
