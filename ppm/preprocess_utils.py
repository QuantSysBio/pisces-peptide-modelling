""" Functions which are used in processing of both canonical/cryptic and spliced peptides.
"""
from math import ceil, floor
import os

import numpy as np
import polars as pl

from ppm.constants import (
    COMMON_FEATURES,
    TRANSCRIPT_FEATURES,
    SPLICE_SPECIFIC_FEATURES,
    STRATUM_SPECIFIC_FEATURES,
)
np.random.seed(42)


def merge_orf_level_data(pep_df, config, stratum, label, id_columns):
    """ Function to merge in data from the reading frame level
    """
    if stratum == 'spliced':
        prot_df = pl.read_parquet(f'{config.antigen_folder}/canonical.parquet')
        final_can_id_df = pl.read_csv(
            f'{config.canonical_results}/unique_peps_scored.csv',
            columns=['peptide', 'proteinID', 'label']
        )
        final_can_id_df = final_can_id_df.filter(pl.col('label').eq(1))
        final_can_id_df = final_can_id_df.group_by('proteinID').agg(
            pl.col('peptide').n_unique().alias('nCanonicalPeptides'),
            pl.col('peptide').alias('canonicalPeptides'),
        )
        prot_df = prot_df.join(final_can_id_df, how='left', on=['proteinID'])
        prot_df = prot_df.with_columns(pl.col('nCanonicalPeptides').fill_null(0))
        pep_df = pep_df.select(['peptide', 'proteinID'] + SPLICE_SPECIFIC_FEATURES)
    else:
        prot_df = pl.read_parquet(f'{config.antigen_folder}/{stratum}.parquet')
        pep_df = pep_df.select(['peptide', 'proteinID'])
    if stratum == 'intergenic': #TODO
        prot_df = prot_df.with_columns(
            *[pl.lit(None).alias(feature) for feature in TRANSCRIPT_FEATURES[config.cell_line]]
        )

    select_features = (
        COMMON_FEATURES +
        TRANSCRIPT_FEATURES[config.cell_line] +
        STRATUM_SPECIFIC_FEATURES[stratum]
    )

    pep_df = pep_df.join(
        prot_df.select(select_features), how='inner', on='proteinID',
    )
    pep_df = pep_df.with_columns(
        pl.lit(label).alias('label')
    )

    final_columns = (
        id_columns + [
            'protSeq', 'label',
            'rnaSeq', 'iupred3_preds', 'proteinHydrophobicity'
        ] + TRANSCRIPT_FEATURES[config.cell_line]
        + STRATUM_SPECIFIC_FEATURES[stratum]
    )

    return pep_df.select(final_columns)

def get_sampled_negative_peps(config, pep_len, is_cryptic, stratum):
    pep_dfs = []
    sample_df = pl.read_csv(f'{config.background_folder}/sample_ratios/ratio_{pep_len}.csv')
    sample_vals = dict(zip(sample_df['dataset'].to_list(), sample_df['fraction'].to_list()))
    print(sample_vals)
    min_count = get_min_counts(config, pep_len, stratum, is_cryptic, sample_vals)

    for dataset in os.listdir(f'{config.background_folder}/remapped/{pep_len}'):
        if _check_dataset(dataset, config.cell_line):
            try:
                pep_df = _get_pep_df(config, is_cryptic, pep_len, dataset, stratum)
            except:
                print(f'Failure for {dataset}, length {pep_len}')
                continue
            sample_goal = min_count*sample_vals[dataset]
            if stratum == 'spliced':
                sample_goal /= 100
            sample_goal_frac = sample_goal - floor(sample_goal)
            if np.random.random() > sample_goal_frac:
                sample_goal = ceil(sample_goal)
            else:
                sample_goal = floor(sample_goal)

            print(f'{dataset} : {min_count*sample_vals[dataset]} {sample_vals[dataset]} {pep_df.shape[0]}')
            if pep_df.shape[0] and pep_df.shape[0] > sample_goal:
                pep_df = pep_df.sample(n=sample_goal, seed=1)

            if pep_df.shape[0]:
                pep_dfs.append(pep_df)

    return pep_dfs

def get_min_counts(config, pep_len, stratum, is_cryptic, sample_vals):
    sampled_min_count = 1_000_000_000
    for dataset in os.listdir(f'{config.background_folder}/remapped/{pep_len}'):
        if _check_dataset(dataset, config.cell_line):
            if not os.path.exists(f'{config.background_folder}/remapped/{pep_len}/{dataset}/peptides.csv'):
                print(f'Missing for {dataset}, length {pep_len}')
                continue

            try:
                pep_df = _get_pep_df(config, is_cryptic, pep_len, dataset, stratum)
            except:
                print(f'Failure for {dataset}, length {pep_len}')
                continue
            pep_count = pep_df.shape[0]/sample_vals[dataset]
            print(f'{dataset} : {pep_count} {sample_vals[dataset]} {pep_df.shape[0]}')
            if pep_count and pep_count < sampled_min_count:
                sampled_min_count = pep_count
    print(sampled_min_count)
    if sampled_min_count == 1_000_000_000:
        return 0
    return sampled_min_count


def _check_dataset(dataset, cell_line):
    """ Helper function to check if a dataset if for the K562 or B721.221 cell line.
    """
    if cell_line == 'K562':
        return dataset.startswith('K562')
    elif cell_line == 'B721.221':
        return dataset.startswith('Sarkizova') or dataset.startswith('Abelin')
    else:
        raise ValueError(f'Unknown cell line {cell_line}')
    

def _get_pep_df(config, is_cryptic, pep_len, dataset, stratum):
    pep_df = pl.read_csv(f'{config.background_folder}/remapped/{pep_len}/{dataset}/peptides.csv')
    if is_cryptic:
        det_df = pl.read_csv(f'{config.background_folder}/remapped/{pep_len}/{dataset}/details/cryptic.csv')
    else:
        det_df = pl.read_csv(f'{config.background_folder}/remapped/{pep_len}/{dataset}/details/{stratum}.csv')

    pep_df = pep_df.unique('peptide')
    det_df = det_df.unique('peptide')

    if stratum == 'spliced':
        pep_df = pep_df.filter(pl.col(f'nCrypticProteins').eq(0) & pl.col(f'nSplicedProteins').gt(0) )
        det_df = det_df.with_columns(
            pl.col('sr1').fill_null('NA')
        )
    else:
        pep_df = pep_df.filter(pl.col(f'{stratum}_nProteins').gt(0) & pl.col(f'nSplicedProteins').eq(0))
    pep_df = pep_df.filter(
        pl.col('fusion_nProteins').eq(0) & pl.col('mutation_nProteins').eq(0) &
        pl.col('TrEMBL_nProteins').eq(0)
    )
    return pep_df.select(['peptide']).join(det_df, how='inner', on='peptide')
