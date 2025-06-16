import multiprocessing as mp
import os
import re
import warnings

import numpy as np
import polars as pl

from ppm.constants import (
    AMINO_ACIDS,
    ID_COLUMNS,
    SPLICE_SPECIFIC_FEATURES,
)
from ppm.preprocess_utils import merge_orf_level_data, get_sampled_negative_peps

warnings.filterwarnings('ignore')


def process_spliced(config):
    """ Function run to process uniquely mapped peptides attributed to each stratum.
    """
    print(f'\tprocessing cell line {config.cell_line}')
    for pep_len in range(9, 13):
        print(f'\t\t\tprocessing length {pep_len}...')
        process_spliced_len(config, pep_len)


def process_spliced_len(config, pep_len):
    """ Function to process data

    stratum : str
        The stratum being processed
    is_cryptic : bool
        Is the stratum cryptic
    cell_line : str
        The cell line being processed
    """
    os.environ['POLARS_MAX_THREADS'] = '1'
    if not os.path.exists(f'{config.output_folder}/trainingDatasets/{pep_len}'):
        os.mkdir(f'{config.output_folder}/trainingDatasets/{pep_len}')

    pos_pep_df = gather_positive_samples_spliced(config, pep_len)
    gather_negative_samples(config, pep_len, pos_pep_df)


def gather_positive_samples_spliced(config, pep_len, is_mm=False):
    """ Function get positive peptides.
    """
    # Read in data, get 9mer, K562 data only:
    pos_pep_df = pl.read_parquet(
        config.peptides_pq,
        columns=[
            'peptide', 'stratum', 'cellLines', 'splicedProteins'
        ] + SPLICE_SPECIFIC_FEATURES
    )
    pos_pep_df = pos_pep_df.filter(
        pl.col('cellLines').list.contains(config.cell_line) &
        pl.col('peptide').str.len_chars().eq(pep_len)
    ).rename({'splicedProteins': 'proteinID'})

    if is_mm:
        pos_pep_df = pos_pep_df.filter(
            pl.col('stratum').eq('multi-mapped')
        )
    else:
        pos_pep_df = pos_pep_df.filter(
            pl.col('stratum').eq('spliced')
        )

    # Explode columns with different possible accessions listed:
    pos_pep_df = pos_pep_df.explode(
        ['proteinID'] + SPLICE_SPECIFIC_FEATURES
    )
    pos_pep_df = pos_pep_df.with_columns(
        pl.col('interveningSeqLengths').cast(pl.Int64)
    )
    pos_pep_df = merge_orf_level_data(
        pos_pep_df, config, 'spliced', 1, ID_COLUMNS + SPLICE_SPECIFIC_FEATURES
    )
    add_features_mp(pos_pep_df, config, 1, pep_len, is_mm)

    return pos_pep_df.select(['peptide'])

def gather_negative_samples(config, pep_len, pos_peps):
    """ Function to get random background peptides for a given stratum and cell line.
    """
    # Collect negative samples from relevant cell line.
    pep_dfs = get_sampled_negative_peps(config, pep_len, False, 'spliced')
    neg_pep_df = pl.concat(pep_dfs)

    neg_pep_df = neg_pep_df.join(
        pos_peps.with_columns(pl.lit('flag').alias('uniqueFlag')), how='left', on='peptide'
    )
    neg_pep_df = neg_pep_df.filter(pl.col('uniqueFlag').is_null()).drop('uniqueFlag')


    neg_pep_df = neg_pep_df.with_columns(
        pl.col('sr1').str.split(' '),
        pl.col('interveningSeqLengths').str.split(' ').cast(pl.List(pl.Int64)),
        pl.col('splicedProteins').str.split(' '),
        pl.col('sr1_Index').str.split(' ').cast(pl.List(pl.Int64)),
        pl.col('sr2_Index').str.split(' ').cast(pl.List(pl.Int64)),
        pl.col('isForward').str.split(' ').cast(pl.List(pl.Int8)),
    )
    print(neg_pep_df.filter(
        pl.col('sr1').list.len().ne(pl.col('interveningSeqLengths').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('splicedProteins').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('sr1_Index').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('sr2_Index').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('isForward').list.len()) 
    ))
    neg_pep_df = neg_pep_df.with_columns(
        pl.col('sr1').map_elements(lambda x : [a for a in x if a], return_dtype=pl.List(pl.String)),
        pl.col('splicedProteins').map_elements(lambda x : [a for a in x if a], return_dtype=pl.List(pl.String)),
    )
    print(neg_pep_df.filter(
        pl.col('sr1').list.len().ne(pl.col('interveningSeqLengths').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('splicedProteins').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('sr1_Index').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('sr2_Index').list.len()) |
        pl.col('sr1').list.len().ne(pl.col('isForward').list.len()) 
    ))
    print(neg_pep_df)
    neg_pep_df = neg_pep_df.explode(
        SPLICE_SPECIFIC_FEATURES + ['splicedProteins']
    )
    neg_pep_df = neg_pep_df.rename({'splicedProteins': 'proteinID'})
    neg_pep_df = neg_pep_df.unique(subset=[
        'peptide', 'proteinID', 'sr1_Index', 'sr2_Index',
    ])

    neg_pep_df = merge_orf_level_data(
        neg_pep_df, config, 'spliced', 0, ID_COLUMNS + SPLICE_SPECIFIC_FEATURES
    )

    add_features_mp(neg_pep_df, config, 0, pep_len, False)


def add_features_mp(pep_df, config, label, pep_len, is_mm):
    pep_df = pep_df.with_row_count('group').with_columns(
        pl.col('group')//1_000
    )
    pos_pep_dfs = pep_df.partition_by('group')
    func_args = []
    for idx, pp_df in enumerate(pos_pep_dfs):
        func_args.append((pp_df, config.output_folder, pep_len, idx, label, is_mm))

    with mp.get_context('spawn').Pool(processes=20) as pool:
        pool.starmap(add_spliced_features, func_args)


def add_spliced_features(total_pep_df, output_folder, pep_len, idx, label, is_mm=False):
    """ Add all required training features to the DataFrame.
    """
    # Apply create_features function to compute all required features:
    total_pep_df = total_pep_df.with_columns(
        pl.struct([
            'peptide', 'sr1', 'protSeq', 'iupred3_preds',
            'sr1_Index', 'sr2_Index', 'canonicalPeptides'
        ]).map_elements(
            lambda x : create_spliced_features(
                x['peptide'], x['sr1'], x['protSeq'], x['iupred3_preds'],
                x['sr1_Index'], x['sr2_Index'], x['canonicalPeptides']
            )
        ).alias('allFeatures')
    )
    total_pep_df = total_pep_df.unnest('allFeatures')

    total_pep_df = total_pep_df.drop(['iupred3_preds', 'protSeq', 'rnaSeq', 'canonicalPeptides'])
    if is_mm:
        total_pep_df.write_parquet(
            f'{output_folder}/mmDatasets/{pep_len}/df_{label}_{idx}.parquet'
        )
    else:
        total_pep_df.write_parquet(
            f'{output_folder}/trainingDatasets/{pep_len}/df_{label}_{idx}.parquet'
        )


def create_spliced_features(peptide, sr1, prot_seq, iupred3_preds, sr1_index, sr2_index, canonical_peptides):
    """ Function to create features that may be relevant for model training.
    """
    sr2 = peptide[len(sr1):]
    sr1 = prot_seq[sr1_index:sr1_index+len(sr1)]
    sr2 = prot_seq[sr2_index:sr2_index+len(sr2)]
    il_prot_seq = prot_seq.replace('I', 'L')
    results = {}
    results['protLength'] = len(prot_seq)

    p1_index = (sr1_index + len(sr1)) - 1
    p2_index = (sr1_index + len(sr1)) - 2
    p_minus_1_index = sr1_index+len(sr1)
    p_minus_2_index = p_minus_1_index + 1

    p_prime_2_index = sr2_index+1
    p_prime_minus_1_index = sr2_index-1
    p_prime_minus_2_index = p_prime_minus_1_index - 1

    p1_aa = prot_seq[p1_index]
    if p_minus_1_index < len(prot_seq):
        p_neg_1_aa = prot_seq[p_minus_1_index]
    else:
        p_neg_1_aa = 'X'

    p1_prime_aa = prot_seq[sr2_index]
    if p_prime_minus_1_index > 0:
        p_neg_1_prime_aa = prot_seq[p_prime_minus_1_index]
    else:
        p_neg_1_prime_aa = 'X'
    
    for amino_acid in AMINO_ACIDS:
        results[f'p1_{amino_acid}'] = int(p1_aa == amino_acid)
        results[f'p_neg_1_{amino_acid}'] = int(p_neg_1_aa == amino_acid)
        results[f'p1_prime_{amino_acid}'] = int(p1_prime_aa == amino_acid)
        results[f'p_neg_1_prime_{amino_acid}'] = int(p_neg_1_prime_aa == amino_acid)

    if p2_index > 0:
        results['p2'] = prot_seq[p2_index]
    else:
        results['p2'] = 'X'

    if p_minus_1_index < len(prot_seq):
        results['p_minus_1'] = prot_seq[p_minus_1_index]
    else:
        results['p_minus_1'] = 'X'

    if p_minus_2_index < len(prot_seq):
        results['p_minus_2'] = prot_seq[p_minus_2_index]
    else:
        results['p_minus_2'] = 'X'


    if p_prime_2_index < len(prot_seq):
        results['p2_prime'] = prot_seq[p_prime_2_index]
    else:
        results['p2_prime'] = 'X'

    if p_prime_minus_1_index >= 0:
        results['p_minus_1_prime'] = prot_seq[p_prime_minus_1_index]
    else:
        results['p_minus_1_prime'] = 'X'

    if p_prime_minus_2_index >= 0:
        results['p_minus_2_prime'] = prot_seq[p_prime_minus_2_index]
    else:
        results['p_minus_2_prime'] = 'X'



    if not canonical_peptides:
        results['sr2_can_dist'] = None
        results['sr1_can_dist'] = None
    else:
        start_of_peps = []
        for peptide in canonical_peptides:
            start_of_peps.extend([
                occ.start() for occ in re.finditer(peptide, il_prot_seq) 
            ])
        results['sr1_can_dist'] = min([
            abs(sr1_index-pep_start) for pep_start in start_of_peps
        ])
        results['sr2_can_dist'] = min([
            abs(sr2_index-pep_start) for pep_start in start_of_peps
        ])

    for a_a in AMINO_ACIDS:
        results[f'{a_a}_p1'] = int(sr1[-1] == a_a)
        results[f'{a_a}_p1_prime'] = int(sr2[0] == a_a)


    results['sr1_localDisorder'] = np.mean(
        iupred3_preds[sr1_index:sr1_index+len(sr1)]
    ) 
    results['sr2_localDisorder'] = np.mean(
        iupred3_preds[sr2_index:sr2_index+len(sr1)]
    )

    return results
