""" Functions for pre-processing canonical or cryptic peptides.
"""
import re
import warnings

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import polars as pl


from ppm.constants import (
    AMINO_ACID_GROUPS,
    CRYPTIC_STRATA,
    ID_COLUMNS,
    START_CODONS,
    STOP_CODONS,
)
from ppm.preprocess_utils import merge_orf_level_data, get_sampled_negative_peps
from ppm.kozak_scoring import kozak_similarity_score
warnings.filterwarnings('ignore')


def preprocess_canonical(config):
    """ Function to preprocess canonical peptides.
    """
    for pep_len in range(9, 13):
        print(f'\t\t\tprocessing length {pep_len}...')
        process_stratum('canonical', False, config, pep_len)


def preprocess_cryptic(config):
    """ Function run to process uniquely mapped peptides attributed to each stratum.
    """
    print(f'\tprocessing cell line {config.cell_line}')
    for stratum in CRYPTIC_STRATA:
        print(f'\t\tprocessing stratum {stratum}...')
        for pep_len in range(9, 13):
            print(f'\t\t\tprocessing length {pep_len}...')
            process_stratum(stratum, True, config, pep_len)


def process_stratum(stratum, is_cryptic, config, pep_len):
    """ Function to process data

    stratum : str
        The stratum being processed
    is_cryptic : bool
        Is the stratum cryptic
    cell_line : str
        The cell line being processed
    """
    # Gather and format positives samples - pisces identified peptides
    # and negative samples - random background peptides with matched amino
    # acid distribution.
    pos_pep_df = gather_positive_samples(stratum, is_cryptic, config, pep_len, False)
    neg_pep_df = gather_negative_samples(stratum, is_cryptic, config, pep_len)

    if pos_pep_df is None and neg_pep_df is None:
        return

    # Concatenate so positive samples before negative, then drop duplicates
    # so no peptides are positive and negative.
    total_pep_df = pl.concat([df for df in (pos_pep_df, neg_pep_df) if df is not None])
    total_pep_df = total_pep_df.unique(
        subset=['peptide', 'proteinID'],
        maintain_order=True,
    )

    # Additional features for cryptic and canonical peptides.
    total_pep_df = add_features(total_pep_df)

    if is_cryptic:
        total_pep_df = total_pep_df.with_columns(
            pl.lit(CRYPTIC_STRATA.index(stratum)).alias('stratum')
        )

    total_pep_df.write_csv(f'{config.output_folder}/trainingDatasets/{stratum}_{pep_len}.csv')


def gather_positive_samples(stratum, is_cryptic, config, pep_len, is_mm):
    """ Function get positive peptides.
    """
    # Read in data, get 9mer, K562 data only:
    pos_pep_df = pl.read_parquet(
        config.peptides_pq,
        columns=[
            'peptide', 'stratum', 'cellLines', 'piscesDiscoverable',
            f'{stratum}_nProteins', f'{stratum}_Proteins',
            'fusion_nProteins', 'mutation_nProteins', 'TrEMBL_nProteins'
        ],
    )
    pos_pep_df = pos_pep_df.filter(
        pl.col('cellLines').list.contains(config.cell_line) &
        pl.col('peptide').str.len_chars().eq(pep_len) &
        pl.col('piscesDiscoverable').eq(1)
    )

    # Filter to unique in that stratum:
    if is_mm:
        pos_pep_df = pos_pep_df.filter(
            pl.col('stratum').eq('multi-mapped') & pl.col(f'{stratum}_nProteins').gt(0)
        )
    elif is_cryptic:
        pos_pep_df = pos_pep_df.filter(
            pl.col('stratum').eq('cryptic') & pl.col(f'{stratum}_nProteins').gt(0) &
            pl.col('fusion_nProteins').eq(0) & pl.col('mutation_nProteins').eq(0) &
            pl.col('TrEMBL_nProteins').eq(0)
        )
    else:
        pos_pep_df = pos_pep_df.filter(
            pl.col('stratum').eq(stratum)
        )

    pos_pep_df = pos_pep_df.drop(['fusion_nProteins', 'mutation_nProteins', 'TrEMBL_nProteins'])
    # Explode columns with different possible accessions listed:
    pos_pep_df = pos_pep_df.rename({f'{stratum}_Proteins': 'proteinID'}).explode('proteinID')

    return merge_orf_level_data(pos_pep_df, config, stratum, 1, ID_COLUMNS)


def gather_negative_samples(stratum, is_cryptic, config, pep_len):
    """ Function to get random background peptides for a given stratum and cell line.
    """
    # Collect negative samples from relevant cell line.
    pep_dfs = get_sampled_negative_peps(config, pep_len, is_cryptic, stratum)
    if not pep_dfs:
        return None
    for pep_df in pep_dfs:
        print(pep_df.columns)
    neg_pep_df = pl.concat(pep_dfs)
    neg_pep_df = neg_pep_df.with_columns(
        pl.col(f'{stratum}_Proteins').map_elements(
            lambda x : [a for a in x.split(' ') if a]
        )
    ).explode(f'{stratum}_Proteins')
    neg_pep_df = neg_pep_df.rename({f'{stratum}_Proteins': 'proteinID'})

    return merge_orf_level_data(neg_pep_df, config, stratum, 0, ID_COLUMNS)


def create_features(peptide, prot_seq, rna_seq, iupred3_preds):
    """ Function to create features that may be relevant for model training.
    """
    pep_analysis = ProteinAnalysis(peptide)
    il_prot_seq = prot_seq.replace('I', 'L')
    pep_len = len(peptide)
    results = {}

    results['peptide_hydrophobicity'] = pep_analysis.gravy()
    results['protLength'] = len(prot_seq)

    for aa_group, aa_list in AMINO_ACID_GROUPS.items():
        results[f'C_term_{aa_group}'] = int(peptide[-1] in aa_list)
    for nucleotide in 'AUGC':
        results[f'{nucleotide}_frac'] = rna_seq.count(nucleotide)/len(rna_seq)

    results['position'] = [
        occ.start() for occ in re.finditer(peptide, il_prot_seq)
    ]
    results['il_peptide'] = [
        prot_seq[pos:pos+pep_len] for pos in results['position']
    ]

    for start_codon in START_CODONS:
        results[f'{start_codon}_upstream'] = [0 for _ in results['position']]
    results['start_dist'] = [None for _ in results['position']]
    results['endCodon'] = [
        rna_seq[(position+pep_len-1)*3:(position+pep_len)*3] for position in results['position']
    ]
    results['postCodon'] = [
        rna_seq[(position+pep_len)*3:(position+pep_len+1)*3] for position in results['position']
    ]
    results['kozakScore'] = None
    for pos_idx, position in enumerate(results['position']):
        for idx in range(position, -1, -1):
            codon = rna_seq[idx*3:(idx+1)*3]
            if codon in STOP_CODONS:
                break
            for start_codon in START_CODONS:
                if codon == start_codon:
                    results[f'{start_codon}_upstream'][pos_idx] = 1
                    distance = position - idx
                    if (
                        results['start_dist'][pos_idx] is None or
                        results['start_dist'][pos_idx] > distance
                    ) and codon == 'AUG':
                        results['start_dist'][pos_idx] = position - idx
                        if idx*3 > 10 and (idx*3) + 13 < len(rna_seq):
                            results['kozakScore'] = kozak_similarity_score(
                                rna_seq[(idx*3)-10:(idx*3)+13]
                            )

    results['relativePosition'] = [
        pos/results['protLength'] for pos in results['position']
    ]
    results['localDisorder'] = [
        np.mean(
            iupred3_preds[pos:pos+len(peptide)]
        ) for pos in results['position']
    ]
    # Mostly relevant for frameshift, look for stop codons
    results['stopDistances'] = [
        il_prot_seq[start_pos+pep_len:].index('*') if '*' in il_prot_seq[start_pos+pep_len:]
        else len(il_prot_seq[start_pos+pep_len:])
        for start_pos in results['position']
    ]
    results['fragmentLength'] = []
    for aa_group in AMINO_ACID_GROUPS:
        results[f'C_term_neg_1_{aa_group}'] = []
    results['C_term_neg_1_end'] = []
    results['C_term_downstream'] = []
    results['C_term_downstream_2'] = []
    results['N_term_upstream'] = []
    results['N_term_upstream_2'] = []
    prot_splits = prot_seq.split('*')
    for prot_frag in prot_splits:
        padded_frag = prot_frag +'XX'
        if peptide in prot_frag.replace('I', 'L'):
            for pos in [
                occ.start() for occ in re.finditer(peptide, prot_frag.replace('I', 'L'))
            ]:
                results['fragmentLength'].append(len(prot_frag))
                results['C_term_downstream'].append(padded_frag[pos+pep_len])
                results['C_term_downstream_2'].append(padded_frag[pos+pep_len+1])
                results['N_term_upstream'].append(('X' + padded_frag)[pos])
                results['N_term_upstream_2'].append(('XX' + padded_frag)[pos])
                for aa_group, aa_list in AMINO_ACID_GROUPS.items():
                    results[f'C_term_neg_1_{aa_group}'].append(
                        int(padded_frag[pos+pep_len] in aa_list)
                    )
                results['C_term_neg_1_end'].append(int(padded_frag[pos+pep_len] == 'X'))

    return results


def add_features(total_pep_df):
    """ Add all required training features to the DataFrame.
    """
    # Apply create_features function to compute all required features:
    total_pep_df = total_pep_df.with_columns(
        pl.struct(['peptide', 'protSeq', 'rnaSeq', 'iupred3_preds']).map_elements(
            lambda x : create_features(x['peptide'], x['protSeq'], x['rnaSeq'], x['iupred3_preds'])
        ).alias('allFeatures')
    )
    total_pep_df = total_pep_df.unnest('allFeatures')

    # Explode in case a peptide mulit-maps within a single ORF.
    total_pep_df = total_pep_df.explode(
        [
            'position', 'relativePosition', 'stopDistances', 'fragmentLength', 'localDisorder',
            'C_term_neg_1_end', 'C_term_downstream', 'N_term_upstream',
            'C_term_downstream_2', 'N_term_upstream_2', 'il_peptide',
            'start_dist', 'endCodon', 'postCodon',
        ] + [
            f'{start_codon}_upstream' for start_codon in START_CODONS
        ] + [
            f'C_term_neg_1_{aa_group}' for aa_group in AMINO_ACID_GROUPS
        ]
    )

    return total_pep_df.drop(['iupred3_preds', 'protSeq', 'rnaSeq'])
