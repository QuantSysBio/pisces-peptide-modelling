
import os

import numpy as np
from numpy.random import choice
import polars as pl
import pandas as pd

from pisces.proteome_utils import distribute_mapping
from pisces.remap import (
    process_fasta_folder, extract_spliced_details, extract_details
)

np.random.seed(42)

N_RANDOM_PEPTIDES = 1_000_000

AMINO_ACIDS = 'ACDEFGHKLMNPQRSTVWY' # Isoleucine ignored for PISCES applications.

def create_bg(config, pep_length):
    get_can_distro(config, pep_length)
    generate_random(config, pep_length)
    remap_random(config, pep_length)

def get_can_distro(config, pep_length):

    if not os.path.exists(f'{config.background_folder}/frequency_dfs/{pep_length}'):
        os.mkdir(f'{config.background_folder}/frequency_dfs/{pep_length}')

    can_df_counts = []
    pep_df = pl.read_parquet('/data/John/ALL_FINAL/shared_dfs/casanovo_241124/peptides.parquet')
    datasets = pep_df['datasets'].explode().unique(maintain_order=True).sort().to_list()
    if config.cell_line == 'K562':
        datasets = [x for x in datasets if x.startswith('K562')]
    else:
        datasets = [x for x in datasets if not x.startswith('K562')]
    print(datasets)
    for dataset in datasets:
        remapped_df = pep_df.filter(pl.col('datasets').list.contains(dataset))
        can_df = remapped_df.filter(
            pl.col('stratum').eq('canonical') & pl.col('piscesDiscoverable') &
            pl.col('peptide').str.len_chars().eq(pep_length) &
            pl.col('cellLines').list.contains(config.cell_line)
        )
        if not can_df.shape[0]:
            continue

        can_df_counts.append(
            {'dataset': dataset, 'count': can_df['peptide'].n_unique()}
        ) 
        counts = np.zeros(shape=(len(AMINO_ACIDS), pep_length))
        for peptide in can_df['peptide'].to_list():
            for pos_idx, a_a in enumerate(peptide):
                aa_idx = AMINO_ACIDS.index(a_a)
                counts[aa_idx, pos_idx] += 1
        counts /= can_df.shape[0]

        df = pd.DataFrame(
            counts, index=list(AMINO_ACIDS), columns=list(range(1,pep_length+1))
        )
        df.to_csv(
            f'{config.background_folder}/frequency_dfs/{pep_length}/{dataset}.csv'
        )

    count_df = pd.DataFrame(can_df_counts)
    count_df['fraction'] = count_df['count']/count_df['count'].max()
    count_df.to_csv(
        f'{config.background_folder}/sample_ratios/ratio_{pep_length}.csv',
        index=False,
    )


def generate_random(config, pep_length):
    if not os.path.exists(f'{config.background_folder}/random_dfs/{pep_length}'):
        os.mkdir(f'{config.background_folder}/random_dfs/{pep_length}')

    for file_name in os.listdir(f'{config.background_folder}/frequency_dfs/{pep_length}'):
        freq_df = pd.read_csv(
            f'{config.background_folder}/frequency_dfs/{pep_length}/{file_name}'
        )
        groups = [
            choice(
                list(AMINO_ACIDS), N_RANDOM_PEPTIDES, p=freq_df[str(pos_idx)]
            ) for pos_idx in range(1,pep_length+1)
        ]
        peptides = []
        for row_idx in range(N_RANDOM_PEPTIDES):
            peptides.append(''.join(
                groups[pos_idx][row_idx] for pos_idx in range(pep_length)
            ))

        random_df = pd.DataFrame({'peptide': peptides})
        random_df.to_csv(
            f'{config.background_folder}/random_dfs/{pep_length}/{file_name}', index=False,
        )

CANONICAL_PROTEOME = 'background_analsis/proteome_CDS_main_ORF.fasta'
EXPANDED_STRATA_FOLDER = 'dataFolder/allStrata'
N_CORES = 70

def remap_random(config, pep_length):
    """ Function to remap peptides to the spliced proteome.
    """
    if not os.path.exists(f'{config.background_folder}/remapped/{pep_length}'):
        os.mkdir(f'{config.background_folder}/remapped/{pep_length}')

    for file_name in os.listdir(
        f'{config.background_folder}/random_dfs/{pep_length}'
        ):
        dataset = file_name.split('.')[0]
        if os.path.exists(f'{config.background_folder}/remapped/{pep_length}/{dataset}/peptides.csv'):
            print(f'Skipping {dataset}...')
            continue
        else:
            print(f'Running {dataset}...')

        unique_pep_df = pl.read_csv(
            f'{config.background_folder}/random_dfs/{pep_length}/{file_name}'
        )
        unique_pep_df = unique_pep_df.unique()
        output_folder = f'{config.background_folder}/remapped/{pep_length}/{dataset}'
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        if not os.path.exists(f'{output_folder}/details'):
            os.mkdir(f'{output_folder}/details')

        # Map to canonical and spliced strata.
        unique_pep_df = distribute_mapping(
            unique_pep_df,
            config.proteome,
            'canonical',
            config.n_cores,
            with_splicing=True,
            max_intervening=None,
        )

        # Extract details for canonical and spliced
        unique_pep_df = extract_details(unique_pep_df, [], output_folder, 'canonical')
        unique_pep_df = extract_spliced_details(unique_pep_df, output_folder)

        print(unique_pep_df)

        # Map to cryptic strata and extract details
        unique_pep_df, _ = process_fasta_folder(
            unique_pep_df,
            config.cryptic_folder,
            config.n_cores,
            output_folder,
            'cryptic',
        )
        print(unique_pep_df)
        print(f'writing to {output_folder}/peptides.csv')
        unique_pep_df.write_csv(f'{output_folder}/peptides.csv')
