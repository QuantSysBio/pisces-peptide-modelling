""" Functions for combining PISCES and t/d dataframes from multiple projects
    into a single dataframe.
"""
import polars as pl
import pandas as pd
from pisces.plot_utils import split_strata, split_nc_strata

MAIN_NAMES = [
    'canonical (t/d)', 'contam (t/d)', 'canonical (pisces)',
    'contaminant', 'non-canonical'
]
NC_NAMES = ['spliced', 'multi-mapped', 'cryptic', 'unmapped',]
ALL_NAMES = [
    'canonical (t/d)', 'contam (t/d)', 'canonical (pisces)', 'contaminant',
    'spliced', 'multi-mapped', 'cryptic', 'unmapped',
]
COLUMNS = [
    'canonical_nProteins', 'canonical_Proteins',
    'nSplicedProteins', 'sr1', 'interveningSeqLengths',
    'splicedProteins', 'sr1_Index', 'sr2_Index', 'isForward',
    'nCrypticProteins', 'fiveUTR_Proteins', 'TrEMBL_Proteins',
    'intergenic_Proteins', 'CDS_frameshift_Proteins', 'threeUTR_Proteins',
    'fusion_Proteins', 'mutation_Proteins',
    'lncRNA_Proteins', 'intronic_Proteins', 'fiveUTR_nProteins',
    'TrEMBL_nProteins', 'intergenic_nProteins', 'CDS_frameshift_nProteins',
    'threeUTR_nProteins', 'lncRNA_nProteins', 'intronic_nProteins',
    'fusion_nProteins', 'mutation_nProteins',
]
COUNT_COLUMNS = {
    'canonical_nProteins': 'canonical_Proteins',
    'nSplicedProteins': 'splicedProteins',
    'fiveUTR_nProteins': 'fiveUTR_Proteins',
    'threeUTR_nProteins': 'threeUTR_Proteins',
    'intergenic_nProteins': 'intergenic_Proteins',
    'intronic_nProteins': 'intronic_Proteins',
    'lncRNA_nProteins': 'lncRNA_Proteins',
    'fusion_nProteins': 'fusion_Proteins',
    'mutation_nProteins': 'mutation_Proteins',
    'CDS_frameshift_nProteins': 'CDS_frameshift_Proteins',
    'TrEMBL_nProteins': 'TrEMBL_Proteins',
}
INT_COLS = [
    'interveningSeqLengths', 'sr1_Index', 'sr2_Index', 'isForward',
]
SPLICED_FEATS = [
    'sr1','interveningSeqLengths','splicedProteins','sr1_Index','sr2_Index','isForward'
]

def create_pisces_db(config):
    """ Function to create a PISCES database from multiple projects.
    """
    all_dfs = {
        'contaminant': [], 'canonical (pisces)': [],
        'spliced': [], 'multi-mapped': [], 'cryptic': [], 'unmapped': [],
        'canonical (t/d)': [], 'contam (t/d)': [],
    }
    meta_df = pd.read_csv(config.meta_df)
    for _, df_row in meta_df.iterrows():
        project_path = df_row['projectPath']
        dataset = df_row['dataset']
        allele = df_row['allele']
        cell_line = df_row['cellLine']

        remapped_df, td_df = get_identified_peptides(project_path)

        can_td_df = td_df[~td_df['mapsToContaminant']]
        contam_td_df = td_df[td_df['mapsToContaminant']]
        can_df, contam_df, nc_df = split_strata(remapped_df)
        spliced_df, mm_df, cryptic_df, unmapped_df = split_nc_strata(nc_df)

        for strat_df, name in zip(
            [
                can_td_df, contam_td_df, can_df, contam_df,
                spliced_df, mm_df, cryptic_df, unmapped_df
            ],
            ALL_NAMES,
        ):
            id_df, details_df = get_details_df(
                strat_df, name, project_path
            )
            id_df = combine_dfs(id_df, details_df, allele, cell_line, dataset, name)
            all_dfs[name].append(id_df)

    combined_df = combine_all_strata(all_dfs)

    unique_df = combined_df.drop_duplicates(subset=['peptide'], keep=False)
    duplicated_df = combined_df[combined_df.duplicated(subset=['peptide'], keep=False)]

    duplicated_df = duplicated_df.groupby('peptide', as_index=False).agg({
        'stratum': lambda x: 'contaminant' if 'contaminant' in x.values else 'error',
        'piscesDiscoverable': 'max',
        'tdDiscoverable': 'max',
        'datasets': _combine_lists,
        'cellLines': _combine_lists,
        'alleles': _combine_lists,
        'maxProbability': 'max',
    })
    duplicated_df['nDatasets'] = duplicated_df['datasets'].apply(len)
    for column in COLUMNS:
        if 'nProteins' in column or column in ['nCrypticProteins', 'nSplicedProteins']:
            duplicated_df[column] = 0
        else:
            duplicated_df[column] = ''

    combined_df = pd.concat([unique_df, duplicated_df[unique_df.columns]])

    combined_df = pl.from_pandas(combined_df)
    for column in COLUMNS:
        if not 'nProteins' in column and column not in ['nCrypticProteins', 'nSplicedProteins']:
            if column in INT_COLS:
                combined_df = combined_df.with_columns(
                    pl.col(column).map_elements(
                        lambda x : [] if (not x or x is None) else [int(a) for a in x.split(' ')],
                        return_dtype=pl.List(pl.Int64)
                    )
                )
                combined_df = combined_df.with_columns(
                    pl.col(column).fill_null([])
                )
            else:
                combined_df = combined_df.with_columns(
                    pl.col(column).map_elements(
                        lambda x : [] if (not x or x is None) else [y for y in x.split(' ') if y and y is not None],
                        return_dtype=pl.List(pl.String)
                    )
                )
                combined_df = combined_df.with_columns(
                    pl.col(column).fill_null([])
                )

    for count_column, list_columns in COUNT_COLUMNS.items():
        combined_df = combined_df.with_columns(
            pl.col(list_columns).list.len().alias(count_column)
        )

    combined_df = combined_df.sort(by=['stratum', 'maxProbability'], descending=[False, True])
    combined_df.write_parquet(config.peptides_pq)

def combine_all_strata(all_dfs):
    """ Function to combine all strata into a single dataframe.
    """
    can_dfs = []
    nc_dfs = []
    contam_dfs = []
    agg_dict = {
            'adjustedProbability': 'max',
            'datasets': lambda x : list(set(x)),
            'alleles': lambda x : list(set(x)),
            'cellLines': lambda x : list(set(x)),
        }
    agg_dict.update({col: 'first' for col in COLUMNS})

    for (stratum_name, stratum_dfs) in all_dfs.items():
        stratum_df = combine_stratum_dfs(stratum_name, stratum_dfs, agg_dict)
        if 'canonical' in stratum_name:
            can_dfs.append(stratum_df)
        elif 'contam' in stratum_name:
            contam_dfs.append(stratum_df)
        else:
            nc_dfs.append(stratum_df)

    combo_df = pd.concat(nc_dfs)
    combo_df = combo_df.drop_duplicates(subset=['peptide'])

    contam_df = combine_td_pisces_dfs(contam_dfs)
    can_df = combine_td_pisces_dfs(can_dfs)
    total_df = pd.concat([contam_df, can_df, combo_df])
    return total_df

def get_identified_peptides(project_path):
    """ Function to get identified peptides via PISCES and t/d.
    """
    remapped_df = pd.read_csv(f'{project_path}/filtered_mapped.csv')
    remapped_df = remapped_df[
        (remapped_df['adjustedProbability'] > 0.85) #&
        # (remapped_df['qValue_PSM'] < 0.05) &
        # (remapped_df['qValue_peptide'] < 0.05)
    ]
    remapped_df = remapped_df.sort_values(by='adjustedProbability', ascending=False)
    remapped_df = remapped_df.drop_duplicates(subset=['peptide'])

    td_df = pd.read_csv(f'{project_path}/canonicalOutput/finalPsmAssignments.csv')
    td_df = td_df[(td_df['qValue'] < 0.01) & (td_df['postErrProb'] < 0.1)]

    return remapped_df, td_df


def get_details_df(strat_df, name, project_path):
    """ Function to get details for each stratum.
    """
    if name in ['cryptic', 'spliced']:
        details_df = pd.read_csv(f'{project_path}/details/{name}.csv')
        details_df = details_df.drop_duplicates(subset=['peptide'])
    elif name == 'multi-mapped':
        details_spliced_df = pd.read_csv(f'{project_path}/details/spliced.csv')
        details_cryptic_df = pd.read_csv(f'{project_path}/details/cryptic.csv')
        details_spliced_df = details_spliced_df.drop_duplicates(subset=['peptide'])
        details_cryptic_df = details_cryptic_df.drop_duplicates(subset=['peptide'])
        details_df = pd.merge(
            details_spliced_df, details_cryptic_df, on=['peptide'], how='outer',
        )
    elif name == 'canonical (pisces)':
        details_df = pd.read_csv(f'{project_path}/details/canonical.csv')
        details_df = details_df.drop_duplicates(subset=['peptide'])
    elif name in ['canonical (t/d)', 'contam (t/d)']:
        code = name.split(' ')[0]
        strat_df = strat_df.rename(columns={'proteins': f'{code}_Proteins'})
        strat_df[f'{code}_nProteins'] = strat_df[f'{code}_Proteins'].apply(
            lambda x : len(x.split(' ')) if x else 0
        )
        details_df = None
    elif name in ['contaminant', 'unmapped']:
        details_df = None
    else:
        raise ValueError(f'Unrecognised name {name}')

    return strat_df, details_df


def combine_dfs(id_df, details_df, allele, cell_line, dataset, name):
    """ Function to combine identification DataFrame with mapping DataFrame.
    """
    id_df['alleles'] = allele
    id_df['cellLines'] = cell_line
    id_df['datasets'] = dataset
    if details_df is not None:
        id_df = pd.merge(id_df, details_df)

    for column in COLUMNS:
        if column not in id_df.columns:
            if 'nProteins' in column or column in ['nCrypticProteins', 'nSplicedProteins']:
                id_df[column] = 0
            else:
                id_df[column] = ''

    if 't/d' in name:
        id_df['adjustedProbability'] = 1 - id_df['postErrProb']

    return id_df[
        ['peptide', 'adjustedProbability', 'datasets', 'alleles', 'cellLines',] +
        COLUMNS
    ].drop_duplicates(subset=[
        'peptide', 'adjustedProbability', 'datasets', 'alleles', 'cellLines',
    ])


def combine_stratum_dfs(stratum_name, stratum_dfs, agg_dict):
    """ Function to combine dataframes for a given stratum.
    """
    combo_df = pd.concat(stratum_dfs)
    combo_df = combo_df.groupby(['peptide'], as_index=False).agg(agg_dict)
    combo_df = combo_df.rename(
        columns={
            'adjustedProbability': 'maxProbability',
        }
    )
    combo_df['nDatasets'] = combo_df['datasets'].apply(len)
    if 'canonical' in stratum_name:
        combo_df['stratum'] = 'canonical'
    elif stratum_name == 'contam (t/d)':
        combo_df['stratum'] = 'contaminant'
    else:
        combo_df['stratum'] = stratum_name

    if 't/d' in stratum_name:
        combo_df['piscesDiscoverable'] = False
        combo_df['tdDiscoverable'] = True
    else:
        combo_df['piscesDiscoverable'] = True
        combo_df['tdDiscoverable'] = False

    combo_df = combo_df[
        [
            'peptide', 'stratum', 'piscesDiscoverable', 'tdDiscoverable', 'nDatasets',
            'datasets', 'cellLines', 'alleles', 'maxProbability'
        ] + COLUMNS
    ]
    return combo_df.sort_values(by=['maxProbability'], ascending=False)

def combine_td_pisces_dfs(dfs):
    """ Function to combine t/d and PISCES dataframes.
    """
    agg_dict = {
        'stratum': 'first',
        'piscesDiscoverable': 'max',
        'tdDiscoverable': 'max',
        'datasets': _combine_lists,
        'cellLines': _combine_lists,
        'alleles': _combine_lists,
        'maxProbability': 'max',
    }
    agg_dict.update({col: 'first' for col in COLUMNS})

    combo_df = pd.concat(dfs)
    combo_df = combo_df.groupby('peptide', as_index=False).agg(agg_dict)
    combo_df['nDatasets'] = combo_df['datasets'].apply(len)
    combo_df = combo_df.sort_values(by=['maxProbability'], ascending=False).reset_index(drop=True)

    return combo_df

def _combine_lists(list_of_lists):
    """ Helper function to combine lists keeping uniqueness."""
    total_list = []
    for entry in list_of_lists:
        total_list.extend(entry)
    return list(set(total_list))
