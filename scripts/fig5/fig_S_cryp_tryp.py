import pandas as pd
import polars as pl
DARK_COLOURS = [
    'dodgerblue',
    'orange',
    'yellow',
    'purple',
    'darkgrey',
    'forestgreen',
    'navy',
    'pink',
    'cyan',
    '#908FB3',
]
CRYPTIC_STRATA = [
    # 'TrEMBL',
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
    # 'mutation',
    # 'fusion',
]

def get_cryptic_stratum_counts(nc_df):
    """ Function to get counts for each of the cryptic strata.
    """
    count_data = []
    total_unique = 0
    for stratum in CRYPTIC_STRATA:
        strat_df = nc_df.filter(
            pl.col(f'{stratum}_nProteins') > 0
        )

        for stratum2 in CRYPTIC_STRATA:
            if stratum == stratum2:
                continue
            strat_df = strat_df.filter(
                pl.col(f'{stratum2}_nProteins') == 0
            )
        unique_count = strat_df['peptide'].n_unique()
        count_data.append({
            'stratum': stratum,
            'uniqueCount': unique_count,
        })
        total_unique += unique_count

    count_data.append({
        'stratum': 'multi-mapped',
        'uniqueCount': nc_df.shape[0] - total_unique,
    })
    return pd.DataFrame(count_data)


TRYPTIC_RESULTS = [
    'projects/pisces_tryptic/K562-tryptic/outputFolder',
    'projects/pisces_tryptic/sarkizova_B721.221/outputFolder',
]
x = pd.concat([pd.concat([
    pd.read_csv(f'{project}/final/cryptic.csv') for project in TRYPTIC_RESULTS
]),pd.concat([
    pd.read_csv(f'{project}/final/multimapped.csv') for project in TRYPTIC_RESULTS
])

])
x = x[x['adjustedProbability'] > 0.85].drop_duplicates(subset=['peptide'])
y = x[(x['mutation_nProteins'] == 0) & (x['fusion_nProteins'] == 0) & (x['TrEMBL_nProteins'] == 0)]
cont_df = get_cryptic_stratum_counts(pl.from_pandas(y))
y[(y['CDS_frameshift_nProteins'] > 0) & (y['lncRNA_nProteins'] == 0) & (y['intronic_nProteins'] == 0) & (y['intergenic_nProteins'] == 0)].to_csv('temp.csv', index=False)
print(cont_df)
# print(y.shape[0], 'cryptic peptides with adjusted probability > 0.85, no mutations or fusions')
# det_df = pd.concat([
#     pd.read_csv(f'{project}/final/cryptic.csv') for project in TRYPTIC_RESULTS
# ]).drop_duplicates(subset=['peptide'])
# print(det_df.columns)

# x = pd.merge(x, det_df, how='inner', on='peptide')
# print(x.shape[0], 'cryptic peptides with adjusted probability > 0.85')
# print(x.columns)