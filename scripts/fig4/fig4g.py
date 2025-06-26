import polars as pl
import os
import pandas as pd
from pisces.plot_utils import create_comparison_logo_plot


pep_df = pl.read_parquet('data/piscesPeptides.parquet')
datasets = pep_df['datasets'].explode().unique(maintain_order=True).to_list()

x = pd.read_csv('scripts/fig4/per_allele_data.csv')
print(x['cisFreq'].quantile(0.75))
print(x['cisFreq'].quantile(0.25))
print(x['crypFreq'].quantile(0.75))
print(x['crypFreq'].quantile(0.25))
cis_alleles = x[x['cisFreq'] > x['cisFreq'].quantile(0.75)]['allele'].tolist()
cis_bg_alleles = x[x['cisFreq'] < x['cisFreq'].quantile(0.75)]['allele'].tolist()
cryp_alleles = x[x['crypFreq'] > x['crypFreq'].quantile(0.75)]['allele'].tolist()
cryp_bg_alleles = x[x['crypFreq'] < x['crypFreq'].quantile(0.75)]['allele'].tolist()
cluster_dfs = {
    'cis': [],
    'cisBg': [],
    'cryp': [],
    'crypBg': [],
}
can_counts = {}
cis_freqs = {}
cryp_freqs = {}
mm_freqs = {}
nc_freqs = {}
mean_lengths = {}
allele_projects = {}
unmapped_freqs = {}
for dataset in datasets:
    remapped_df = pep_df.filter(pl.col('datasets').list.contains(dataset))
    can_df = remapped_df.filter(
        pl.col('stratum').eq('canonical')
    ).select(['peptide']).unique(maintain_order=True).to_pandas()
    can_df = can_df[can_df['peptide'].apply(len) == 9]
    if dataset in cis_alleles:
        cluster_dfs['cis'].append(can_df)
    elif dataset in cis_bg_alleles:
        cluster_dfs['cisBg'].append(can_df)
    if dataset in cryp_alleles:
        cluster_dfs['cryp'].append(can_df)
    elif dataset in cryp_bg_alleles:
        cluster_dfs['crypBg'].append(can_df)
NAMES = {'cis': 'spliced', 'cryp': 'cryptic'}
for cluster in 'cis', 'cryp':
    # if not os.path.exists(f'scripts/dataset_level/{cluster}_comp'):
    #     os.mkdir(f'scripts/dataset_level/{cluster}_comp')
    create_comparison_logo_plot(
        [
            pd.concat(cluster_dfs[cluster]),
            pd.concat(cluster_dfs[f'{cluster}Bg']),
        ],
        [
            f'top quartile',
            'remainder',
        ],
        9,
        f'manuscript_figs/fig4',
        file_name=f'fig4g_{cluster}.svg',
    )