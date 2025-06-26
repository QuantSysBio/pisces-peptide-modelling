""" Function to show distribution of frequency across 
"""
from statistics import variance

import polars as pl
import plotly.graph_objects as go
from pisces.plot_utils import clean_plotly_fig, split_strata, split_nc_strata
import plotly.io as pio
from scipy.stats import wilcoxon, mannwhitneyu

pep_df = pl.read_parquet('data/piscesPeptides.parquet')

datasets = pep_df['datasets'].explode().unique().to_list()
cis_freqs = {}
cryp_freqs = {}
mean_lengths = {}
dataset_projects = {}
for dataset in datasets:
    remapped_df = pep_df.filter(pl.col('datasets').list.contains(dataset))
    can_df = remapped_df.filter(
        pl.col('stratum').eq('canonical') & pl.col('piscesDiscoverable')
    )
    contam_df = remapped_df.filter(
        pl.col('stratum').eq('contaminant')
    )
    spliced_df = remapped_df.filter(
        pl.col('stratum').eq('spliced')
    )
    cryptic_df = remapped_df.filter(
        pl.col('stratum').eq('cryptic')
    )
    mm_df = remapped_df.filter(
        pl.col('stratum').eq('multi-mapped')
    )

    mean_lengths[dataset] = can_df['peptide'].map_elements(len).mean()
    cis_freqs[dataset] = 100*spliced_df['peptide'].n_unique()/(
        can_df['peptide'].n_unique()+spliced_df['peptide'].n_unique()
    )
    cryp_freqs[dataset] = 100*cryptic_df['peptide'].n_unique()/(
        can_df['peptide'].n_unique()+cryptic_df['peptide'].n_unique()
    )

fig = go.Figure()
fig.add_trace(
    go.Box(
        x=['spliced']*len(cis_freqs),
        y=list(cis_freqs.values()),
        boxpoints=False,
        line_color='black',
        line_width=1,
        fillcolor='#9BBFE5',
        opacity=0.8,
    )
)
print(wilcoxon(list(cis_freqs.values()), list(cryp_freqs.values())))
print(mannwhitneyu(list(cis_freqs.values()), list(cryp_freqs.values())))
print(variance(list(cis_freqs.values())), variance(list(cryp_freqs.values())))
fig.add_trace(
    go.Box(
        x=['cryptic']*len(cryp_freqs),
        y=list(cryp_freqs.values()),
        boxpoints=False,
        line_color='black',
        line_width=1,
        fillcolor='#BA69BE',
        opacity=0.8,
    )
)
fig = clean_plotly_fig(fig)
fig = fig.update_layout(height=500, width=200)
fig.update_yaxes(range=[0,12])
fig.show()
pio.write_image(fig, 'manuscript_figs/fig4/fig4e.svg')