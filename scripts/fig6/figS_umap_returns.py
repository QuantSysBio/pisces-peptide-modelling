
import sys

import numpy as np
import pandas as pd
import polars as pl
# import plotly.express as px
from plotly.subplots import make_subplots
from pisces.constants import AMINO_ACIDS
import plotly.graph_objects as go
from pisces.plot_utils import clean_plotly_fig, get_count_df, split_strata, split_nc_strata
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import umap.umap_ as umap
import plotly.io as pio
import plotly.express as px
from ppm.constants import CRYPTIC_STRATA

SPLICED_FOLDERS = ['output/splicedK562_241129', 'output/spliced721_241206']
CRYPTIC_FOLDERS = ['output/crypticK562_241129', 'output/cryptic721_241205']


pep_df = pl.read_parquet('/data/John/ALL_FINAL/shared_dfs/casanovo_241124/peptides.parquet')
datasets = pep_df['datasets'].explode().unique(maintain_order=True).sort().to_list()
ys = []
for s_f, c_f in zip(SPLICED_FOLDERS, CRYPTIC_FOLDERS):
    spliced_df = pd.read_csv(
        f'{s_f}/mm_scored.csv', usecols=['peptide', 'meanScore']
    )
    spliced_df = spliced_df.rename(columns={'meanScore': 'splicedScore'})

    cryptic_df = pd.read_csv(
        f'{c_f}/mm_scored.csv',
        usecols=['peptide', 'meanScore'],
    )
    cryptic_df = cryptic_df.rename(columns={'meanScore': 'crypticScore'})

    mm_df = pd.merge(spliced_df, cryptic_df, how='inner', on=['peptide'])
    print(cryptic_df.shape)
    print(mm_df.shape)
    ys.append(mm_df)

mm_df_scored = pd.concat(ys)
mm_df_scored = mm_df_scored.groupby('peptide', as_index=False).agg({'splicedScore': 'max', 'crypticScore': 'max'})
mm_df_scored['totalScore'] = mm_df_scored['splicedScore'] + mm_df_scored['crypticScore']
print(mm_df_scored)
mm_df_scored['splicedScore'] /= mm_df_scored['totalScore']
mm_df_scored['crypticScore'] /= mm_df_scored['totalScore']
mm_df_scored = pl.from_pandas(mm_df_scored)
print(mm_df_scored)

can_counts = {}
cis_freqs = {}
cryp_freqs = {}

mean_lengths = {}
allele_projects = {}

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
    mm_df_raw = remapped_df.filter(
        pl.col('stratum').eq('multi-mapped')
    )
    print('raw', mm_df_raw.shape)
    mm_df = mm_df_scored.join(mm_df_raw, how='inner', on='peptide')
    print('merged', mm_df.shape)

    mean_lengths[dataset] = can_df['peptide'].apply(len).mean()
    cis_freqs[dataset] = 100*(spliced_df['peptide'].n_unique()+mm_df['splicedScore'].sum())/(
        can_df['peptide'].n_unique()+spliced_df['peptide'].n_unique()
    )
    cryp_freqs[dataset] = 100*(cryptic_df['peptide'].n_unique()+mm_df['crypticScore'].sum())/(
        can_df['peptide'].n_unique()+cryptic_df['peptide'].n_unique()
    )

print(cryp_freqs)
print(cis_freqs)
alleles = []
arrays = []
for dataset in datasets:
    alleles.append(dataset)
    remapped_df = pep_df.filter(pl.col('datasets').list.contains(dataset))
    can_df = remapped_df.filter(
        pl.col('stratum').eq('canonical')
    )

    can_df = can_df.filter(pl.col('peptide').str.len_chars() == 9)
    can_df = can_df.select(['peptide']).unique(maintain_order=True).to_pandas()
    can_count2 = can_df['peptide'].nunique()
    count_df = get_count_df(can_df, 9, AMINO_ACIDS)
    # count_df = count_df[(count_df.index == 2) | (count_df.index == 3) | (count_df.index == 9)]
    norm_df = count_df.div(count_df.sum(axis=1), axis=0)
    norm_array = norm_df.values.flatten()
    # norm_array = np.append(norm_array, mean_lengths[dataset])

    # norm_array = np.append(norm_array, can_counts[allele]/can_count2)


    arrays.append(norm_array)



total_df = pd.DataFrame(arrays, index=alleles)


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
fig.add_trace(
    go.Box(
        x=['cryptic']*len(cis_freqs),
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
fig.update_yaxes(range=[0,20])
fig.show()
pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_A.svg')


pca = umap.UMAP(random_state=42)
# pca = PCA(n_components=2, random_state=42)
pca_df = pca.fit_transform(total_df)
# print(pca.explained_variance_)
print(pca_df)
print(alleles)
kmeans_df = pd.DataFrame(pca_df, columns=['pca1', 'pca2',])
kmeans = KMeans(n_clusters=6, random_state=1)
kmeans_df['kmeans_res'] = kmeans.fit_predict(kmeans_df)
kmeans_df['allele'] = alleles
kmeans_df['cisFreq'] = cis_freqs.values()
kmeans_df['crypFreq'] = cryp_freqs.values()
kmeans_df = kmeans_df.sort_values('crypFreq').reset_index(drop=True)
kmeans_df.to_csv('fig6/umap_returns/per_allele_data.csv', index=False)
print(kmeans_df)
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=kmeans_df['pca1'],
        y=kmeans_df['pca2'],
        # text=alleles,
        customdata=kmeans_df[['allele', 'crypFreq']], hovertemplate='%{customdata[0]}<br>'+'%{customdata[1]}',
        marker_color=kmeans_df['crypFreq'],
        marker_opacity=0.85,
        mode='markers+text',marker_colorscale='Magenta',textposition='bottom center', marker_cmin=2, marker_cmax=12,
        marker_size=9,
        marker_line_color='black', marker_line_width=0.3,
        marker_colorbar=dict(
            title="frequency %"
        ),
    ),
)
fig = clean_plotly_fig(fig)
fig.update_layout(height=400, width=420)
fig.update_xaxes(range=[6,15], dtick=3)
fig.update_yaxes(range=[3,15], dtick=3)
pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_B.svg')
fig.show()

fig = go.Figure()
kmeans_df = kmeans_df.sort_values('cisFreq').reset_index(drop=True)
fig.add_trace(
    go.Scatter(
        x=kmeans_df['pca1'],
        y=kmeans_df['pca2'],
        # text=alleles,
        customdata=kmeans_df[['allele', 'cisFreq']], hovertemplate='%{customdata[0]}<br>'+'%{customdata[1]}',
        marker_color=kmeans_df['cisFreq'],
        marker_opacity=0.85,
        mode='markers+text',marker_colorscale='deep',textposition='bottom center', marker_cmin=2, marker_cmax=8,
        marker_line_color='black', marker_line_width=0.3,
        marker_size=9,
        marker_colorbar=dict(
            title="frequency %"
        ),
    ),
)
fig = clean_plotly_fig(fig)
fig.update_layout(height=400, width=420)
fig.update_xaxes(range=[6,15], dtick=3)
fig.update_yaxes(range=[3,15], dtick=3)
pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_C.svg')
fig.show()


COL_MAP = {
    0: 'red', 1: 'orange', 2: 'green', 3:'blue', 4:'cyan', 5: 'pink', 6: 'brown', 7:'cyan'
}

marker_colors = []
for allele in alleles:
    if cis_freqs.get(allele, 0) > 3.5:
        marker_colors.append('blue')
    elif cryp_freqs.get(allele, 0) > 6:
        marker_colors.append('purple')
    else:
        marker_colors.append('gray')

# fig = go.Figure()
# fig.add_trace(
#     go.Scatter(
#         x=kmeans_df['pca1'],
#         y=kmeans_df['pca2'],
#         customdata=kmeans_df[['allele', 'kmeans_res']], hovertemplate='%{customdata[0]}<br>'+'%{customdata[1]}',
#         # text=[allel if allel in cis_freqs else '' for allel in alleles ],
#         marker_color=[COL_MAP.get(x, 'black') for x in kmeans_df['kmeans_res'].tolist()],
#         marker_line_color='black', marker_line_width=0.5,
#         # marker_opacity=0.8,
#         # marker_color_discrete_sequence=[COL_MAP[x] for x in kmeans_res.tolist()],
#         # marker_opacity=[1 if allele in cis_freqs else 0.2 for allele in alleles],
#         mode='markers+text',textposition='bottom center',
#         marker_size=9,
#     ),
# )
# fig = clean_plotly_fig(fig)
# fig.update_layout(height=400, width=390)
# fig.update_xaxes(range=[6,15], dtick=3)
# fig.update_yaxes(range=[3,15], dtick=3)
# pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_C.svg')
# fig.show()
g_df = kmeans_df.groupby('kmeans_res', as_index=False).agg({'cisFreq': 'mean', 'crypFreq': 'mean'})
print(g_df)
fig = go.Figure()
fig.add_trace(
    go.Box(
        x=kmeans_df['kmeans_res'], y=kmeans_df['cisFreq'], fillcolor='#9BBFE5', line_color='black', line_width=0.5,
        opacity=0.8,
        boxpoints=False,
    )
)
fig.add_trace(
    go.Box(
        x=kmeans_df['kmeans_res'], y=kmeans_df['crypFreq'], fillcolor='#BA69BE', line_color='black', line_width=0.5,
        opacity=0.8,
        boxpoints=False,
    )
)
fig = clean_plotly_fig(fig)
# fig = fig.update_layout(bargap=0)
fig.update_yaxes(range=[0,20], title='frequency %')
fig.update_xaxes(title='cluster')
fig.update_layout(boxmode='group')
pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_D.svg')
fig.show()

COL_MAP_2 = {'A': 'blue', 'B': 'green', 'C': 'orange'}
# fig = go.Figure()
kmeans_df['alleleGroup'] = kmeans_df['allele'].apply(lambda x : x.split('-')[1][0])
# fig.add_trace(
#     go.Scatter(
#         x=kmeans_df['pca1'],
#         y=kmeans_df['pca2'],
#         customdata=kmeans_df[['allele', 'alleleGroup']], hovertemplate='%{customdata[0]}<br>'+'%{customdata[1]}',
#         # text=[allel if allel in cis_freqs else '' for allel in alleles ],
#         marker_color=[COL_MAP_2.get(x, 'black') for x in kmeans_df['alleleGroup'].tolist()],
#         marker_line_color='black', marker_line_width=0.5,
#         # marker_opacity=0.8,
#         # marker_color_discrete_sequence=[COL_MAP[x] for x in kmeans_res.tolist()],
#         # marker_opacity=[1 if allele in cis_freqs else 0.2 for allele in alleles],
#         mode='markers+text',textposition='bottom center',
#         marker_size=9,
#     ),
# )
# fig = clean_plotly_fig(fig)
# fig.update_layout(height=400, width=390, showlegend=True)
# fig.update_xaxes(range=[6,15], dtick=3)
# fig.update_yaxes(range=[3,15], dtick=3)
# pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_D.svg')
# fig.show()

g_df = kmeans_df.groupby('alleleGroup', as_index=False).agg({'cisFreq': 'mean', 'crypFreq': 'mean'})
fig = go.Figure()
fig.add_trace(
    go.Box(
        x=kmeans_df['alleleGroup'], y=kmeans_df['cisFreq'], fillcolor='#9BBFE5', line_color='black', line_width=0.5,
        opacity=0.8,
        boxpoints=False,
    )
)
fig.add_trace(
    go.Box(
        x=kmeans_df['alleleGroup'], y=kmeans_df['crypFreq'], fillcolor='#BA69BE', line_color='black', line_width=0.5,
        opacity=0.8,
        boxpoints=False,
    )
)
fig = clean_plotly_fig(fig)
# fig = fig.update_layout(bargap=0)
fig.update_yaxes(range=[0,20], dtick=4, title='frequency %')
fig.update_xaxes(title='cluster')
fig.update_layout(boxmode='group')
pio.write_image(fig, 'fig6/umap_returns/figS_umapReturns_E.svg')
fig.show()
# fig = go.Figure()
# fig.add_trace(
#     go.Scatter(
#         x=[pos[1] for pos in pca_df],
#         y=[pos[0] for pos in pca_df],
#         text=alleles,
#         # text=[allel if allel in cis_freqs else '' for allel in alleles ],
#         marker_color=[unmapped_freqs.get(allele, 'black') for allele in alleles],
#         marker_opacity=[1 if allele in cis_freqs else 0.2 for allele in alleles],
#         mode='markers+text',marker_colorscale='deep',textposition='bottom center', #marker_cmin=2, marker_cmax=6,
#         marker_line_color='black', marker_line_width=0.5,
#         marker_size=15,
#         marker_colorbar=dict(
#             title="unmapped %"
#         ),
#     ),
# )
# fig = clean_plotly_fig(fig)
# fig.update_layout(height=800, width=800)
# fig.show()