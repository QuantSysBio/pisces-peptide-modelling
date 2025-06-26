from math import ceil
import numpy as np
import polars as pl
from pisces.plot_utils import clean_plotly_fig, DARK_COLOURS, LIGHT_COLORS
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import pandas as pd

MAIN_NAMES = ['canonical (t/d)', 'canonical (pisces)', 'contaminant', 'non-canonical']
NAMES = ['spliced', 'multi-mapped', 'cryptic', 'unmapped',]
ALL_NAMES = MAIN_NAMES + NAMES
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
    'TrEMBL',
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
    'mutation',
    'fusion',
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

def main():
    remapped_dfs = []
    pep_df = pl.read_parquet('data/piscesPeptides.parquet')
    # pep_df = pep_df.filter(pl.col('cellLines').list.contains('K562') & pl.col('peptide').str.len_chars().is_in([9,10,11,12]))

    # print(cryptic_df.columns)
    fig = plot_cryptic_breakdown(
        pep_df, 
    )
    fig.show()
    pio.write_image(fig, 'manuscript_figs/figS/cryptic_counts/figS_cryptic_AB.svg')


def plot_cryptic_breakdown(pep_df, ):
    """ Create bar plots of the different cryptic strata.
    """
    mm_df = pep_df.filter(pl.col('stratum') == 'multi-mapped')
    count_df_spliced = get_cryptic_stratum_counts(mm_df)

    cryptic_df = pep_df.filter(pl.col('stratum') == 'cryptic')
    count_df_not_spliced = get_cryptic_stratum_counts(cryptic_df)

    fig = make_subplots(
        rows=1, cols=2,
    #     subplot_titles=[
    #         'Unique to Cryptic Strata', 'Multi-mapped to Spliced',
    #     ],
    )
    # fig = go.Figure()

    for col_idx, c_df in enumerate([count_df_not_spliced, count_df_spliced]):
        for idx, df_row in c_df.iterrows():
            fig.add_trace(
                go.Bar(
                    x=[df_row['stratum']],
                    y=[df_row['uniqueCount']],
                    marker_color=DARK_COLOURS[idx],
                    marker_line_color='black',
                    marker_line_width=0.25,
                ),
                row=1,
                col=col_idx+1,
            )


    fig.update_layout(
        barmode='overlay',
    )
    fig = clean_plotly_fig(fig)
    fig.update_yaxes(title='', range=[0, 3000], dtick=1000)
    fig['layout']['yaxis']['title'] = '# Peptides'

    fig.update_layout(
        height=300, width=600, bargap=0,
    )

    return fig

if __name__ == '__main__':
    main()