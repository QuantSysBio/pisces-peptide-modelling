
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import linregress, pearsonr
import numpy as np
from ppm.constants import NUCLEOTIDE_COLOUR_SCHEME, BGD_COLOURS
from ppm.analysis_utils import plot_distros, get_expr_freq
from scipy.stats import fisher_exact, pearsonr, spearmanr


STRATUM_COLOUR_SCHEME = {
    'fiveUTR': 'orange',
    'threeUTR': 'yellow',
    'CDS_frameshift': 'purple',
    'lncRNA': 'darkgrey',
    'intronic': 'forestgreen',
    'intergenic': 'navy',
}
CRYPTIC_STRATA = [
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
]
COLOURS = ['#00BFBF', 'orange']


def plot_proteomics(unique_pep_df, cell_line, stratum):
    if 'stratum' in unique_pep_df.columns:
        unique_pep_df = unique_pep_df[
            (unique_pep_df['stratum'].isin([
                CRYPTIC_STRATA.index('fiveUTR'),
                CRYPTIC_STRATA.index('threeUTR'),
                CRYPTIC_STRATA.index('CDS_frameshift'),
                CRYPTIC_STRATA.index('intronic'),
            ]))
        ]
    if cell_line == 'K562':
        proteomics_column = 'proteomics_K562'
    else:
        proteomics_column = 'proteomics_B721'
    print(cell_line, stratum, unique_pep_df[unique_pep_df['label'] == 1][proteomics_column].mean(), fisher_exact([
        [
            unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df[proteomics_column] == 1)].shape[0],
            unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df[proteomics_column] == 0)].shape[0],
        ],
        [
            unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df[proteomics_column] == 1)].shape[0],
            unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df[proteomics_column] == 0)].shape[0],
        ],
    ]))
    fig = go.Figure()
    traces = plot_distros(
        unique_pep_df, proteomics_column, 'bar', display_range=['x']
    )
    for trace in traces:
        fig.add_trace(trace)


    fig.update_layout({
            'yaxis1': {'range': [0,100]},
        },
    )

    fig = clean_plotly_fig(fig)

    fig.update_layout(
        # {
        #     'yaxis1': {'range': [0, 100], 'title_text': 'percentage expressed'},
        # },
        width=100,height=300,bargap=0,
    )
    # fig.show()
    # fig.update_xaxes(linecolor='white')
    # fig.show()
    pio.write_image(fig, f'manuscript_figs/figS/proteomics/proteomics_{cell_line}_{stratum}.svg')
    # if cell_line == 'K562':
    #     return
    # else:
    #     fig.show()
PROJECTS = {
    'K562': [
        'output2/final/canonicalK562',
        'output2/final/crypticK562',
        'output2/final/splicedK562',
    ],
    '721.221': [
        'output2/final/canonical721',
        'output2/final/cryptic721',
        'output2/final/spliced721',
    ],
}
STRATA = ['canonical', 'cryptic', 'spliced']
if __name__ == '__main__':
    for cell_line in PROJECTS:
        for idx, df_loc in enumerate(PROJECTS[cell_line]):
            pep_df = pd.read_csv(f'{df_loc}/unique_peps_scored.csv')
            plot_proteomics(pep_df, cell_line, STRATA[idx])