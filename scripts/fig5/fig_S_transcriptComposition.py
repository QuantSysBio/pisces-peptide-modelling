
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import linregress, pearsonr
import numpy as np
from ppm.constants import NUCLEOTIDE_COLOUR_SCHEME
from ppm.analysis_utils import plot_distros
from scipy.stats import mannwhitneyu#, fisher_exact, pearsonr, spearmanr

PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']

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
def plot_rna_fracs():
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])

    fig = make_subplots(rows=2, cols=4, subplot_titles=list('ACGU'), vertical_spacing=0.1)
    nucleotide_list = list('ACGU')


    for idx, nuc in enumerate(nucleotide_list):
        print(nuc, mannwhitneyu(
            unique_pep_df[unique_pep_df['label'] == 1][f'{nuc}_frac'],
            unique_pep_df[unique_pep_df['label'] == 0][f'{nuc}_frac']
        ))
        cor_res = linregress(
            unique_pep_df[f'{nuc}_frac'],
            unique_pep_df[f'{nuc}_frac_shap']
        )
        print(nuc, cor_res)
        vals = np.linspace(0, 0.6)
        res = vals*cor_res.slope + cor_res.intercept
        plot_df = unique_pep_df.sample(n=2000, random_state=42)
        fig.add_trace(
            go.Scatter(
                x=plot_df[f'{nuc}_frac'], y=plot_df[f'{nuc}_frac_shap'],
                mode='markers', line_color=NUCLEOTIDE_COLOUR_SCHEME[nuc], name=nuc,
                # marker_line_color='black', marker_line_width=0.5,
                opacity=0.8,
            ),
            row=2, col=1+idx,
        )
        fig.add_trace(
            go.Scatter(
                x=vals, y=res, mode='lines', line_color='black', line_width=0.5, name=nuc,
                # marker_line_color='black', marker_line_width=0.5,
                opacity=0.8,
            ),
            row=2, col=1+idx,
        )

    fig = clean_plotly_fig(fig)
    for idx, nucleotide in enumerate(nucleotide_list):
        traces = plot_distros(unique_pep_df, f'{nucleotide}_frac', 'violin')
        for trace in traces:
            fig.add_trace(trace, row=1, col=1+idx)

    fig = clean_plotly_fig(fig)
    fig.update_layout(
        width=800, height=500,
    )

    for idx in range(1,5):
        fig.update_layout({f'yaxis{idx}': {'range': [0, 0.6]}})
        fig.add_hline(y=0, line_width=0.5, row=2, col=idx)
    for idx in range(5,9):
        fig.update_layout({f'yaxis{idx}': {'range': [-2, 2]}, f'xaxis{idx}': {'title_text': 'fraction', 'linecolor': 'white', 'range':[0, 0.6]}})
    fig.update_layout({f'yaxis5': {'title_text': 'impact on model'}})
    fig.update_layout({f'yaxis1': {'title_text': 'fraction'}})
    fig.show()

    pio.write_image(fig, 'manuscript_figs/figS/transcriptComposition/transcriptComposition.svg')

if __name__ == '__main__':
    plot_rna_fracs()