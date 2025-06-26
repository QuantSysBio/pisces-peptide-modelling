
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import linregress, pearsonr
import numpy as np
from ppm.constants import NUCLEOTIDE_COLOUR_SCHEME,BGD_COLOURS, START_CODONS
from ppm.analysis_utils import plot_distros
from scipy.stats import mannwhitneyu, fisher_exact, pearsonr, spearmanr
import warnings
warnings.simplefilter(action='ignore')
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

def plot_codons():
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
    # unique_pep_df = unique_pep_df[unique_pep_df['peptide'].apply(lambda x : x[-1] == 'R')]
    fig = make_subplots(rows=2, cols=len(START_CODONS), subplot_titles=START_CODONS,vertical_spacing=0.1)
    group_labels = ['background', 'detected']
    for idx, codon in enumerate(START_CODONS):
        has_codon_pos = unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df[f'{codon}_upstream'] > 0)].shape[0]
        no_codon_pos = unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df[f'{codon}_upstream'] == 0)].shape[0]
        has_codon_neg = unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df[f'{codon}_upstream'] > 0)].shape[0]
        no_codon_neg = unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df[f'{codon}_upstream'] == 0)].shape[0]
        print(codon, fisher_exact([
            [
                has_codon_pos,
                no_codon_pos,
            ],
            [
                has_codon_neg,
                no_codon_neg,
            ],
        ]))
        for label in [0, 1]:
            label_samples = unique_pep_df[unique_pep_df['label'] == label]
            label_samples['group'] = group_labels[label]

            fig.add_trace(
                go.Bar(
                    x=[group_labels[label]], y=[100*label_samples[label_samples[f'{codon}_upstream'] > 0].shape[0]/label_samples.shape[0]],
                    marker_color=BGD_COLOURS[label], marker_line_color='black',
                ),
                row=1, col=1+idx,
            )

    colors = ['#EC008C', '#99813C']
    for idx, codon in enumerate(START_CODONS):
        print(
            codon,
            mannwhitneyu(
                unique_pep_df[unique_pep_df[f'{codon}_upstream'] > 0][f'{codon}_upstream_shap'].dropna(),
                unique_pep_df[unique_pep_df[f'{codon}_upstream'] == 0][f'{codon}_upstream_shap'].dropna()
            )
        )
        for label in [0, 1]:
            if label == 1:
               mini_col = unique_pep_df[unique_pep_df[f'{codon}_upstream'] > 0]
            else:
               mini_col = unique_pep_df[unique_pep_df[f'{codon}_upstream'] == 0]
            fig.add_trace(
                go.Violin(
                    x=mini_col[f'{codon}_upstream'], y=mini_col[f'{codon}_upstream_shap'],
                    fillcolor=colors[label], opacity=0.8,
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
                row=2, col=idx+1
            )
        fig.add_hline(y=0, line_width=0.5, row=2, col=idx+1)
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        width=750, height=600, bargap=0,
    )
    for idx in range(1,6):
        fig.update_layout({f'yaxis{idx}': {'range': [0, 100]}})
    for idx in range(6,11):
        fig.update_layout({f'yaxis{idx}': {'range': [-0.5, 1]}, f'xaxis{idx}': {'title_text': 'codon present', 'linecolor': 'white'}})
    fig.update_layout({f'yaxis6': {'title_text': 'impact on model'}})
    fig.update_layout({f'yaxis1': {'title_text': 'percentage with codon upstream'}})
    fig.show()
    pio.write_image(fig, 'manuscript_figs/figS/startCodon/startCodon.svg')

if __name__ == '__main__':
    plot_codons()