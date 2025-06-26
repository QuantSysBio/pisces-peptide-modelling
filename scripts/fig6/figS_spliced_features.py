import numpy as np
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import plotly.figure_factory as ff
from ppm.constants import AMINO_ACIDS

PROJECTS = ['output/splicedK562_241129', 'output/spliced721_241206']


COLOURS = ['#00BFBF', 'orange']
AA_COLOUR_SCHEME = {
    'P': 'deeppink',
    'M': 'orange',
    'A': 'orange',
    'V': 'orange',
    'I': 'orange',
    'L': 'orange',
    'F': 'orange',
    'Y': 'orange',
    'W': 'orange',
    'H': 'seagreen',
    'R': 'seagreen',
    'K': 'seagreen',
    'D': 'firebrick',
    'E': 'firebrick',
    'N': 'dodgerblue',
    'Q': 'dodgerblue',
    'S': 'dodgerblue',
    'T': 'dodgerblue',
    'G': 'dodgerblue',
    'C': 'dodgerblue',
    'X': 'black',
}

P1_FEATURES = [f'p1_{a_a}' for a_a in AMINO_ACIDS]
P_NEG_1_FEATURES = [f'p_neg_1_{a_a}' for a_a in AMINO_ACIDS]
P1_PRIME_FEATURES = [f'p1_prime_{a_a}' for a_a in AMINO_ACIDS]
P_NEG_1_PRIME_FEATURES = [f'p_neg_1_prime_{a_a}' for a_a in AMINO_ACIDS]
def plot_shap_splice_site():
    """
    """
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
    unique_pep_df['p1'] = unique_pep_df.apply(
        lambda x : [a_a for a_a in AMINO_ACIDS if x[f'p1_{a_a}'] == 1][0],
        axis=1,
    )
    unique_pep_df['p1 shap'] = unique_pep_df[[f'{a}_shap' for a in P1_FEATURES]].sum(axis=1)
    unique_pep_df['p1 prime'] = unique_pep_df.apply(
        lambda x : [a_a for a_a in AMINO_ACIDS if x[f'p1_prime_{a_a}'] == 1][0],
        axis=1,
    )
    unique_pep_df['p1 prime shap'] = unique_pep_df[[f'{x}_shap' for x in P1_PRIME_FEATURES]].sum(axis=1)
    positions = ['p1', "p1 prime"]

    for idx, position in enumerate(positions):
        fig = go.Figure()
        for a_a in AMINO_ACIDS:
            mini_col = unique_pep_df[unique_pep_df[position] == a_a]
            fig.add_trace(
                go.Violin(
                    x=mini_col[position], y=mini_col[f'{position} shap'],
                    fillcolor=AA_COLOUR_SCHEME[a_a], opacity=0.8,
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
            )
        fig = clean_plotly_fig(fig)
        fig.add_hline(y=0, line_width=0.5)
        fig.update_layout(
            {
                'yaxis1': {'range': [-2, 2]},
            },
            width=500,height=200,
        )
        fig.update_xaxes(linecolor='white')
        fig.update_yaxes(title='impact on model', dtick=1)
        if position == 'p1':
            pio.write_image(fig, 'fig6/img/FigS_splicedFeats_A.svg')
        else:
            pio.write_image(fig, 'fig6/img/FigS_splicedFeats_B.svg')
        fig.show()
        # else:
        #     pio.write_image(fig, f'{config.output_folder}/imgs/FigS_slicing_features_A.svg')

if __name__ == '__main__':
    plot_shap_splice_site()
