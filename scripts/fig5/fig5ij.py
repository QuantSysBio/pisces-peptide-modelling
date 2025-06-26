
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from scipy.stats import mannwhitneyu
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
GROUP_COLOUR_SCHEME = {
    'acidic': 'firebrick',
    'basic': 'seagreen',
    'hydrophobic': 'orange',
    'polar': 'dodgerblue',
    'special-case': 'deeppink',
    'end': 'black',
}
C_TERM_FEATS = [
    'C_term_basic', 'C_term_special-case', 'C_term_hydrophobic',
    'C_term_acidic', 'C_term_polar', 
]
C_TERM_POST_FEATS = [
    'C_term_neg_1_basic', 'C_term_neg_1_special-case', 'C_term_neg_1_hydrophobic',
    'C_term_neg_1_acidic', 'C_term_neg_1_polar', 'C_term_neg_1_end',
]
COLOURS = ['#00BFBF', 'orange']
def plot_aug():
    """
    """
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
    unique_pep_df['C-term'] = unique_pep_df.apply(
        lambda x : [m.split('_')[-1] for m in C_TERM_FEATS if x[m] == 1][0],
        axis=1,
    )
    unique_pep_df['C-term shap'] = unique_pep_df[[f'{x}_shap' for x in C_TERM_FEATS]].sum(axis=1)
    unique_pep_df['post C-term'] = unique_pep_df.apply(
        lambda x : [m.split('_')[-1] for m in C_TERM_POST_FEATS if x[m] == 1][0],
        axis=1,
    )
    unique_pep_df['post C-term shap'] = unique_pep_df[[f'{x}_shap' for x in C_TERM_POST_FEATS]].sum(axis=1)

    positions = ['C-term', 'post C-term']
    for idx, position in enumerate(positions):
        if position == 'C-term':
            for aa_grp in ['basic', 'acidic', 'hydrophobic', 'polar', 'special-case']:
                print(f'{aa_grp} C-term')
                print(mannwhitneyu(
                    unique_pep_df[unique_pep_df['C_term_basic'] == 1]['C-term shap'].dropna(),
                    unique_pep_df[unique_pep_df[f'C_term_{aa_grp}'] == 1]['C-term shap'].dropna(),
                ))
        if position == 'post C-term':
            for aa_grp in ['basic', 'acidic', 'hydrophobic', 'polar', 'special-case']:
                print(f'{aa_grp} C-term')
                print(mannwhitneyu(
                    unique_pep_df[unique_pep_df['C_term_neg_1_end'] == 1]['post C-term shap'].dropna(),
                    unique_pep_df[unique_pep_df[f'C_term_neg_1_{aa_grp}'] == 1]['post C-term shap'].dropna(),
                ))
        fig = go.Figure()
        for group, colour in GROUP_COLOUR_SCHEME.items():
            mini_col = unique_pep_df[unique_pep_df[position] == group]
            fig.add_trace(
                go.Violin(
                    x=mini_col[position], y=mini_col[f'{position} shap'],
                    fillcolor=colour, opacity=0.8,
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
            )
        fig = clean_plotly_fig(fig)
        fig.add_hline(y=0, line_width=0.5)

        fig.update_xaxes(linecolor='white')
        fig.update_yaxes(title='impact on model', dtick=1)
        if idx == 0:
            fig.update_layout(
                {
                    'yaxis1': {'range': [-2, 2], 'dtick': 1,},
                },
                width=300,height=300,
            )
            fig.show()
            fig.write_image('manuscript_figs/fig5/fig5i.svg')
        else:
            fig.update_layout(
                {
                    'yaxis1': {'range': [-3, 3]},
                },
                width=330,height=300,
            )
            fig.show()
            fig.write_image('manuscript_figs/fig5/fig5j.svg')

if __name__ == '__main__':
    plot_aug()