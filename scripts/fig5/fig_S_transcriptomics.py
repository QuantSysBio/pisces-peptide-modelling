from functools import reduce
from operator import or_
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from plotly.subplots import make_subplots
from scipy.stats import linregress, pearsonr
import numpy as np
from ppm.constants import NUCLEOTIDE_COLOUR_SCHEME, BGD_COLOURS
from ppm.analysis_utils import plot_distros, get_expr_freq
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

def get_expressed_count(label_peps, fractions):
    return (
        label_peps[reduce(
            or_,
            [label_peps[f'tr_TPM_K562_{fraction}'] > 0 for fraction in fractions]
        )].shape[0]
    )

def plot_transcriptomics(unique_pep_df, cell_line, stratum):
    if 'stratum' in unique_pep_df.columns:
        unique_pep_df = unique_pep_df[unique_pep_df['stratum'] != CRYPTIC_STRATA.index('intergenic')]
    if cell_line == 'K562':
        fractions = ['poly', 'S80', 'free', 'bulk']
        fig = make_subplots(rows=1, cols=5, subplot_titles=[
            'any fraction > 0',
            'poly abundance',
            'S80 abundance',
            'free abundance',
            'bulk abundance',
        ])
        group_labels = ['background', 'detected']
        expressed_count_pos = get_expressed_count(unique_pep_df[unique_pep_df['label'] == 1], fractions)
        not_expressed_count_pos = unique_pep_df[unique_pep_df['label'] == 1].shape[0] - expressed_count_pos
        expressed_count_neg = get_expressed_count(unique_pep_df[unique_pep_df['label'] == 0], fractions)
        not_expressed_count_neg = unique_pep_df[unique_pep_df['label'] == 0].shape[0] - expressed_count_neg
        print(expressed_count_pos, not_expressed_count_pos, expressed_count_neg, not_expressed_count_neg)
        print(cell_line, stratum, fisher_exact([
            [
                expressed_count_pos,
                not_expressed_count_pos,
            ],
            [
                expressed_count_neg,
                not_expressed_count_neg,
            ],
        ]))
        for label in [0, 1]:
            label_peps = unique_pep_df[unique_pep_df['label'] == label]
            pct_expr = get_expr_freq(label_peps, fractions)
            if 'stratum' in unique_pep_df.columns and label == 1:
                stratum_counts = []
                for stratum in CRYPTIC_STRATA:
                    if stratum != 'intergenic':
                        print(stratum)
                        strat_df = label_peps[label_peps['stratum'] == CRYPTIC_STRATA.index(stratum)]
                        strat_pct_expr = get_expr_freq(strat_df, fractions)
                        stratum_counts.append(strat_pct_expr)
            print(label, pct_expr)


        

            fig.add_trace(
                go.Bar(
                    x=[group_labels[label]], y=[pct_expr],
                    marker_color=BGD_COLOURS[label], marker_line_color='black',
                ),
                row=1, col=1,
            )
        for idx, fraction in enumerate(fractions):
            pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
            neg_peps = unique_pep_df[unique_pep_df['label'] == 0]
            pos_peps['group'] = 'detected'
            neg_peps['group'] = 'background'
            pos_tr_df = pos_peps[pos_peps[f'tr_TPM_K562_{fraction}'] > 0]
            neg_tr_df = neg_peps[neg_peps[f'tr_TPM_K562_{fraction}'] > 0]
            print(
                stratum, fraction,
                mannwhitneyu(
                    pos_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10).dropna(),
                    neg_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10).dropna(),
            ))
            fig.add_trace(
                go.Violin(
                    x=neg_tr_df['group'], y=neg_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10),
                    fillcolor=BGD_COLOURS[0],
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
                row = 1, col=2+idx,
            )
            fig.add_trace(
                go.Violin(
                    x=pos_tr_df['group'], y=pos_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10),
                    fillcolor=BGD_COLOURS[1],
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
                row = 1, col=2+idx,
            )


        fig = clean_plotly_fig(fig)
        fig.update_layout({
                'yaxis1': {'range': [0,100]},
                'yaxis2': {'range': [-6,6]},
                'yaxis3': {'range': [-6,6]},
                'yaxis4': {'range': [-6,6]},
                'yaxis5': {'range': [-6,6]},
            },
            width=1000, height=300,
        )
    else:
        fig = make_subplots(rows=1, cols=2, subplot_titles=[
            'bulk > 0',
            'bulk abundance',
        ])
        traces = plot_distros(
            unique_pep_df, 'tr_TPM_721', 'bar', display_range=['x']
        )
        for trace in traces:
            fig.add_trace(trace, row=1, col=1)

        
        print(cell_line, stratum, fisher_exact([
            [
                unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df['tr_TPM_721'] > 0)].shape[0],
                unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df['tr_TPM_721'] == 0)].shape[0],
            ],
            [
                unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df['tr_TPM_721'] > 0)].shape[0],
                unique_pep_df[(unique_pep_df['label'] == 0) & (unique_pep_df['tr_TPM_721'] == 0)].shape[0],
            ],
        ]))


        unique_pep_df['tr_TPM_721'] = unique_pep_df['tr_TPM_721'].apply(np.log10)

        print(
            stratum, cell_line,
            mannwhitneyu(
                unique_pep_df[unique_pep_df['label'] == 1]['tr_TPM_721'].dropna(),
                unique_pep_df[unique_pep_df['label'] == 0]['tr_TPM_721'].dropna(),
        ))
        traces = plot_distros(
            unique_pep_df, 'tr_TPM_721', 'violin'
        )
        for trace in traces:
            fig.add_trace(trace, row=1, col=2)

        fig.update_layout({
                'yaxis1': {'range': [0,100]},
                'yaxis2': {'range': [-6,6]},
                # 'yaxis3': {'range': [-6,6]},
                # 'yaxis4': {'range': [-6,6]},
                # 'yaxis5': {'range': [-6,6]},
            },
            # width=1000, height=300,
        )


    # if 'stratum' in unique_pep_df.columns:
    #     stratum_counts = []
    #     fig = go.Figure()
    #     for idx, stratum in enumerate(CRYPTIC_STRATA):
    #         if stratum == 'intergenic':
    #             continue
    #         s_df = unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df['stratum'] == idx)]
    #         s_c = 100*s_df[s_df['tr_TPM_721'] > 0].shape[0]/s_df.shape[0]

    #         fig.add_trace(
    #             go.Bar(
    #                 x=[stratum], y=[s_c],
    #                 marker_color=STRATUM_COLOUR_SCHEME[stratum], opacity=0.8,
    #                 marker_line_color='black',
    #                 marker_line_width=0.5,
    #             )
    #         )
    fig = clean_plotly_fig(fig)

    width=400
    if cell_line == 'K562':
        width=1000
    fig.update_layout(
        # {
        #     'yaxis1': {'range': [0, 100], 'title_text': 'percentage expressed'},
        # },
        width=width,height=300,bargap=0,
    )
    # fig.update_xaxes(linecolor='white')
    # fig.show()
    pio.write_image(fig, f'manuscript_figs/figS/transcriptomics/transcriptomic_distro_{cell_line}_{stratum}.svg')
PROJECTS = {
    'K562': [
        'output/final/canonicalK562',
        'output/final/crypticK562',
        'output/final/splicedK562',
    ],
    '721.221': [
        'output/final/canonical721',
        'output/final/cryptic721',
        'output/final/spliced721',
    ],
}
STRATA = ['canonical', 'cryptic', 'spliced']
if __name__ == '__main__':
    for cell_line in PROJECTS:
        for idx, df_loc in enumerate(PROJECTS[cell_line]):
            pep_df = pd.read_csv(f'{df_loc}/unique_peps_scored.csv')
            plot_transcriptomics(pep_df, cell_line, STRATA[idx])