""" Functions for plotting results of the combined model.
"""
import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, create_comparison_logo_plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from sklearn.metrics import (
    precision_recall_curve, roc_curve, log_loss, precision_score, recall_score, auc, 
)
from scipy.stats import linregress

from ppm.analysis_utils import (
    plot_feature_importances, plot_transcriptomics, plot_ubi_distro, plot_go_distro,
)
from ppm.constants import (
    AA_COLOUR_SCHEME,
    AMINO_ACIDS,
    BGD_COLOURS,
    NUCLEOTIDE_COLOUR_SCHEME,
)
from ppm.report_template_spliced import create_spliced_report


COLOUR_DICT = {
    'canonical': '#EC9A56',
    'spliced': '#9BBFE5',
    'cryptic': '#BA69BE',
}
P1_FEATURES = [f'p1_{a_a}' for a_a in AMINO_ACIDS]
P_NEG_1_FEATURES = [f'p_neg_1_{a_a}' for a_a in AMINO_ACIDS]
P1_PRIME_FEATURES = [f'p1_prime_{a_a}' for a_a in AMINO_ACIDS]
P_NEG_1_PRIME_FEATURES = [f'p_neg_1_prime_{a_a}' for a_a in AMINO_ACIDS]


def analyse_spliced_model(config):
    print(f'{config.output_folder}/unique_peps_scored.csv')
    unique_pep_df = pd.read_csv(
        f'{config.output_folder}/unique_peps_scored.csv'
    )

    # print(ok)
    # plot_ubi_distro(unique_pep_df)
    plot_transcriptomics(unique_pep_df, config)
    plot_hotspots(unique_pep_df, config)

    plot_intv_seq_len(unique_pep_df, config)

    # plot_rna_fracs(unique_pep_df, config)
    model_stats = plot_aucs(unique_pep_df, config)

    plot_feature_importances(config)

    plot_shap_splice_site(unique_pep_df, config)
    # plot_go_distro(unique_pep_df)

    create_spliced_report(config, model_stats)


def plot_intv_seq_len(unique_pep_df, config):
    fig = make_subplots(rows=1, cols=2, subplot_titles=['', ''])
    groups = ['background', 'detected']
    pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
    print(pos_peps['interveningSeqLengths'].median())
    print(pos_peps['interveningSeqLengths'].quantile(0.25))
    print(pos_peps['interveningSeqLengths'].quantile(0.75))
    fig.add_trace(
        go.Scatter(
            x=pos_peps['interveningSeqLengths'].apply(np.log10), y=pos_peps['interveningSeqLengths_shap'],
            mode='markers', line_color=COLOUR_DICT[config.model]
        ),
        row=1, col=2,
    )
    for label in [0,1]:
        label_peps = unique_pep_df[unique_pep_df['label'] == label]
        label_peps['group'] = groups[label]

        fig.add_trace(
            go.Box(
                x=label_peps['group'], y=label_peps['interveningSeqLengths']/label_peps['protLength'], #meanline_visible=True,
                fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, boxpoints=False
            ),
            row=1, col=1, 
        )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 1], 'title_text': 'intervening sequence div. protein length'},
            'yaxis2': {'range': [-3, 3], 'title_text': 'impact on model'},
            'xaxis2': {'range': [0, 6], 'title_text': 'log10(intervening sequence)', 'linecolor': 'white'},
        },
        width=500,
    )
    fig.add_hline(y=0, line_width=0.5, row=1, col=2)

def plot_hotspots(unique_pep_df, config):
    fig = go.Figure()
    print(unique_pep_df.columns)

    can_peps = []
    mean_shaps = []
    for idx in range(25):
        mini_col = unique_pep_df[unique_pep_df['nCanonicalPeptides'] == idx]
        if mini_col.shape[0]:
            mean_shaps.append(
                mini_col['nCanonicalPeptides_shap'].mean()
            )
            can_peps.append(idx)
        # fig.add_trace(
        #     go.Violin(
        #         x=mini_col['nCanonicalPeptides'], y=mini_col[f'nCanonicalPeptides_shap'],
        #         fillcolor='pink', opacity=0.9,
        #         line_color='black',
        #         line_width=0.5,
        #         meanline_visible=True,
        #         points=False,
        #     ),
        # )
    fig.add_trace(
        go.Scatter(
            x=can_peps, y=mean_shaps,
            mode='lines', opacity=0.9,
            line_color='firebrick',
            # line_width=0.5,
        ),
    )

    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [-0.5, 1]},
        },
        width=800,height=300,
    )
    fig.update_xaxes(linecolor='white', title='canonical peptides from protein')
    fig.update_yaxes(title='mean impact on model')

    fig = go.Figure()
    print(unique_pep_df.columns)
    colours = ['wheat', 'dodgerblue']
    for idx in range(2):
        mini_col = unique_pep_df[unique_pep_df['label'] == idx]
        a_df = mini_col.groupby('nCanonicalPeptides', as_index=False)['peptide'].count()
        a_df['peptide'] /= a_df['peptide'].sum()
        fig.add_trace(
            go.Bar(
                x=a_df['nCanonicalPeptides'], y=a_df[f'peptide'],
                marker_color=colours[idx], opacity=0.9,
                marker_line_color='black',
                marker_line_width=0.5,
                # meanline_visible=True,
                # points=False,
            ),
        )

    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 1]},
            'xaxis1': {'range': [-0.5, 25.5], 'tick0':0, 'dtick':5,},
        },
        width=800,height=300,bargap=0,
    )
    fig.update_xaxes(linecolor='white', title='canonical peptides from protein')
    fig.update_yaxes(title='fraction of peptides')
    # pio.write_image(fig, f'{config.output_folder}/imgs/stratum_shap.svg')

    fig = make_subplots(rows=2, cols=2, subplot_titles=['SR1', 'SR2'], vertical_spacing=0.1)
    groups = ['background', 'detected']
    for sr in [1,2]:
        filt_peps = unique_pep_df[unique_pep_df['nCanonicalPeptides'] > 0]
        for label in [0, 1]:
            lab_peps = filt_peps[filt_peps['label'] == label]
            fig.add_trace(
                go.Violin(
                    x=[groups[label]]*len(lab_peps), y=lab_peps[f'sr{sr}_can_dist']/lab_peps['protLength'],
                    fillcolor=BGD_COLOURS[label], opacity=0.8,
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
                row=1, col=sr
            )
            if label == 1:
                fig.add_trace(
                    go.Scatter(
                        x=lab_peps[f'sr{sr}_can_dist'].apply(np.log10), y=lab_peps[f'sr{sr}_can_dist_shap'],
                        fillcolor=BGD_COLOURS[label], opacity=0.8, mode='markers',
                        marker_color=COLOUR_DICT[config.model]
                        # line_color='black',
                        # line_width=0.5,
                        # meanline_visible=True,
                        # boxpoints=False,
                    ),
                    row=2, col=sr
                )

    fig = clean_plotly_fig(fig)
    # fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 1], 'title_text': 'sr distance/protein length'},
            'yaxis2': {'range': [0, 1]},
            'xaxis3': {'title_text': 'log10(sr distance)', 'range': [0,4]},
            'yaxis3': {'range': [-1, 2], 'title_text': 'impact on model', 'dtick': 1},
            'xaxis4': {'title_text': 'log10(sr distance)', 'range': [0,4]},
            'yaxis4': {'range': [-1, 4]},
        },
        width=500,height=600,bargap=0,
    )
    # fig.update_xaxes(linecolor='white', title='canonical peptides from protein')
    # fig.update_yaxes(title='fraction of peptides')

def plot_aucs(unique_pep_df, config):

    fig = make_subplots(rows=1, cols=2, subplot_titles=['ROC', 'PR'])
    fprs, tprs, _ = roc_curve(unique_pep_df['label'], unique_pep_df['prediction_xgb2'])
    precisions, recalls, _ = precision_recall_curve(unique_pep_df['label'], unique_pep_df['prediction_xgb2'])

    fig.add_trace(
        go.Scatter(x=fprs, y=tprs, mode='lines', line_color=COLOUR_DICT[config.model], name='spliced'),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(x=recalls, y=precisions, mode='lines', line_color=COLOUR_DICT[config.model], name='spliced'),
        row=1, col=2,
    )

    fig = clean_plotly_fig(fig)
    fig.update_layout(width=600)
    fig.update_xaxes(range=[0,1])
    fig.update_yaxes(range=[0,1])
    pio.write_image(fig, f'{config.output_folder}/imgs/model_performance.svg')
    return {
        'log_loss': log_loss(unique_pep_df['label'], unique_pep_df['prediction_xgb2']),
        'auc-pr': auc(recalls, precisions),
        'auc-roc': auc(fprs, tprs),
        'precision': precision_score(unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)),
        'recall': recall_score(unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)),
    }


def plot_shap_splice_site(unique_pep_df, config):
    """
    """
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
    fig = make_subplots(
        rows=2, cols=1, subplot_titles=positions
    )
    for idx, position in enumerate(positions):
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
                row=idx+1, col=1
            )
    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [-2, 2]},
            'yaxis2': {'range': [-2, 2]},
        },
        width=800,height=600,
    )
    fig.update_xaxes(linecolor='white')
    fig.update_yaxes(title='impact on model', dtick=1)
    pio.write_image(fig, f'{config.output_folder}/imgs/spliced_site-shap.svg')
    get_splice_site_entropy(unique_pep_df, config)

def plot_shap_stratum(unique_pep_df, config):
    """
    """
    fig = go.Figure()
    for stratum, colour in STRATUM_COLOUR_SCHEME.items():
        mini_col = unique_pep_df[unique_pep_df['stratum'] == stratum]
        fig.add_trace(
            go.Violin(
                x=mini_col['stratum'], y=mini_col[f'stratum_shap'],
                fillcolor=colour, opacity=0.9,
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
            'yaxis1': {'range': [-2, 6]},
        },
        width=800,height=300,
    )
    fig.update_xaxes(linecolor='white')
    fig.update_yaxes(title='impact on model')
    pio.write_image(fig, f'{config.output_folder}/imgs/stratum_shap.svg')


def get_splice_site_entropy(unique_pep_df, config):
    unique_pep_df = unique_pep_df.drop_duplicates(subset='peptide')
    unique_pep_df = unique_pep_df[[
        'peptide', 'label',
        'p2', 'p1', 'p_minus_1', 'p_minus_2',
        'p_minus_2_prime', 'p_minus_1_prime', 'p1 prime', 'p2_prime'
    ]].drop_duplicates('peptide')
    unique_pep_df['peptide'] = unique_pep_df['p2'] + unique_pep_df['p1'] + unique_pep_df['p_minus_1'] + unique_pep_df['p_minus_2']

    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'sr1 splice site: detected', 'background'],
        4, f'{config.output_folder}/imgs',
        file_name=f'sr1_splice_site.svg', amino_acids=AMINO_ACIDS+'X',
        y_lim=0.04, plot_size=5,
        x_ticks=['', 'p1', '', ''],
        vline=2.5,
    )

    unique_pep_df['peptide'] = unique_pep_df['p_minus_2_prime'] + unique_pep_df['p_minus_1_prime'] + unique_pep_df['p1 prime'] + unique_pep_df['p2_prime']

    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'sr 2 splice site: detected', 'background'],
        4, f'{config.output_folder}/imgs',
        file_name=f'sr2_splice_site.svg', amino_acids=AMINO_ACIDS+'X',
        y_lim=0.04, plot_size=5,
        x_ticks=['', '', "p1'", ''],
        vline=2.5,
    )


def plot_rna_fracs(unique_pep_df, config):
    fig = make_subplots(rows=2, cols=2, subplot_titles=list('ACGU'))
    nucleotide_list = list('ACGU')
    for idx, nuc in enumerate(nucleotide_list):
        cor_res = linregress(
            unique_pep_df[f'{nuc}_frac'],
            unique_pep_df[f'{nuc}_frac_shap']
        )

        vals = np.linspace(0, 0.6)
        res = vals*cor_res.slope + cor_res.intercept
        plot_df = unique_pep_df.sample(n=2000, random_state=42)
        fig.add_trace(
            go.Scatter(x=plot_df[f'{nuc}_frac'], y=plot_df[f'{nuc}_frac_shap'], mode='markers', line_color=NUCLEOTIDE_COLOUR_SCHEME[nuc], name=nuc),
            row=1+(idx//2), col=1+(idx%2),
        )
        fig.add_trace(
            go.Scatter(x=vals, y=res, mode='lines', line_color='black', line_width=0.5, name=nuc),
            row=1+(idx//2), col=1+(idx%2),
        )


    fig = clean_plotly_fig(fig)
    fig.update_layout({
            'yaxis1': {'title': 'impact on model'},
            'yaxis3': {'title': 'impact on model'}
        },
        width=600, height=600,
    )
    fig.update_xaxes(range=[0,1], title='nucleotide fraction')
    fig.update_yaxes(title='impact on model')

    fig.add_hline(y=0, line_width=0.5)
    fig.update_xaxes(linecolor='white')
    pio.write_image(fig, f'{config.output_folder}/imgs/rna_frac_shap.svg')

    fig = make_subplots(rows=1, cols=4, subplot_titles=nucleotide_list)
    pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
    neg_peps = unique_pep_df[unique_pep_df['label'] == 0]
    pos_peps['group'] = 'detected'
    neg_peps['group'] = 'background'
    for idx, nucleotide in enumerate(nucleotide_list):
        fig.add_trace(
            go.Violin(
                x=pos_peps['group'], y=pos_peps[f'{nucleotide}_frac'], meanline_visible=True,
                fillcolor='#ADD8E6', line_color='black', line_width=0.5, points=False),
            row=1, col=1+idx, 
        )
        fig.add_trace(
            go.Violin(x=neg_peps['group'], y=neg_peps[f'{nucleotide}_frac'], fillcolor='wheat', line_color='black',
                      line_width=0.5, points=False, meanline_visible=True),
            row=1, col=1+idx
        )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        width=800, height=400,
    )
    fig.update_yaxes(
        range=[0,0.5]
    )
    pio.write_image(fig, f'{config.output_folder}/imgs/rna_frac_distro.svg')

