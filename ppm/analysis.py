""" Functions for plotting results of the combined model.
"""

import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, create_comparison_logo_plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from scipy.stats import linregress, pearsonr

from ppm.analysis_utils import (
    plot_feature_importances, plot_aucs, plot_transcriptomics, plot_distros, plot_shap, plot_ubi_distro, plot_go_distro
)
from ppm.constants import (
    AMINO_ACIDS,
    BGD_COLOURS,
    NUCLEOTIDE_COLOUR_SCHEME,
    START_CODONS,
    STRATUM_COLOUR_SCHEME,
)
from ppm.report_template_cryptic import create_cryptic_report
from ppm.report_template_canonical import create_canonical_report

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
GROUP_COLOUR_SCHEME = {
    'acidic': 'firebrick',
    'basic': 'seagreen',
    'hydrophobic': 'orange',
    'polar': 'dodgerblue',
    'special-case': 'deeppink',
    'end': 'black',
}
CELL_LINE = 'K562'
C_TERM_FEATS = [
    'C_term_basic', 'C_term_special-case', 'C_term_hydrophobic',
    'C_term_acidic', 'C_term_polar', 
]
C_TERM_POST_FEATS = [
    'C_term_neg_1_basic', 'C_term_neg_1_special-case', 'C_term_neg_1_hydrophobic',
    'C_term_neg_1_acidic', 'C_term_neg_1_polar', 'C_term_neg_1_end',
]
CRYPTIC_STRATA = [
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
]



def plot_hydrophobicity(unique_pep_df, config):
    fig = make_subplots(rows=1, cols=2)
    print(unique_pep_df['peptide_hydrophobicity'])
    violin_traces = plot_distros(unique_pep_df, 'peptide_hydrophobicity', 'violin')
    print(len(violin_traces))
    for trace in violin_traces:
        print('\ttrace', type(trace))
        fig.add_trace(trace, row=1, col=1)

    scatter_traces = plot_shap(unique_pep_df, 'peptide_hydrophobicity', 'scatter')
    # print(scatter_traces)
    for trace in scatter_traces:
        fig.add_trace(trace, row=1, col=2)

    fig = clean_plotly_fig(fig)
    fig.update_layout(
        {'yaxis1': dict(range=[-4, 4], title_text='peptide hydrophobicity'),
        'yaxis2': dict(range=[-8, 4], title_text='impact on model'),
        'xaxis2': dict(range=[-4, 4], title_text='peptide hydrophobicity')},
    )
    unique_pep_df[['peptide', 'peptide_hydrophobicity']].to_csv('ok.csv', index=False)

def plot_fragment_length(pep_df, config):
    pep_df['fragmentLength'] = pep_df['start_dist'] + pep_df['stopDistances']
    titles = ['all', '']
    for stratum in CRYPTIC_STRATA:
        titles.append(stratum)
        titles.append('')
    fig = go.Figure()

    traces = plot_distros(pep_df, 'fragmentLength', 'scatter', display_range=[0, 51])
    for trace in traces:
        fig.add_trace(trace)



    fig = clean_plotly_fig(fig)


def plot_start_dist(pep_df, config):
    print(pep_df['AUG_upstream'].mean())
    print(pep_df[['AUG_upstream', 'start_dist']])
    titles = ['all', '']
    for stratum in CRYPTIC_STRATA:
        titles.append(stratum)
        titles.append('')
    fig = make_subplots(rows=7, cols=2, subplot_titles=titles)

    traces = plot_distros(pep_df, 'start_dist', 'scatter', display_range=[0, 10])
    for trace in traces:
        fig.add_trace(trace, row=1, col=1)

    traces = plot_shap(pep_df, 'start_dist', 'line', display_range=[0,10])
    for trace in traces:
        fig.add_trace(trace, row=1, col=2)

    if config.model == 'cryptic':
        for strat_idx, stratum in enumerate(CRYPTIC_STRATA):
            strat_df = pep_df[pep_df['stratum'] == strat_idx]
            traces = plot_distros(strat_df, 'start_dist', 'scatter', display_range=[0, 10])
            for trace in traces:
                fig.add_trace(trace, row=2+strat_idx, col=1)

            traces = plot_shap(strat_df, 'start_dist', 'line', display_range=[0,10])
            for trace in traces:
                fig.add_trace(trace, row=2+strat_idx, col=2)


    fig = clean_plotly_fig(fig)
    fig.update_layout({
        'yaxis1': dict(range=[0, 1], title_text='fraction of peptides'),
        'yaxis2': dict(range=[-1, 2], title_text='impact on model'),
        'xaxis1': dict(range=[0, 10], title_text='start_dist'),
        'xaxis2': dict(range=[0, 10], title_text='start_dist')
    }, width=1000, height=1200,)


def plot_stop_dist(pep_df, config):
    print(pep_df['stopDistances'].mean())
    titles = ['all', '']
    for stratum in CRYPTIC_STRATA:
        titles.append(stratum)
        titles.append('')
    fig = make_subplots(rows=7, cols=2, subplot_titles=titles)

    traces = plot_distros(pep_df, 'stopDistances', 'scatter', display_range=[0, 10])
    for trace in traces:
        fig.add_trace(trace, row=1, col=1)


    if config.model == 'cryptic':

        traces = plot_shap(pep_df, 'stopDistances', 'line', display_range=[0,10])
        for trace in traces:
            fig.add_trace(trace, row=1, col=2)

        for strat_idx, stratum in enumerate(CRYPTIC_STRATA):
            strat_df = pep_df[pep_df['stratum'] == strat_idx]
            traces = plot_distros(strat_df, 'stopDistances', 'scatter', display_range=[0, 10])
            for trace in traces:
                fig.add_trace(trace, row=2+strat_idx, col=1)

            traces = plot_shap(strat_df, 'stopDistances', 'line', display_range=[0,10])
            for trace in traces:
                fig.add_trace(trace, row=2+strat_idx, col=2)


    fig = clean_plotly_fig(fig)
    fig.update_layout({
        'yaxis1': dict(range=[0, 1], title_text='fraction of peptides'),
        'yaxis2': dict(range=[-1, 2], title_text='impact on model'),
        'xaxis1': dict(range=[0, 10], title_text='stopDistances'),
        'xaxis2': dict(range=[0, 10], title_text='stopDistances')
    }, width=1000, height=1200,)

def plot_g_frac(unique_pep_df):
    groups = ['background', 'detected']
    fig = make_subplots(rows=1, cols=3)
    frame_df = unique_pep_df[unique_pep_df['stratum'] == 1]
    other_df = unique_pep_df[unique_pep_df['stratum'] != 1]

    for g_idx, st_df in enumerate([frame_df, other_df]):
        traces = plot_distros(st_df, 'G_frac', 'box')
        for trace in traces:
            fig.add_trace(trace, row=1, col=1+g_idx)


    pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
    for idx, stratum in enumerate(CRYPTIC_STRATA):
        label_peps = pos_peps[pos_peps['stratum'] == idx]
        fig.add_trace(
            go.Box(
                x=[stratum]*len(label_peps), y=label_peps['G_frac'], #meanline_visible=True,
                fillcolor=STRATUM_COLOUR_SCHEME[stratum], line_color='black', line_width=0.5, boxpoints=False
            ),
            row=1, col=3
        )
    fig.update_yaxes(range=[0,1])
    fig.update_layout(
        {'yaxis1': {'title_text': 'Position in protein'}},
    )
    fig = clean_plotly_fig(fig)


def plot_kozak(unique_pep_df, config):
    groups = ['background', 'detected']
    fig = go.Figure()

    traces = plot_distros(unique_pep_df, 'kozakScore', 'box')
    for trace in traces:
        fig.add_trace(trace)
    # for label in [0,1]:
    #     label_peps = st_df[st_df['label'] == label]
    #     label_peps['group'] = groups[label]

    #     fig.add_trace(
    #         go.Violin(
    #             x=label_peps['group'], y=label_peps['relativePosition'], meanline_visible=True,
    #             fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, points=False, spanmode = 'hard'
    #         ),
    #         row=1, col=1
    #     )
    #     if g_idx == 0:
    #         m_c = STRATUM_COLOUR_SCHEME['CDS_frameshift']
    #     else:
    #         m_c = 'green'
    #     fig.add_trace(
    #         go.Scatter(
    #             x=st_df['relativePosition'], y=st_df['relativePosition_shap'], mode='markers',
    #             marker_color=m_c
    #             #meanline_visible=True,
    #             # fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, boxpoints=False
    #         ),
    #         row=2, col=1+g_idx
    #     )

    # pos_peps = unique_pep_df[unique_pep_df['label'] == label]
    # for idx, stratum in enumerate(CRYPTIC_STRATA):
    #     label_peps = pos_peps[pos_peps['stratum'] == idx]

    #     fig.add_trace(
    #         go.Violin(
    #             x=[stratum]*len(label_peps), y=label_peps['relativePosition'], meanline_visible=True,
    #             fillcolor=STRATUM_COLOUR_SCHEME[stratum], line_color='black', line_width=0.5, points=False, spanmode = 'hard'
    #         ),
    #         row=1, col=3
    #     )
    # fig.update_yaxes(range=[0,1])
    # fig.update_layout(
    #     {'yaxis1': {'title_text': 'Position in protein'}},
    # )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        {
            'yaxis': {'range': [0, 1], 'title_text': 'Kozak Similarity Score'},
            # 'xaxis4': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            # 'xaxis5': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            # 'yaxis5': {'range': [-2, 2], 'dtick': 1,},
            # 'yaxis1': {'title_text': 'relative position'},
        },
        width=200,height=300
    )
    # fig.add_hline(y=0, line_width=0.5, row=2, col=1)
    # fig.add_hline(y=0, line_width=0.5, row=2, col=2)


def plot_rel_pos(unique_pep_df):
    groups = ['background', 'detected']
    fig = make_subplots(rows=2, cols=4, specs=[
        [{}, {},{'colspan': 2}, None],
        [{}, {},{'colspan': 2}, None],
        ]
    )
    frame_df = unique_pep_df[unique_pep_df['stratum'] == 2]
    other_df = unique_pep_df[unique_pep_df['stratum'] != 2]

    for g_idx, st_df in enumerate([frame_df, other_df]):
        for label in [0,1]:
            label_peps = st_df[st_df['label'] == label]
            label_peps['group'] = groups[label]

            fig.add_trace(
                go.Violin(
                    x=label_peps['group'], y=label_peps['relativePosition'], meanline_visible=True,
                    fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, points=False, spanmode = 'hard'
                ),
                row=1, col=1+g_idx
            )
        if g_idx == 0:
            m_c = STRATUM_COLOUR_SCHEME['CDS_frameshift']
        else:
            m_c = 'green'
        fig.add_trace(
            go.Scatter(
                x=st_df['relativePosition'], y=st_df['relativePosition_shap'], mode='markers',
                marker_color=m_c
                #meanline_visible=True,
                # fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, boxpoints=False
            ),
            row=2, col=1+g_idx
        )

    pos_peps = unique_pep_df[unique_pep_df['label'] == label]
    for idx, stratum in enumerate(CRYPTIC_STRATA):
        label_peps = pos_peps[pos_peps['stratum'] == idx]

        fig.add_trace(
            go.Violin(
                x=[stratum]*len(label_peps), y=label_peps['relativePosition'], meanline_visible=True,
                fillcolor=STRATUM_COLOUR_SCHEME[stratum], line_color='black', line_width=0.5, points=False, spanmode = 'hard'
            ),
            row=1, col=3
        )
    fig.update_yaxes(range=[0,1])
    fig.update_layout(
        {'yaxis1': {'title_text': 'Position in protein'}},
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        {
            'yaxis4': {'range': [-2, 2], 'dtick': 1, 'title_text': 'impact on model'},
            'xaxis4': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            'xaxis5': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            'yaxis5': {'range': [-2, 2], 'dtick': 1,},
            'yaxis1': {'title_text': 'relative position'},
        },
        width=800,height=600, bargap=0,
    )
    fig.add_hline(y=0, line_width=0.5, row=2, col=1)
    fig.add_hline(y=0, line_width=0.5, row=2, col=2)

def plot_codons(unique_pep_df):
    # unique_pep_df = unique_pep_df[unique_pep_df['peptide'].apply(lambda x : x[-1] == 'R')]
    fig = make_subplots(rows=3, cols=1)
    for pos_idx, position in enumerate(['end', 'post', 'combined']):
        if position == 'C-term':
            arg_df = unique_pep_df[unique_pep_df['peptide'].apply(lambda x : x[-1] == 'R')]
        elif position == 'post':
            arg_df = unique_pep_df[unique_pep_df['C_term_downstream'] == 'R']
        else:
            arg_df = unique_pep_df[unique_pep_df['C_term_downstream'] == 'R']
            arg_df = arg_df[arg_df['peptide'].apply(lambda x : x[-1] == 'R')]
            arg_df['combinedCodon'] = arg_df['endCodon'] +arg_df['postCodon']
        print(arg_df.groupby('label')['peptide'].count())

        for label in [0,1]:
            lab_df = arg_df[arg_df['label'] == label]
            x = lab_df.groupby(f'{position}Codon', as_index=False)['peptide'].count()
            x['peptide'] /= lab_df.shape[0]
            fig.add_trace(
                go.Bar(x=x[f'{position}Codon'], y=x['peptide'],
                    marker_color=BGD_COLOURS[label], marker_line_color='black',
                ), row=pos_idx+1, col=1,
            )
    fig = clean_plotly_fig(fig)
    fig.update_layout(height=600, bargap=0, width=800)

def plot_groups(unique_pep_df):
    # unique_pep_df = unique_pep_df[unique_pep_df['peptide'].apply(lambda x : x[-1] == 'R')]
    unique_pep_df['endPost'] = unique_pep_df['peptide'].apply(lambda x : x[-1]) + unique_pep_df['C_term_downstream']
    fig = go.Figure()
    # print(arg_df.groupby('label')['peptide'].count())
    alls = []
    for label in [0,1]:
        lab_df = unique_pep_df[unique_pep_df['label'] == label]
        x = lab_df.groupby(f'endPost', as_index=False)['peptide'].count()
        x['peptide'] /= lab_df.shape[0]
        x= x.sort_values(by='endPost')
        # x = x[x['peptide'] > 0.001]
        alls.append(x)
    x = alls[1]
    x['peptide']/=alls[0]['peptide']
    print(alls)
    fig.add_trace(
        go.Bar(x=x[f'endPost'], y=x['peptide'],
            marker_color='cyan', marker_line_color='black',
        ),# row=pos_idx+1, col=1,
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(height=600, bargap=0, width=800)


def analyse_model(config):
    unique_pep_df = pd.read_csv(
        f'{config.output_folder}/unique_peps_scored.csv'
    )
    plot_codons(unique_pep_df)
    plot_groups(unique_pep_df)
    plot_feature_importances(config)

    plot_stop_dist(unique_pep_df, config)
    plot_start_dist(unique_pep_df, config)
    plot_fragment_length(unique_pep_df, config)
    plot_transcriptomics(unique_pep_df, config)

    plot_shap_termini(unique_pep_df, config)
    get_entropy_plots(unique_pep_df, config)

    model_stats = plot_aucs(unique_pep_df, config)

    if config.model == 'cryptic':
        plot_rna_fracs(unique_pep_df, config)
        plot_rel_pos(unique_pep_df)
        plot_kozak(unique_pep_df, config)
        plot_start_codons(unique_pep_df, config)
        plot_shap_stratum(unique_pep_df, config)
        # coverage_plot(unique_pep_df, config)
        create_cryptic_report(config, model_stats)
    else:
        plot_hydrophobicity(unique_pep_df, config)
        plot_rna_fracs(unique_pep_df, config, prefix='UTR_')
        plot_utr(unique_pep_df)
        # plot_ubi_distro(unique_pep_df)
        # plot_go_distro(unique_pep_df)

        create_canonical_report(config, model_stats)

def plot_utr(pep_df):
    fig = go.Figure()
    traces = plot_distros(pep_df, 'UTR_length', 'violin')
    fig.add_traces(traces)
    fig = clean_plotly_fig(fig)

def coverage_plot(unique_pep_df, config):

    fig = make_subplots(rows=2, cols=1)
    fig = clean_plotly_fig(fig)
    for stratum in CRYPTIC_STRATA:
        strat_df = unique_pep_df[unique_pep_df['stratum'] == CRYPTIC_STRATA.index(stratum)]
        for label in [1, 0]:
            lab_df = strat_df[strat_df['label'] == label]
            print(lab_df['relativePosition'].max())
            lab_df['relPepEnd'] = lab_df['relativePosition'] + (lab_df['peptide'].apply(len)/lab_df['protLength'])
            print(lab_df['relPepEnd'].max())
            distros = [0]*10_000
            for _, df_row in lab_df.iterrows():
                for idx in range(round(df_row['relativePosition']*10_000), round(df_row['relPepEnd']*10_000)):
                    distros[idx] += 1
            distros = [x/max(distros) for x in distros]
            fig.add_trace(
                go.Scatter(
                    x=[x/10_000 for x in range(10_000)], y=distros, marker_color=STRATUM_COLOUR_SCHEME[stratum],
                    opacity=0.8, marker_line_color=BGD_COLOURS[label], marker_line_width=0.5, mode='lines',
                ), row=2-label, col=1,
            )


def plot_shap_termini(unique_pep_df, config):
    """
    """
    positions = ['C-term', 'post C-term']

    for pos_name, pos_feats in zip(positions, [C_TERM_FEATS, C_TERM_POST_FEATS]):
        unique_pep_df[pos_name] = unique_pep_df.apply(
            lambda x : [m.split('_')[-1] for m in pos_feats if x[m] == 1][0],
            axis=1,
        )
        unique_pep_df[f'{pos_name} shap'] = unique_pep_df[
            [f'{x}_shap' for x in pos_feats]
        ].sum(axis=1)

    fig = make_subplots(
        rows=2, cols=2, subplot_titles=positions, vertical_spacing=0.1,
    )
    pos_counts = {}
    neg_counts = {}
    for idx, position in enumerate(positions):
        pos_counts[position] = []
        neg_counts[position] = []
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
                row=2, col=1+idx
            )
            pos_counts[position].append(mini_col[mini_col['label'] == 1].shape[0])
            neg_counts[position].append(mini_col[mini_col['label'] == 0].shape[0])
    for idx, position in enumerate(positions):
        pos_sum = sum(pos_counts[position])
        pos_fracs = [100*x/pos_sum for x in pos_counts[position]]
        neg_sum = sum(neg_counts[position])
        neg_fracs = [100*x/neg_sum for x in neg_counts[position]]

        for group_idx, (group, colour) in enumerate(GROUP_COLOUR_SCHEME.items()):
            fig.add_trace(
                go.Bar(
                    x=[group], y=[pos_fracs[group_idx] - neg_fracs[group_idx]],
                    marker_color=colour, opacity=0.8,
                    marker_line_color='black',
                    marker_line_width=0.5,
                ),
                row=1, col=1+idx
            )
    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis3': {'range': [-2, 1], 'dtick': 1, 'title_text': 'impact on model'},
            'yaxis4': {'range': [-3, 3], 'dtick': 1},
            'yaxis1': {'range': [-40, 40], 'dtick': 20, 'title_text': 'change in percentage'},
            'yaxis2': {'range': [-40, 40], 'dtick': 20},
        },
        width=800,height=600, bargap=0,
    )
    fig.update_xaxes(linecolor='white')

    pio.write_image(fig, f'{config.output_folder}/imgs/C-term-shap.svg')

def plot_shap_stratum(unique_pep_df, config):
    """
    """
    unique_pep_df['stratum'] = unique_pep_df['stratum'].apply(lambda x : CRYPTIC_STRATA[x])

    fig = make_subplots(rows=2, cols=1)
    pos_counts = []
    neg_counts = []
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
            ), row=2, col=1
        )
        pos_counts.append(mini_col[mini_col['label'] == 1].shape[0])
        neg_counts.append(mini_col[mini_col['label'] == 0].shape[0])

    pos_sum = sum(pos_counts)
    pos_fracs = [100*x/pos_sum for x in pos_counts]
    neg_sum = sum(neg_counts)
    neg_fracs = [100*x/neg_sum for x in neg_counts]

    for strat_idx, stratum in enumerate(CRYPTIC_STRATA):
        fig.add_trace(
            go.Bar(
                x=[stratum], y=[pos_fracs[strat_idx] - neg_fracs[strat_idx]],
                marker_color=STRATUM_COLOUR_SCHEME[stratum], opacity=0.8,
                marker_line_color='black',
                marker_line_width=0.5,
            ),
            row=1, col=1
        )

    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [-60, 60], 'title_text': 'percentage difference'},
            'yaxis2': {'range': [-2, 6], 'title_text': 'impact on model'},
        },
        width=800,height=600,bargap=0,
    )
    fig.update_xaxes(linecolor='white')
    pio.write_image(fig, f'{config.output_folder}/imgs/stratum_shap.svg')


def get_entropy_plots(unique_pep_df, config):
    unique_pep_df = unique_pep_df.drop_duplicates(subset='il_peptide')
    unique_pep_df = unique_pep_df[['il_peptide', 'label', 'C_term_downstream', 'N_term_upstream', 'C_term_downstream_2', 'N_term_upstream_2']].drop_duplicates('il_peptide')
    unique_pep_df['C_term_surround'] = unique_pep_df.apply(
        lambda x : x['il_peptide'][-3:] + x['C_term_downstream'] + x['C_term_downstream_2'],
        axis=1
    )
    unique_pep_df['N_term_surround'] = unique_pep_df.apply(
        lambda x : x['N_term_upstream_2'] + x['N_term_upstream'] + x['il_peptide'][:3],
        axis=1
    )
    unique_pep_df = unique_pep_df.drop('il_peptide', axis=1).rename(columns={
        'C_term_surround': 'peptide',
    })
    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'C-term: detected', 'background'],
        5, f'{config.output_folder}/imgs',
        file_name=f'c_term_js.svg', amino_acids=AMINO_ACIDS+'X',
        y_lim=0.1, plot_size=5,
        x_ticks=['', '', 'C-term', '', ''], #axis='right'
    )

    unique_pep_df = unique_pep_df.drop('peptide', axis=1).rename(columns={
        'N_term_surround': 'peptide',
    })
    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)

    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'N-term: detected', 'background'],
        5, f'{config.output_folder}/imgs',
        file_name=f'n_term_js.svg', amino_acids=AMINO_ACIDS+'X',
        y_lim=0.1, plot_size=5,
        x_ticks=['', '', 'N-term', '', ''],
    )

def plot_rna_fracs(unique_pep_df, config, prefix=''):
    fig = make_subplots(rows=2, cols=4, subplot_titles=list('ACGU'), vertical_spacing=0.1)
    nucleotide_list = list('ACGU')


    for idx, nuc in enumerate(nucleotide_list):
        cor_res = linregress(
            unique_pep_df[f'{prefix}{nuc}_frac'],
            unique_pep_df[f'{prefix}{nuc}_frac_shap']
        )
        print(nuc, cor_res)
        vals = np.linspace(0, 0.6)
        res = vals*cor_res.slope + cor_res.intercept
        plot_df = unique_pep_df.sample(n=2000, random_state=42)
        fig.add_trace(
            go.Scatter(
                x=plot_df[f'{prefix}{nuc}_frac'], y=plot_df[f'{prefix}{nuc}_frac_shap'],
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
        traces = plot_distros(unique_pep_df, f'{prefix}{nucleotide}_frac', 'violin')
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

    pio.write_image(fig, f'{config.output_folder}/imgs/rna_frac_distro.svg')

def plot_start_codons(unique_pep_df, config):
    fig = make_subplots(rows=2, cols=len(START_CODONS), subplot_titles=START_CODONS,vertical_spacing=0.1)
    group_labels = ['background', 'detected']
    for idx, codon in enumerate(START_CODONS):
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

    colors = ['magenta', 'cyan']
    for idx, codon in enumerate(START_CODONS):
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
        width=1000, height=600, bargap=0,
    )
    for idx in range(1,6):
        fig.update_layout({f'yaxis{idx}': {'range': [0, 100]}})
    for idx in range(6,11):
        fig.update_layout({f'yaxis{idx}': {'range': [-2, 1]}, f'xaxis{idx}': {'title_text': 'codon present', 'linecolor': 'white'}})
    fig.update_layout({f'yaxis6': {'title_text': 'impact on model'}})
    fig.update_layout({f'yaxis1': {'title_text': 'percentage'}})

    # fig.update_yaxes(
    #     range=[0,100]
    # )


        # fig = make_subplots(rows=1, cols=4, subplot_titles=fraction)
        # pos_peps['group'] = 'detected'
        # neg_peps['group'] = 'background'
        # for idx, nucleotide in enumerate(nucleotide_list):
        #     fig.add_trace(
        #         go.Violin(
        #             x=pos_peps['group'], y=pos_peps[f'{nucleotide}_frac'], meanline_visible=True,
        #             fillcolor='#ADD8E6', line_color='black', line_width=0.5, points=False),
        #         row=1, col=1+idx, 
        #     )
        #     fig.add_trace(
        #         go.Violin(x=neg_peps['group'], y=neg_peps[f'{nucleotide}_frac'], fillcolor='wheat', line_color='black',
        #                 line_width=0.5, points=False, meanline_visible=True),
        #         row=1, col=1+idx
        #     )
        # fig = clean_plotly_fig(fig)
        # fig.update_layout(
        #     width=800, height=400,
        # )
        # fig.update_yaxes(
        #     range=[0,0.5]
        # )
        # pio.write_image(fig, f'{config.output_folder}/imgs/rna_frac_distro.svg')
