
from functools import reduce
from operator import or_

import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from sklearn.metrics import (
    precision_recall_curve, roc_curve, log_loss, precision_score, recall_score, auc,
)

from ppm.constants import (
    BGD_COLOURS, CRYPTIC_STRATA, N_CV_GROUPS, COLOUR_DICT, STRATUM_COLOUR_SCHEME
)

def plot_ubi_distro(pep_df):
    fig = go.Figure()
    traces = plot_distros(pep_df, 'ubiCounts', 'bar', display_range=[0,50])
    for trace in traces:
        fig.add_trace(trace)
    fig = clean_plotly_fig(fig)
    fig.show()

def plot_go_distro(pep_df):
    fig = make_subplots(cols=3, rows=2)
    for idx, feat in enumerate(
        ['Cytoplasmic translation', 'RNA processing', 'Nucleic acid metabolic process']
    ):
        print(pep_df.groupby('label')[feat].mean())
        print(pep_df.groupby(feat)['label'].mean())
        print(pep_df[feat].isna().sum())
        traces = plot_distros(pep_df, feat, 'bar', display_range=[0,2])
        for trace in traces:
            fig.add_trace(trace, row=1, col=1+idx)

        traces = plot_shap(pep_df, feat, mode='violin', display_range=[0,2])
        for trace in traces:
            fig.add_trace(trace, row=2, col=1+idx)
    fig = clean_plotly_fig(fig)
    fig.show()

def plot_feature_importances(config):
    import_df = pd.read_csv(f'{config.output_folder}/models/importances.csv')
    import_df['meanImportance'] = import_df[
        [f'model_{idx}' for idx in range(N_CV_GROUPS)]
    ].mean(axis=1)
    import_df['stdImportance'] = import_df[
        [f'model_{idx}' for idx in range(N_CV_GROUPS)]
    ].std(axis=1)
    print(import_df)
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=import_df['meanImportance'],
            y=import_df['feature'],
            marker_color='cyan',
            marker_line_color='black',
            marker_line_width=0.5,
            orientation='h',
            error_x={
                'type': 'data',
                'array': import_df['stdImportance'].values,
                'color': 'black',
                'thickness': 1,
            },
            opacity=0.8,
        )
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(height=600, width=500, bargap=0)
    fig.update_xaxes(range=[0,0.3])
    fig.update_yaxes(dtick=1)
    pio.write_image(fig, f'{config.output_folder}/imgs/feature_importance.svg')


def plot_aucs(unique_pep_df, config):
    """ Plot area under receiver operator and PR curves
    """
    fig = make_subplots(rows=1, cols=2, subplot_titles=['ROC', 'PR'])
    fprs, tprs, _ = roc_curve(
        unique_pep_df['label'], unique_pep_df['prediction_xgb2']
    )
    precisions, recalls, _ = precision_recall_curve(
        unique_pep_df['label'], unique_pep_df['prediction_xgb2']
    )

    fig.add_trace(
        go.Scatter(
            x=fprs, y=tprs, name=config.model,
            line_color=COLOUR_DICT[config.model], mode='lines',
        ), row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=recalls, y=precisions, name=config.model,
            line_color=COLOUR_DICT[config.model], mode='lines',
        ), row=1, col=2,
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
        'precision': precision_score(
            unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)
        ),
        'recall': recall_score(
            unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)
        ),
    }

def get_expr_freq(label_peps, fractions):
    return (
        100*label_peps[reduce(
            or_,
            [label_peps[f'tr_TPM_K562_{fraction}'] > 0 for fraction in fractions]
        )].shape[0]/label_peps.shape[0]
    )

def plot_transcriptomics(unique_pep_df, config):
    if 'stratum' in unique_pep_df.columns:
        unique_pep_df = unique_pep_df[unique_pep_df['stratum'] != CRYPTIC_STRATA.index('intergenic')]
    if config.cell_line == 'K562':
        fractions = ['poly', 'S80', 'free', 'bulk']
        fig = make_subplots(rows=1, cols=5, subplot_titles=[
            'any fraction > 0',
            'poly abundance',
            'S80 abundance',
            'free abundance',
            'bulk abundance',
        ])
        group_labels = ['background', 'detected']
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
            fig.add_trace(
                go.Violin(
                    x=pos_tr_df['group'], y=pos_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10),
                    fillcolor='#ADD8E6',
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                ),
                row = 1, col=2+idx,
            )
            fig.add_trace(
                go.Violin(
                    x=neg_tr_df['group'], y=neg_tr_df[f'tr_TPM_K562_{fraction}'].apply(np.log10),
                    fillcolor='wheat',
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
        pio.write_image(fig, f'{config.output_folder}/imgs/transcriptomic_distro.svg')
    else:
        fig = make_subplots(rows=1, cols=2, subplot_titles=[
            'bulk > 0',
            'bulk abundance',
        ])
        print(unique_pep_df.groupby('label')['tr_TPM_721'].mean())
        traces = plot_distros(
            unique_pep_df, 'tr_TPM_721', 'bar', display_range=['x']
        )
        for trace in traces:
            fig.add_trace(trace, row=1, col=1)


        unique_pep_df['tr_TPM_721'] = unique_pep_df['tr_TPM_721'].apply(np.log10)
        traces = plot_distros(
            unique_pep_df, 'tr_TPM_721', 'violin'
        )
        for trace in traces:
            fig.add_trace(trace, row=1, col=2)


    if 'stratum' in unique_pep_df.columns:
        stratum_counts = []
        fig = go.Figure()
        for idx, stratum in enumerate(CRYPTIC_STRATA):
            if stratum == 'intergenic':
                continue
            s_df = unique_pep_df[(unique_pep_df['label'] == 1) & (unique_pep_df['stratum'] == idx)]
            if config.cell_line == 'K562':
                s_c = get_expr_freq(s_df, fractions)
            else:
                s_c = 100*s_df[s_df['tr_TPM_721'] > 0].shape[0]/s_df.shape[0]

            fig.add_trace(
                go.Bar(
                    x=[stratum], y=[s_c],
                    marker_color=STRATUM_COLOUR_SCHEME[stratum], opacity=0.8,
                    marker_line_color='black',
                    marker_line_width=0.5,
                )
            )
    fig = clean_plotly_fig(fig)

    fig.update_layout(
        {
            'yaxis1': {'range': [0, 100], 'title_text': 'percentage expressed'},
        },
        width=300,height=300,bargap=0,
    )
    fig.update_xaxes(linecolor='white')
    fig.show()
    pio.write_image(fig, f'{config.output_folder}/imgs/transcriptomic_distro.svg')


def plot_distros(unique_pep_df, feature_name, mode, display_range=None, color_scheme=BGD_COLOURS, sep_column='label', sep_col_range=[0,1]):
    traces = []
    group_labels = 'background', 'detected'
    for label in sep_col_range:
        label_peps = unique_pep_df[unique_pep_df[sep_column] == label]
        if mode == 'violin':
            traces.append(
                go.Violin(
                    x=[group_labels[label]]*len(label_peps),
                    y=label_peps[feature_name],
                    fillcolor=color_scheme[label],
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                    spanmode='hard',
                )
            )
        elif mode == 'box':
            traces.append(
                go.Box(
                    x=[group_labels[label]]*len(label_peps),
                    y=label_peps[feature_name],
                    fillcolor=color_scheme[label],
                    line_color='black',
                    line_width=0.5,
                    boxpoints=False,
                )
            )
        elif mode == 'boxdot':
            traces.append(
                go.Box(
                    x=[group_labels[label]]*len(label_peps),
                    y=label_peps[feature_name],
                    fillcolor=color_scheme[label],
                    line_color='black',
                    line_width=0.5,
                    boxpoints='outliers',
                )
            )
        elif mode == 'scatter':
            fractions = []
            for distance in range(*display_range):
                pep_count = label_peps[
                    (label_peps[feature_name] <= distance) &
                    (label_peps[feature_name].notna())
                ].shape[0]
                fractions.append(pep_count/label_peps.shape[0])

            traces.append(
                go.Scatter(
                    x=list(range(*display_range)),
                    y=fractions,
                    mode='lines',
                    line_color=color_scheme[label],
                    marker_line_color='black',
                    marker_line_width=0.25,
                    opacity=0.8,
                    name=f'{sep_column}={label}'
                )
            )
        elif mode == 'bar':
            if 'x' not in display_range:
                feat_values = list(range(*display_range))
            else:
                feat_values = display_range
            fractions = []
            for val in feat_values:
                if val == 'x':
                    fractions.append(
                        100*label_peps[label_peps[feature_name] != 0].shape[0]/label_peps.shape[0]
                    )
                else:
                    fractions.append(
                        100*label_peps[label_peps[feature_name] == val].shape[0]/label_peps.shape[0]
                    )
                    
            traces.append(
                go.Bar(
                    x=feat_values,
                    y=fractions,
                    marker_color=color_scheme[label],
                    marker_line_color='black',
                    marker_line_width=0.5,
                    name=f'{sep_column}={label}'
                )
            )


    return traces


def plot_shap(unique_pep_df, feature_name, mode, display_range=None):
    traces = []
    if mode == 'scatter':
        traces.append(go.Scatter(
            x=unique_pep_df[feature_name],
            y=unique_pep_df[f'{feature_name}_shap'],
            marker_color=unique_pep_df['label'].apply(lambda x : BGD_COLOURS[x]),
            mode='markers',
        ))
    elif mode == 'line':
        mean_impacts = []
        for f_val in range(*display_range):
            impact = unique_pep_df[
                unique_pep_df[feature_name] == f_val
            ][f'{feature_name}_shap'].mean()
            mean_impacts.append(impact)
        traces.append(
            go.Scatter(
                x=list(range(*display_range)),
                y=mean_impacts,
                mode='lines',
                line_color='firebrick',
            )
        )
    elif mode == 'violin':
        colors = ['magenta', 'cyan']
        for f_val in range(*display_range):
            sub_df = unique_pep_df[
                unique_pep_df[feature_name] == f_val
            ]
            traces.append(
                go.Violin(
                    x=[f_val]*len(sub_df),
                    y=sub_df[f'{feature_name}_shap'],
                    fillcolor=colors[f_val],
                    line_color='black',
                    line_width=0.5,
                    meanline_visible=True,
                    points=False,
                )
            )


    return traces
