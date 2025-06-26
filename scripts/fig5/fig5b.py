import pandas as pd
from pisces.plot_utils import clean_plotly_fig, create_comparison_logo_plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from sklearn.metrics import (
    precision_recall_curve, roc_curve, log_loss, precision_score, recall_score, auc,
)

STRATA = [
    'canonical',
    'cryptic',
    'spliced',
]
CRYPTIC_STRATA = [
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
]
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
COLOUR_DICT = {
    'canonical': '#EC9A56',
    'spliced': '#9BBFE5',
    'cryptic': '#BA69BE',
    'fiveUTR': 'orange',
    'threeUTR': 'goldenrod',
    'CDS_frameshift': 'purple',
    'lncRNA': 'darkgrey',
    'intronic': 'forestgreen',
    'intergenic': 'navy',
}

def get_curves(unique_pep_df, stratum, dash):
    fprs, tprs, _ = roc_curve(unique_pep_df['label'], unique_pep_df['prediction_xgb2'])
    precisions, recalls, _ = precision_recall_curve(unique_pep_df['label'], unique_pep_df['prediction_xgb2'])
    roc_cur = go.Scatter(x=fprs, y=tprs, mode='lines', line_color=COLOUR_DICT[stratum], name=stratum, line_dash=dash)
    pr_curve = go.Scatter(x=recalls, y=precisions, mode='lines', line_color=COLOUR_DICT[stratum], name=stratum, line_dash=dash)
    return roc_cur, pr_curve

def clean_curve(fig):
    fig.add_trace(
        go.Scatter(x=[a*0.1 for a in range(11)], y=[a*0.1 for a in range(11)], mode='lines', line={'color': 'black', 'dash': 'dot'}, name='baseline'),
    )
    fig.update_layout(width=400, showlegend=True)
    fig.update_xaxes(range=[0,1], dtick=0.2)
    fig.update_yaxes(range=[0,1], dtick=0.2)
    return fig

def plot_aucs():
    fig = go.Figure()
    fig_pr = go.Figure()
    for value in PROJECTS:
        for (stratum, project) in zip(STRATA, PROJECTS[value]):
            unique_pep_df = pd.read_csv(f'{project}/unique_peps_scored.csv')
            if value == 'K562':
                dash = 'solid'
            else:
                dash = 'dash'
            roc_cur, pr_curve = get_curves(unique_pep_df, stratum, dash)
            print(stratum, value, auc(roc_cur.x, roc_cur.y), auc(pr_curve.x, pr_curve.y))
            fig.add_trace(roc_cur)
            fig_pr.add_trace(pr_curve)
  
    fig = clean_plotly_fig(fig)
    fig = clean_curve(fig)
    fig.show()
    fig.write_image('manuscript_figs/fig5/fig5b.svg')

    for value in PROJECTS:
        fig2 = go.Figure()
        cryptic_df = pd.read_csv(f'{PROJECTS[value][1]}/unique_peps_scored.csv')
        for idx, stratum in enumerate(CRYPTIC_STRATA):
            # if value == 'K562':
            #     dash = 'solid'
            # else:
            #     dash = 'dash'
            sub_df = cryptic_df[cryptic_df['stratum'] == idx]
            roc_cur, _ = get_curves(sub_df, stratum, 'solid')
            fig2.add_trace(roc_cur)

        fig2 = clean_plotly_fig(fig2)
        fig2 = clean_curve(fig2)
        # fig2.show()
        if value == 'K562':
            fig2.write_image('manuscript_figs/figS/rocCurves/figROC_a.svg')
        else:
            fig2.write_image('manuscript_figs/figS/rocCurves/figROC_b.svg')

    # fig_pr = clean_plotly_fig(fig_pr)
    # fig_pr.update_layout(width=300)
    # fig_pr.update_xaxes(range=[0,1], dtick=0.2)
    # fig_pr.update_yaxes(range=[0,1], dtick=0.2)
    # fig_pr.show()
    # pio.write_image(fig, 'fig6/img/FigS_model_A.svg')
    # return {
    #     'log_loss': log_loss(unique_pep_df['label'], unique_pep_df['prediction_xgb2']),
    #     'auc-pr': auc(recalls, precisions),
    #     'auc-roc': auc(fprs, tprs),
    #     'precision': precision_score(unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)),
    #     'recall': recall_score(unique_pep_df['label'], unique_pep_df['prediction_xgb2'].apply(round)),
    # }

if __name__ == '__main__':
    plot_aucs()
