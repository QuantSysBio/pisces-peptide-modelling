
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
from ppm.analysis_utils import plot_distros, plot_shap
from ppm.constants import (
    BGD_COLOURS, CRYPTIC_STRATA, N_CV_GROUPS, COLOUR_DICT, STRATUM_COLOUR_SCHEME
)
from scipy.stats import mannwhitneyu
from scipy.stats import linregress, pearsonr


PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']
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

def plot_kozak():
    """
    """
    pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])

    # print(pep_df['AUG_upstream'].mean())
    # print(pep_df[['AUG_upstream', 'start_dist']])
    # titles = ['all', '']
    
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.15)

    traces = plot_distros(pep_df, 'kozakScore', 'boxdot')
    for trace in traces:
        fig.add_trace(trace, row=1, col=1)
    pos_peps = pep_df[pep_df['label'] ==1].dropna(subset=['kozakScore'])
    traces = plot_shap(pos_peps, 'kozakScore', 'scatter', display_range=[0,1])
    for trace in traces:
        trace.update(marker_color=COLOUR_DICT['cryptic'])
        fig.add_trace(trace, row=1, col=2)


    cor_res = linregress(
        pos_peps['kozakScore'],
        pos_peps['kozakScore_shap'],
    )
    print('linregrss', cor_res)
    vals = np.linspace(0.2, 1.0)
    res = vals*cor_res.slope + cor_res.intercept
    fig.add_trace(
        go.Scatter(
            x=vals, y=res, mode='lines', line_color='black', line_width=0.5, name='linregrss',
            # marker_line_color='black', marker_line_width=0.5,
            opacity=0.8,
        ),
        row=1, col=2,
    )


    print(mannwhitneyu(pep_df[pep_df['label'] == 1]['kozakScore'].dropna(), pep_df[pep_df['label'] == 0]['kozakScore'].dropna()))

    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5, row=1, col=2)
    fig.update_layout({
        'yaxis1': dict(range=[0, 1], title_text='Kozak similarity score'),
        'yaxis2': dict(range=[-1.5, 1.5], title_text='impact on model'),
        # 'xaxis1': dict(range=[0, 50], title_text='distance to start codon'),
        'xaxis2': dict(range=[0, 1], title_text='Kozak similarity score')
    }, width=500, height=300,)
    fig.write_image('manuscript_figs/figS/kozak/figS_kozak.svg')
    fig.show()
if __name__ == '__main__':
    plot_kozak()