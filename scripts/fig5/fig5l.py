
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from ppm.constants import BGD_COLOURS
import numpy as np
from ppm.analysis_utils import plot_distros, plot_shap
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu
from scipy.stats import linregress, pearsonr
import warnings
warnings.simplefilter(action='ignore')
PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']
CAN_PROJECTS = ['output/final/canonicalK562', 'output/final/canonical721']

STRATUM_COLOUR_SCHEME = {
    'fiveUTR': '#FFA756',
    'threeUTR': '#F2EEB8',
    'CDS_frameshift': '#B790D4',
    'lncRNA': 'lightgrey',
    'intronic': '#4ece78',
    'intergenic': '#528AAE',
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

def plot_rel_pos():
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
    can_unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in CAN_PROJECTS])
    label_can_peps = can_unique_pep_df[can_unique_pep_df['label'] == 1]
    groups = ['background', 'detected']
    fig = make_subplots(rows=2, cols=2)
    fig_5i = go.Figure()
    frame_df = unique_pep_df[unique_pep_df['stratum'] == 2]
    other_df = unique_pep_df[unique_pep_df['stratum'] != 2]


    for g_idx, st_df in enumerate([frame_df, other_df]):

        print(g_idx, mannwhitneyu(
            st_df[st_df['label'] == 0]['relativePosition'].dropna(),
            st_df[st_df['label'] == 1]['relativePosition'].dropna(),
        ))
        for label in [0,1]:
            label_peps = st_df[st_df['label'] == label]
            label_peps['group'] = groups[label]

            fig.add_trace(
                go.Box(
                    x=label_peps['group'], y=label_peps['relativePosition'], #meanline_visible=True,
                    fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, 
                    boxpoints='outliers',
                    notched=False, #spanmode = 'hard'
                ),
                row=1, col=1+g_idx
            )
        if g_idx == 0:
            m_c = STRATUM_COLOUR_SCHEME['CDS_frameshift']
        else:
            m_c = 'green'
        st_df = st_df[st_df['label'] == 1]

        cor_res = linregress(
            st_df[f'relativePosition'],
            st_df[f'relativePosition_shap']
        )
        print(g_idx, cor_res)
        vals = np.linspace(0, 1.0)
        res = vals*cor_res.slope + cor_res.intercept
        fig.add_trace(
            go.Scatter(
                x=st_df['relativePosition'], y=st_df['relativePosition_shap'], mode='markers',
                marker_color=m_c
                #meanline_visible=True,
                # fillcolor=BGD_COLOURS[label], line_color='black', line_width=0.5, boxpoints=False
            ),
            row=2, col=1+g_idx
        )
        fig.add_trace(
            go.Scatter(
                x=vals, y=res, mode='lines', line_color='black', line_width=0.5, name=str(g_idx),
                # marker_line_color='black', marker_line_width=0.5,
                opacity=0.8,
            ),
            row=2, col=1+g_idx
        )

    pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
    for idx, stratum in enumerate(CRYPTIC_STRATA):
        label_peps = pos_peps[pos_peps['stratum'] == idx]

        fig_5i.add_trace(
            go.Box(
                x=[stratum]*len(label_peps), y=label_peps['relativePosition'],# meanline_visible=True,
                fillcolor=STRATUM_COLOUR_SCHEME[stratum], line_color='black', line_width=0.5, 
                    boxpoints='outliers',
                    notched=False,# spanmode = 'hard',
                # opacity=0.8,
            ),
        )
        print(stratum, mannwhitneyu(
            frame_df[frame_df['label'] == 1]['relativePosition'].dropna(),
            label_peps['relativePosition'].dropna(),
        ))
    fig_5i.add_trace(
        go.Box(
            x=['canonical']*len(label_can_peps), y=label_can_peps['relativePosition'],# meanline_visible=True,
            fillcolor='#EC9A56', line_color='black', line_width=0.5, boxpoints=False,# spanmode = 'hard',
            # opacity=0.8,
        ),
    )
    print('canonical', mannwhitneyu(
        frame_df[frame_df['label'] == 1]['relativePosition'].dropna(),
        label_can_peps['relativePosition'].dropna(),
    ))
    fig.update_yaxes(range=[0,1])
    fig.update_layout(
        {'yaxis1': {'title_text': 'Position in protein'}},
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        {
            'yaxis4': {'range': [-2, 2], 'dtick': 1, 'title_text': 'impact on model'},
            'xaxis3': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            'xaxis4': {'range': [0,1],'title_text': 'relaitve position', 'linecolor': 'white'},
            'yaxis3': {'range': [-2,2],'title_text': 'impact on model',},
            'yaxis4': {'range': [-2,2],'title_text': 'impact on model',},
            'yaxis5': {'range': [-2, 2], 'dtick': 1,},
            'yaxis1': {'title_text': 'relative position'},
        },
        width=700,height=600, bargap=0,
    )
    fig.add_hline(y=0, line_width=0.5, row=2, col=1)
    fig.add_hline(y=0, line_width=0.5, row=2, col=2)
    fig.show()
    fig.write_image('manuscript_figs/figS/relativePosition/fig_relativePosition.svg')

    fig_5i = clean_plotly_fig(fig_5i)
    fig_5i.update_layout(yaxis=dict(range=[0,1]))
    fig_5i.show()
    fig_5i.write_image('manuscript_figs/fig5/fig5l.svg')
# 
if __name__ == '__main__':
    plot_rel_pos()