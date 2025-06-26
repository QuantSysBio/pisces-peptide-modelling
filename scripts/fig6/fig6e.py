import numpy as np
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from ppm.constants import BGD_COLOURS
import plotly.figure_factory as ff

PROJECTS = ['output/splicedK562_241129', 'output/spliced721_241205']


COLOURS = ['#00BFBF', 'orange']
def plot_aug():
    """
    """
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])

    # unique_pep_df = unique_pep_df[unique_pep_df['interveningSeqLengths'] < 100]
    fig = go.Figure()
    # data=go.Scatter(

    #     x=unique_pep_df['interveningSeqLengths'],
    #     y=unique_pep_df['interveningSeqLengths_shap'],
    #     mode='markers', marker_color='#9BBFE5', marker_line={'width': 0.5, 'color': 'black'},
    #     opacity=0.8,
    # ))

    # pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    # neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)
    pos_peps = unique_pep_df[unique_pep_df['label'] == 1]
    neg_peps = unique_pep_df[unique_pep_df['label'] == 0]
    pos_peps['group'] = 'detected'
    neg_peps['group'] = 'background'
    for (label, name, df) in zip([1,0], ['detected', 'background'], [pos_peps, neg_peps]):
        fig.add_trace(
            go.Bar(
                x=[name], y=[100*df[df['nCanonicalPeptides'] > 0].shape[0]/df.shape[0]],
                marker_color=BGD_COLOURS[label], opacity=0.8,
                marker_line_color='black',
                marker_line_width=0.5,
                # meanline_visible=True,
                # points=False,
            ),
        )
    # fig = ff.create_2d_density(
    #     x=unique_pep_df['interveningSeqLengths'].apply(np.log10),
    #     y=unique_pep_df['interveningSeqLengths_shap'],
    #     hist_color='rgb(255, 237, 222)', point_size=0, ncontours=20,
    # )


    fig = clean_plotly_fig(fig)
    # fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 100]},
        },
        width=150,height=300,bargap=0,
    )
    # fig.update_xaxes(title='intervening seq length', linecolor='white', tickmode='array')
    fig.update_yaxes(title='antigens with 1+ canonical peptides')
    pio.write_image(fig, 'fig6/img/Fig6_E.svg')
    fig.show()

if __name__ == '__main__':
    plot_aug()