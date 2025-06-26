import numpy as np
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import plotly.figure_factory as ff

PROJECT = 'splicedK562_241010'


COLOURS = ['#00BFBF', 'orange']
def plot_aug():
    """
    """
    unique_pep_df = pd.read_csv(f'{PROJECT}/unique_peps_scored.csv')
    print([x for x in unique_pep_df.columns if 'int' in x])
    unique_pep_df = unique_pep_df[unique_pep_df['interveningSeqLengths'] < 100]
    fig = go.Figure(data=go.Scatter(

        x=unique_pep_df['interveningSeqLengths'],
        y=unique_pep_df['interveningSeqLengths_shap'],
        mode='markers', marker_color='#9BBFE5', marker_line={'width': 0.5, 'color': 'black'},
        opacity=0.8,
    ))
    # fig = ff.create_2d_density(
    #     x=unique_pep_df['interveningSeqLengths'].apply(np.log10),
    #     y=unique_pep_df['interveningSeqLengths_shap'],
    #     hist_color='rgb(255, 237, 222)', point_size=0, ncontours=20,
    # )


    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [-2, 4]},
            'xaxis1': {'range': [0, 100]},
        },
        width=300,height=300,
    )
    fig.update_xaxes(title='intervening seq length', linecolor='white', tickmode='array')
    fig.update_yaxes(title='impact on model')
    pio.write_image(fig, 'fig6/img/Fig6_B.svg')
    fig.show()

if __name__ == '__main__':
    plot_aug()