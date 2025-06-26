import numpy as np
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import plotly.figure_factory as ff
from ppm.analysis_utils import plot_distros, plot_shap
from plotly.subplots import make_subplots

PROJECTS = ['output/splicedK562_241129', 'output/spliced721_241206']


COLOURS = ['#00BFBF', 'orange']
def plot_aug():
    """
    """
    pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])

    # titles = ['all', '']
    
    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.15)
    print(pep_df[pep_df['label'] == 1]['interveningSeqLengths'].sort_values().tail(20))
    fig1 = go.Figure()
    traces = plot_distros(pep_df, 'interveningSeqLengths', 'scatter', display_range=[0, 101])
    for trace in traces:
        fig.add_trace(trace, row=1, col=1)

    
    traces = plot_distros(pep_df, 'interveningSeqLengths', 'violin')
    for trace in traces:
        fig.add_trace(trace, row=1, col=2)

    
    traces = plot_distros(pep_df, 'interveningSeqLengths', 'boxdot')
    for trace in traces:
        fig.add_trace(trace, row=1, col=3)

    # fig2 = go.Figure()
    # traces = plot_shap(pep_df, 'interveningSeqLengths', 'line', display_range=[0,101])
    # for trace in traces:
    #     fig2.add_trace(trace, row=1, col=2)
    # fig = ff.create_2d_density(
    #     x=unique_pep_df['interveningSeqLengths'].apply(np.log10),
    #     y=unique_pep_df['interveningSeqLengths_shap'],
    #     hist_color='rgb(255, 237, 222)', point_size=0, ncontours=20,
    # )


    fig = clean_plotly_fig(fig)
    fig.add_hline(y=35, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 1]},
            'yaxis2': {'range': [0, 2_000], 'dtick': 500},
            'yaxis3': {'range': [0, 2_000], 'dtick': 500},
            'xaxis1': dict(title='intervening seq length', range=[0,1000]),
        },
        width=600,height=300,
    )
    # fig.update_xaxes(title='intervening seq length', linecolor='white', tickmode='array')
    fig.update_xaxes
    pio.write_image(fig, 'fig6/img/Fig6_B_alt.svg')
    fig.show()

if __name__ == '__main__':
    plot_aug()