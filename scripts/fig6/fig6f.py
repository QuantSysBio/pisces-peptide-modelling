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
    print(pep_df.columns.tolist())
    pep_df = pep_df[pep_df['nCanonicalPeptides'] > 0]
    
    fig1 = go.Figure()
        
    traces = plot_distros(pep_df, 'sr1_can_dist', 'boxdot')
    for trace in traces:
        fig1.add_trace(trace)

    # traces = plot_distros(pep_df, 'sr1_can_dist', 'scatter', display_range=[0, 101])
    # for trace in traces:
    #     fig1.add_trace(trace)

    fig2 = go.Figure()
        
    traces = plot_distros(pep_df, 'sr2_can_dist', 'boxdot')
    for trace in traces:
        fig2.add_trace(trace,)

    # traces = plot_distros(pep_df, 'sr2_can_dist', 'scatter', display_range=[0,101])
    # for trace in traces:
    #     fig2.add_trace(trace)
    # fig = ff.create_2d_density(
    #     x=unique_pep_df['interveningSeqLengths'].apply(np.log10),
    #     y=unique_pep_df['interveningSeqLengths_shap'],
    #     hist_color='rgb(255, 237, 222)', point_size=0, ncontours=20,
    # )

    for fig in [fig1, fig2]:
        fig = clean_plotly_fig(fig)
        # fig.add_hline(y=0, line_width=0.5)
        fig.update_layout(
            {
                'yaxis1': {'range': [0, 500]},
                'yaxis2': {'range': [0, 500]},
                # 'yaxis2': {'range': [-1, 2], 'dtick': 1},
            },
            width=200,height=300,
        )
        # fig.update_xaxes(title='intervening seq length', linecolor='white', tickmode='array')
        # fig.update_xaxes(title='distance', range=[0,100])
        fig.show()
    pio.write_image(fig1, 'fig6/img/Fig6_F_sr1.svg')
    pio.write_image(fig2, 'fig6/img/Fig6_F_sr2.svg')
    # fig.show()

if __name__ == '__main__':
    plot_aug()