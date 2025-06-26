
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
PROJECT = 'crypticK562_241011'

STRATUM_COLOUR_SCHEME = {
    'fiveUTR': 'orange',
    'threeUTR': 'yellow',
    'CDS_frameshift': 'purple',
    'lncRNA': 'darkgrey',
    'intronic': 'forestgreen',
    'intergenic': 'navy',
}
CRYPTIC_STRATA = [
    'fiveUTR',
    'threeUTR',
    'CDS_frameshift',
    'lncRNA',
    'intronic',
    'intergenic',
]
def plot_cryptic():
    """
    """
    unique_pep_df = pd.read_csv(f'{PROJECT}/unique_peps_scored.csv')
    unique_pep_df['stratum'] = unique_pep_df['stratum'].apply(lambda x : CRYPTIC_STRATA[x])
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
    pio.write_image(fig, 'fig6/img/stratum_shap.svg')
    fig.show()
    fig = go.Figure()
    for stratum, colour in STRATUM_COLOUR_SCHEME.items():
        mini_col = unique_pep_df[unique_pep_df['stratum'] == stratum]
        print(stratum, mini_col.shape[0])
        fig.add_trace(
            go.Bar(
                x=[stratum],
                y=[100*mini_col[mini_col['label']==1].shape[0]/mini_col.shape[0]],
                marker_color=colour, opacity=0.9,
                marker_line_color='black',
                marker_line_width=0.5,
            ),
        )

    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5)
    fig.update_layout(
        {
            'yaxis1': {'range': [0, 100]},
        },
        width=800,height=300,
    )
    # fig.update_xaxes(linecolor='white')
    fig.update_yaxes(title='fraction detected')
    pio.write_image(fig, 'fig6/img/stratum_detected.svg')
    fig.show()

if __name__ == '__main__':
    plot_cryptic()