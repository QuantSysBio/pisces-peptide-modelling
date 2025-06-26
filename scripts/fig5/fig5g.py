
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from ppm.analysis_utils import plot_distros, plot_shap
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu

PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']

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
COLOURS = ['#00BFBF', 'orange']
def plot_start_dist():
    """
    """
    pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])

    print(pep_df.columns)
    print(pep_df[['C_term_downstream', 'C_term_downstream_2', 'stopDistances']])
    print(pep_df[pep_df['stopDistances'] == 1][['C_term_downstream', 'C_term_downstream_2', 'stopDistances']])
    # titles = ['all', '']
    print(mannwhitneyu(pep_df[pep_df['label'] == 1]['stopDistances'].dropna(), pep_df[pep_df['label'] == 0]['stopDistances'].dropna()))

    fig = go.Figure()

    # traces = plot_distros(pep_df, 'stopDistances', 'scatter', display_range=[0, 21])
    # for trace in traces:
    #     fig.add_trace(trace)
    traces = plot_distros(pep_df, 'stopDistances', 'boxdot')
    for trace in traces:
        fig.add_trace(trace)

    # traces = plot_shap(pep_df, 'stopDistances', 'line', display_range=[0,51])
    # for trace in traces:
    #     fig.add_trace(trace, row=1, col=2)

    # if config.model == 'cryptic':
    #     for strat_idx, stratum in enumerate(CRYPTIC_STRATA):
    #         strat_df = pep_df[pep_df['stratum'] == strat_idx]
    #         traces = plot_distros(strat_df, 'start_dist', 'scatter', display_range=[0, 10])
    #         for trace in traces:
    #             fig.add_trace(trace, row=2+strat_idx, col=1)

    #         traces = plot_shap(strat_df, 'start_dist', 'line', display_range=[0,10])
    #         for trace in traces:
    #             fig.add_trace(trace, row=2+strat_idx, col=2)


    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5, row=1, col=2)
    fig.update_layout({
        'yaxis1': dict(range=[0, 100], title_text='distance to stop codon'),
        # 'yaxis2': dict(range=[-1, 1], title_text='impact on model'),
        # 'xaxis1': dict(range=[0, 20], title_text='distance to stop codon'),
        # 'xaxis2': dict(range=[0, 20], title_text='distance to start codon', linecolor='white')
    }, width=100, height=300,)
    fig.write_image('manuscript_figs/fig5/fig5g.svg')
    fig.show()


if __name__ == '__main__':
    plot_start_dist()