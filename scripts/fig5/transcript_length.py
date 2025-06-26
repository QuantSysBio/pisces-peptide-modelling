
from pisces.plot_utils import clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from ppm.analysis_utils import plot_distros, plot_shap
from plotly.subplots import make_subplots

PROJECTS = ['output/crypticK562_241129', 'output/cryptic721_241205']

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
    pep_df['transcriptLength'] = pep_df['start_dist'] + pep_df['peptide'].apply(len) + pep_df['stopDistances']

    print(pep_df['AUG_upstream'].mean())
    print(pep_df[['AUG_upstream', 'start_dist']])
    # titles = ['all', '']
    
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.15)

    traces = plot_distros(pep_df, 'start_dist', 'violin')
    for trace in traces:
        fig.add_trace(trace, row=1, col=1)

    traces = plot_distros(pep_df, 'stopDistances', 'violin')
    for trace in traces:
        fig.add_trace(trace, row=1, col=2)


    # traces = plot_shap(pep_df, 'stopDistances', 'boxdot')
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
        'yaxis1': dict(range=[0, 60], title_text='distance to start codon'),
        'yaxis2': dict(range=[0, 60], title_text='distance to stop codon'),
        # 'yaxis2': dict(range=[-1, 1], title_text='impact on model'),
        # 'xaxis1': dict(range=[0, 20], title_text='distance to start codon'),
        # 'xaxis2': dict(range=[0, 20], title_text='distance to start codon', linecolor='white')
    }, width=500, height=300,)
    fig.write_image('fig5/img/Fig_transcript_length.svg')
    fig.show()


if __name__ == '__main__':
    plot_start_dist()