
from pisces.plot_utils import create_comparison_logo_plot
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']

from ppm.constants import AMINO_ACIDS
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
GROUP_COLOUR_SCHEME = {
    'acidic': 'firebrick',
    'basic': 'seagreen',
    'hydrophobic': 'orange',
    'polar': 'dodgerblue',
    'special-case': 'deeppink',
    'end': 'black',
}
C_TERM_FEATS = [
    'C_term_basic', 'C_term_special-case', 'C_term_hydrophobic',
    'C_term_acidic', 'C_term_polar', 
]
C_TERM_POST_FEATS = [
    'C_term_neg_1_basic', 'C_term_neg_1_special-case', 'C_term_neg_1_hydrophobic',
    'C_term_neg_1_acidic', 'C_term_neg_1_polar', 'C_term_neg_1_end',
]
COLOURS = ['#00BFBF', 'orange']

def plot_entro():
    unique_pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
    unique_pep_df = unique_pep_df.drop_duplicates(subset='il_peptide')
    unique_pep_df = unique_pep_df[['il_peptide', 'label', 'C_term_downstream', 'N_term_upstream', 'C_term_downstream_2', 'N_term_upstream_2']].drop_duplicates('il_peptide')
    unique_pep_df['peptide'] = unique_pep_df.apply(
        lambda x : x['N_term_upstream_2'] + x['N_term_upstream'] + x['il_peptide'][:1] + 'AAA' + x['il_peptide'][-1:] + x['C_term_downstream'] + x['C_term_downstream_2'],
        axis=1
    )

    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'N-term: detected', 'background'],
        9, 'manuscript_figs/fig5',
        file_name=f'fig5h.svg', amino_acids=AMINO_ACIDS+'X',
        y_lim=0.1, plot_size=5,
        x_ticks=['', '', 'N-term', '', '...', '', 'C-term', '', ''],
    )

if __name__ == '__main__':
    plot_entro()