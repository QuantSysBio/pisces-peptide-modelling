
from pisces.plot_utils import create_comparison_logo_plot
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd

from ppm.constants import AMINO_ACIDS

PROJECTS = ['output/splicedK562_241129', 'output/spliced721_241205']

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
    # unique_pep_df = pd.read_csv('output/splicedK562_241129/unique_peps_scored.csv')

    unique_pep_df['p1'] = unique_pep_df.apply(
        lambda x : [a_a for a_a in AMINO_ACIDS if x[f'p1_{a_a}'] == 1][0],
        axis=1,
    )
    unique_pep_df['p1 prime'] = unique_pep_df.apply(
        lambda x : [a_a for a_a in AMINO_ACIDS if x[f'p1_prime_{a_a}'] == 1][0],
        axis=1,
    )
    # unique_pep_df = unique_pep_df[[
    #     'peptide', 'label',
    #     'p2', 'p1', 'p_minus_1', 'p_minus_2'
    # ]].drop_duplicates('peptide')
    unique_pep_df['peptide'] = unique_pep_df['p2'] + unique_pep_df['p1'] + unique_pep_df['p_minus_1'] + unique_pep_df['p_minus_2']

    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'splice site: detected', 'background'],
        4, f'fig6/img',
        file_name=f'Fig6_C_sr1.svg', amino_acids=AMINO_ACIDS,
        y_lim=0.04, plot_size=5,
        x_ticks=['', 'p1', 'p1"', ''],
        vline=2.5,
    )

    unique_pep_df = unique_pep_df[[
        'peptide', 'label',
        'p2', 'p1', 'p_minus_1', 'p_minus_2',
        'p_minus_2_prime', 'p_minus_1_prime', 'p1 prime', 'p2_prime'
    ]].drop_duplicates('peptide')
    unique_pep_df['peptide'] = unique_pep_df['p_minus_2_prime'] + unique_pep_df['p_minus_1_prime'] + unique_pep_df['p1 prime'] + unique_pep_df['p2_prime']

    pos_pep_df = unique_pep_df[unique_pep_df['label'] == 1].drop('label', axis=1)
    neg_pep_df = unique_pep_df[unique_pep_df['label'] == 0].drop('label', axis=1)


    create_comparison_logo_plot(
        [pos_pep_df, neg_pep_df],
        [f'splice site: detected', 'background'],
        4, f'fig6/img',
        file_name=f'Fig6_C_sr2.svg', amino_acids=AMINO_ACIDS,
        y_lim=0.04, plot_size=5,
        x_ticks=['', '', "p1'", ''],
        vline=2.5,
    )

if __name__ == '__main__':
    plot_entro()