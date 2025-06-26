
from pisces.plot_utils import create_comparison_logo_plot, clean_plotly_fig
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
STRATA = [
    'canonical', 'cryptic', 'spliced'
]
COLOUR_DICT = {
    'canonical': '#EC9A56',
    'spliced': '#9BBFE5',
    'cryptic': '#BA69BE',
}
PROJECTS = [
    [
        'output2/final/canonicalB721',
        'output2/final/crypticB721',
        'output2/final/spliced721',
    ],
    [
        'output2/final/canonical721',
        'output2/final/cryptic721',
        'output2/final/spliced721',
    ]
]
TRYPTIC_RESULTS = [
    # 'projects/pisces_tryptic/K562-tryptic/outputFolder',
    'projects/pisces_tryptic/sarkizova_B721.221/outputFolder',
]
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

    unique_pep_df = pd.read_csv('output2/final/canonical721/unique_peps_scored.csv')
    unique_pep_df = unique_pep_df[unique_pep_df['label'] == 1]

    proteins = unique_pep_df[['proteinID', 'proteomics_B721']].drop_duplicates(subset=['proteinID'])
    print(proteins.shape[0])
    tryptic_df = pd.read_csv(f'{TRYPTIC_RESULTS[0]}/final/canonical.csv')
    tryptic_df = tryptic_df[tryptic_df['piscesDiscoverable'] == 1]

    tryptic_df = tryptic_df[tryptic_df['postErrProb'] < 0.15]
    peptide_df = tryptic_df[['peptide']].drop_duplicates()
    tryptic_det_df = pd.read_csv(f'{TRYPTIC_RESULTS[0]}/details/canonical.csv')
    tryptic_det_df = tryptic_det_df.drop_duplicates(subset=['peptide'])

    total_df = pd.merge(peptide_df, tryptic_det_df, how='inner', on='peptide')
    print(total_df.shape[0])

    total_df['proteinID'] = total_df['canonical_Proteins'].apply(lambda x : x.split(' '))
    tryp_prot_df = total_df[['proteinID', 'peptide']].explode('proteinID').groupby('proteinID', as_index=False)['peptide'].apply(list)

    combo_prot_df = pd.merge(proteins, tryp_prot_df, how='left', on='proteinID', indicator=True)
    combo_prot_df = combo_prot_df[(combo_prot_df['_merge'] == 'both') & (combo_prot_df['proteomics_B721'] == 0)]
    combo_prot_df.to_csv('output2/final/canonical721/missed.csv', index=False)

if __name__ == '__main__':
    plot_entro()