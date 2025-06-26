
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
        'output/final/canonicalK562',
        'output/final/crypticK562',
        'output/final/splicedK562',
    ],
    [
        'output/final/canonical721',
        'output/final/cryptic721',
        'output/final/spliced721',
    ]
]
TRYPTIC_RESULTS = [
    'projects/pisces_tryptic/K562-tryptic/outputFolder',
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
    fracs = {
        'K562': [],
        '721.221': [],
    }
    for projects, tryptic_folder, cell_line in zip(PROJECTS, TRYPTIC_RESULTS, ['K562', '721.221']):
        for stratum, project in zip(STRATA, projects):
            unique_pep_df = pd.read_csv(f'{project}/unique_peps_scored.csv')
            unique_pep_df = unique_pep_df[unique_pep_df['label'] == 1]
            proteins = unique_pep_df[['proteinID']].drop_duplicates()
            print(proteins.shape[0], stratum, cell_line)
            if stratum == 'cryptic':
                tryptic_df = pd.concat([
                    pd.read_csv(f'{tryptic_folder}/final/{sub_group}.csv') for sub_group in ['cryptic', 'multimapped']
                ])
                tryptic_df = tryptic_df[tryptic_df['adjustedProbability'] > 0.85]
                peptide_df = tryptic_df[['peptide']].drop_duplicates()
                tryptic_det_df = pd.concat([
                    pd.read_csv(f'{tryptic_folder}/details/{sub_group}.csv') for sub_group in ['cryptic',]
                ])
                tryptic_det_df = tryptic_det_df.drop_duplicates(subset=['peptide'])
                total_df = pd.merge(peptide_df, tryptic_det_df, how='inner', on='peptide')
                print(total_df.shape[0])
                for stratum in CRYPTIC_STRATA:
                    total_df['crypticProteins'] = total_df.apply(
                        lambda x : ' '.join([x[f'{stratum}_Proteins'] for stratum in CRYPTIC_STRATA if x[f'{stratum}_Proteins'] is not None and isinstance(x[f'{stratum}_Proteins'], str)]), axis=1,
                    )
                total_df['proteinID'] = total_df['crypticProteins'].apply(lambda x : x.split(' '))
                tryp_prot_df = total_df[['proteinID']].explode('proteinID').drop_duplicates()

            else:
                tryptic_df = pd.read_csv(f'{tryptic_folder}/final/canonical.csv')
                tryptic_df = tryptic_df[tryptic_df['piscesDiscoverable'] == 1]

                tryptic_df = tryptic_df[tryptic_df['postErrProb'] < 0.15]
                peptide_df = tryptic_df[['peptide']].drop_duplicates()
                tryptic_det_df = pd.read_csv(f'{tryptic_folder}/details/canonical.csv')
                tryptic_det_df = tryptic_det_df.drop_duplicates(subset=['peptide'])

                total_df = pd.merge(peptide_df, tryptic_det_df, how='inner', on='peptide')
                print(total_df.shape[0])

                total_df['proteinID'] = total_df['canonical_Proteins'].apply(lambda x : x.split(' '))
                tryp_prot_df = total_df[['proteinID']].explode('proteinID').drop_duplicates()


            combo_prot_df = pd.merge(proteins, tryp_prot_df, how='left', on='proteinID', indicator=True)
            fracs[cell_line].append(100*combo_prot_df[combo_prot_df['_merge'] == 'both'].shape[0]/combo_prot_df.shape[0])
            print(combo_prot_df.shape[0], stratum, cell_line)

    print(fracs)
    fig = go.Figure()
    # mean_fracs = [fracs[cell_line][idx] + fracs[cell_line][idx] for idx in range(3)]
    for stratum, frac in zip(STRATA, fracs[cell_line]):
        fig.add_trace(
            go.Bar(
                x=[stratum], y=[frac],
                marker_color=COLOUR_DICT[stratum], opacity=0.9,
                marker_line_color='black',
                marker_line_width=0.5,
                # meanline_visible=True,
                # points=False,
            ),
        )
    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0, line_width=0.5, row=1, col=2)
    fig.update_layout({
        'yaxis': dict(range=[0, 100], title_text='antigens found in tryptic proteomics'),
        # 'yaxis2': dict(range=[-1, 1], title_text='impact on model'),
        'xaxis': dict(title_text='stratum'),
        # 'xaxis2': dict(range=[0, 20], title_text='distance to start codon', linecolor='white')
    }, width=200, height=300, bargap=0)
    # fig.write_image('fig5/img/Fig5_H.svg')
    fig.show()


if __name__ == '__main__':
    plot_entro()