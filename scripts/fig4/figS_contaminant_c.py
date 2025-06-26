from math import ceil
import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, get_cryptic_stratum_counts, DARK_COLOURS, LIGHT_COLORS
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.express as px
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
PROJECTS = {
    'pisces_ip': [
        'Abelin-HF-A2902',
        'Abelin-HF-B3501',
        'Abelin-HF-B5101',
        'Abelin-HF-B5401',
        'Abelin-HF-B5701',
        'Abelin-QE-A0101',
        'Abelin-QE-A0201',
        'Abelin-QE-A0203',
        'Abelin-QE-A0204',
        'Abelin-QE-A0207',
        'Abelin-QE-A0301',
        'Abelin-QE-A2402',
        'Abelin-QE-A3101',
        'Abelin-QE-A6802',
        'Abelin-QE-B4402',
        'Abelin-QE-B4403',
        'Abelin-QE-B5101',
        'Abelin-QE-B5401',
        'Abelin-QE-B5701',
        'K562_MMlab-Exploris-A0101',
        'K562_MMlab-Exploris-A0201',
        'K562_MMlab-Exploris-A0301',
        'K562_MMlab-Exploris-A1101',
        'K562_MMlab-Exploris-A3001',
        'K562_MMlab-Exploris-A3101',
        'K562_MMlab-Exploris-A6801',
        'K562_MMlab-Exploris-B0702',
        'K562_MMlab-Exploris-B0801',
        'K562_MMlab-Exploris-B1501',
        'K562_MMlab-Exploris-B2705',
        'K562_MMlab-Exploris-B4001',
        'K562_MMlab-Lumos-A0201',
        'K562_MMlab-Lumos-B0702',
        'Sarkizova-Lumos-A0202',
        'Sarkizova-Lumos-A0205',
        'Sarkizova-Lumos-A0206',
        'Sarkizova-Lumos-A0211',
        'Sarkizova-Lumos-A1101',
        'Sarkizova-Lumos-A2407',
        'Sarkizova-Lumos-A2501',
        'Sarkizova-Lumos-A2601',
        'Sarkizova-Lumos-A3201',
        'Sarkizova-Lumos-A3301',
        'Sarkizova-Lumos-A3402',
        'Sarkizova-Lumos-A3601',
        'Sarkizova-Lumos-A6601',
        'Sarkizova-Lumos-A6801',
        'Sarkizova-Lumos-B1510',
        'Sarkizova-Lumos-B1517',
        'Sarkizova-Lumos-B3503',
        'Sarkizova-Lumos-B3802',
        'Sarkizova-Lumos-B4001',
        'Sarkizova-Lumos-B4201',
        'Sarkizova-Lumos-B4501',
        'Sarkizova-Lumos-B5703',
        'Sarkizova-Lumos-B5801',
        'Sarkizova-Lumos-B5802',
        'Sarkizova-Lumos-C0701',
        'Sarkizova-Lumos-C0704',
        'Sarkizova-Lumos-C0801',
        'Sarkizova-Lumos-C0802',
        'Sarkizova-QEplus-A2301',
        'Sarkizova-QEplus-A3001',
        'Sarkizova-QEplus-A3002',
        'Sarkizova-QEplus-A3303',
        'Sarkizova-QEplus-A3401',
        'Sarkizova-QEplus-A7401',
        'Sarkizova-QEplus-B0702',
        'Sarkizova-QEplus-B0704',
        'Sarkizova-QEplus-B0801',
        'Sarkizova-QEplus-B1301',
        'Sarkizova-QEplus-B1302',
        'Sarkizova-QEplus-B1402',
        'Sarkizova-QEplus-B1501',
        'Sarkizova-QEplus-B1502',
        'Sarkizova-QEplus-B1503',
        'Sarkizova-QEplus-B1801',
        'Sarkizova-QEplus-B2705',
        'Sarkizova-QEplus-B3507',
        'Sarkizova-QEplus-B3701',
        'Sarkizova-QEplus-B3801',
        'Sarkizova-QEplus-B4002',
        'Sarkizova-QEplus-B4006',
        'Sarkizova-QEplus-B4201',
        'Sarkizova-QEplus-B4601',
        'Sarkizova-QEplus-B4901',
        'Sarkizova-QEplus-B5001',
        'Sarkizova-QEplus-B5201',
        'Sarkizova-QEplus-B5301',
        'Sarkizova-QEplus-B5501',
        'Sarkizova-QEplus-B5502',
        'Sarkizova-QEplus-B5601',
        'Sarkizova-QEplus-C0102',
        'Sarkizova-QEplus-C0202',
        'Sarkizova-QEplus-C0302',
        'Sarkizova-QEplus-C0303',
        'Sarkizova-QEplus-C0304',
        'Sarkizova-QEplus-C0401',
        'Sarkizova-QEplus-C0403',
        'Sarkizova-QEplus-C0501',
        'Sarkizova-QEplus-C0602',
        'Sarkizova-QEplus-C0701',
        'Sarkizova-QEplus-C1202',
        'Sarkizova-QEplus-C1203',
        'Sarkizova-QEplus-C1402',
        'Sarkizova-QEplus-C1403',
        'Sarkizova-QEplus-C1502',
        'Sarkizova-QEplus-C1601',
        'Sarkizova-QEplus-C1701',
        'Sarkizova-QEplus-C0702',
        'Sarkizova-QEplus-C0704',
        'Sarkizova-QEplus-C0801',
        'Sarkizova-QEplus-C0802',
    ]
}

MAIN_NAMES = ['canonical (t/d)', 'canonical (pisces)', 'contaminant', 'non-canonical']
NAMES = ['spliced', 'multi-mapped', 'cryptic', 'unmapped',]
ALL_NAMES = MAIN_NAMES + NAMES

DARK_COLOURS = [
    '#23362b',
    '#1bb28c',
    '#e86a58',
    '#fed45b',
    '#9bc7c5',
    'black',
]
LIGHT_COLORS = [
    '#23362b',
    '#1bb28c',
    '#e86a58',
    '#fed45b',
    '#9bc7c5',
    'black',
]
def main():
    remapped_counts = []
    # cryptic_dfs = []
    for user in PROJECTS:
        for project in PROJECTS[user]:
            print(project)
            remapped_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/filtered_mapped.csv')
            remapped_df = remapped_df[remapped_df['adjustedProbability'] > 0.85]
            remapped_df['dataset'] = project
            remapped_df = remapped_df.drop([
                'source', 'scan', 'modifiedSequence', 'piscesScore1',
                'secondScore', 'deltaScore', 'delta3Score', 'isoDeltaScore',
                'isoDelta2Score', 'isoDelta3Score', 'mutated', 'piscesScore',
                'adjustedProbability', 'qValue_PSM', 'qValue_peptide',
                'canonical_nProteins', 'nSplicedProteins','fiveUTR_nProteins',
                'TrEMBL_nProteins', 'intergenic_nProteins', 'CDS_frameshift_nProteins',
                'threeUTR_nProteins', 'lncRNA_nProteins', 'intronic_nProteins',
                'nCrypticProteins', 'charge',], axis=1
            )
            remapped_df = remapped_df[(remapped_df['nSpecific_ContamsProteins'] > 0) | (remapped_df['nContamProteins'] > 0)]
            remapped_counts.append({
                'genericCount': remapped_df[remapped_df['nContamProteins'] > 0]['peptide'].nunique(),
                'specificCount': remapped_df[remapped_df['nSpecific_ContamsProteins'] > 0]['peptide'].nunique(),
            })
            # cryptic_dfs.append(cryptic_df)

    remapped_df = pd.DataFrame(remapped_counts)
    
    fig = go.Figure()
    fig = px.strip(
        x=['generic']*remapped_df.shape[0] + ['specific']*remapped_df.shape[0],
        y=remapped_df['genericCount'].tolist() + remapped_df['specificCount'].tolist(),
        color = ['generic']*remapped_df.shape[0] + ['specific']*remapped_df.shape[0],
        color_discrete_map={'generic': 'turquoise', 'specific': 'tomato'},
    )

    fig.update_layout(
        barmode='overlay',
    )
    fig = clean_plotly_fig(fig)
    fig.update_yaxes(title='# Peptides', range=[0, 800], dtick=200)
    fig.update_xaxes(title='')

    fig.update_layout(
        height=400, width=200, bargap=0, boxmode="overlay",
    )
    fig.show()
    pio.write_image(fig, f'manuscript_figs/figS/contaminant/figS_contaminant_c.svg')


if __name__ == '__main__':
    main()