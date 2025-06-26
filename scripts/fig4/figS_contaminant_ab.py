from math import ceil
import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, get_cryptic_stratum_counts, DARK_COLOURS, LIGHT_COLORS
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

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
    'blue',
    '#e86a58',
    '#fed45b',
    '#9bc7c5',
    'cyan',
    '#908FB3',
]
LIGHT_COLORS = [
    '#23362b',
    'lightblue',
    '#e86a58',
    '#fed45b',
    '#9bc7c5',
    'black',
]
def main():
    remapped_dfs = []
    # cryptic_dfs = []
    for user in PROJECTS:
        for project in PROJECTS[user]:
            print(project)
            remapped_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/filtered_mapped.csv')
            remapped_df = remapped_df[remapped_df['adjustedProbability'] > 0.85]
            remapped_df = remapped_df.drop([
                'source', 'scan', 'modifiedSequence', 'piscesScore1',
                'secondScore', 'deltaScore', 'delta3Score', 'isoDeltaScore',
                'isoDelta2Score', 'isoDelta3Score', 'mutated', 'piscesScore',
                'adjustedProbability', 'qValue_PSM', 'qValue_peptide',
                'canonical_nProteins', 'nSplicedProteins','fiveUTR_nProteins',
                'TrEMBL_nProteins', 'intergenic_nProteins', 'CDS_frameshift_nProteins',
                'threeUTR_nProteins', 'lncRNA_nProteins', 'intronic_nProteins',
                'nCrypticProteins', 'charge', 'perPositionScores'], axis=1
            )

            remapped_df = remapped_df[(remapped_df['nSpecific_ContamsProteins'] > 0) | (remapped_df['nContamProteins'] > 0)]
            remapped_dfs.append(remapped_df)
            # cryptic_dfs.append(cryptic_df)

    remapped_df = pd.concat(remapped_dfs)
    print(remapped_df.columns)
    remapped_df = remapped_df.fillna(0)
    remapped_df = remapped_df.groupby('peptide', as_index=False).sum()
    pep_df = pd.read_parquet('data/piscesPeptides.parquet', columns=['peptide', 'stratum'])
    pep_df = pep_df[pep_df['stratum'] == 'contaminant']
    print(remapped_df.shape)
    print(pep_df.shape)
    remapped_df = pd.merge(pep_df, remapped_df, on=['peptide'], how='inner') 
    print(pep_df.shape[0] - remapped_df.shape[0])
    remapped_df['nContamProteins'] = remapped_df['nContamProteins'].fillna(1)
    print(remapped_df.shape)
    remapped_df = remapped_df.drop(['stratum'], axis=1)
    
    for col in remapped_df.columns:
        if col not in ('dataset', 'peptide'):
            print(col)
            remapped_df[col] = remapped_df[col].apply(lambda x : 1 if x > 0 else 0)
    print(remapped_df)
    print(remapped_df.sum())
    # cryptic_df = pd.concat(cryptic_dfs)
    generic_peps = set(remapped_df[remapped_df['nContamProteins'] > 0]['peptide'].tolist())
    specific_peps = set(remapped_df[remapped_df['nSpecific_ContamsProteins'] > 0]['peptide'].tolist())

    v = venn2([generic_peps, specific_peps], ('Generic Contaminants', 'Specific Contaminants'))


    v.get_patch_by_id('10').set_color('purple')
    v.get_patch_by_id('01').set_color('cyan')
    v.get_patch_by_id('11').set_color('#729efd')
    plt.tight_layout()
    plt.savefig('manuscript_figs/figS/contaminant/figS_contaminant_a.svg', format='svg')
    # print(cryptic_df.columns)
    # fig = plot_cryptic_breakdown(
    #     remapped_df, 
    # )
    # fig.show()
    e_coli_strains = [col for col in remapped_df.columns if col.startswith('Escherichia coli')]
    remapped_df['Escherichia coli (various strains)'] = remapped_df[e_coli_strains].apply(
        lambda x : 1 if sum([x[col] for col in e_coli_strains]) else 0, axis=1
    )
    remapped_df = remapped_df.drop(e_coli_strains, axis=1)
    print(remapped_df.sum())
    strep_cols = [col for col in remapped_df.columns if col.startswith('Streptococcus ')]
    remapped_df['Streptococcus_nProteins'] = remapped_df[strep_cols].apply(
        lambda x : 1 if sum([x[col] for col in strep_cols]) else 0, axis=1
    )
    remapped_df = remapped_df.drop(strep_cols, axis=1)
    print(remapped_df.sum(numeric_only=True).sort_values())
    non_grouped_cols = [
        'Streptococcus_nProteins',
        'Bos taurus_nProteins',
        'maxquant_nProteins',
        'MaConDa_nProteins',
        'nContamProteins',
        'nSpecific_ContamsProteins',
        'Mus musculus_nProteins',
        'Rattus norvegicus_nProteins',
    ]
    group_cols = [col for col in remapped_df.columns if col not in non_grouped_cols + ['peptide', 'dataset']]
    print(group_cols)
    remapped_df['Other'] = remapped_df[group_cols].apply(
        lambda x : 1 if sum([x[col] for col in group_cols]) else 0, axis=1
    )
    remapped_df = remapped_df.drop(group_cols, axis=1)
    print(remapped_df.sum(numeric_only=True).sort_values())

    out_cols = [
        'nContamProteins',
        'Rattus norvegicus_nProteins',
        'Mus musculus_nProteins',
        'Bos taurus_nProteins',
        'Streptococcus_nProteins',
        'Other',
    ]
    results = []
    total_unique = 0
    for col in out_cols:
        sub_df = remapped_df[remapped_df[col] == 1]
        total_count = sub_df.shape[0]
        for col2 in out_cols:
            if col != col2:
                sub_df = sub_df[sub_df[col2] == 0]
        unique_count = sub_df.shape[0]
        if col == 'nContamProteins':
            name = 'MaxQuant/MaConDa'
        else:
            name = col.split('_')[0].lower()
        results.append({
            'proteome': name,
            # 'totalCount': total_count,
            'uniqueCount': unique_count
        })
        total_unique += unique_count
    
    results.append({
        'proteome': 'multi-mapped',
        'uniqueCount': remapped_df.shape[0] - total_unique,
    })
    res_df = pd.DataFrame(results)
    print(res_df)
    fig = go.Figure()
    for idx, df_row in res_df.iterrows():
        # fig.add_trace(
        #     go.Bar(
        #         x=[df_row['proteome']],
        #         y=[df_row['totalCount']],
        #         marker_color=LIGHT_COLORS[idx],
        #         marker_opacity=0.5,
        #         marker_line_color='black',
        #     )
        # )
        fig.add_trace(
            go.Bar(
                x=[df_row['proteome']],
                y=[df_row['uniqueCount']],
                marker_color=DARK_COLOURS[idx],
                marker_line_color='black',
            )
        )

    fig.update_layout(
        barmode='overlay',
    )
    fig = clean_plotly_fig(fig)
    fig.update_yaxes(title='# Peptides', range=[0, 3000])

    fig.update_layout(
        height=400, bargap=0,
    )
    fig.show()
    pio.write_image(fig, f'manuscript_figs/figS/contaminant/figS_contaminant_b.svg')


if __name__ == '__main__':
    main()