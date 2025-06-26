
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, split_strata, split_nc_strata
import plotly.graph_objects as go
import plotly.io as pio

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

COLOURS = {
    'canonical (t/d)': '#F4BF91',
    'canonical (pisces)': '#EC9A56',
    'non-canonical': 'pink',
    'spliced': '#9BBFE5',
    'multi-mapped': '#8AA53D',
    'cryptic': '#BA69BE',
    'contaminant': '#B5AE70',
    'unmapped': '#B9B2C2',
}

def main():
    all_dfs = {
        'spliced': [], 'multi-mapped': [], 'cryptic': [], 'unmapped': [],
        'canonical (t/d)': [], 'canonical (pisces)': [], 'contaminant': [], 'non-canonical': [],
    }
    for user in PROJECTS:
        for project in PROJECTS[user]:
            remapped_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/filtered_mapped.csv')
            remapped_df = remapped_df[
                remapped_df['adjustedProbability'] > 0.85
            ]
            can_td_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/canonicalOutput/finalPsmAssignments.csv')
            can_td_df = can_td_df[can_td_df['qValue'] < 0.01]
            can_df, contam_df, nc_df = split_strata(remapped_df)
            spliced_df, mm_df, cryptic_df, unmapped_df = split_nc_strata(nc_df)
            print(project, nc_df['peptide'].nunique())
            for df, name in zip(
                [can_td_df, can_df, contam_df, nc_df, spliced_df, mm_df, cryptic_df, unmapped_df],
                ALL_NAMES,
            ):
                all_dfs[name].append(df['peptide'].drop_duplicates())

    all_counts = {
        a_k: pd.concat([x for x in a_v]).nunique() for (a_k, a_v) in all_dfs.items()
    }
        

    fig1 = go.Figure()
    fig_4a_data = []
    for dataset_name in [
        'canonical (t/d)', 'canonical (pisces)', 'contaminant', 'non-canonical'
    ]:
        fig_4a_data.append({'name': dataset_name, 'count': all_counts[dataset_name]})
        fig1.add_trace(
            go.Bar(
                x=[dataset_name],
                y=[all_counts[dataset_name]],
                marker_color=COLOURS[dataset_name],
                marker_line_color='black',
                opacity=0.8,
            )
        )
    fig_4a_df = pd.DataFrame(fig_4a_data)
    fig_4a_df.to_excel('manuscript_figs/fig4/fig4a.xlsx', index=False)
    fig1 = clean_plotly_fig(fig1)
    fig1.update_layout(height=300, width=300)
    fig1['layout']['yaxis']['range'] = [0,250_000]
    fig1['layout']['yaxis']['title'] = '# Peptides'


    fig2 = go.Figure()
    fig_4b_data = []
    for dataset_name in [
        'spliced', 'multi-mapped', 'cryptic', 'unmapped'
    ]:
        fig_4b_data.append({'name': dataset_name, 'count': all_counts[dataset_name]})
        fig2.add_trace(
            go.Bar(
                x=[dataset_name],
                y=[all_counts[dataset_name]],
                marker_color=COLOURS[dataset_name],
                marker_line_color='black',
                opacity=0.8,
            )
        )
    fig_4b_df = pd.DataFrame(fig_4b_data)
    fig_4b_df.to_excel('manuscript_figs/fig4/fig4b.xlsx', index=False)
    fig2 = clean_plotly_fig(fig2)
    fig2.update_layout(height=300, width=300)
    fig2['layout']['yaxis']['range'] = [0, 8_000]
    fig2['layout']['yaxis']['title'] = '# Peptides'


    fig1.show()
    fig2.show()
    pio.write_image(fig1, 'manuscript_figs/fig4/fig4a.svg')
    pio.write_image(fig2, 'manuscript_figs/fig4/fig4b.svg')

if __name__ == '__main__':
    main()
