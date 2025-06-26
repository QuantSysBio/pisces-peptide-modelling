
import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig, split_strata, split_nc_strata
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import mannwhitneyu
PEP_DF_LOC = 'data/piscesPeptides.parquet'


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
    'canonical (PISCES)': '#EC9A56',
    'non-canonical': 'pink',
    'spliced': '#9BBFE5',
    'multi-mapped': '#8AA53D',
    'cryptic': '#BA69BE',
    'contaminant': '#B5AE70',
    'unmapped': '#B9B2C2',
}
PLOT_NAMES = {
    'spectralAngle': 'Spectral Angle',
    'spearmanR': 'Spearman Correlation',
    'deltaRT': 'iRT Prediction Error',
    'mhcpanPrediction': 'Predicted BA (log10)',
    'nuggetsPrediction': 'Predicted BA (log10)',
}

def main():

    all_dfs =  []
    binding_affinities = []
    for user in PROJECTS:
        for project in PROJECTS[user]:
            qc_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/nonCanonicalMetrics.csv')
            qc_df['mhcpanPrediction'] = qc_df['mhcpanPrediction'].apply(np.log10)
            qc_df['nuggetsPrediction'] = qc_df['nuggetsPrediction'].apply(np.log10)
            qc_df = qc_df.rename(columns={'adjustedProbability': 'probability'})
        
            can_df = pd.read_csv(f'projects/{user}/{project}/outputFolder/canonicalOutput/finalPsmAssignments.csv')
            can_df = can_df[can_df['postErrProb'] < 0.01]

            ba_cols = [col for col in can_df.columns if col.endswith('BindingAffinity')]
            can_df['mhcpanPrediction'] = can_df[ba_cols].apply(
                lambda df_row : np.log10(np.min([df_row[col] for col in ba_cols])), axis=1,
            )
            can_df['probability'] = 1 - can_df['postErrProb']
            all_dfs.append(qc_df[[
                'source', 'scan', 'peptide', 'modifiedSequence', 'spectralAngle', 'spearmanR',
                'mhcpanPrediction', 'deltaRT', 'probability'
            ]])
            all_dfs.append(can_df[[
                'source', 'scan', 'peptide', 'modifiedSequence', 'spectralAngle', 'spearmanR',
                'mhcpanPrediction', 'deltaRT', 'probability'
            ]])
            binding_affinities.append({
                'dataset': project,
                'mhcpanPrediction_canonical_mean': can_df['mhcpanPrediction'].apply(lambda x : 10**x).mean(),
                'mhcpanPrediction_canonical_median': can_df['mhcpanPrediction'].apply(lambda x : 10**x).median(),
                'mhcpanPrediction_noncanonical_mean': qc_df['mhcpanPrediction'].apply(lambda x : 10**x).mean(),
                'mhcpanPrediction_noncanonical_median': qc_df['mhcpanPrediction'].apply(lambda x : 10**x).median(),
            })
    ba_df = pd.DataFrame(binding_affinities)
    ba_df.to_csv('binding_affinities.csv', index=False)
    all_metric_df = pd.concat(all_dfs)
    all_metric_df = all_metric_df.sort_values('probability', ascending=False)
    print(all_metric_df.shape)
    all_metric_df = all_metric_df.drop_duplicates(subset=['peptide'])
    pep_df = pd.read_parquet(PEP_DF_LOC, columns=[
        'peptide', 'stratum', 'piscesDiscoverable', 'tdDiscoverable',
    ])
    pep_df['sequenceLength'] = pep_df['peptide'].apply(len)
    total_df = pd.merge(all_metric_df, pep_df, how='inner', on=['peptide'], indicator=True)
    baseline_df = s_df = total_df[
        (total_df['stratum'] == 'canonical') &
        (total_df['tdDiscoverable'] == True)
    ]
    signif_data = []
    for feature in ['spectralAngle', 'spearmanR', 'deltaRT', 'mhcpanPrediction']:
        fig = go.Figure()
        for stratum in ['canonical (t/d)', 'canonical (PISCES)']:
            if stratum == 'canonical (t/d)':
                s_df = total_df[
                    (total_df['stratum'] == 'canonical') &
                    (total_df['tdDiscoverable'] == True)
                ]
            elif stratum == 'canonical (PISCES)':
                s_df = total_df[
                    (total_df['stratum'] == 'canonical') &
                    (total_df['piscesDiscoverable'] == True)
                ]

            signif_data.append({
                'stratum': stratum,
                'feature': feature,
                'p-value': mannwhitneyu(s_df[feature].dropna(), baseline_df[feature].dropna()).pvalue
            })
            if feature == 'mhcpanPrediction':
                fig.add_trace(
                    go.Box(
                        x=[stratum]*len(s_df),
                        y=s_df[feature],
                        boxpoints=False,
                        line_color='black',
                        line_width=0.5,
                        fillcolor=COLOURS[stratum],
                    )
                )
            else:
                fig.add_trace(
                    go.Violin(
                        x=[stratum]*len(s_df),
                        y=s_df[feature],
                        points=False,
                        line_color='black',
                        line_width=0.5,
                        fillcolor=COLOURS[stratum],
                        meanline_visible=True,
                    )
                )
        for stratum in ['spliced', 'multi-mapped', 'cryptic', 'unmapped', 'contaminant']:
            s_df = total_df[total_df['stratum'] == stratum]

            if feature in ('mhcpanPrediction', 'nuggetsPrediction'):
                s_df = s_df[s_df['sequenceLength'] > 7]

            signif_data.append({
                'stratum': stratum,
                'feature': feature,
                'p-value': mannwhitneyu(s_df[feature].dropna(), baseline_df[feature].dropna()).pvalue
            })
            if feature == 'mhcpanPrediction':
                fig.add_trace(
                    go.Box(
                        x=s_df['stratum'],
                        y=s_df[feature],
                        boxpoints=False,
                        line_color='black',
                        line_width=0.5,
                        fillcolor=COLOURS[stratum],
                    )
                )
            else:
                fig.add_trace(
                    go.Violin(
                        x=s_df['stratum'],
                        y=s_df[feature],
                        points=False,
                        line_color='black',
                        line_width=0.5,
                        fillcolor=COLOURS[stratum],
                        meanline_visible=True,
                    )
                )
            fig.update_yaxes(title=PLOT_NAMES[feature])
            if feature in ('spectralAngle', 'spearmanR'):
                fig.update_yaxes(range=[0,1])
            elif feature in ('mhcpanPrediction', 'nuggetsPrediction'):
                fig.update_yaxes(range=[0,5])
            else:
                fig.update_yaxes(range=[0,50])

        fig = clean_plotly_fig(fig)
        fig.show()
        if feature == 'spectralAngle':
            pio.write_image(fig, f'manuscript_figs/fig4/fig4c.svg')
        if feature == 'spearmanR':
            pio.write_image(fig, f'manuscript_figs/figS/quality/figS_Quality_a.svg')
        if feature == 'deltaRT':
            pio.write_image(fig, f'manuscript_figs/figS/quality/figS_Quality_b.svg')
        if feature == 'mhcpanPrediction':
            pio.write_image(fig, f'manuscript_figs/figS/quality/figS_Quality_c.svg')

    signif_data_df = pd.DataFrame(signif_data)
    signif_data_df.to_csv('manuscript_figs/fig4/fig4c_significance.csv', index=False)

if __name__ == '__main__':
    main()
