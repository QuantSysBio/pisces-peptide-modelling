import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pisces.plot_utils import clean_plotly_fig
from scipy.stats import pearsonr, linregress, spearmanr

def main():
    cell_line_df = pd.read_csv('scripts/fig4/per_allele_data.csv')
    cell_line_df['HLA_allele'] = cell_line_df['allele'].apply(lambda x : x.split('-')[1])
    cell_line_df['HLA_allele'] = cell_line_df['HLA_allele'].apply(lambda x : x[0] + '*' + x[1:3] + ':' + x[3:])
    cell_line_df['cellLine'] = cell_line_df['allele'].apply(
        lambda x : 'K562' if x.split('-')[0] == 'K562_MMlab' else 'B721.221'
    )
    b721_df = cell_line_df[cell_line_df['cellLine'] == 'B721.221'][
        ['HLA_allele', 'cisFreq', 'crypFreq']
    ].rename(columns={'cisFreq': 'b721cisFreq', 'crypFreq': 'b721crypFreq'})
    k562_df = cell_line_df[cell_line_df['cellLine'] == 'K562'][
        ['HLA_allele', 'cisFreq', 'crypFreq']
    ].rename(columns={'cisFreq': 'k562cisFreq', 'crypFreq': 'k562crypFreq'})
    tot_df = pd.merge(b721_df, k562_df, on='HLA_allele', how='inner')
    tot_df = tot_df.sort_values(by='k562cisFreq').reset_index(drop=True)
    print(linregress(tot_df['k562cisFreq'], tot_df['b721cisFreq']))
    print(spearmanr(tot_df['k562cisFreq'], tot_df['b721cisFreq']))
    print(pearsonr(tot_df['k562cisFreq'], tot_df['b721cisFreq']))
    lm_res = linregress(tot_df['k562cisFreq'], tot_df['b721cisFreq'])
    x_range = np.linspace(0, 3)
    y_range = (lm_res.slope*x_range) + lm_res.intercept
    print(tot_df)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=tot_df['k562cisFreq'],
            y=tot_df['b721cisFreq'],
            # text=tot_df['HLA_allele'],
            marker_size=9,
        marker_line_color='black', marker_line_width=0.3,
            mode='markers',
            # textposition='bottom left',
            marker_color='#9BBFE5',
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=x_range, y=y_range, mode='lines',
            line={'width':0.5, 'dash': 'solid', 'color': 'black'},
        ),
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        xaxis_title='K562 spliced frequency',
        yaxis_title='B721.221 spliced frequency',
        width=300,
        height=300,
    )
    fig.update_xaxes(range=[0,3])
    fig.update_yaxes(range=[0,3])
    fig.show()
    fig.write_image('manuscript_figs/figS/cellLines/figS_CellLines_b.svg')
    fig = go.Figure()
    tot_df = tot_df.sort_values(by='b721crypFreq').reset_index(drop=True)
    print(linregress(tot_df['k562crypFreq'], tot_df['b721crypFreq']))
    print(spearmanr(tot_df['k562crypFreq'], tot_df['b721crypFreq']))
    print(pearsonr(tot_df['k562crypFreq'], tot_df['b721crypFreq']))
    lm_res = linregress(tot_df['k562crypFreq'], tot_df['b721crypFreq'])
    x_range = np.linspace(0, 10)
    y_range = (lm_res.slope*x_range) + lm_res.intercept
    print(tot_df)
    tot_df.to_excel('manuscript_figs/figS/cellLines/cellLines_data.xlsx', index=False)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_range, y=y_range, mode='lines',
            line={'width':0.5, 'dash': 'solid', 'color': 'black'},
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=tot_df['k562crypFreq'],
            y=tot_df['b721crypFreq'],
            marker_color='#BA69BE',        marker_size=9,
        marker_line_color='black', marker_line_width=0.3,
            # text=tot_df['HLA_allele'],
            mode='markers',
            # textposition='bottom right',
        ),
    )
    fig = clean_plotly_fig(fig)
    fig.update_layout(
        xaxis_title='K562 cryptic frequency',
        yaxis_title='B721.221 cryptic frequency',
        width=300,
        height=300,
    )
    fig.update_xaxes(range=[0,10], dtick=2)
    fig.update_yaxes(range=[0,10], dtick=2)
    fig.show()
    fig.write_image('manuscript_figs/figS/cellLines/figS_CellLines_a.svg')

if __name__ == '__main__':
    main()