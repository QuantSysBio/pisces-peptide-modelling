import pandas as pd
from plotly.colors import n_colors
import plotly.express as px
from pisces.plot_utils import clean_plotly_fig

SPLICED_FOLDER = 'output/splicedK562_241129'
CRYPTIC_FOLDER = 'output/crypticK562_241129'

PROJECTS = ['output/splicedK562_241129', 'output/spliced721_241205']


COLOURS = ['#00BFBF', 'orange']
def plot_aug():
    cryp_cols = n_colors('rgb(34,139,34)', 'rgb(128, 0, 128)', 128, colortype='rgb')
    # print(spliced)
    splice_cols = n_colors('rgb(34,139,34)', 'rgb(0, 0, 255)', 128, colortype='rgb')
    spliced_df = pd.read_csv(
        f'{SPLICED_FOLDER}/mm_scored.csv', usecols=['peptide', 'meanScore']
    )
    spliced_df = spliced_df.rename(
        columns={'meanScore': 'splicedScore'}, 
    )

    cryptic_df = pd.read_csv(
        f'{CRYPTIC_FOLDER}/mm_scored.csv',
        usecols=['peptide', 'meanScore','stratum'],
    )
    cryptic_df = cryptic_df.rename(
        columns={'meanScore': 'crypticScore', 'stratum': 'crypticStratum'}, 
    )

    mm_df = pd.merge(spliced_df, cryptic_df, how='inner', on='peptide')
    mm_df['scoreDiff'] = mm_df['splicedScore'] - mm_df['crypticScore']
    mm_df['color'] = mm_df['scoreDiff'].apply(lambda x : splice_cols[round(128*x)] if x > 0 else cryp_cols[-round(128*x)])
    mm_df.to_csv('tmep.csv', index=False)
    mm_df = mm_df.sort_values('scoreDiff', ascending=False)
    mm_df.to_csv('temp.csv', index=False)
    print(mm_df)
    print(mm_df[(mm_df['splicedScore'] < 0.5) & (mm_df['crypticScore'] < 0.5)].shape[0])
    print(mm_df[(mm_df['splicedScore'] >= 0.5) & (mm_df['crypticScore'] >= 0.5)].shape[0])
    print(mm_df[(mm_df['splicedScore'] < 0.5) & (mm_df['crypticScore'] > 0.5)].shape[0])
    print(mm_df[(mm_df['splicedScore'] > 0.5) & (mm_df['crypticScore'] < 0.5)].shape[0])
    fig = px.scatter(mm_df, x='crypticScore', y='splicedScore', color='color', color_discrete_map=dict(zip(cryp_cols+splice_cols,cryp_cols+splice_cols)))
    fig = clean_plotly_fig(fig)
    fig.add_hline(y=0.5, line_dash='dash', line_color='black', line_width=0.5)
    fig.add_vline(x=0.5, line_dash='dash', line_color='black', line_width=0.5)
    fig.update_traces(marker_line_color='black',marker_line_width=0.25, marker_opacity=0.8)
    fig.update_xaxes(range=[0,1])
    fig.update_yaxes(range=[0,1])
    fig.update_layout(height=300, width=300)
    fig.show()
    fig.write_image('fig6/img/Fig6_G.svg')

    # for stratum in mm_df['stratum'].unique().tolist():
    #     mm_df = mm_df[mm_df['stratum'] == stratum]
    # fig = px.scatter(mm_df, x='crypticScore', y='splicedScore', color='crypticStratum')#, color_discrete_map={'dodgerblue': 'dodgerblue', 'purple': 'purple'})
    # fig = clean_plotly_fig(fig)
    # fig.update_traces(marker_line_color='black',marker_line_width=0.25, marker_opacity=0.8)
    # fig.update_xaxes(range=[0,1])
    # fig.update_yaxes(range=[0,1])
    # # fig.update_layout(title=stratum, title_x=0.5)
    # fig.update_layout(showlegend=True)
    # fig.show()

if __name__ == '__main__':
    plot_aug()
