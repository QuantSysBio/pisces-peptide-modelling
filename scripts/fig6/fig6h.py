import plotly.graph_objects as go
from pisces.plot_utils import clean_plotly_fig
from ppm.constants import CRYPTIC_STRATA

import pandas as pd

STRATUM_COLOUR_SCHEME = {
    'fiveUTR': '#FFA756',
    'threeUTR': '#F2EEB8',
    'CDS_frameshift': '#B790D4',
    'lncRNA': 'lightgrey',
    'intronic': '#4ece78',
    'intergenic': '#528AAE',
    'spliced': '#9BBFE5',
}
SPLICED_FOLDERS = ['output/splicedK562_241129', 'output/spliced721_241206']
CRYPTIC_FOLDERS = ['output/crypticK562_241129', 'output/cryptic721_241205']
ys = []
for s_f, c_f in zip(SPLICED_FOLDERS, CRYPTIC_FOLDERS):
    unique_spliced_df = pd.read_csv(
        f'{s_f}/unique_peps_scored.csv', usecols=['peptide', 'label']
    )
    unique_spliced_df = unique_spliced_df[unique_spliced_df['label'] == 1]
    unique_spliced_df['stratum'] = 'spliced'
    unique_spliced_df['meanScore'] = 1.0

    spliced_df = pd.read_csv(
        f'{s_f}/mm_scored.csv', usecols=['peptide', 'meanScore']
    )
    spliced_df['stratum'] = 'spliced'
    # spliced_df = spliced_df.rename(
    #     columns={'meanScore': 'splicedScore'}, 
    # )


    unique_cryptic_df= pd.read_csv(
        f'{c_f}/unique_peps_scored.csv', usecols=['peptide', 'stratum', 'label']
    )
    unique_cryptic_df['stratum'] = unique_cryptic_df['stratum'].apply(lambda x : CRYPTIC_STRATA[x])
    unique_cryptic_df['meanScore'] = 1.0
    unique_cryptic_df = unique_cryptic_df[unique_cryptic_df['label'] == 1]
    print(unique_cryptic_df.shape)
    print(unique_cryptic_df)

    cryptic_df = pd.read_csv(
        f'{c_f}/mm_scored.csv',
        usecols=['peptide', 'meanScore','stratum'],
    )
    print(cryptic_df)
    # cryptic_df['stratum'] = cryptic_df['stratum'].apply(lambda x : CRYPTIC_STRATA[x])

    mm_df = pd.concat([spliced_df, cryptic_df, unique_cryptic_df, unique_spliced_df])
    print(mm_df.shape)
    mm_df['totalScore'] = mm_df.groupby('peptide')['meanScore'].transform(sum)
    mm_df['softmaxScore'] = mm_df['meanScore']/mm_df['totalScore']
    ys.append(mm_df.groupby('stratum')['softmaxScore'].sum())
    print(mm_df.shape)
print(ys[0])
print(ys[1])

y = ys[0] + ys[1]
print(y)
y = y.to_dict()
fig = go.Figure()
for stratum in ['spliced'] + CRYPTIC_STRATA:
    fig.add_trace(
        go.Bar(x=[stratum], y=[y[stratum] ],
                    marker_color=STRATUM_COLOUR_SCHEME.get(stratum, '#9BBFE5'),
                    marker_line_color='black',)
    )
fig = clean_plotly_fig(fig)
fig.update_yaxes(range=[0,8_000], title='peptides')
fig.update_layout(bargap=0)
fig.show()
fig.write_image('fig6/img/Fig6_H.svg')
