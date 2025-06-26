import pandas as pd
import plotly.graph_objects as go
from ppm.analysis_utils import plot_distros
from pisces.plot_utils import clean_plotly_fig
from scipy.stats import linregress, fisher_exact
from ppm.constants import AMINO_ACIDS
from statsmodels.stats.multitest import multipletests
PROJECTS = ['output/final/crypticK562', 'output/final/cryptic721']

pep_df = pd.concat([pd.read_csv(f'{project}/unique_peps_scored.csv') for project in PROJECTS])
pep_df['KR'] = pep_df[['peptide', 'C_term_downstream', 'C_term_downstream_2']].apply(
    lambda x : int(
        x['peptide'][-1] in 'KRH'
    ), axis=1,
)
print(pep_df[['KR', 'peptide', 'C_term_basic']])
count = pep_df[pep_df['label'] == 1].shape[0]
for label in [0,1]:
    pos_pep_df = pep_df[pep_df['label'] == label]
    if label == 0:
        pos_pep_df = pos_pep_df.sample(n=count, random_state=42)
    color_scheme=['#908FB3', '#2E8B58']
    fig = go.Figure()
    count_res = []
    for kr in [0,1]:
        sub_df = pos_pep_df[pos_pep_df['C_term_basic'] == kr]
        counts = []
        for a_a in AMINO_ACIDS + 'X':
            count_res.append({
                'aa': a_a,
                'kr': kr,
                'count': sub_df[sub_df['C_term_downstream'] == a_a].shape[0],
                'count_not': sub_df[sub_df['C_term_downstream'] != a_a].shape[0],
            })
            counts.append(100*sub_df[sub_df['C_term_downstream'] == a_a].shape[0]/sub_df.shape[0])
        fig.add_trace(go.Bar(x=list(AMINO_ACIDS+'X'), y=counts, marker_color=color_scheme[kr], marker_line_color='black', name=f'KR at C-term={kr}'))
    count_df = pd.DataFrame(count_res)
    no_kr_df = count_df[count_df['kr'] == 0].rename(columns={
        'count': 'no_KR_count',
        'count_not': 'no_KR_count_not',
    }).drop('kr', axis=1)
    kr_df = count_df[count_df['kr'] == 1].rename(columns={
        'count': 'KR_count',
        'count_not': 'KR_count_not',
    }).drop('kr', axis=1)
    tot_df = pd.merge(no_kr_df, kr_df, on='aa')
    tot_df['p_val'] = tot_df.apply(
        lambda x : (
            fisher_exact([
                [x['KR_count'], x['KR_count_not']],
                [x['no_KR_count'], x['no_KR_count_not']]
            ]).pvalue
        ), axis=1
    )
    tot_df['adjusted_pValue'] = multipletests(tot_df['p_val'], method='fdr_bh')[1]
    tot_df.to_csv(f'manuscript_figs/figS/KR/kr_counts_{label}.csv', index=False)
    fig = clean_plotly_fig(fig)
    fig.update_layout(bargap=0, )
    fig.update_yaxes(range=[0,30])

    fig.write_image(f'manuscript_figs/figS/KR/figS_KR_{label}.svg')
    fig.show()