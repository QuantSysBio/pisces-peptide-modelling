import numpy as np
import pandas as pd
from pisces.plot_utils import clean_plotly_fig
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.stats import pearsonr, linregress, spearmanr
SR1_PATH = 'fig6/invitro_comp/SR1_diffLogo.csv'
SR2_PATH = 'fig6/invitro_comp/SR2_diffLogo.csv'
PISCES_SR1_ENTROPY_PATH = 'fig6/img/Fig6_C_sr1.csv'
PISCES_SR2_ENTROPY_PATH = 'fig6/img/Fig6_C_sr2.csv'

COMPARISON_POSITIONS_SR1 = [ 'p1', 'p_neg_1',]
COMPARISON_POSITIONS_SR2 = [ 'p_neg_1_prime', 'p1_prime',]
AA_COLOUR_SCHEME = {
    'P': 'deeppink',
    'M': 'orange',
    'A': 'orange',
    'V': 'orange',
    'I': 'orange',
    'L': 'orange',
    'F': 'orange',
    'Y': 'orange',
    'W': 'orange',
    'H': 'seagreen',
    'R': 'seagreen',
    'K': 'seagreen',
    'D': 'firebrick',
    'E': 'firebrick',
    'N': 'dodgerblue',
    'Q': 'dodgerblue',
    'S': 'dodgerblue',
    'T': 'dodgerblue',
    'G': 'dodgerblue',
    'C': 'dodgerblue',
    'X': 'black',
}
def main():
    entropy_dfs = []
    sr_dfs = []
    for ent_path, comp_pos, sr_path in zip(
        [PISCES_SR1_ENTROPY_PATH, PISCES_SR2_ENTROPY_PATH],
        [COMPARISON_POSITIONS_SR1, COMPARISON_POSITIONS_SR2],
        [SR1_PATH, SR2_PATH],
    ):
        entropy_df = pd.read_csv(ent_path)
        entropy_df = entropy_df.transpose()
        entropy_df = entropy_df.drop([0, 3], axis=1)

        entropy_df.columns = comp_pos
        entropy_df = entropy_df[entropy_df.index != 'X']
        entropy_df['colour'] = entropy_df.index
        entropy_df['colour'] = entropy_df['colour'].apply(lambda x : AA_COLOUR_SCHEME[x] )
        entropy_dfs.append(entropy_df)
        # sr_dfs.append(pd.read_csv(sr_path, index_col='aa'))
        sr_df = pd.read_csv(sr_path)
        sr_df.index = sr_df['aa']
        if 'sr1' in ent_path:
            sr_df = sr_df[['p1', 'p_neg_1']]
        else:
            sr_df = sr_df[['p_neg_1_prime','p1_prime']]
        # sr_df = sr_df.transpose()

        sr_df.columns = comp_pos
        sr_df = sr_df[sr_df.index != 'X']
        sr_df['colour'] = sr_df.index
        sr_df['colour'] = sr_df['colour'].apply(lambda x : AA_COLOUR_SCHEME[x] )
        sr_dfs.append(sr_df)

    row_idx = 0
    for sr_df, ent_df, comp_pos in zip(sr_dfs, entropy_dfs, [COMPARISON_POSITIONS_SR1, COMPARISON_POSITIONS_SR2]):
        row_idx += 1
        for idx, pos in enumerate(comp_pos):
            fig = go.Figure()
            # sr_df = sr_df[sr_df.index!='R']
            # sr_df = sr_df[sr_df.index!='D']
            # sr_df = sr_df[sr_df.index!='E']
            # ent_df = ent_df[ent_df.index!='R']
            # ent_df = ent_df[ent_df.index!='D']
            # ent_df = ent_df[ent_df.index!='E']
            print(pearsonr(sr_df[pos], ent_df[pos]))
            print(spearmanr(sr_df[pos], ent_df[pos]))
            print(linregress(sr_df[pos], ent_df[pos]))
            lm_res = linregress(sr_df[pos], ent_df[pos])
            x_range = np.linspace(sr_df[pos].min(), sr_df[pos].max())
            y_range = (lm_res.slope*x_range) + lm_res.intercept

            col_idx = 1+idx%2
            fig.add_trace(
                go.Scatter(x=sr_df[pos], y=ent_df[pos], mode='markers+text', text=ent_df.index,
                        textposition="bottom center", marker_color=ent_df['colour'],
                        marker_line_width=0.25, marker_size=10, marker_line_color='black',
                        ),
                # row=row_idx, col=col_idx
            )
            fig.add_trace(
                go.Scatter(
                    x=x_range, y=y_range, mode='lines',
                    line={'width':1, 'dash': 'solid', 'color': 'red'},
                ),
                # row=row_idx, col=col_idx,
            )
            fig.add_vline(x=0, line={'width':0.5, 'dash': 'dash', 'color': 'black'})
            fig.add_hline(y=0, line={'width':0.5, 'dash': 'dash', 'color': 'black'})
            fig = clean_plotly_fig(fig)
            fig.update_layout(height=400, width=400)
            if row_idx == 1 and col_idx == 1:
                fig.update_layout({
                    # 'title': 'sr1', 'title_x': 0.5,
                    'xaxis': {'range': [-0.01, 0.015]},
                    'yaxis': {'range': [-0.01, 0.01]},
                })
                fig.write_image('fig6/img/Fig6_D.svg')
            elif row_idx == 1:
                fig.update_layout({
                    # 'title': 'sr-1', 'title_x': 0.5,
                    'xaxis': {'range': [-0.004, 0.004]},
                    'yaxis': {'range': [-0.002, 0.004], 'dtick': 0.002},
                })
                fig.write_image('fig6/img/FigS_invitro_A.svg')
            elif row_idx == 2 and col_idx == 1:
                fig.update_layout({
                    # 'title': 'sr1"', 'title_x': 0.5,
                    'xaxis': {'range': [-0.005, 0.01]},
                    'yaxis': {'range': [-0.002, 0.006], 'dtick': 0.002},
                })
                fig.write_image('fig6/img/FigS_invitro_B.svg')
            else:
                fig.update_layout({
                    # 'title': 'sr-1"', 'title_x': 0.5,
                    'xaxis': {'range': [-0.004, 0.004]},
                    'yaxis': {'range': [-0.006, 0.008]},
                })
                fig.write_image('fig6/img/FigS_invitro_C.svg')

            fig.show()

if __name__ == '__main__':
    main()
