import pandas as pd
from pisces.plot_utils import clean_plotly_fig, create_comparison_logo_plot
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

SR_COLS = {
    1: ['p2','p1','p_neg_1','p_neg_2'],
    2: ['p_neg_2_prime', 'p_neg_1_prime', 'p1_prime', 'p2_prime'],
}
def main():

    for data_set in ['polypeptide','protein']:
        count_dfs = []
        for label in ['_bg', '']:
            sr1_df = pd.read_csv(
                f'fig6/prot_vs_poly/{data_set}{label}.csv',
                # usecols=['aa'] + SR_COLS[sr]
            )
            # for col in SR_COLS[sr]:
            #     sr1_df[col] /= sr1_df[col].sum()
            # sr1_df.to_csv(f'fig6/invitro_comp/fracs_{sr}_{label}.csv', index=False)
            sr1_df = sr1_df.set_index('aa')
            print(sr1_df)
            sr1_df = sr1_df.transpose()

            sr1_df = sr1_df.reset_index(drop=True)
            count_dfs.append(sr1_df)
        print(count_dfs[0])
        print(count_dfs[1])
        create_comparison_logo_plot(
            [count_dfs[1], count_dfs[0]],
            [f'{data_set}_d', f'{data_set}_bg'],
            16, f'fig6/prot_vs_poly',
            file_name=f'{data_set}_d_vs_bg.svg', amino_acids='ACDEFGHIKLMNPQRSTVWY' + 'X',
            # y_lim=0.03,
            plot_size=5,
            # x_ticks=['', 'p1', '', ''],
            vline=2.5, is_count_df=True
        )

if __name__ == '__main__':
    main()
