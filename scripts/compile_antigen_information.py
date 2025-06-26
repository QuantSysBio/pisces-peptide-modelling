""" Functions for compiling antigen information.
"""

import polars as pl

PER_EXON = '/data/John/ALL_FINAL/background_analsis/K562_expr_exon.csv'
B721_TRANSCRIPTOMICS_PATH = '/home/sina.garazhian01/spi_art/rna_seq/trx_721_seq.parquet'
K562_TRANSCRIPTOMICS_PATH = '/home/sina.garazhian01/spi_art/rna_seq/frac_seq.parquet'

def main():
    per_exon_df = pl.read_csv(PER_EXON)
    per_exon_df = per_exon_df.select('GENEID', 'TXNAME', 'exon_id', 'EnsemblProt').rename({
        'GENEID': 'geneID',
        'TXNAME': 'transcriptID',
        'exon_id': 'exonID',
        'EnsemblProt': 'proteinID',
    })
    print(per_exon_df)

    b721_df = pl.read_parquet(B721_TRANSCRIPTOMICS_PATH)
    b721_df = b721_df.with_columns(
        pl.mean_horizontal(
            'A5701', 'E01.C0701', 'A5401', 'F01.C0701', 'D01.C0401',
            'B01.C0401', 'C01.C0401', 'A2902', 'A01.C0401', 'G01.C0701',
            'A5101', 'H01.C0701',
        ).alias('tr_TPM_721'),
    )
    b721_df = b721_df.select('TXNAME', 'tr_TPM_721').rename({'TXNAME': 'transcriptID'})
    b721_df = b721_df.unique()

    k562_df = pl.read_parquet('/home/sina.garazhian01/spi_art/rna_seq/frac_seq.parquet')
    k562_df = k562_df.select(
        'TXNAME', 'tot', 'free', 'S80', 'poly'
    ).rename({
        'TXNAME': 'transcriptID',
        'tot': 'tr_TPM_K562_bulk',
        'free': 'tr_TPM_K562_free',
        'S80': 'tr_TPM_K562_S80',
        'poly': 'tr_TPM_K562_poly',
    })
    k562_df = k562_df.unique()

if __name__ == "__main__":
    main()