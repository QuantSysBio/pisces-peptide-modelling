import polars as pl

TEST_PATH = 'datasets/250515_piscesDB.parquet'

def test_pisces_db():
    pisces_df = pl.read_parquet(TEST_PATH)
    print(pisces_df.columns)
    # Check that peptide is not duplicated
    assert pisces_df['peptide'].n_unique() == pisces_df.shape[0]
    assert pisces_df.filter(pl.col('stratum').eq('error')).shape[0] == 0


    for stratum in ['canonical', 'fiveUTR', 'threeUTR', 'intronic', 'intergenic',]:
        pisces_df = pisces_df.with_columns(
            pl.struct([f'{stratum}_nProteins', f'{stratum}_Proteins']).map_elements(
                lambda x: x[f'{stratum}_nProteins'] == len(x[f'{stratum}_Proteins']),
                return_dtype=pl.Boolean,
            ).alias(f'{stratum}TestEquality'),
            pl.col(f'{stratum}_Proteins').map_elements(
                lambda x: min([a != '' and a is not None for a in x]) if len(x) > 0 else True,
                return_dtype=pl.Boolean,
            ).alias(f'{stratum}TestNoNulls'),
        )
        assert pisces_df[f'{stratum}TestEquality'].min()
        assert pisces_df[f'{stratum}TestNoNulls'].min()

if __name__ == "__main__":
    test_pisces_db()
