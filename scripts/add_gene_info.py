import polars as pl

def main():
    protein_df = pl.read_parquet(
        '/data/John/pisces_peptide_modelling/antigen/prots/canonical.parquet',
    )
    gene_df = pl.read_parquet(
        '/data/John/pisces_peptide_modelling/antigen/transcript/gene.parquet',
    )
    # print(protein_df['proteomics_K562'].value_counts())
    # print(protein_df['proteomics_B721'].value_counts())
    protein_df = protein_df.drop(['proteomics_K562', 'proteomics_B721'])
    for name, project_name in zip(['K562', 'B721'], ['K562-tryptic-250625', 'B721-fragger-250724']):
        tryp_df = pl.read_csv(
            f'/data/John/ALL_FINAL/projects/pisces_tryptic/{project_name}/outputFolder/canonicalOutput/final.percolatorSeparate.proteins.txt',
            separator='\t',
        )
        tryp_df = tryp_df.filter(pl.col('q-value').lt(0.01) & pl.col('posterior_error_prob').lt(0.1))
        tryp_df = tryp_df.with_columns(pl.col('ProteinId').str.split(' '))
        tryp_df = tryp_df.explode('ProteinId').unique()
        tryp_df = tryp_df.rename({'ProteinId': 'proteinID'}).select(['proteinID']).with_columns(
            pl.lit(1).alias(f'proteomics_{name}')
        )
        # print(tryp_df.filter(pl.col('proteinID').eq('ENSP00000462667.1')))
        protein_df = protein_df.join(
            tryp_df, how='left', on='proteinID',
        )

    protein_df = protein_df.with_columns(
        pl.col('proteomics_K562').fill_null(0).cast(pl.Int32),
        pl.col('proteomics_B721').fill_null(0).cast(pl.Int32),
    )
    print(protein_df.shape, protein_df['geneID'].n_unique())
    gene_proteomics_df = protein_df.group_by('geneID').agg(
        [
            pl.col('proteomics_K562').max().alias('proteomics_K562'),
            pl.col('proteomics_B721').max().alias('proteomics_B721'),
        ]
    )
    print(gene_proteomics_df.filter(pl.col('geneID').eq('ENSG00000270882.3')))
    gene_df = gene_df.join(gene_proteomics_df, how='left', on='geneID')
    gene_df = gene_df.with_columns(
        pl.col('proteomics_K562').fill_null(0).cast(pl.Int32),
        pl.col('proteomics_B721').fill_null(0).cast(pl.Int32),
    )

    print(gene_df['proteomics_K562'].value_counts())
    print(gene_df['proteomics_B721'].value_counts())
    gene_df.write_parquet(
        '/data/John/pisces_peptide_modelling/antigen/prots/gene.parquet',
    )

if __name__ == '__main__':
    main()
