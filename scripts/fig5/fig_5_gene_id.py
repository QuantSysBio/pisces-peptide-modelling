import polars as pl

ALL_CRYPTIC_STRATA = [
    'fiveUTR', 'threeUTR', 'CDS_frameshift',
    'lncRNA', 'intronic', 'intergenic',
]
STRATA_WITH_GENES = [
    'fiveUTR', 'threeUTR', 'CDS_frameshift', 'intronic'
]

def main():
    cryp_df = pl.read_csv("output/crypticK562_241129/unique_peps_scored.csv")
    can_df = pl.read_csv("output/canonicalK562_241129/unique_peps_scored.csv")

    cryp_df = cryp_df.filter(pl.col('label').eq(1))
    can_df = can_df.filter(pl.col('label').eq(1))
    print(can_df.columns)
    print(cryp_df.columns)

    cryp_df = cryp_df.select(['peptide', 'proteinID', 'stratum'])
    can_df = can_df.select(['peptide', 'proteinID'])

    gene_df = pl.read_parquet('antigen/transcript/gene.parquet')
    prot_df = pl.read_parquet('antigen/transcript/protein.parquet')
    prot_df = prot_df.select(['geneID', 'proteinID'])
    per_exon_df = gene_df.select(['geneID', 'exonID']).explode('exonID')
    per_transcript_df = gene_df.select(['geneID', 'transcriptID',]).explode('transcriptID')
    sub_dfs = []

    for stratum in STRATA_WITH_GENES:
        sub_df = cryp_df.filter(pl.col('stratum') == ALL_CRYPTIC_STRATA.index(stratum))
        if stratum in ['fiveUTR', 'threeUTR']:
            sub_df = sub_df.with_columns(
                pl.col('proteinID').str.split('_').list.first().alias('exonID')
            )
            print(stratum, sub_df['peptide'].n_unique())
            sub_df = sub_df.join(per_exon_df, how='inner', on='exonID')
            print(stratum, sub_df['peptide'].n_unique())
            sub_df = sub_df.select(['peptide', 'geneID'])
        else:
            if stratum in ['CDS_frameshift']:
                sub_df = sub_df.with_columns(
                    pl.col('proteinID').str.split('|').list.first().alias('transcriptID')
                )
            else:
                sub_df = sub_df.with_columns(
                    pl.col('proteinID').str.split('_').list.first().alias('transcriptID')
                )
            print(stratum, sub_df['peptide'].n_unique())
            sub_df = sub_df.join(per_transcript_df, how='inner', on='transcriptID')
            print(stratum, sub_df['peptide'].n_unique())
            sub_df = sub_df.select(['peptide', 'geneID'])
        sub_dfs.append(sub_df)
    print(cryp_df['peptide'].n_unique())
    cryp_df = pl.concat(sub_dfs)
    cryp_df = cryp_df.unique()
    print(cryp_df['peptide'].n_unique())
    print(cryp_df)
    cryp_df = cryp_df.group_by('geneID').agg(pl.col('peptide').alias('crypticPeptides'))

    can_df = can_df.join(prot_df, how='inner', on='proteinID')
    can_df = can_df.group_by('geneID').agg(pl.col('peptide').alias('canonicalPeptides'))

    total_df = cryp_df.join(can_df, how='inner', on='geneID').with_columns(
        pl.col('crypticPeptides').list.lengths().alias('crypticPeptideCount'),
        pl.col('canonicalPeptides').list.lengths().alias('canonicalPeptideCount')
    )
    total_df = total_df.sort(['crypticPeptideCount', 'canonicalPeptideCount'], descending=True)
    print(total_df.head(10))

if __name__ == '__main__':
    main()