import pandas as pd

x = pd.read_csv('fig6/umap_returns/per_allele_data.csv')
x = x.rename(columns={'cisFreq': 'resolved spliced frequency', 'crypFreq': 'resolved cryptic frequency'})
y = pd.read_csv('/data/John/ALL_FINAL/scripts/fig4/per_allele_data.csv')
y = y.rename(columns={'cisFreq': 'unique spliced frequency', 'crypFreq': 'unique cryptic frequency'})
print(x.columns)
print(y.columns)

z = pd.merge(y, x[['allele', 'resolved spliced frequency', 'resolved cryptic frequency']])
z = z[[
    'allele', 'pca1', 'pca2', 'kmeans_res',
    'unique spliced frequency',  'unique cryptic frequency',
    'resolved spliced frequency', 'resolved cryptic frequency'
]]
z = z.rename(columns={'allele': 'dataset'})
z['dataset'] = z['dataset'].apply(
    lambda x : x.replace('Abelin', 'B721.221').replace('Sarkizova', 'B721.221').replace('K562_MMlab', 'K562')
)
for col in [
    'pca1', 'pca2',
    'unique spliced frequency',  'unique cryptic frequency',
    'resolved spliced frequency', 'resolved cryptic frequency'
]:
    z[col] = z[col].apply(lambda x : round(x, 2))
print(z)
z.to_csv('PISCES_Table_S3.csv', index=False)
