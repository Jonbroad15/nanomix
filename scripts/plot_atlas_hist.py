import pandas as pd
import matplotlib.pyplot as plt

plt.figure(figsize=(16,8), dpi=100)
df = pd.read_csv('bermanAtlas.tsv', sep='\t')
cell_types = ['CD4T-cells EPIC','CD8T-cells EPIC','B-cells EPIC','Erythrocyte progenitors','Cortical neurons','Monocytes EPIC','Prostate']
df[cell_types].hist()
plt.savefig('bermanAtlasHist.png')
