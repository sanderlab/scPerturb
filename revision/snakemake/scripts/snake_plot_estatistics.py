import matplotlib.pyplot as pl
import seaborn as sns
import numpy as np
import pandas as pd

# export table
df = pd.read_csv(snakemake.input[0], index_col=0)

# plot result
with sns.axes_style('whitegrid'):
    fig, ax = pl.subplots(1,1, figsize=[6,4], dpi=100)
sns.scatterplot(data=df[df.index!='control'], x='log10_edist', y='pvalue_adj', hue='significant_adj', palette={True: 'tab:green', False: 'tab:red'})
sig = np.sum(df['significant_adj'])
total = len(df)-1  # (removes control)
ax.set_xticks([0,1,2])
# ax.set_xticklabels([r'$10^0$', r'$10^1$', r'$10^2$'])
ax.set_xticklabels([1, 10, 100])
ax.set_title(f'E-tests for {snakemake.wildcards.dataset}\n{sig}/{total} are significant (pv<0.05)')
ax.set_ylabel('Adjusted p-value')
ax.set_xlabel('E-distance to unperturbed')
ax.legend(title='Significant')
ax.axhline(0.05, c='r', linestyle='--')
small = df[(df['significant_adj']) & (df.index!='control')].sort_values('edist').iloc[0]
big = df[(~df['significant_adj']) & (df.index!='control')].sort_values('edist').iloc[-1]
ax.annotate(f'E-distance\n{np.round(small.edist,2)}', xy=(small.log10_edist, small.pvalue_adj),  xycoords='data',
            xytext=(0.5, small.pvalue_adj), textcoords='data', fontsize=10,
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='center',
            )
ax.annotate(f'E-distance\n{np.round(big.edist,2)}', xy=(big.log10_edist, big.pvalue_adj),  xycoords='data',
            xytext=(2, big.pvalue_adj+0.2), textcoords='data', fontsize=10,
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='center',
            )
pl.savefig(snakemake.output[0], bbox_inches='tight', dpi=120)
pl.close()