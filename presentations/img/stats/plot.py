import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('stats.csv')
print(df)
print(df.columns)
df['tau'] = 100. * df['tau']

dfc = df[['tau', '#visited_cells', '#filtered_cells', '#visited_mcells', '#visited_ccells']].copy()
dfv = df[['tau', '#filtered_vectors', '#checked_pvectors', '#checked_mvectors']].copy()

df_list = [('cells_stats', dfc.melt('tau', var_name='metric', value_name='val')), ('vector_stats', dfv.melt('tau', var_name='metric', value_name='val'))]

print(dfc)
print(dfv)

for (tt, dfm) in df_list:
    g = sns.catplot(x="tau", y="val", hue='metric', data=dfm, kind='point', legend=False,
    palette=sns.color_palette('icefire', 10), scale = 0.5, markers=['o', 'v', '*', 'X', 's', 'D', 'p'])
    plt.xlabel(r'$\mathrm{tau}$', fontsize = 11)
    plt.ylabel(r'$\mathrm{value}$', fontsize = 11)
    plt.legend(loc='best')
    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    plt.grid()  #just add this
    plt.yscale('log')
    # plt.ylim([0,1e5])
    plt.savefig(f"{tt}.png")
    plt.close()
