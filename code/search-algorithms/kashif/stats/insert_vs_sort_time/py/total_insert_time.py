import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = '../img/'
csv_file = '../csv/nn_insert_time.csv'

plt.rcParams["figure.figsize"] = (8,7)
df = pd.read_csv(csv_file)

df['k'] = df['k'].astype('str')
df['total_insert_time'] = (df['total_insert_sort_time'] - df['total_sort_time']).abs()
df = df.drop(['total_insert_sort_time'], axis=1)

dfm = df.melt(['k', '#th', 'type'], var_name='time', value_name='val')
print(dfm)

flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=dfm, x='k', y='val', hue='time',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v'])

# g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.15)
plt.xlabel(r'$\mathrm{k}$', fontsize = 14)
plt.ylabel(r'$\mathrm{time\ (sec)}$', fontsize = 14)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.savefig(f"{outdir}nn_insert_vs_sort_time.png")
plt.legend(loc='best', title="-")
plt.close()



