import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = '../img/'
csv_file = '../csv/sort.csv'

plt.rcParams["figure.figsize"] = (8,7)
df = pd.read_csv(csv_file)

df['k'] = df['k'].astype('str')
# df = df.drop(['total_insert_sort_time'], axis=1)

dfm = df.melt(['k', '#th', 'Kashif-version'], var_name='time', value_name='val')
print(dfm)


dfsa = dfm.copy()
dfsa = dfsa.drop(dfsa[dfsa.time == 'total_insert_time'].index)
print(dfsa)

# exit(1)


flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=dfsa, x='k', y='val', hue='Kashif-version',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v'])

# g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.15, top=0.9, left=0.22)
plt.xlabel(r'$\mathrm{k}$', fontsize = 14)
plt.ylabel(r'$\mathrm{time\ (sec)}$', fontsize = 14)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.title("Kashif Total Query time \n (1 query, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors)");
plt.savefig(f"{outdir}sort_running_vs sort_all_nns(1).png")
plt.legend(loc='best', title="-")
plt.close()



