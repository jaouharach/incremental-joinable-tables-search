import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = '../img/'
csv_file = '../csv/insert_and_get_nn_count.csv'

plt.rcParams["figure.figsize"] = (8,7)
df_ref = pd.read_csv(csv_file)
print(df_ref)
# exit(1)
# df_ref['k'] = df_ref['k'].astype('str')
df = df_ref.groupby(
    ['k',]
).agg(
    insert_nn_count = ('#insert_nn','mean'),
    get_ith_nn_count = ('#get_ith_nn','mean'),
    get_kth_nn_count = ('#get_kth_nn','mean'),
).reset_index()

print(df)
dfm = df.melt(['k'], var_name='metric', value_name='avg. value')
# dfm['k'] = dfm['k'].astype('str')


print(dfm)

flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=dfm, x='k', y='avg. value', hue='metric',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v', '*'])

# g.set_xticklabels(rotation=30)
# plt.subplots_adjust(bottom=0.1, top=0.85, left=0.26)
plt.xlabel(r'$\mathrm{k}$', fontsize = 12)
plt.ylabel(r'$\mathrm{avg. value}$', fontsize = 12)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.yscale("log")
# plt.title("Kashif Average insert NN count \n (10 queries, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors) \n")
plt.savefig(f"{outdir}insert_get_nn_avg_count.png")
plt.legend(loc='best', title="-")
plt.close()



