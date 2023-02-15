# measuring kashif memory usage when using diff data structures to store knns
# 10 queries column of size 100


# amount of memory required for running knn search for a column of size |Q| =
# (24 * |Q| * k ) + (#worker_threads * k * 72)
# size of the results stored in the coordinator 24 bytes per nn




import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def format_nb(nb):
    if(nb < 1000):
        return str(nb)
    elif (nb < 1000000):
        return str(nb//1000) + "k"
    else:
        return str(nb//1000000) + "m"

def kb_to_gb(nb):
    return nb / (1024 * 1024)

outdir = '../img/'
csv_file = '../csv/memory_usage.csv'

plt.rcParams["figure.figsize"] = (8,7)
df = pd.read_csv(csv_file)
print(df)

pmu = df.groupby(
    ['k','nb_threads','knn data structure']
).agg(
    peak_momory_usage_kb = ('peak_momory_usage_kb','mean'),
).reset_index()


pmu['k'] = pmu['k'].apply(lambda x: format_nb(x))
pmu['peak_momory_usage_kb'] = pmu['peak_momory_usage_kb'].apply(lambda x: kb_to_gb(x))


print(pmu)


flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=pmu, x='k', y='peak_momory_usage_kb', hue='knn data structure',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v', '*'])

g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.16)
plt.xlabel(r'$\mathrm{k}$', fontsize = 12)
plt.ylabel(r'$\mathrm{peak\ memory\ usage\ (gb)}$', fontsize = 12)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
# plt.yscale("log")
# plt.title("Kashif Average insert NN count \n (10 queries, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors) \n")
plt.savefig(f"{outdir}memory_usage.png")
plt.legend(loc='best', title="-")
plt.close()



