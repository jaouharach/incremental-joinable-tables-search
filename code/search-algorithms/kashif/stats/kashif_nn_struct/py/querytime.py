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

outdir = '../img/'
csv_file = '../csv/querytime_3struct.csv'

plt.rcParams["figure.figsize"] = (8,7)
df = pd.read_csv(csv_file)
print(df)

qt = df.groupby(
    ['k','nb_threads','knn data structure']
).agg(
    querytime = ('querytime','mean'),
).reset_index()


qt['k'] = qt['k'].apply(lambda x: format_nb(x))

print(qt)

flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=qt, x='k', y='querytime', hue='knn data structure',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v', '*'])

g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.16)
plt.xlabel(r'$\mathrm{k}$', fontsize = 12)
plt.ylabel(r'$\mathrm{mean\ query\ time\ (sec)}$', fontsize = 12)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.yscale("log")
# plt.title("Kashif Average insert NN count \n (10 queries, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors) \n")
plt.savefig(f"{outdir}querytime_3struct(log).png")
plt.legend(loc='best', title="-")
plt.close()



