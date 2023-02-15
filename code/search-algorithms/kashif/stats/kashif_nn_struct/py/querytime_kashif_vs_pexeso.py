# measuring mean query time at kmax
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
kcsv_file = '../csv/querytime_3struct_kmax.csv'
pcsv_file = '../../performance/csv/pexeso/pexeso_querytime_all.csv'


plt.rcParams["figure.figsize"] = (8,7)
kdf = pd.read_csv(kcsv_file)
pdf = pd.read_csv(pcsv_file)

# kashif mean query time
kqt = kdf.groupby(
    ['k','nb_threads','knn data structure']
).agg(
    querytime = ('querytime','mean'),
).reset_index()

# pexeso mean query time
pqt = pdf.groupby(
    ['tau', 'join_threashold']
).agg(
    mean_querytime = ('querytime','mean'),
).reset_index()

kqt['k'] = kqt['k'].apply(lambda x: format_nb(x))
print(kqt)
print(pqt)


flatui = ["#3498db", "#e74c3c", "#34495e", "#2ecc71"]
g = sns.catplot(data=kqt, x='k', y='querytime', hue='knn data structure',  kind='point', scale = 0.5,
    palette=sns.color_palette(flatui),  markers=['o', 'v', '*'], legend=False)

colors =  ['violet', 'magenta', 'indigo', 'black']

for i, tau in enumerate(pqt['tau'].values):
    g.map(plt.axhline, y=pqt['mean_querytime'][i], ls='--', c=f'{colors[i]}', label=f'pexeso (tau = {tau*100}%)')

g.set_xticklabels(rotation=30)
plt.subplots_adjust(bottom=0.16)
plt.xlabel(r'$\mathrm{k_{max}}$', fontsize = 12)
plt.ylabel(r'$\mathrm{mean\ query\ time\ (sec)}$', fontsize = 12)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this
plt.yscale("log")
plt.legend(loc='best', title="-")

# plt.title("Kashif Average insert NN count \n (10 queries, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors) \n")
plt.savefig(f"{outdir}querytime_kashif_3struct_vs_pexeso(log).png")
plt.close()



