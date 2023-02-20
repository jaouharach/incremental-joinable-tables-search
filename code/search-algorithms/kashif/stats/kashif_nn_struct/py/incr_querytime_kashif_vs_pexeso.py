# measuring mean query time at diferent k positions (time measure while returning incremental results)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def format_nb(nb):
    nb = int(nb)
    if(nb < 1000):
        return str(nb)
    elif (nb < 1000000):
        return str(nb//1000) + "k"
    else:
        return str(nb//1000000) + "m"

outdir = '../img/'
kcsv_file = '../csv/querytime_mmheap_incr.csv'
pcsv_file = '../csv/pexeso_querytime_all.csv'
krecall_file = "../../recall/recall.csv"

plt.rcParams["figure.figsize"] = (8,7)
kdf = pd.read_csv(kcsv_file)
krdf = pd.read_csv(krecall_file)
pdf = pd.read_csv(pcsv_file)

# kashif mean query time
kqt = kdf.groupby(
    ['k','nb_threads','knn data structure']
).agg(
    querytime = ('querytime','mean'),
).reset_index()


# kashif mean query recall
krecall = krdf.groupby(
    ['k']
).agg(
    recall = ('recall','mean'),
).reset_index()


# kashif mean query precision
kprecision = krdf.groupby(
    ['k']
).agg(
    precision = ('precision','mean'),
).reset_index()



# pexeso mean query time
pqt = pdf.groupby(
    ['tau', 'join_threashold']
).agg(
    mean_querytime = ('querytime','mean'),
).reset_index()

kdf['k'] = kdf['k'].astype(str)
krecall['k'] = krecall['k'].astype(str)


krecall['k'] = krecall['k'].apply(lambda x: format_nb(x))
kprecision['k'] = kprecision['k'].apply(lambda x: format_nb(x))
kdf['k'] = kdf['k'].apply(lambda x: format_nb(x))
kqt['k'] = kqt['k'].apply(lambda x: format_nb(x))

# only keep results of kashif using min max heap
# kqt = kqt.drop(kqt[kqt['knn data structure'] == 'kashif (sorted-arr)'].index)
pqt = pqt.drop(pqt[pqt['tau'] == 0.08].index)

# print(krecall)
# print(kprecision)
# print(kqt)
# print(pqt)

fig, ax1 = plt.subplots()
ax1.yaxis.grid(True, color="lightgray")
ax2 = ax1.twinx()
ax2.yaxis.grid(True, linestyle='dashed', color="lightgray")



print(krecall)
# plot kashif recall
sns.barplot(data=krecall, x=krecall['k'], y="recall", ax = ax1, color="lightgreen")
# plot kashif  querytime
g = sns.lineplot(data=kqt["querytime"], sort=False, ax=ax2, marker='o', label="Kashif (minmax-heap)", linewidth=2, color="royalblue")


# plot all pexeso lines
# ax2.axhline(pqt['mean_querytime'][2], c=f'b', label=f'pexeso (tau = 6%)', ls='--')


# plot all pexeso lines
colors =  ['pink', 'tomato', 'red', 'black']
for i, tau in enumerate(pqt['tau'].values):
    ax2.axhline(pqt['mean_querytime'][i], c=f'{colors[i]}', label=f'pexeso (tau = {tau*100}%)', ls='--')

# set x and y lables
ax1.set_ylabel("Average recall", fontsize=13)
ax1.set_ylim([0, 1])
ax2.set_ylabel("Mean query time", fontsize=13)
ax1.set_xlabel(r'$\mathrm{k_{i}\ (k_{max}\ =\ 1m)}$', fontsize = 15)

plt.grid()  #just add this
plt.legend(loc='best')
plt.tight_layout(pad=2.3)

# plt.yscale("log")
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10),
#           ncol=2, fancybox=True, shadow=True)

# plt.title("Kashif Average insert NN count \n (10 queries, query size = 100, dataset = 100k tables, 490k cols, ~5M vectors) \n")
plt.savefig(f"{outdir}querytime_recall_kashif_mmheap_incr_vs_pexeso.png")
plt.close()



