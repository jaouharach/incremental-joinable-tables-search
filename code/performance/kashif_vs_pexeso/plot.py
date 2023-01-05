# evaluate query time for k = {10, 100, 1.000, 10.000}

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# set text font
plt.rcParams.update({
                    "text.usetex": True,
                    "font.family": "serif",
                    "font.serif": "Computer Modern",
                    "savefig.dpi": 130})
plt.rc('axes', axisbelow=True)
plt.rcParams["figure.figsize"] = (6,5)

k_settings = ["1", "3", "5", "10", "30", "50", "100", "300", "500", "1000"]
fout = "recall_diff_settings.csv"
gt_file = "pexeso_results_tau=6%_T=60%.csv"
results_file = f"kashif_querytime.csv"

dfp = pd.read_csv(gt_file, delimiter=",") 
dfk = pd.read_csv(results_file, delimiter=",") 
dft = dfk.groupby(['k']).agg(kashif = ('time','mean')).reset_index()
dft.insert(2, "pexeso", [dfp['time'].mean() for k in k_settings], True)
dft = dft.melt('k', var_name='algo', value_name='avg. querytime')
print(dft)

fig, ax = plt.subplots()
ax.yaxis.grid(True, color="lightgray")
g = sns.lineplot(data=dft, x='k', y='avg. querytime', hue='algo', legend=False, palette=['cornflowerblue', 'coral'], marker='v')

plt.legend(labels=['kashif', 'pexeso'], fontsize = 13)
plt.ylabel('avg. querytime', fontsize = 13)
plt.xlabel('k', fontsize = 13)
plt.savefig(f"querytime_kashif_vs_pexeso.png")
plt.yscale('log')
plt.show()
