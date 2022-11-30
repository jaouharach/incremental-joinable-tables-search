# evaluate query time for k = {10, 100, 1.000, 10.000}
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import math 


output_dir = "./img/"
csv = [("seq.csv", "seq"), ("th4.csv", "th4"), ("th8.csv", "th8"), ("th16.csv", "th16"), ("th32.csv", "th32"), ("th64.csv", "th64")]

dtt = pd.DataFrame()
dtt['k'] = ["1", "3", "5", "10", "30", "50", "100", "300", "500", "1000"]
ktt = ["1", "3", "5", "10", "30", "50", "100", "300", "500", "1000"]
for c, t in csv:
    df = pd.read_csv(c, delimiter=',')
    plt.rc('axes', axisbelow=True)
    df['k'] = df['k'].astype(int)
    df['time'] = df['time'].astype(float)

    querytime = df.groupby(['k']).agg(time = ('time','mean'),).reset_index()
    querytime.rename(columns={'time': t}, inplace=True)
    querytime_seq = querytime.sort_values('k')
    querytime['k'] = querytime['k'].astype(str)

    dtt = pd.merge(dtt, querytime, on='k')


# plt.style.use('classic')
print(dtt.head(13))
dtt = dtt.melt('k', var_name='algo', value_name='avg. querytime (sec)')
print(dtt.head(13))
sns.color_palette("Set2")
g = sns.catplot(x="k", y="avg. querytime (sec)", hue='algo', data=dtt, kind='point', order=ktt, legend=False, 
palette=sns.color_palette('icefire', 10), markers=['o', 'v', '*', 'x', '|', '3'], scale = 0.5)

plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.yscale("log")
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
plt.grid()  #just add this

plt.savefig(f"th_ylog.png")
plt.close()



