# evaluate query time for k = {10, 100, 1.000, 10.000}
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import math 


output_dir = "./img/"
csv_file = "./csv/rec4.csv"


df = pd.read_csv(csv_file, delimiter=',')
plt.rc('axes', axisbelow=True)
df['k'] = df['k'].astype(str)
df['recall'] = df['recall'].astype(float)
# df['precision'] = df['precision'].astype(float)
df['time'] = df['time'].astype(float)

# # df['time'] = df['time'].apply(lambda x: x*(10**3))

# querytime
querytime_df = df.groupby(
    ['k']
).agg(
    time = ('time','mean'),
).reset_index()

# recall
recall_df = df.groupby(
    ['k']
).agg(
    recall = ('recall','mean'),
).reset_index()

# precision
precision_df = df.groupby(
    ['k']
).agg(
    precision = ('precision','mean'),
).reset_index()
precision_df = precision_df.round(2)

querytime_df['k'] = querytime_df['k'].astype(int)
querytime_df = querytime_df.sort_values('k')
querytime_df['k'] = querytime_df['k'].astype(str)
print(querytime_df.head(13))

recall_df['k'] = recall_df['k'].astype(int)
recall_df = recall_df.sort_values('k')
recall_df['k'] = recall_df['k'].astype(str)
print(recall_df.head(13))

df = pd.merge(recall_df, querytime_df, on='k')
print(df.head(13))


# BAR PLOT
fig, ax1 = plt.subplots()
ax1.yaxis.grid(True, color="lightgray")
ax2 = ax1.twinx()
ax2.yaxis.grid(True, linestyle='dashed', color="lightgray")
sns.barplot(data = df,x="k", y="recall", ax=ax1, ci=None, order=df['k'], color="cornflowerblue")
sns.lineplot(data = df["time"], sort = False, ax=ax2, color="red")
# ax1.bar_label(ax1.containers[0], labels=df['precision'], color="slategray")
# ax1.margins(y=0.1)

ax1.set_ylabel(r'$\mathrm{avg.\ recall}$')
ax2.set_ylabel(r'$\mathrm{avg. time \ (sec)}$')
# ax1.set_ylabel(r'$\mathrm{recall}$')
# ax2.set_ylabel(r'$\mathrm{time}$')
ax1.set_ylim([0, 1.5])
ax2.set_yscale('linear')
# ax2.legend(loc=0)
# ax2.set_xlabel(r'$\mathrm{k}$')

    # # Set these based on your column counts
    # columncounts = [25 for i in range(k_count)]

    # # Maximum bar width is 1. Normalise counts to be in the interval 0-1. Need to supply a maximum possible count here as maxwidth
    # def normaliseCounts(widths,maxwidth):
    #     widths = np.array(widths)/float(maxwidth)
    #     return widths

    # widthbars = normaliseCounts(columncounts,100)

    # # Loop over the bars, and adjust the width (and position, to keep the bar centred)
    # for bar,newwidth in zip(ax2.patches,widthbars):
    #     x = bar.get_x()
    #     width = bar.get_width()
    #     centre = x+width/2.

    #     bar.set_x(centre-newwidth/2.)
    #     bar.set_width(newwidth)

plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)
# plt.legend(loc='upper left')
# plt.title("Kashif: mean recall (10 query columns of size [50 - 100])")
plt.savefig(f"{output_dir}/rec4.png")
plt.close()



