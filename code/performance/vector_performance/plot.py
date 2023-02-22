import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

output_dir = "./img/"
# csv = "vec_70_0_0_1000k.csv"
csv = "vec_70_0_0_10000k.csv"

df = pd.read_csv(csv, delimiter=',')
df['k'] = df['k'].astype(str)
k = list( df['k'])
df['d'] = df['d'].astype(float)
df['recall'] = df['recall'].astype(float)
df['time'] = df['time'].astype(float)
# df = df.drop('time', axis=1)


# plt.style.use('classic')
# print(df.head(10))
dtt = df.melt('k', var_name='metric', value_name='value')
# print(dtt.head(10))

g = sns.catplot(x="k", y="value", hue='metric', data=dtt, kind='point', order=k, 
legend=False, markers=['o', 'v', 's'], scale = 0.5, aspect = 2, palette=sns.color_palette('bright', 50))

plt.xticks([1, 500, 1000, 5000, 10000], ["1", "500", "1000", "5000", "10000"])
plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
plt.legend(loc='best')
plt.xticks(fontsize = 11)
plt.yticks(fontsize = 11)

plt.savefig(f"qvec_perf_k=10000.png")
plt.close()



