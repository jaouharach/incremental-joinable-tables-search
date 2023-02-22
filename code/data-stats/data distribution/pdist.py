import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv("pdist.csv")

print(df.head(3))
df["pair"] = df["src"] + "," + df[" dest"]

counts, bin_edges = np.histogram(df[" d"], bins=10, density=True)
print("counts : edges")

for c, e in zip(counts, bin_edges):
    print(f"c = {c} \t:\t e = {e}")

pdf = counts / sum(counts)
cdf = np.cumsum(pdf)

plt.plot(bin_edges[1:], pdf)
plt.plot(bin_edges[1:], cdf)
plt.xlabel("")
plt.ylabel("cdf and pdf")
plt.title("Pairwise distance")
plt.show()


plt.bar(df.index, df[" d"])
 
plt.xlabel("")
plt.ylabel("d(.,.)")
plt.title("Pairwise distance")
plt.show()