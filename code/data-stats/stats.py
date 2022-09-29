import pandas as pd
import re
import numpy as np

csv_file = "ncols_horizontal_tables.csv"
csv_file = "ncols_vertical_tables.csv"

cols = pd.read_csv(csv_file, delimiter=":") 
cols.columns=['idx', 'Label', 'value']

cols.drop('idx',  axis=1, inplace=True)
for index, row in cols.iterrows():
    cols.at[index,'Label'] =  int(re.search('[0-9]+', row['Label']).group())


pd.to_numeric(cols['Label'])

print(np.dot(cols['Label'],  cols['value']))