### Objective
# The objective of this script is to read data from NOAA ISD stations into useful data formats.

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', None)
data_dir = "data/"      # Directory where raw data is stored
data_cols = ['STATION', 'DATE', 'TMP', 'SLP', 'QUALITY_CONTROL']
data_date = [20200601, 20200605]
idx = [-48, -1]

data_list = []
for file in os.listdir(data_dir):
    if file.endswith(".csv"):
        df = pd.read_csv(os.path.join(data_dir, file), low_memory=False)
        # data_list.append(pd.read_csv(os.path.join(data_dir, file), usecols=data_cols, low_memory=False))
        data_list.append(pd.read_csv(os.path.join(data_dir, file), low_memory=False))

# print(df.columns)

data = pd.concat(data_list)
data[data_cols[2]] = data[data_cols[2]].str.replace(',', '.')
data[data_cols[2]] = data[data_cols[2]].str.replace('+', '')
data[data_cols[2]] = data[data_cols[2]].str.replace('A', '0')
data[data_cols[2]] = data[data_cols[2]].astype("float")
data[data_cols[2]] = data[data_cols[2]].astype("float")

data.loc[data[data_cols[2]] >= 200.0, [data_cols[2]]] = np.nan

print(data.iloc[-2:-1])

# fig, ax = plt.subplots(figsize=[8,6])
# plt.plot(data[data_cols[1]][0:200], data[data_cols[2]][0:200])
# plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor")
# ax.set_xticks(ax.get_xticks()[::12])
# plt.show()