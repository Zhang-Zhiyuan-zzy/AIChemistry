'''
-*- coding: utf-8 -*-
@Time    : 2021/4/17 10:04
@Author  : Zhang Zhiyuan
@File    : visualization of chelating.py
@Software: PyCharm
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Setting parameters

path_root = os.path.join(os.path.dirname(__file__), '..')
path_data_reader = os.path.join(path_root, 'results', 'chelating_ring.xlsx')
sheet_names = ['Na2', 'K2', 'Rb2', 'Cs2', 'Mg2', 'Ca2', 'Sr2', 'Ba2']
#sheet_names = ['Cs2', 'Sr2']
df_list = []
for sheet_name in sheet_names:
    locals()['df_{}'.format(sheet_name)] = pd.read_excel(path_data_reader, sheet_name=sheet_name, index_col=0)
    if sheet_name != 'K2':
        locals()['df_{}'.format(sheet_name)]['metal'] = sheet_name[0:2]
    else:
        locals()['df_{}'.format(sheet_name)]['metal'] = 'K'
    df_list.append(locals()['df_{}'.format(sheet_name)])

df1 =pd.concat(df_list, axis=0)

df1['rbl0'] = df1['bl0'] / df1['ibl0']
df1['rbl1'] = df1['bl1'] / df1['ibl1']

df2 = df1[df1['shortest_path01'] == 2]
sample = []
for name, df_g in df2.groupby('shortest_path01'):
    sample.append(df_g.sample(n=266))

df = pd.concat(sample, axis=0)

df['shortest_path01'] = df['shortest_path01'] + 1

sns.set()
sns.scatterplot(x='rbl0', y='rbl1', hue='shortest_path01', data=df,)
plt.show()


