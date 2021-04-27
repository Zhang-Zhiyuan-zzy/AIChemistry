'''
-*- coding: utf-8 -*-
@Time    : 2021/4/23 16:32
@Author  : Zhang Zhiyuan
@File    : Sr_system.py
@Software: PyCharm
'''
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Setting parameters
path_root = os.path.join(os.path.dirname(__file__), '..')
path_excel_reader = os.path.join(path_root, 'results', 'crystal_system_results.xlsx')
list_bl = ['Orbl', 'Nrbl', 'Srbl', 'Prbl']
df = pd.read_excel(path_excel_reader, sheet_name='Sr', engine='openpyxl')
df['bl'] = 0
for i, row in df.iterrows():
    list_uzero = []
    for n_bl in list_bl:
        if not not row[n_bl]:
            list_uzero.append(row[n_bl])
    if not not list_uzero:
        df.loc[i, 'bl'] = np.mean(list_uzero)
df.drop(df[df['bl'] == 0].index, inplace=True)
print(len(df))
sns.set()
sns.scatterplot(x='labels_', y='rcm', hue='bl', data=df)
plt.show()

