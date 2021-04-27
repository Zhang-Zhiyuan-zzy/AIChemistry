'''
-*- coding: utf-8 -*-
@Time    : 2021/4/21 17:33
@Author  : Zhang Zhiyuan
@File    : clusters.py
@Software: PyCharm
'''
import os

import numpy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.datasets as data
from sklearn.manifold import TSNE
import pandas as pd
import hdbscan

# Setting parameter
path_root = os.path.join(os.path.dirname(__file__), '..')
path_excel_read = os.path.join(path_root, 'results', 'crystal_system1.xlsx')
path_cluster_result = os.path.join(path_root, 'results', 'crystal_system_results.xlsx')
elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']
c_atoms = ['O', 'N', 'S', 'P']






# Getting data

with pd.ExcelWriter(path_cluster_result) as writer:
    for ele in elements:
        print(ele)
        df = pd.read_excel(path_excel_read, sheet_name=ele, engine='openpyxl', index_col=0).iloc[:, 1:14]
        df_norm = (df-df.mean()) / df.std()
        print(df_norm)
        df_results = pd.read_excel(path_excel_read, sheet_name=ele, engine='openpyxl', index_col=0).iloc[:, 0:14]
        projection = TSNE().fit_transform(df_norm)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=int(df.shape[0] ** 0.5), min_samples=1)
        print(ele)
        clusterer.fit(df)
        df_results['labels_'] = clusterer.labels_
        df_results['probabilities_'] = clusterer.probabilities_
        df_results.to_excel(writer, sheet_name=ele, engine='openpyxl')
        print(ele)
        color_palette = sns.color_palette('Paired', clusterer.labels_.max() + 1)
        cluster_colors = [color_palette[x] if x >= 0
                          else (0.5, 0.5, 0.5)
                          for x in clusterer.labels_]
        cluster_member_colors = [sns.desaturate(x, p) for x, p in
                                 zip(cluster_colors, clusterer.probabilities_)]
        plt.scatter(*projection.T, s=60, linewidth=0, c=cluster_member_colors, alpha=0.25)
        plt.plot()
        plt.show()
        plt.close()


'''
color_palette = sns.color_palette('Paired', clusterer.labels_.max() + 1)
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in clusterer.labels_]
cluster_member_colors = [sns.desaturate(x, p) for x, p in
                         zip(cluster_colors, clusterer.probabilities_)]
plt.scatter(*projection.T, s=60, linewidth=0, c=cluster_member_colors, alpha=0.25)
plt.plot()
plt.show()
'''






