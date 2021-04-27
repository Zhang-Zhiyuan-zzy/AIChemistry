#!/usr/bin/env python
# _*_ coding:utf-8 _*_
"""
Time    :   2021/3/15 11:43
Author  :   Cheng Min
File    :   statistics_atoms.py
Software:   PyCharm
"""

"""获取包含下列元素的晶体名称及其化学式"""

import os

import pandas as pd

path_root = os.path.join(os.path.dirname(__file__), '..')
path_module = os.path.join(path_root, 'public_module')
from public_module import statistics_from_crystal

# 参数设置
list_elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']  # 考虑的元素
path_mol2_file_dir = os.path.join(path_root, 'results/mol2_file')  # 储存*.mol2文件的文件夹
path_save_result = os.path.join(path_root, 'results/files_includ_atoms_neighbors.xlsx')

# 根据identifier判断其neighbors
# 判断保存结果的文件夹是否存在，若不存在，则创建
if not os.path.exists(os.path.join(path_save_result, '..')):
    os.mkdir(os.path.dirname(path_save_result))

# 保存结果
with pd.ExcelWriter(path_save_result) as writer:
    for element in list_elements:
        path_current_mol_file_dir = os.path.join(path_mol2_file_dir, element)
        test = statistics_from_crystal.StatisticsFromCrystal(list_elements)
        locals()['dict_result_{}'.format(element)] = test.get_neighbor_atoms(element, path_current_mol_file_dir)
        Ser = pd.Series(locals()['dict_result_{}'.format(element)])
        Ser.name = 'count'
        Ser.index.name = 'neighbors'
        Ser.to_excel(writer, sheet_name=element, engine='openpyxl')

