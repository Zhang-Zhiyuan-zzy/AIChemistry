#!/usr/bin/env python
# _*_ coding:utf-8 _*_
"""
Time    :   2021/3/27 14:12
Author  :   Cheng Min
File    :   statistics_neighbor_groups.py
Software:   PyCharm
"""

import os

import pandas as pd

path_root = os.path.join(os.path.dirname(__file__), '..')
path_module = os.path.join(path_root, 'public_module')


# 参数设置
list_elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']  # 考虑的元素
path_con = os.path.join(path_root, 'input_data/function_group')
path_mols_dir = os.path.join(path_root, 'results/mol2_file')
path_save_data = os.path.join(path_root, 'results/statistics_neighbour_groups/statistics_neighbour_groups.xlsx')

# 统计
for element in list_elements:
    path_current_mol = os.path.join(path_mols_dir, element)
    test = statistics_from_crystal.StatisticsFromCrystal(list_elements)
    locals()['dict_result_{}'.format(element)], list_con_names = test.get_neighbor_function_groups(path_current_mol, path_con, element)

# 判断结果保存文件夹是否存在，若不存在，则创建
path_save_data_dir = os.path.dirname(path_save_data)
if not os.path.exists(path_save_data_dir):
    os.makedirs(path_save_data_dir)

# 保存结果
with pd.ExcelWriter(path_save_data) as writer:
    for element in list_elements:
        df = pd.DataFrame(locals()['dict_result_{}'.format(element)]).T
        df.columns = list_con_names
        df.index.name = 'Refcodes'
        df.to_excel(writer, sheet_name=element)