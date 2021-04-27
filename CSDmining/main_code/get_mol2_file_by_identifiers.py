#!/usr/bin/env python
# _*_ coding:utf-8 _*_
"""
Time    :   2021/3/27 12:33
Author  :   Cheng Min
File    :   get_mol2_file_by_identifiers.py
Software:   PyCharm
"""

"""获取去除溶剂之后的*.mol2文件，并将其保存到指定文件夹中"""

import os

import pandas as pd

path_root = os.path.join(os.path.dirname(__file__), '..')
path_module = os.path.join(path_root, 'public_module')
from public_module import statistics_from_crystal

# 参数设置
list_elements  = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']  # 考虑的元素
path_file_excel = os.path.join(path_root, 'results/files_include_specific_atoms.xlsx')  # 包含指定原子的文件
path_save_result_dir = os.path.join(path_root, 'results/mol2_file')  # 储存*.mol2文件的文件夹

# 取出包含各个原子的identifier，并赋值给变量list_identifier_Na，list_identifier_K，...，list_identifier_Ba
for element in list_elements:
    arr_identifier_temp = pd.read_excel(path_file_excel, sheet_name=element, engine='openpyxl').iloc[:, 2].values.flatten()
    locals()['list_identifier_{}'.format(element)] = arr_identifier_temp

for element in list_elements:
    test = statistics_from_crystal.StatisticsFromCrystal(list_elements)
    path_save_current_element_result_dir = os.path.join(path_save_result_dir, element)
    test.get_mol_files_by_indentifiers(locals()['list_identifier_{}'.format(element)], path_save_current_element_result_dir)
