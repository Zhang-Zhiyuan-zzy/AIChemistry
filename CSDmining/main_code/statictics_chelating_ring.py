'''
-*- coding: utf-8 -*-
@Time    : 2021/4/16 14:25
@Author  : Zhang Zhiyuan
@File    : statictics_chelating_ring.py
@Software: PyCharm
'''
import copy
import os
import pandas as pd
from public_module.statistics_from_crystal import StatisticsFromCrystal
from itertools import combinations

# Setting parameters
path_root = os.path.join(os.path.dirname(__file__), '..')
path_mol2_file = os.path.join(path_root, 'results', 'mol2', 'Remove_solvent_and_anions')
path_results_save = os.path.join(path_root, 'results', 'chelating_ring.xlsx')
path_anion = os.path.join(path_root, 'input_data', 'anion')
path_id_reader = os.path.join(path_root, 'input_data', 'files_include_specific_atoms.xlsx')
elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']
basic_head = ['chlt_num', 'chlt_atom', 'id']
bl_module = ['bond type', 'bl', 'ibl']
chlt_module = ['start_end_atom', 'shortest_path', 'chlt_ring_chain']


# Analyzing chelating ring
for ele in elements:
    list_id = list(pd.read_excel(path_id_reader, sheet_name=ele, engine='openpyxl').iloc[:, 0])
    mining = StatisticsFromCrystal(elements)
    mining.get_entry(list_id)
    mining.delete_solvents()
    mining.delete_anion(path_anion)
    locals()['{}'.format(ele)] = mining.get_chelating_ring_count(elements)

# Save result to Excel
with pd.ExcelWriter(path_results_save) as writer:
    for ele in elements:
        for chelating_num in locals()['{}'.format(ele)]:
            column = copy.deepcopy(basic_head)
            list_num = range(chelating_num)
            for i in list_num:
                column.extend([name + str(i) for name in bl_module])
            for i, j in combinations(list_num, 2):
                column.extend([name + str(i) + str(j) for name in chlt_module])

            df = pd.DataFrame(locals()['{}'.format(ele)][chelating_num], columns=column)
            df.to_excel(writer, sheet_name=ele + str(chelating_num), engine='openpyxl')