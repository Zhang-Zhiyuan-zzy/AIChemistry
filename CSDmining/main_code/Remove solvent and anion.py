'''
-*- coding: utf-8 -*-
@Time    : 2021/4/13 13:21
@Author  : Zhang Zhiyuan
@File    : Remove solvent and anion.py
@Software: PyCharm
'''

import os
from public_module import statistics_from_crystal
import pandas as pd


#Setting parameters

path_root = os.path.join(os.path.dirname(__file__), '..')
path_anion = os.path.join(path_root, 'input_data', 'anion')
path_id_reader = os.path.join(path_root, 'input_data', 'files_include_specific_atoms.xlsx')
path_mol2_save = os.path.join(path_root, 'results', 'mol2', 'Remove_solvent_and_anions')
elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']

for ele in elements:
    path_mol2_file = os.path.join(path_mol2_save, ele)
    list_id = list(pd.read_excel(path_id_reader, sheet_name=ele, engine='openpyxl').iloc[:, 0])
    mining = statistics_from_crystal.StatisticsFromCrystal(elements)
    mining.get_entry(list_id)
    mining.delete_solvents()
    mining.delete_anion(path_anion)
    mining.csd_writer(path_mol2_file)

with pd.ExcelWriter(path_results_save) as writer:
    for ele in elements:
        for chelating_num in locals()['{}'.format(ele)]:
            df = pd.DataFrame(locals()['{}'.format(ele)][chelating_num], columns=sheet_head)
            df.to_excel(writer, sheet_name=ele, engine='openpyxl')
