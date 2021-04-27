'''
-*- coding: utf-8 -*-
@Time    : 2021/4/19 14:46
@Author  : Zhang Zhiyuan
@File    : statictics_bond_length.py
@Software: PyCharm
'''

import os
import pandas as pd
from public_module.statistics_from_crystal import StatisticsFromCrystal

# Setting paramaters
path_root = os.path.join(os.path.dirname(__file__), '..')
path_id = os.path.join(path_root, 'input_data', 'files_include_specific_atoms.xlsx')
path_sub = os.path.join(path_root, 'input_data', 'Function group')
path_anion = os.path.join(path_root, 'input_data', 'anion')
path_excel_save = os.path.join(path_root, 'results', 'bond_length.xlsx')
elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']

# Statistic
for ele in elements:
    list_id = list(pd.read_excel(path_id, sheet_name=ele, engine='openpyxl').iloc[:, 0])
    mining = StatisticsFromCrystal(elements)
    mining.get_entry(list_id)
    mining.delete_solvents()
    mining.delete_anion(path_anion)
    mining.define_substructures(con_path=path_sub)
    locals()['dict_bond_length_{}'.format(ele)] = mining.coordination_bond_length(ele)

# Save results
with pd.ExcelWriter(path_excel_save) as writer:
    for ele in elements:
        list_result = []
        for bond_type in locals()['dict_bond_length_{}'.format(ele)]:
            list_result.extend(locals()['dict_bond_length_{}'.format(ele)][bond_type])
        df1 = pd.DataFrame(list_result, columns=['bond type', 'atomic pair', 'iBL', 'BL', 'sub', 'main', 'coo_num',
                                                 'metal', 'c_atom', 'identifier'])
        df1.to_excel(writer, sheet_name=ele, engine='openpyxl')
