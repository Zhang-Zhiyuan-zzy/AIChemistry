'''
-*- coding: utf-8 -*-
@Time    : 2021/4/15 14:49
@Author  : Zhang Zhiyuan
@File    : crystal_system.py
@Software: PyCharm
'''

import clusters
import os
import pandas as pd
from public_module.statistics_from_crystal import StatisticsFromCrystal

# Setting parameters
path_root = os.path.join(os.path.dirname(__file__), '..')
path_mol2_file = os.path.join(path_root, 'results', 'mol2', 'Remove_solvent_and_anions')
path_results_save = os.path.join(path_root, 'results', 'crystal_system.xlsx')
path_anion = os.path.join(path_root, 'input_data', 'anion')
path_id_reader = os.path.join(path_root, 'input_data', 'files_include_specific_atoms.xlsx')
elements = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']
sheet_head = ['Id', 'ra', 'rb', 'rc', 'm_num_ml', 'rcm',
             'Oc', 'Ombl', 'Oibl',' Orbl',
             'Nc', 'Nmbl', 'Nibl', 'Nrbl',
             'Sc', 'Smbl', 'Sibl', 'Srbl',
             'Pc', 'Pmbl', 'Pibl', 'Prbl',
             'mfsmam', 'mfsmak', 'mfsmwc']


for ele in elements:
    list_id = list(pd.read_excel(path_id_reader, sheet_name=ele, engine='openpyxl').iloc[:, 0])
    mining = StatisticsFromCrystal(elements)
    mining.get_entry(list_id)
    mining.delete_solvents()
    mining.delete_anion(path_anion)
    locals()['{}'.format(ele)] = mining.get_system_coordinated_information(ele)

# Save results
with pd.ExcelWriter(path_results_save) as writer:
    for ele in elements:
        df = pd.DataFrame(locals()['{}'.format(ele)], columns=sheet_head)
        df.to_excel(writer, sheet_name=ele, engine='openpyxl')


