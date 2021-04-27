#!/usr/bin/env python
# _*_ coding:utf-8 _*_
"""
Time    :   2021/3/15 16:18
Author  :   Cheng Min
File    :   statistics_identifiers_include_specific_elements.py
Software:   PyCharm
"""

"""获取包含下列元素的晶体名称及其化学式
注：再获取包含特定元素的晶体时，未进行溶剂去除。未考虑去除溶剂不影响得到的结果；相反，考虑溶剂会增加计算过程。
"""

if __name__=='__main__':

    import os
    import pandas as pd

    path_root = os.path.join(os.path.dirname(__file__), '..')
    path_module = os.path.join(path_root, 'public_module')
    from public_module import statistics_from_crystal

    # 参数设置
    crystal_number = None
    list_elements  = ['Na', 'Si', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']  # 考虑的元素
    path_save_result = os.path.join(path_root, 'results/files_include_specific_atoms.xlsx')

    # 统计主程序
    test = statistics_from_crystal.StatisticsFromCrystal(list_elements=list_elements, crystal_number=crystal_number)
    list_result = test.get_all_files_include_specific_elements()

    # 保存结果
    if not os.path.exists(path_save_result):
        os.mkdir(os.path.dirname(path_save_result))

    with pd.ExcelWriter(path_save_result) as writer:
        for i in range(len(list_elements)):
            df = pd.DataFrame(list_result[i], columns=['element', 'formula', 'identifier']).to_excel(writer, sheet_name=list_elements[i], index=None, engine='openpyxl')


