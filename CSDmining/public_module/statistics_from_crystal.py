#!/usr/bin/env python
# _*_ coding:utf-8 _*_
"""
Time    :   2021/3/15 11:06
Author  :   Cheng Min
File    :   statistics_from_crystal.py
Software:   PyCharm
"""

""""从CCDC数据库中统计相关信息
"""

import os
import glob
from tqdm import tqdm
from ccdc.search import QueryAtom
from ccdc import search
from ccdc import io
from itertools import combinations


class StatisticsFromCrystal(object):

    def __init__(self, list_elements, crystal_number=None):
        """初始化

        :param list_elements:考虑的配位原子构成的列表，type：list
        :param crystal_number:筛选的晶体的数量，type：int。若为None，则表示对整个CCDC数据库筛选
        """
        self.list_elements = list_elements
        self.crystal_number = crystal_number
        self.dict_substructure = dict()
        self.entry_reader = self.__get_crystal()
        self.dict_crystal = dict()


    def __get_crystal(self):
        count = 0
        list_entry = []
        entry_reader = io.EntryReader('CSD')

        if self.crystal_number is None:
            return entry_reader
        else:
            for i in entry_reader:
                list_entry.append(i)
                count = count + 1
                if count > self.crystal_number:
                    break

        return list_entry


    def get_entry(self, list_identifier):

        list_entry = list()
        if isinstance(list_identifier, list):
            p_bar = tqdm(list_identifier)
            for Id in p_bar:
                try:
                    list_entry.append(self.entry_reader.entry(Id))
                except:
                    BaseException()
                p_bar.set_description('正在加载entry：')
            print('%s entry load successful' % len(list_entry))
        else:
            raise TypeError('The input should be a list of str!')
        self.entry_reader = list_entry


    def is_solvent(self, c, solvent_smiles):
        """检查该组分是否是溶剂"""
        return c.smiles == 'O' or c.smiles in solvent_smiles


    def has_metal(self, c):
        """检查该组分是否含有金属"""
        return any(a.is_metal for a in c.atoms)


    def is_multidentate(self, c, mol):
        """Check for components bonded to metals more than once.

        If monodentate is not specified in the arguments, skip this test.
        """
        monodentate = False
        if not monodentate:
            return True
        got_one = False
        for a in c.atoms:
            orig_a = mol.atom(a.label)
            if any(x.is_metal for b in orig_a.bonds for x in b.atoms):
                if got_one:
                    return True
                got_one = True
        return False


    def csd_writer(self, path_save, *, file_type='.mol2'):
        '''

        :param path_save:
        :param file_type:
        :return:
        '''
        if not os.path.isdir(path_save):
            os.makedirs(path_save, exist_ok=True)
        os.chdir(path_save)
        p_bar = tqdm(self.entry_reader)
        for entry in p_bar:
            try:
                with io.CrystalWriter(entry.identifier + '.' + file_type) as writer:
                    writer.write(entry.crystal)
            except:
                BaseException()
            p_bar.set_description('Crystal' + '.' + file_type + 'writing:')


    def csd_read(self, path_read, *, file_type='mol2'):
        '''
        从指定文件夹读取其中的cif文件，文件夹中只能有cif，然后以dict{cif文件名：cif中的分子}保存到：
        self.__Dict_mol_from_cif中
        :param path_read: <class:str>用于专门储存Cifs的文件夹路径
        :return: None
        '''
        # 确定每个已经去除了溶剂的cif文件的名称和绝对路径
        list_cif_names = os.listdir(path_read)
        list_path_cifs = glob.glob(os.path.join(path_read, '*.' + file_type))
        count = -1
        nub = 0
        # 读取cif文件中的官能团
        dict_crys_temp = dict()
        p_bar = tqdm(list_path_cifs)
        for path_cif_temp in p_bar:
            count += 1
            try:
                mol_temp = io.CrystalReader(path_cif_temp)[0]
                dict_crys_temp[list_cif_names[count]] = mol_temp
            except:
                nub += 1
                print(list_cif_names[count])
                continue
            p_bar.set_description('Crystal.' + file_type + 'reading:')
        self.dict_crystal = dict_crys_temp
        print(nub)


    def coordination_bond_length(self, element=None):
        '''

        :param dict_coo_sub:
        :param element:
        :return:
        '''
        p_bar = tqdm(self.entry_reader)
        # {bond type: [bond type, ideal bond length, bond length, substructure,
        #              The main of substructure connection, Connection with ligands
        #              metal, coordinated atom, identifier]
        dict_bond = dict()
        for entry in p_bar:
            # Getting molecule
            mol = entry.molecule
            # Remove all of hydrogen
            mol.remove_hydrogens()
            # Ensure labels are unique
            mol.normalise_labels()
            # remove single_atoms
            mol.remove_atoms([single_atom for single_atom in mol.atoms if not single_atom.bonds])
            # dict to save coordinated bonds information
            dict_mol_bonds = dict()
            for comp in mol.components:
                dict_coo_atoms = self.__coo_atom(comp, element)
                # Getting the bond length
                dict_bond_length = self.__measurement_bond_length(comp, element)
                # Get coordinate substructures which the coordinated atoms belong to
                for metal in dict_coo_atoms:
                    # Set of atoms in the matched substructure in the molecule
                    set_sub_atoms = set()

                    set_coo_atoms_label = dict_coo_atoms[metal]
                    # Set of coordinated atoms in defined substructures
                    set_coo_atoms_label_in_sub = set()
                    '''
                    Get information of coordinated atoms about which substructure belongs to
                    '''
                    list_sub_name = sorted(self.dict_substructure, reverse=True)
                    for sub_name in list_sub_name:
                        # Defining method of substructure searching
                        substructure_search = search.SubstructureSearch()
                        substructure_search.add_substructure(self.dict_substructure[sub_name])
                        # Searching
                        hits = substructure_search.search(comp)
                        if not not hits:
                            for hit in hits:
                                hit_atoms_label = set(atom.label for atom in hit.match_atoms())
                                if not any(atom_label in set_sub_atoms for atom_label in hit_atoms_label):
                                    set_sub_atoms.update([atom.label for atom in hit.match_atoms()
                                                          if atom.atomic_symbol != 'C'])
                                    coo_atoms_in_sub_label = set_coo_atoms_label & hit_atoms_label
                                    if not not coo_atoms_in_sub_label:
                                        set_coo_atoms_label_in_sub.update(coo_atoms_in_sub_label)
                                        # Add the substructures type to bond information
                                        # Add the type of the main that the substructure connects with
                                        for coo_atom_in_sub_label in coo_atoms_in_sub_label:
                                            if len(dict_bond_length[metal + '-' + coo_atom_in_sub_label]) < 4:
                                                dict_bond_length[metal + '-' + coo_atom_in_sub_label]. \
                                                    extend([sub_name, self.__main_type_sub_connect(hit)])

                    # Set of coordinated atoms out defined substructures
                    set_coo_atoms_label_out_sub = set_coo_atoms_label - set_coo_atoms_label_in_sub
                    # Filling coordination atoms with Nan that do not belong to any defined substructure
                    for coo_atom_out_sub_label in set_coo_atoms_label_out_sub:
                        dict_bond_length[metal + '-' + coo_atom_out_sub_label]. \
                            extend(['NaN', 'NaN'])

                    '''
                    Get information of Connection of metal with ligands
                    '''
                    c_comp = comp.copy()
                    c_comp.remove_atoms(atom for atom in c_comp.atoms if atom.atomic_symbol == element)
                    set_coo_atoms_label_in_comp = set()
                    for sub_c_comp in c_comp.components:
                        set_sub_c_comp_atoms = set(atom.label for atom in sub_c_comp.atoms)
                        common_coo_atoms = set_sub_c_comp_atoms & set_coo_atoms_label
                        set_coo_atoms_label_in_comp.update(common_coo_atoms)
                        len_common_coo_atoms = len(common_coo_atoms)
                        if len_common_coo_atoms != 0:
                            for common_coo_atom in common_coo_atoms:
                                dict_bond_length[metal + '-' + common_coo_atom].append(len_common_coo_atoms)
                    set_coo_atoms_label_out_comp = set_coo_atoms_label - set_coo_atoms_label_in_comp
                    for common_coo_atom in set_coo_atoms_label_out_comp:
                        dict_bond_length[metal + '-' + common_coo_atom].append(0)
                # Setting basic information for bonds

                for metal in dict_coo_atoms:
                    for coo_atom in dict_coo_atoms[metal]:
                        dict_bond_length[metal + '-' + coo_atom].extend([metal, coo_atom, entry])
                        dict_bond_length[metal + '-' + coo_atom].insert(0, metal + '-' + coo_atom)
                        dict_bond_length[metal + '-' + coo_atom].insert(0, element + '-' + mol.atom(
                            coo_atom).atomic_symbol)
                # Save information of bond length
                dict_mol_bonds.update(dict_bond_length)

            # Save results
            for bond in dict_mol_bonds:
                if dict_mol_bonds[bond][0] in dict_bond:
                    dict_bond[dict_mol_bonds[bond][0]].append(dict_mol_bonds[bond])
                else:
                    dict_bond[dict_mol_bonds[bond][0]] = []
                    dict_bond[dict_mol_bonds[bond][0]].append(dict_mol_bonds[bond])

        return dict_bond


    def delete_solvents(self, list_solvent_names=None):
        """删除晶体中的溶剂，若没有指定溶剂列表，则默认为CCDC数据库自带的溶剂列表

        :param list_solvent_names: 溶剂名称构成的列表，type：list or tuple
        :return: None
        """

        # CSD数据库的溶剂所在的路径
        solvent_file = os.path.join(
            os.path.dirname(io.csd_directory()),
            'Mercury',
            'molecular_libraries',
            'ccdc_solvents'
        )

        # 若没指定需要去除的溶剂列表，则会将CSD数据库中指定的74个溶剂都考虑进去。以下代码得到溶剂的smiles字符串
        if not list_solvent_names:
            if os.path.isdir(solvent_file):
                solvent_smiles = [
                    io.MoleculeReader(f)[0].smiles
                    for f in glob.glob(os.path.join(solvent_file, '*.mol2'))
                ]
            else:
                raise FileExistsError('路径不存在！')
        else:
            if os.path.isdir(solvent_file):
                solvent_smiles = [
                    io.MoleculeReader(os.path.join(solvent_file, solvent + '.mol2')[0].smiles
                    for solvent in list_solvent_names)
                ]
            else:
                raise FileExistsError('路径不存在！')

        # 去除溶剂
        list_crystals_remove_solvents = []
        p_bar = tqdm(self.entry_reader)
        for entry in p_bar:
            try:
                if entry.has_3d_structure:
                    # Ensure labels are unique
                    mol = entry.molecule
                    mol.normalise_labels()
                    # Use a copy
                    clone = mol.copy()
                    # Remove all bonds containing a metal atom
                    clone.remove_bonds(b for b in clone.bonds if any(a.is_metal for a in b.atoms))
                    # Work out which components to remove
                    to_remove = [
                        c
                        for c in clone.components
                        if not self.has_metal(c) and (not self.is_multidentate(c, mol) or self.is_solvent(c, solvent_smiles))
                    ]
                    # Remove the atoms of selected components
                    mol.remove_atoms(
                        mol.atom(a.label) for c in to_remove for a in c.atoms
                    )
                    # Write the CIF
                    entry.crystal.molecule = mol
                    list_crystals_remove_solvents.append(entry)
                else:
                    list_crystals_remove_solvents.append(entry)
            except BaseException:
                list_crystals_remove_solvents.append(entry)
            p_bar.set_description('正在去除溶剂：')
        self.entry_reader = list_crystals_remove_solvents
        return None


    def delete_anion(self, path_anion):
        '''
        removing anions which are defined by mol2 file in a entry
        :param path_anion: the defined anions files
        :return: None
        '''

        if os.path.isdir(path_anion):
            anion_list = [
                    search.MoleculeSubstructure(io.MoleculeReader(f)[0].components[0])
                    for f in glob.glob(os.path.join(path_anion, '*.mol2'))
                     ]
        else:
            raise FileExistsError('do not find the path!')

        list_crystals_remove_anion = []
        p_bar = tqdm(self.entry_reader)

        for entry in p_bar:
            if entry.has_3d_structure:
                # Ensure labels are unique
                mol = entry.molecule
                mol.normalise_labels()
                # Use a copy
                clone = mol.copy()
                # Remove all metal atoms
                clone.remove_atoms(a for a in clone.atoms
                                   if a.is_metal or not a.bonds)
                for c in clone.components:
                    for anion in anion_list:
                        ani_search = search.SubstructureSearch()
                        ani_search.add_substructure(anion)
                        hits = ani_search.search(c)
                        for hit in hits:
                            hit_atoms = hit.match_atoms()
                            if len(hit_atoms) == len(c.atoms):
                                mol.remove_atoms(mol.atom(a.label) for a in hit_atoms)
                entry.crystal.molecule = self.__delete_isolated_atoms(mol)
                list_crystals_remove_anion.append(entry)
            p_bar.set_description('Anions removing...')

        self.entry_reader = list_crystals_remove_anion


    def __delete_isolated_atoms(self, mol):
        '''
        remove all of isolated atoms in the molecule
        :param mol: the molecule
        :return: the molecule after remove isolated atoms
        '''
        mol.remove_atoms(a for a in mol.atoms if not a.bonds)
        return mol


    def get_all_files_include_specific_elements(self):
        """获取包含指定元素的晶体的相关信息

        :return: 包含晶体信息的列表
        """

        # 定义需要统计的晶体
        p_bar = tqdm(self.entry_reader)

        # 定义多个列表，每个元素单独一个列表用于存放数据
        list_result = []
        for _ in self.list_elements:
            list_result.append([])

        # 开始统计
        for crystal in p_bar:
            crystal = crystal.molecule
            for i in range(len(self.list_elements)):
                settings = search.Search.Settings()
                settings.must_have_elements = [self.list_elements[i]]
                judge = settings.test(crystal)
                if judge:
                    list_result[i].append([self.list_elements[i], crystal.formula, crystal.identifier])
            p_bar.set_description('查找进程：')

        return list_result


    def get_mol_files_by_indentifiers(self, list_identifiers, path_save_mol_dir, list_solvent_names=None):
        """保存*.mol2到指定文件夹中

        :param list_identifiers: 晶体的id,type：list or tuple or array
        :param path_save_mol_dir: 保存晶体mol2文件的文件夹绝对路径，type：string
        :param list_solvent_names: 溶剂的名称，type：list or None
        :return: None
        """
        # 重新定义self.entry_reader
        self.entry_reader = [self.entry_reader.entry(id_name) for id_name in list_identifiers]
        # 删除溶剂
        self.delete_solvents(list_solvent_names=list_solvent_names)  # 去除溶剂，并且重新生成去除溶剂过后的self.entry_reader
        # 判断储存mol文件的路径是否存在，若不存在，则创建
        if not os.path.exists(path_save_mol_dir):
            os.makedirs(path_save_mol_dir)
        # 保存mol文件到指定文件夹
        p_bar = tqdm(self.entry_reader)
        for entry in p_bar:
            mol_file = entry.molecule
            os.chdir(path_save_mol_dir)
            with io.MoleculeWriter('%s.mol2' % (entry.identifier)) as writer:
                writer.write(mol_file)
            p_bar.set_description('正在保存文件：')

        return None


    def get_neighbor_atoms(self, query_atom, path_mol2_file_dir):
        """获取配位原子

        :param query_atom:金属元素，Sr,K,Na等，type:string
        :param path_mol2_file_dir: *.mol2文件所在的绝对路径，type：string
        :return: 配位原子名称及数量，type：dict
        """

        # 重新定义self.entry_reader,从mol2文件中读取文件，保存为entry
        list_entry = []
        list_mol2_file_path = glob.glob(os.path.join(path_mol2_file_dir, '*.mol2'))
        for path_temp in list_mol2_file_path:
            entry_temp = io.EntryReader(path_temp)[0]
            list_entry.append(entry_temp)
        self.entry_reader = list_entry.copy()

        # 创建搜索的原子，此处为Na/Mg/Cr/.../Sr
        s = search.QuerySubstructure()
        q = QueryAtom(query_atom)
        s.add_atom(q)

        # 搜索neighbors
        entry_reader = self.entry_reader
        dict_ligating_atom_statistics = dict()  # 用于储存金属原子的配位原子类型及数量
        pbar = tqdm(range(len(entry_reader)))
        for count in pbar:
            set_ligating_atom = set()  # 空集合，用于存放该分子中所涉及的配位原子类型
            mol = entry_reader[count].crystal.molecule
            # 找出该晶体文件中金属原子的配位原子
            for atom in mol.atoms:
                bool_judge = s.match_atom(atom)  # 判断该原子是否是金属原子
                # 若和金属原子匹配
                if bool_judge:
                    neighbors = atom.neighbours  # 寻找其配位原子
                    # 若配位原子不存在，则跳过
                    if len(neighbors) == 0:
                        pass
                    # 配位原子不为0，则将配位原子增加到对应的集合当中
                    else:
                        for neighbour in neighbors:
                            set_ligating_atom.add(neighbour.atomic_symbol)

            # 对金属原子的配位原子数统计
            for element in set_ligating_atom:
                # 字典中若存在该element，则计数增加1
                try:
                    dict_ligating_atom_statistics[element] += 1
                # 字典中若不存在该element，则新建该element，并计数1
                except BaseException:
                    dict_ligating_atom_statistics[element] = 1
            pbar.set_description('正在统计配位原子：')

        return dict_ligating_atom_statistics


    @staticmethod
    def get_all_function_groups(path_mols, path_con):
        """从*.mol2文件中找到指定基团的类型及数量

        :param path_cifs:
        :param path_con:
        :return:
        """

        # 确定每个已经去除了溶剂的*.mol2文件的名称和绝对路径
        list_mol_names = os.listdir(path_mols)
        list_path_mols = glob.glob(os.path.join(path_mols, '*.mol2'))

        # 通过con定义功能基团
        list_con_names = os.listdir(path_con)
        path_conner_list = glob.glob(os.path.join(path_con, '*.con'))
        list_connser_substructure = []
        for path in path_conner_list:
            connser_substructure = search.ConnserSubstructure(path)
            list_connser_substructure.append(connser_substructure)

        # 读取mol2文件中
        dict_result = dict()
        count = 0
        pbar = tqdm(list_path_mols)
        for path_cif_temp in pbar:
            list_temp = []  # 维度为len(list_connser_substructure)，即维度为定义的官能团个数；该列表用于储存当前cif文件中包含基团的数目
            mol_temp = io.MoleculeReader(path_cif_temp)[0]  # 读取cif文件
            for func_group in list_connser_substructure:
                substructure_search = search.SubstructureSearch()
                _ = substructure_search.add_substructure(func_group)
                hits = substructure_search.search(mol_temp)
                list_temp.append(len(hits))
            dict_result[list_mol_names[count]] = list_temp
            count += 1
        pbar.set_description('正在统计所有的指定基团：')

        return dict_result, list_con_names


    @staticmethod
    def get_neighbor_function_groups(path_mols, path_con, query_atom):

        # 确定每个已经去除了溶剂的*.mol2文件的名称和绝对路径
        list_mol_names = os.listdir(path_mols)
        list_path_mols = glob.glob(os.path.join(path_mols, '*.mol2'))

        # 通过con定义功能基团
        list_con_names = os.listdir(path_con)
        path_conner_list = glob.glob(os.path.join(path_con, '*.con'))
        list_connser_substructure = []
        for path in path_conner_list:
            connser_substructure = search.ConnserSubstructure(path)
            list_connser_substructure.append(connser_substructure)

        # 统计配位基团的类型及数量
        dict_result = dict()
        pbar = tqdm(range(len(list_path_mols)))

        for i in pbar:
            # 读取分子，并且读取出其中的components
            path_mol = list_path_mols[i]
            mol = io.MoleculeReader(path_mol)[0]
            list_components = mol.components
            mol.normalise_labels()
            # 统计每个基团在分子中出现的次数
            list_temp = []  # 储存每个mol2文件中匹配到的配位基团的数量
            for con in list_connser_substructure:
                count_temp = 0  # 基团出现数量
                for component in list_components:
                    set_temp = set()  # 用于存放出现的基团的字符串
                    # 查询金属原子
                    m = QueryAtom(query_atom)
                    s = search.QuerySubstructure()
                    s.add_atom(m)
                    sub_search = search.SubstructureSearch()
                    sub_search.add_substructure(s)
                    mol_metals = sub_search.search(component)

                    if len(mol_metals) > 0:
                        substructure_search = search.SubstructureSearch()
                        substructure_search.add_substructure(con)
                        hits = substructure_search.search(component)

                        if len(hits) > 0:
                            for hit in hits:
                                temp_hit_atoms = hit.match_atoms() # 匹配到的基团的原子
                                for temp_metal in mol_metals:
                                    temp_metal = temp_metal.match_atoms()[0]
                                    common_elements = set(temp_metal.neighbours) & set(temp_hit_atoms)
                                    if len(common_elements) > 0:
                                        set_temp.add(str(temp_hit_atoms))
                                # for num in range(len(mol_metals)):
                                #     metal_label = query_atom + str(num + 1)
                                #     temp_metal = component.atom(metal_label)
                                #     common_elements = set(temp_metal.neighbours) & set(temp_hit_atoms)
                                #     if len(common_elements) > 0:
                                #         set_temp.add(str(temp_hit_atoms))
                    count_temp += len(set_temp)
                list_temp.append(count_temp)
            dict_result[list_mol_names[i]] = list_temp

        return dict_result, list_con_names


    def get_system_coordinated_information(self, m_atom):
        '''
        Determining coordinated information for a system include:
            1) Whether the system is a extensible frameworks and
            How much role do metal atoms play in the ductility of the framework?
            ra, rb, rc
            2) Binding strength between specific metal atoms and ligands and correlative coordinated atoms
            m_num_ml, rcm,
            Oc, Ombl, Oibl, Orbl,
            Nc, Nmbl, Nibl, Nrbl,
            Sc, Smbl, Sibl, Srbl,
            Pc, Pmbl, Pibl, Prbl
            3) Atomic selectivity and atomic economy
            mfsmam, mfsmak, mfsmwc
        :param m_atom: <str> specific metal atomic_symbol
        :return: <list> list_result:
        [
            mol.identifier,
            ra: ratio of specific metals which aside polymeric framework,
            rb: ratio of specific metals which build polymeric framework,
            rc: ratio of specific metals connect to polymeric framework,
            m_num_ml:
            the mean number of bonds between specific element and the extended ligands in extensible framework
            or the ligands whose molecular mass more than 50 or the maximal ligands in inextensible framework,
            rcm: the ratio of chelating metals in all of bonds with special element
            Oc: O_count, Ombl: O_mean_bond_length, Oibl: O_ideal_bond_length, Orbl: O_ralative_bond_length,
            Nc: N_count, Nmbl: N_mean_bond_length, Nibl: N_ideal_bond_length, Nrbl: N_ralative_bond_length,
            Sc: S_count, Smbl: S_mean_bond_length, Sibl: S_ideal_bond_length, Srbl: S_ralative_bond_length,
            Pc: P_count, Pmbl: P_mean_bond_length, Pibl: P_ideal_bond_length, Prbl: P_ralative_bond_length,
            mfsmam: molar fraction of specific metals in all matal,
            mfsmak: molar fraction of specific metals in Alkali and alkaline earth metals,
            mfsmwc: mass fraction of the specified element in the whole crystal
        ]
        '''
        list_coo_atoms = ['O', 'N', 'S', 'P']  # list of analyzed coordinated atoms
        # list of considered Alkali and alkaline earth metals
        list_ak_matals = ['Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']
        list_results = list()  # Return
        p_bar = tqdm(self.entry_reader)  # Crystal information loaded by func: csd_read()
        error_num = 0
        all_num = 0
        for entry in p_bar:

            # getting the molecule in the crystal
            mol = entry.molecule
            # Remove all of hydrogen
            mol.remove_hydrogens()
            # Ensure labels are unique
            mol.normalise_labels()
            # List to stored results for a crystal
            list_crys_results = list()
            # sum number of bonds connect metals to the extensible framework or main molecule in inextensible framework
            sum_b_cf = 0
            sum_c_a = 0        # sum number of chelating atoms in a crystal
            num_m_all = 0      # sun number of special atoms in a crystal
            num_all_metal = 0  # count for all metals
            num_alk_metal = 0  # count for alkali and alkaline earth metals
            num_spc_metal = 0  # count for special metals
            '''
            determine ra, rb, rc, m_num_ml, rcm
            '''
            # for extensible frameworks
            if mol.is_polymeric:
                num_m_build_frame = 0
                num_m_aside_frame = 0
                for comp in mol.components:
                    all_num += 1
                    try:
                        dict_coo_atoms = self.__coo_atom(comp, m_atom)
                    except:
                        error_num += 1
                        break
                    num_m_all += len(dict_coo_atoms)
                    # skipping inextensible components in a framework when analyze ra, rb, rc
                    if comp.is_polymeric:
                        # Analyze each metal one by one for ra, rb, rc, m_mum_ml
                        for m_label in dict_coo_atoms:
                            bonds = comp.atom(m_label).bonds
                            # if the metal have polymeric bonds
                            if any(b._bond.polymeric() for b in bonds):
                                num_m_build_frame += 1
                                sum_b_cf += len(bonds)
                            else:
                                clone = comp.copy()
                                clone.remove_atom(clone.atom(m_label))
                                sub_clone_comps = clone.components
                                if len(sub_clone_comps) == 1:
                                    num_m_aside_frame += 1
                                    sum_b_cf += len(bonds)
                                else:
                                    num_poly = 0
                                    # Set of metal label connect to framework
                                    set_mc = set()
                                    for sub_clone_comp in sub_clone_comps:
                                        if sub_clone_comp.is_polymeric:
                                            num_poly += 1
                                            set_scc_as = set(a.label for a in sub_clone_comp.atoms)
                                            set_mc.update(set_scc_as & dict_coo_atoms[m_label])
                                    if num_poly >= 2:
                                        num_m_build_frame += 1
                                    else:
                                        num_m_aside_frame += 1
                                    # Count the number of coordinated atoms in the framework
                                    sum_b_cf += len(set_mc)
                        # Determining rcm
                        clone = comp.copy()
                        clone.remove_atoms(a for a in clone.atoms if a.is_metal)
                        # List of set of atoms in each main(M.W >= 50) sub_clone_comp(ligand)
                        list_a_lab_ml = list()
                        for sub_clone_comp in clone.components:
                            set_a_lab_pml = set(a.label for a in sub_clone_comp.atoms)
                            list_a_lab_ml.append(set_a_lab_pml)
                        for m_label in dict_coo_atoms:

                            for set_a_lab_pml in list_a_lab_ml:
                                # Coordinated atoms between metal and a ligand
                                set_ca_ml = dict_coo_atoms[m_label] & set_a_lab_pml
                                if len(set_ca_ml) >= 2:
                                    sum_c_a += 1
                                    break
                # results
                if num_m_all == 0:
                    ra = 0
                    rb = 0
                    rc = 0
                    m_num_ml = 0
                    rcm = 0
                else:
                    ra = round(num_m_aside_frame / num_m_all, 5)
                    rb = round(num_m_build_frame / num_m_all, 5)
                    rc = ra + rb
                    m_num_ml = round(sum_b_cf / num_m_all, 5)
                    rcm = round(sum_c_a / num_m_all, 5)

            # for inextensible frameworks
            else:
                for comp in mol.components:
                    all_num += 1
                    try:
                        dict_coo_atoms = self.__coo_atom(comp, m_atom)
                    except:
                        error_num += 1
                        break
                    # Total of specific
                    num_m_all += len(dict_coo_atoms)
                    # Set to store atoms in the ligand whose molecular weight more than 50
                    set_a_lab_mL = set()
                    clone = comp.copy()
                    clone.remove_atoms(a for a in clone.atoms if a.is_metal)
                    # List of set of atoms in each main(M.W >= 50) sub_clone_comp(ligand)
                    list_a_lab_ml = list()
                    for sub_clone_comp in clone.components:
                        if sub_clone_comp.molecular_weight >= 50:
                            # set of atoms per main(M.W >= 50) sub_clone_comp(ligand)
                            set_a_lab_pml = set(a.label for a in sub_clone_comp.atoms)
                            set_a_lab_mL.update(set_a_lab_pml)
                            list_a_lab_ml.append(set_a_lab_pml)
                    for m_label in dict_coo_atoms:
                        # Set of metal label connect to framework
                        set_mc = dict_coo_atoms[m_label] & set_a_lab_mL
                        sum_b_cf += len(set_mc)
                        for set_a_lab_pml in list_a_lab_ml:
                            # Coordinated atoms between metal and a ligand
                            set_ca_ml = dict_coo_atoms[m_label] & set_a_lab_pml
                            if len(set_ca_ml) >=2:
                                sum_c_a += 1
                                break
                # Results
                ra = 0
                rb = 0
                rc = 0
                if num_m_all == 0:
                    m_num_ml = 0
                    rcm = 0
                else:
                    m_num_ml = round(sum_b_cf / num_m_all, 5)
                    rcm = round(sum_c_a / num_m_all)


            '''
            Determining Oc, Nc, Sc, Pc, Ombl, Nmbl, Smbl, Psbl, Oibl, Nibl, Sibl, Pibl, Orbl, Nrbl, Srbl, Prbl
                         mfsmam, mfsmak, mfsmwc
            '''
            # initializing bond length information
            dict_bl = {
                'Oc': 0, 'Nc': 0, 'Sc': 0, 'Pc': 0,
                'Osbl': 0, 'Nsbl': 0, 'Ssbl': 0, 'Psbl': 0,
                'Oibl': 0, 'Nibl': 0, 'Sibl': 0, 'Pibl': 0
            }
            # mining
            for comp in mol.components:
                all_num += 1
                try:
                    dict_coo_atoms = self.__coo_atom(comp, m_atom)
                    num_spc_metal += len(dict_coo_atoms)
                except:
                    error_num += 1
                    break
                # mfsmam, mfsmak, mfsmwc
                for atom in comp.atoms:
                    if atom.atomic_symbol in list_ak_matals:
                        num_alk_metal += 1
                        num_all_metal += 1
                    elif atom.is_metal:
                        num_all_metal += 1
                for atom in comp.atoms:
                    if atom.atomic_symbol == m_atom:
                        aaw = atom.atomic_weight
                        break
                # Oc, Nc, Sc, Pc, Ombl, Nmbl, Smbl, Psbl, Oibl, Nibl, Sibl, Pibl, Orbl, Nrbl, Srbl, Prbl
                for m_label in dict_coo_atoms:
                    for a_label in dict_coo_atoms[m_label]:
                        a_symbol = comp.atom(a_label).atomic_symbol
                        bond = comp.bond(m_label, a_label)
                        if a_symbol in list_coo_atoms:
                            dict_bl[a_symbol + 'c'] += 1
                            dict_bl[a_symbol + 'sbl'] += round(bond.length, 5)
                            dict_bl[a_symbol + 'ibl'] += round(bond.ideal_bond_length, 5)
            for ele in list_coo_atoms:
                if dict_bl[ele + 'c'] == 0:
                    locals()['{}c'.format(ele)] = 0
                    locals()['{}mbl'.format(ele)] = 0
                    locals()['{}ibl'.format(ele)] = 0
                    locals()['{}rbl'.format(ele)] = 0
                else:
                    locals()['{}c'.format(ele)] = dict_bl[ele + 'c']
                    locals()['{}mbl'.format(ele)] = dict_bl[ele + 'sbl'] / dict_bl[ele + 'c']
                    locals()['{}ibl'.format(ele)] = dict_bl[ele + 'ibl'] / dict_bl[ele + 'c']
                    locals()['{}rbl'.format(ele)] = dict_bl[ele + 'sbl'] / (dict_bl[ele + 'c'] * dict_bl[ele + 'ibl'])
            try:
                mfsmam = num_spc_metal / num_all_metal
                mfsmak = num_spc_metal / num_alk_metal
                mfsmwc = (num_spc_metal * aaw) / mol.molarcular_weight
            except:
                mfsmam = 0
                mfsmak = 0
                mfsmwc = 0
            '''
            Save results
            '''
            list_crys_results.extend([mol.identifier, ra, rb, rc, m_num_ml, rcm])
            for ele in list_coo_atoms:
                list_crys_results.extend(
                    [locals()['{}c'.format(ele)], locals()['{}mbl'.format(ele)],
                     locals()['{}ibl'.format(ele)], locals()['{}rbl'.format(ele)]]
                )
            list_crys_results.extend([mfsmam, mfsmak, mfsmwc])
            list_results.append(list_crys_results)
        print(error_num, all_num, error_num / all_num)

        return list_results


    def get_chelating_ring_count(self, m_ele):

        p_bar = tqdm(self.entry_reader)
        # Return5
        dict_results = dict()

        for entry in p_bar:
            # Get molecule
            mol = entry.molecule
            # Remove all of hydrogen
            mol.remove_hydrogens()
            # Ensure labels are unique
            mol.normalise_labels()
            for comp in mol.components:

                dict_coo_atoms = self.__coo_atom(comp, m_ele)
                clone = comp.copy()
                clone.remove_atoms(a for a in clone.atoms if a.is_metal)
                # Get set of atoms' labels in each ligand in the molecular component
                for sub_clone_lg in clone.components:
                    set_a_lab_pml = set(a.label for a in sub_clone_lg.atoms)
                    for m_lbl in dict_coo_atoms:
                        set_cdn_a_lbl = dict_coo_atoms[m_lbl]
                        set_cnxt_a_lbl = set_cdn_a_lbl & set_a_lab_pml
                        len_cnxt_a_lbl = len(set_cnxt_a_lbl)
                        if len_cnxt_a_lbl >= 2:
                            list_chlt_inf = []
                            chlt_a_symbol = ','.join([sub_clone_lg.atom(chlt_a_lbl).atomic_symbol
                                                      for chlt_a_lbl in set_cnxt_a_lbl])
                            chlt_num = len_cnxt_a_lbl
                            list_chlt_inf.extend([chlt_num, chlt_a_symbol])
                            list_chlt_inf.append(mol.identifier)
                            # Bond length of chelating bonds
                            for a_lbl in set_cnxt_a_lbl:
                                b = comp.bond(m_lbl, a_lbl)
                                list_chlt_inf.extend([m_lbl[0:1] + '-' + comp.atom(a_lbl).atomic_symbol,
                                                      b.length, b.ideal_bond_length])
                            # chelating rings
                            for two_atom_lbl in combinations(set_cnxt_a_lbl, 2):
                                start_end_piont = 'to'.join(two_atom_lbl)
                                shortest_path_atoms = sub_clone_lg.shortest_path_atoms(
                                    sub_clone_lg.atom(two_atom_lbl[0]), sub_clone_lg.atom(two_atom_lbl[1]))
                                shortest_path = len(shortest_path_atoms)
                                shortest_path_atoms_symbol = '>>'.join([a.atomic_symbol for a in shortest_path_atoms])
                                list_chlt_inf.extend([start_end_piont, shortest_path, shortest_path_atoms_symbol])


                            # Save result
                            if list_chlt_inf[0] in dict_results.keys():
                                dict_results[list_chlt_inf[0]].append(list_chlt_inf)
                            else:
                                dict_results[list_chlt_inf[0]] = []
                                dict_results[list_chlt_inf[0]].append(list_chlt_inf)
            p_bar.set_description('Chelating ring analyzing:')

        return dict_results


    def __coo_atom(self, comp, m_ele):
        '''
        Given a molecule<ccdc.molecule.Molecule> and a metal symbol<str>
        Return a dict about the metal with its coordinated atoms
        :param comp: <ccdc.molecule.Molecule>molecule
        :param m_ele: <str>metal elements symbol
        :param label: <bool>return atom label or atom entity
        :return:<dict>{<str>metal label: <set>{<str> coordinated atoms label}
        '''
        dict_coo_atoms = dict()  # dict for storing the coordination atom around the target metal ie.atoms one by one
        qm = QueryAtom(m_ele)
        qs = search.QuerySubstructure()
        qs.add_atom(qm)
        sub_search = search.SubstructureSearch()
        sub_search.add_substructure(qs)
        mol_metals = sub_search.search(comp)

        for mol_metal in mol_metals:
            set_atom = set()
            metal = mol_metal.match_atoms()
            set_atom.update([N_atom.label for N_atom in metal[0].neighbours])
            metal_label = metal[0].label
            dict_coo_atoms[metal_label] = set_atom

        return dict_coo_atoms


    def __measurement_bond_length(self, mol, m_atom):
        '''
        Given a molecule<ccdc.molecule.Molecule> and a metal symbol<str>
        to calculate the ideal bond length and bond length
        that the bond is between the Given metals and its coordinated atoms
        Return a dict contained the information about the bond ideal bond length and length

        Notice: the mol should be a component of <ccdc.molecule.Molecule> which atom labels are normalized
                or the process might suffer a expected break
        :param mol: <ccdc.molecule.Molecule>molecule
        :param m_atom: <str>metal symbol
        :return: <dict> {the bond<str>: [the ideal bond length<float>, the bond length<float>]} XXX
        '''

        # to stored {'metal_label - coordination atom label': [the ideal bond length, the bond length]}
        dict_bond_length = dict()
        dict_coo_atoms = self.__coo_atom(mol, m_atom)
        for m_label in dict_coo_atoms:
            for coo_a_label in dict_coo_atoms[m_label]:
                bond = mol.bond(m_label, coo_a_label)  # get bond from atoms' label between the bond
                dict_bond_length[m_label + '-' + coo_a_label] =\
                    [round(bond.ideal_bond_length, 4), round(bond.length, 4)]

        return dict_bond_length


    def __main_type_sub_connect(self, hit):
        '''
        To determine which type of the molecule main that the hit substructure connect with
        :param hit:
        :return:
        '''
        mol = hit.molecule
        set_hit_atoms = set(atom.label for atom in hit.match_atoms())
        set_sub_neigh = set()
        # get all atom in the hit substructure and around it
        for atom in hit.match_atoms():
            set_sub_neigh.update([na.label for na in atom.neighbours])
        # judge if aromatic
        if any(r.is_aromatic for na_label in set_sub_neigh for r in mol.atom(na_label).rings):
            return 'aromatic'
        # judge if conjugate
        elif any(r.is_fully_conjugated for na_label in set_sub_neigh for r in mol.atom(na_label).rings):
            return 'conjugate'
        else:
            set_hit_edge_atoms = set_sub_neigh - set_hit_atoms
            if any(eb.is_conjugated for ea_label in set_hit_edge_atoms for eb in mol.atom(ea_label).bonds):
                return 'conjugate'
            else:
                for ea_label in set_hit_edge_atoms:
                    for eb in mol.atom(ea_label).bonds:
                        if eb.bond_type == 'Double':
                            for a_eb in eb.atoms:
                                for na_a_eb in a_eb.neighbours:
                                    if (mol.bond(a_eb.label, na_a_eb.label).bond_type == 'Single'
                                            and any(b_na_a_eb.bond_type == 'Double' for b_na_a_eb in na_a_eb.bonds)):
                                        return 'conjugate'
            # judge if unsaturated
            for ea_label in set_hit_edge_atoms:
                for na_ea in mol.atom(ea_label).neighbours:
                    if (not (na_ea in set_hit_atoms)) and (na_ea.atomic_symbol == 'C'):
                        if mol.bond(ea_label, na_ea.label).bond_type in ['Double', 'Triple']:
                            return 'unsaturated'
            # judge if saturated
            for ea_label in set_hit_edge_atoms:
                for na_ea in mol.atom(ea_label).neighbours:
                    if (not (na_ea in set_hit_atoms)) and (na_ea.atomic_symbol == 'C'):
                        if mol.bond(ea_label, na_ea.label).bond_type == 'Single':
                            return 'saturated'
            # return other type
            return 'other'


    def define_substructures(self, *, define_method='con', con_path=None, SMILES_list=None):

        if define_method == 'con' and isinstance(con_path, str):
            Conner_list = os.listdir(con_path)
            for Conner in Conner_list:
                Conner_path = os.path.join(con_path, Conner)
                connser_substructure = search.ConnserSubstructure(Conner_path)
                self.dict_substructure[Conner[0:-4]] = connser_substructure
        elif define_method == 'SMILES' and isinstance(SMILES_list,list) and any(s for s in SMILES_list if s is str):
            for s in SMILES_list:
                SMARTS_substructure = search.SMARTSSubstructure(s)
                self.dict_substructure[s] = SMARTS_substructure
        else:
            raise ImportError('define_method include "con" and "SMILES",'
                              'and then You should give the corresponding path or SMILES list')
        print('%s substructure define successful!' % len(self.dict_substructure))








