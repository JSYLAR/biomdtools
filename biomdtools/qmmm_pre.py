from Bio.PDB import *
from scipy.spatial import cKDTree

class QMMM_PRE:
    def __init__(self,targetfile='',select_dis='',z=1,path='./',save_path='./',save=False):
        self.io = PDBIO()
        self.select_dis=select_dis
        self.z=z
        self.save=save
        self.path=path
        self.save_path=save_path
        parser = PDBParser(PERMISSIVE=True)
        self.structure_name = (str(targetfile).strip('pdb'))
        self.pdb = parser.get_structure(self.structure_name,targetfile)
        self.resList = Selection.unfold_entities(self.pdb, 'R')
    def get_select_atom(self,qm_dict):
        select_atom_list = []
        select_idx_list = []
        for res_idx,atom_list in qm_dict.items():
            if 'all' in atom_list:
                for atom in self.resList[res_idx-1]:
                    select_atom_list.append(atom)
            elif 'side_chain' in atom_list:
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in ['CA','C','N','O','H','HA','H1','H2','H3']:
                        continue
                    select_atom_list.append(atom)
            else:
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in atom_list:
                        select_atom_list.append(atom)

        select_atom_list=list(set(select_atom_list))
        
        self.io.set_structure(self.pdb)
        class Select_pdb(Select):
            def accept_atom(self, atom):
                if atom in select_atom_list:
                    return 1
                else:
                    return 0
        if self.save==True:
            self.io.save(f"{self.save_path}/qm_region.pdb",Select_pdb())
        set(select_atom_list)
        for i in select_atom_list:
            select_idx_list.append(i.get_serial_number()-self.z)
        select_idx_list = list(set(select_idx_list))
        select_idx_list.sort()
        return select_idx_list
    
    def get_active_atom(self,qm_atom_list):
        res_list = []
        active_atom_list = []
        atomList = Selection.unfold_entities(self.pdb, 'A')
        coords = [atom.get_coord() for atom in atomList]
        kdtree = cKDTree(coords)
        for atom in atomList:
            if atom.get_serial_number()-self.z in qm_atom_list:
                coor = atom.get_coord()
                indices = kdtree.query_ball_point(coor, self.select_dis)
                ns = [atomList[i].get_parent() for i in indices]
                res_list += ns
        re_list = list(set(res_list))

        self.io.set_structure(self.pdb)
        class Select_pdb(Select):
            def accept_residue(self, residue):
                if residue in re_list:
                    return 1
                else:
                    return 0
        if self.save == True:
            self.io.save(f"{self.save_path}/active_region.pdb", Select_pdb())
        for res in re_list:
            for atom in res.get_atoms():
                active_atom_list.append(atom.get_serial_number()-self.z)
        active_atom_list.sort()
        return active_atom_list

    def get_reaction_atom(self,reac_dict):
        select_atom_list = []
        select_idx_list = []
        for res_idx,atom_list in reac_dict.items():
            if 'all' in atom_list:
                for atom in self.resList[res_idx-1]:
                    select_atom_list.append(atom)
            elif 'side_chain' in atom_list:
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in ['CA','C','N','O','H','HA','H1','H2','H3']:
                        continue
                    select_atom_list.append(atom)
            else:
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in atom_list:
                        select_atom_list.append(atom)

        select_atom_list=list(set(select_atom_list))
        
        self.io.set_structure(self.pdb)
        class Select_pdb(Select):
            def accept_atom(self, atom):
                if atom in select_atom_list:
                    return 1
                else:
                    return 0
        set(select_atom_list)
        for i in select_atom_list:
            select_idx_list.append(i.get_serial_number()-self.z)
        select_idx_list = list(set(select_idx_list))
        select_idx_list.sort()
        return select_idx_list

    def get_constraint(self,cos_list):
        constraint_block = '  Constraints \n'
        for cos_pair in cos_list:
            cos_dis = str(cos_pair[2])
            cos_atom1 = cos_pair[0]
            for res_idx,atom_tar in cos_atom1.items():
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in atom_tar:
                        atom_idx_1 = str(atom.get_serial_number()-self.z)
            cos_atom2 = cos_pair[1]
            for res_idx,atom_tar in cos_atom2.items():
                for atom in self.resList[res_idx-1]:
                    if atom.get_name() in atom_tar:
                        atom_idx_2 = str(atom.get_serial_number()-self.z)
            constraint = '    {B '+atom_idx_1+' '+atom_idx_2+' '+cos_dis+' C}\n'
            constraint_block += constraint
        constraint_block += '    end\n'
        return constraint_block
