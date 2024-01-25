from Bio.PDB import *
import numpy as np
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.filterwarnings("ignore", category=PDBConstructionWarning)
import numpy as np

bohr_to_angstrom = 1.8897161646321
def get_ele_coord(xyzfile):
    ele_list=[]
    coord_list=[]
    with open(xyzfile, 'r') as file:
        lines = file.readlines()[2:]
        for line in lines:
            if line.strip() != '':
                ele=line.split()[0]
                coordx=float(line.split()[1])
                coordy=float(line.split()[2])
                coordz=float(line.split()[3])
                coords=[coordx,coordy,coordz]
                ele_list.append([ele])
                coord_list.append(coords)
    return np.array(ele_list), np.array(coord_list)*bohr_to_angstrom

def write_xyz(coordinates, ele, ene, file, type = 'a'):
    """
    Write the coordinates into a xyz file.
    """
    header = f'{len(coordinates)}\n{ene}\n'
    xyz = np.concatenate((ele, coordinates), axis=1)
    with open(file, type) as f:
        f.write(header)
        for e, coord in zip(ele, coordinates):
            line = '{:<2}     {: .10}     {: .10}     {: .10}\n'.format(e[0], *coord)
            f.write(line)
            
def get_dimer2(max,dimer1,dimer2,atom1,atom2,dis):

    with open(max, 'r') as file1:
        lines = file1.readlines()
        atom1_line = lines[atom1+2]
        atom2_line = lines[atom2+2]
    #get vector
    atom1_coord = np.array([float(i) for i in atom1_line.split()[1:4]])
    atom2_coord = np.array([float(i) for i in atom2_line.split()[1:4]])
    vector = atom2_coord - atom1_coord
    #Get the vector required to move dist along the vector direction
    vector_unit = vector/np.linalg.norm(vector)
    vector_unit = vector_unit*dis
    #plus vector_unit for xyz1 to get xyz2 and save
    dimer1_atom1_coord = atom1_coord + vector_unit
    dimer2_atom1_coord = atom1_coord - vector_unit
    dimer1_atom2_coord = atom2_coord - vector_unit
    dimer2_atom2_coord = atom2_coord + vector_unit
    with open(max, 'r') as x0, open(dimer1, 'w') as x1, open(dimer2, 'w') as x2:
        lines = x0.readlines()
        for counter,line in enumerate(lines):
            if counter == atom1+2:
                line1 = '{:<2}     {: .10}     {: .10}     {: .10}\n'.format(line.split()[0], *dimer1_atom1_coord)
                x1.write(line1)
                line2 = '{:<2}     {: .10}     {: .10}     {: .10}\n'.format(line.split()[0], *dimer2_atom1_coord)
                x2.write(line2)
            elif counter == atom2+2:
                line1 = '{:<2}     {: .10}     {: .10}     {: .10}\n'.format(line.split()[0], *dimer1_atom2_coord)
                x1.write(line1)
                line2 = '{:<2}     {: .10}     {: .10}     {: .10}\n'.format(line.split()[0], *dimer2_atom2_coord)
                x2.write(line2)
            else:
                line = line
                x1.write(line)
                x2.write(line)

def get_nuclear_charges(ele):
    nc_dict = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,
               'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
               'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36}
    nuclear_charges = np.array([nc_dict[i[0]] for i in ele])
    return nuclear_charges