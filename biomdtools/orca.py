import subprocess
from Bio.PDB import *
import numpy as np
import os
import glob
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from biomdtools.tools import write_xyz
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

class ORCA:
    """
    Class for ORCA calculations.
    """
    def __init__(self, path = '/opt/software/Orca/orca_5_0_4_linux_x86-64_shared_openmpi411/', method = 'b3lyp', basis = 'def2-SVP', D = 'D4', charge = 0, mult = 1, nprocs = 10):
        """
        Initialize the ORCA class.
        """
        self.path = path
        self.method = method
        self.basis = basis
        self.D = D
        self.charge = charge
        self.mult = mult
        self.nprocs = nprocs
        if self.mult > 1:
            self.HFTyp = 'UHF'
            self.check_stb = (f'  STABPerform true\n'
                              f'  STABRestartUHFifUnstable true\n')
        else:
            self.HFTyp = 'RHF'
            self.check_stb = ''

    def engrad(self, coordinates, ele, qm_atom):
        if os.path.exists('engrad.xyz'):
            #remove
            os.remove('engrad.xyz')
        write_xyz(coordinates, ele, '','engrad.xyz','w')
        open(f'engrad.inp','w').write(f'! QMMM\n'
            f'! {self.method} {self.basis} {self.D} ENGRAD BOHRS\n'
            f'! TightSCF \n'
            f'! RIJCOSX \n'
            f'! defgrid2 \n'
            f'%pal\n'
            f'  nprocs {self.nprocs}\n'
            f'end\n'
            f'%maxcore 2000\n'
            f'%scf\n'
            f'  MOInp "engrad.gbw"\n'
            f'  HFTyp {self.HFTyp}\n'
            f'  DirectResetFreq 1\n'
            f'  DIIS Start 0.1 MaxIt 5 MaxEq 20 BFac 1.2 MaxC 15.0 end \n'
            f'  AutoStart true\n'
            f'  SCFMode Direct\n'
            f'  MaxIter 500\n'
            f'end\n'                            
            f'%qmmm\n'
            'QMAtoms {'+f'{qm_atom}'+'}end\n'
            f'ORCAFFFilename "prmtop.ORCAFF.prms" end\n'
            f'*xyzfile {self.charge} {self.mult} engrad.xyz\n')
        subprocess.run(f'{self.path}/orca engrad.inp "--oversubscribe" > engrad.out',stdin=None, input=None, stdout=None, stderr=None, timeout=None, check=False, universal_newlines=False,shell=True)

        if os.path.exists('engrad.engrad'):
            print('Orca run successful')
        else:
            print('Orca run failed')
            self.rm(keep = ['engrad.out','engrad.inp'])
            os._exit(1)

    def result(self):
        with open('engrad.engrad', 'r') as file:
            gradient_data = []
            reading_gradient = False
            for line in file:
                if line.startswith("# Number of atoms"):
                    next(file)
                    atom_num = int(next(file).split()[-1])
                    grad_num = atom_num*3
                elif line.startswith("# The current total energy in Eh"):
                    next(file)
                    total_energy = float(next(file).split()[-1])
                elif line.startswith("# The current gradient in Eh/bohr"):
                    reading_gradient = True
                    gradient_block = []
                    # Skip the line with the number of atoms
                    next(file)
                elif reading_gradient:
                    # Read the gradient data for each atom
                    gradient_block.append(float(line))
                    # Check if we have read the required number of lines
                    if len(gradient_block) == grad_num:
                        reading_gradient = False
                        gradient_data.append(gradient_block)
        self.rm()
        gradient_data = np.array(gradient_data)
        return total_energy, gradient_data
    def rm(self,keep = []):
        files_to_keep = ['engrad.gbw','engrad.out']+keep
        files_to_delete = [file for file in glob.glob('engrad*') if file not in files_to_keep]
        for file in files_to_delete:
                os.remove(file)
        os.rename('engrad.out', 'engrad.out_old')
    def IRC(self,xyz,qm_atom):
        open(f'IRC.inp','w').write(f'! QMMM numfreq IRC\n'
            f'! {self.method} {self.basis} {self.D}  \n'
            f'! TightSCF \n'
            f'! RIJCOSX \n'
            f'! defgrid2 \n'
            f'%pal\n'
            f'  nprocs {self.nprocs}\n'
            f'end\n'
            f'%maxcore 2000\n'
            f'%scf\n'
            f'  MOInp "engrad.gbw"\n'
            f'  HFTyp {self.HFTyp}\n'
            f'  DirectResetFreq 1\n'
            f'  DIIS Start 0.1 MaxIt 5 MaxEq 20 BFac 1.2 MaxC 15.0 end \n'
            f'  AutoStart true\n'
            f'  SCFMode Direct\n'
            f'  MaxIter 500\n'
            f'end\n'                            
            f'%qmmm\n'
            'QMAtoms {'+f'{qm_atom}'+'}end\n'
            'ActiveAtoms {'+f'{qm_atom}'+'}end\n'
            f'ORCAFFFilename "prmtop.ORCAFF.prms" end\n'
            f'*xyzfile {self.charge} {self.mult} {xyz}\n')
        subprocess.run(f'{self.path}/orca IRC.inp "--oversubscribe" > IRC.out',stdin=None, input=None, stdout=None, stderr=None, timeout=None, check=False, universal_newlines=False,shell=True)
    def LOPT(self,xyz,qm_atom,activ_atom):
        xyzname = xyz.strip(".xyz")
        open(f'LOPT.inp','w').write(f'! QMMM l-opt\n'
            f'! {self.method} {self.basis} {self.D}  \n'
            f'! TightSCF \n'
            f'! RIJCOSX \n'
            f'! defgrid2 \n'
            f'%pal\n'
            f'  nprocs {self.nprocs}\n'
            f'end\n'
            f'%maxcore 2000\n'
            f'%scf\n'
            f'  MOInp "{xyzname}.gbw"\n'
            f'  HFTyp {self.HFTyp}\n'
            f'  DirectResetFreq 1\n'
            f'  DIIS Start 0.1 MaxIt 5 MaxEq 20 BFac 1.2 MaxC 15.0 end \n'
            f'  AutoStart true\n'
            f'  SCFMode Direct\n'
            f'  MaxIter 500\n'
            f'end\n'         
            f'%geom\n' 
            f'  MaxIter 500\n' 
            f'end\n' 
            f'%qmmm\n'
            'QMAtoms {'+f'{qm_atom}'+'}end\n'
            'ActiveAtoms {'+f'{activ_atom}'+'}end\n'
            f'ORCAFFFilename "prmtop.ORCAFF.prms" end\n'
            f'*xyzfile {self.charge} {self.mult} {xyz}\n')
        subprocess.run(f'{self.path}/orca LOPT.inp "--oversubscribe" > LOPT.out',stdin=None, input=None, stdout=None, stderr=None, timeout=None, check=False, universal_newlines=False,shell=True)

    def ENE(self,xyz,qm_atom):
            xyzname = xyz.strip(".xyz")
            open(f'ENE.inp','w').write(f'! QMMM\n'
                f'! {self.method} {self.basis} {self.D}  \n'
                f'! TightSCF \n'
                f'! RIJCOSX \n'
                f'! defgrid2 \n'
                f'%pal\n'
                f'  nprocs {self.nprocs}\n'
                f'end\n'
                f'%maxcore 2000\n'
                f'%scf\n'
                f'  MOInp "{xyzname}.gbw"\n'
                f'  HFTyp {self.HFTyp}\n'
                f'  DirectResetFreq 1\n'
                f'  DIIS Start 0.1 MaxIt 5 MaxEq 20 BFac 1.2 MaxC 15.0 end \n'
                f'  {self.check_stb}'
                f'  AutoStart true\n'
                f'  SCFMode Direct\n'
                f'  MaxIter 500\n'
                f'end\n'         
                f'%geom\n' 
                f'  MaxIter 500\n' 
                f'end\n' 
                f'%qmmm\n'
                'QMAtoms {'+f'{qm_atom}'+'}end\n'
                f'ORCAFFFilename "prmtop.ORCAFF.prms" end\n'
                f'*xyzfile {self.charge} {self.mult} {xyz}\n')
            subprocess.run(f'{self.path}/orca ENE.inp "--oversubscribe" > ENE.out',stdin=None, input=None, stdout=None, stderr=None, timeout=None, check=False, universal_newlines=False,shell=True)

