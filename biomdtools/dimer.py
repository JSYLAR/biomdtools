from Bio.PDB import *
from libdlfind.callback import (
    dlf_get_gradient_wrapper,
    dlf_put_coords_wrapper,
    make_dlf_get_params,)
from libdlfind import dl_find
from typing import List
import functools
import numpy as np
from numpy.typing import NDArray
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.filterwarnings("ignore", category=PDBConstructionWarning)
from pebble import concurrent
from biomdtools.tools import write_xyz,get_nuclear_charges

bohr_to_angstrom = 1.8897161646321
@dlf_get_gradient_wrapper
def engrad_func(coordinates: NDArray[np.float_], iimage: int, kiter: int, orca, ele, qm_atom: str):
    orca.engrad(coordinates, ele, qm_atom)
    energy, gradient = orca.result()
    return energy, gradient

@dlf_put_coords_wrapper
def store_results(
    switch: int,
    energy: float,
    coordinates: NDArray[np.float_],
    iam: int,
    energy_final,
    ele,
    qm_list: List[int],
    active_list: List[int],
) -> None:
    energy_final = energy
    coord = coordinates/bohr_to_angstrom
    """Store results from optimization."""
    if switch == 2:
        write_xyz(coord, ele, f'energy: {energy}','last_dimer_full.xyz','w')
        write_xyz(coord[qm_list], ele[qm_list], f'energy: {energy}','trj_qm.xyz')
        write_xyz(coord[qm_list], ele[qm_list], f'energy: {energy}','last_qm.xyz','w')
        write_xyz(coord[active_list], ele[active_list], f'energy: {energy}','trj_active.xyz')
        write_xyz(coord[active_list], ele[active_list], f'energy: {energy}','last_active.xyz','w')

@concurrent.process
def DIMER(orca, coord, coord2, qm_list, active_list, reac_list, ele, restart = 0, maxcycle = 1000):
    atom_num = len(coord)
    qm_atom = (' '.join('%s' %id for id in qm_list))
    active = np.full((atom_num), -1) #Frozen
    active[active_list] = 0 #Active region
    inner = np.full((atom_num), 0) #outer region
    #inner[qm_list] = 1 #inner region of Microiterative optimisation 
    nuclear_charges = get_nuclear_charges(ele) #nuclear charge
    spec = np.concatenate((active,nuclear_charges,inner), axis=0)#
    weight = np.full((atom_num), 1)
    weight[reac_list] = 2
    coord2 = coord2.reshape(-1)
    coord2 = np.concatenate((coord2,weight), axis=-1)
    energy_final = 0
    orca = orca
    dlf_get_params = make_dlf_get_params(coords=coord, 
                                         coords2 = coord2, 
                                         spec = spec, 
                                         nz = atom_num, 
                                         nweight = atom_num, 
                                         icoord = 220,         #Dimer 
                                         iopt = 3,             #LBFGS
                                         nframe = 1,
                                         printl = 2, 
                                         iline = 1,            #trustradius 0 const 1 energy 2 gradient
                                         delta = 0.01,
                                         tolerance = 0.0003, 
                                         maxstep = 0.5, 
                                         maxcycle = maxcycle,
                                         restart = restart,    #0 False 1 True
                                         dump = 1,)            #save ckpt every step
    dlf_get_gradient = functools.partial(engrad_func, orca = orca,qm_atom = qm_atom, ele = ele)
    dlf_put_coords = functools.partial(store_results, 
                                    energy_final = energy_final,
                                    ele = ele,
                                    qm_list=qm_list,
                                    active_list=active_list,)
    dl_find(nvarin=atom_num * 3,
            nvarin2=atom_num * 4,
            nspec = atom_num*3,
            dlf_get_gradient = dlf_get_gradient,
            dlf_get_params = dlf_get_params,
            dlf_put_coords = dlf_put_coords,)
    return energy_final