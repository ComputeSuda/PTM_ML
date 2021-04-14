# -*- coding: utf-8 -*-
"""
Created on Fri May 17 15:32:31 2019

@author: Jixiao he
"""

import sys,os
from prody import *
from pylab import *
from matplotlib.pylab import *

# ZERO = 1e-6
ioff()


list_file = 'pdb_list.txt'
pdb_path = './'
out_path = './output/'
slow_num = 20

print ("""Usage:
./enm_batch [list_file] [pdb_path] [out_path]

This program read pdb file list from a text file as input. 
Outputs:
GNM defined eigenvectors and squared fluctuations for slow modes;
GNM defined PRS and Cross-Correlations;
ANM defined PRS;
Tittint and commute time.""")

if len(sys.argv) == 2:
    list_file = sys.argv[1]
elif len(sys.argv) == 3:
    list_file = sys.argv[1]
    pdb_path = sys.argv[2]
elif len(sys.argv) == 4:
    list_file = sys.argv[1]
    pdb_path = sys.argv[2]
    out_path = sys.argv[3]
elif len(sys.argv) == 5:
    list_file = sys.argv[1]
    pdb_path = sys.argv[2]
    out_path = sys.argv[3]
    slow_num = sys.argv[4]
else:
    print ("No arguments, using default")


if not os.path.exists(out_path):
    try:
        os.makedirs(out_path)
    except:
        raise IOError('Error: failed to create output folder ' + out_path)

try:
    f_list = open(list_file, 'r')
except IOError:
    raise IOError('failed to open ' + list_file + '.')
for line in f_list:
    pdb_name = line.strip()
    # pdb_name = '3HHR.pdb'
    print ('\n' + pdb_name)
    pdb_file = pdb_path + pdb_name
    try:
        pdb_name = pdb_name[:pdb_name.rindex('.')]
    except:
        pass
    structure = parsePDB(pdb_file)
    # calphas = structure.select('calpha')
    calphas = structure.select('not ion and name CA')
    print (str(calphas.numAtoms()) + ' nodes are selected.')
    gnm = GNM(pdb_name)
    gnm.buildKirchhoff(calphas, cutoff=10.0, gamma=1)
    
    resnum = gnm.numAtoms()
    
    # n_modes_n = 20
    n_modes_n = None
    
    gnm.calcModes(n_modes=n_modes_n)
    n_modes = gnm._n_modes
    
    gnm_Eigvecs_slow = gnm[0:slow_num].getEigvecs()
    gnm_Eigvecs_all = gnm.getEigvecs()
    gnm_Eigvecs_1 = gnm[0].getEigvec()
    gnm_Eigvecs_2 = gnm[1].getEigvec()
    gnm_Eigvecs_3 = gnm[2].getEigvec()
    gnm_Eigvecs_1_to_2 = gnm[0:2].getEigvecs()
    gnm_Eigvecs_1_to_3 = gnm[0:3].getEigvecs()
    # msFlucts_slow=sqFlucts_10/gnm[0:slow_num].getVariances().sum()
    sqFlucts_all = calcSqFlucts(gnm)
    sqFlucts_1 = calcSqFlucts(gnm[0])
    sqFlucts_2 = calcSqFlucts(gnm[1])
    sqFlucts_3 = calcSqFlucts(gnm[2])
    sqFlucts_1_to_2 = calcSqFlucts(gnm[0:2])
    sqFlucts_1_to_3 = calcSqFlucts(gnm[0:3])
    
    CC_gnm = calcCrossCorr(gnm[0:slow_num])
    
    #prs_slow, effectiveness, sensitivity = calcPerturbResponse(model=gnm[0:slow_num], atoms=calphas)
    prs_all, effectiveness, sensitivity = calcPerturbResponse(model=gnm, atoms=calphas)
    
    resname_array = calphas.getResnames()
    resid_array = calphas.getResnums()    
     
    
    #out_array = concatenate((reshape(resname_array, (resnum, 1)), reshape(resid_array, (resnum, 1))), axis=1)
    #f_out_name = out_path + os.sep + pdb_name + '_resinfo.txt'
    #header_str = 'resname	resid	'
    #savetxt(f_out_name, out_array, fmt='%s', delimiter=' ', newline='\n', header=header_str, footer='', comments='')
    
    #out_array = concatenate((U_slow, reshape(sqFlucts_slow.round(8), (resnum, 1)), reshape(sqFlucts_all.round(8), (resnum, 1))), axis=1)
    f_out_name = out_path + os.sep + pdb_name + '_gnm_slow_mode.txt'
    #header_str = '	'.join(['slow_'+str(x+1) for x in range(slow_num)]) + '	sqFlucts_slow'+str(slow_num)+'	sqFlucts_all'
    savetxt(f_out_name, gnm_Eigvecs_slow, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_all_mode.txt'
    savetxt(f_out_name, gnm_Eigvecs_all, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_mode_1.txt'
    savetxt(f_out_name, gnm_Eigvecs_1, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_mode_2.txt'
    savetxt(f_out_name, gnm_Eigvecs_2, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_mode_3.txt'
    savetxt(f_out_name, gnm_Eigvecs_3, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_mode_1_to_2.txt'
    savetxt(f_out_name, gnm_Eigvecs_1_to_2, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_mode_1_to_3.txt'
    savetxt(f_out_name, gnm_Eigvecs_1_to_3, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_all_sq.txt'
    savetxt(f_out_name, sqFlucts_all, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_sq_1.txt'
    savetxt(f_out_name, sqFlucts_1, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='') 
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_sq_2.txt'
    savetxt(f_out_name, sqFlucts_2, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')     
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_sq_3.txt'
    savetxt(f_out_name, sqFlucts_3, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')     
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_sq_1_to_2.txt'
    savetxt(f_out_name, sqFlucts_1_to_2, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')    
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_sq_1_to_3.txt'
    savetxt(f_out_name, sqFlucts_1_to_3, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    #f_out_name = out_path + os.sep + pdb_name + '_gnm_slow_prs.txt'
    #savetxt(f_out_name, prs_slow[0], fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_all_prs_effector.txt'
    savetxt(f_out_name, effectiveness, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_all_prs_sensor.txt'
    savetxt(f_out_name, sensitivity, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_gnm_cc.txt'
    savetxt(f_out_name, CC_gnm, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    #f_out_name = out_path + os.sep + pdb_name + '_gnm_modes.nmd'
    #writeNMD(f_out_name, gnm[:20], calphas)
    #f_out_name = out_path + os.sep + pdb_name + '_gnm_prs_eff.txt'
    #savetxt(f_out_name, prs_gnm[1], fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    #f_out_name = out_path + os.sep + pdb_name + '_gnm_prs_sen.txt'
    #savetxt(f_out_name, prs_gnm[2], fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    anm = ANM(pdb_name)
    anm.buildHessian(calphas, cutoff=15.0, gamma=1)
    anm.calcModes(n_modes=n_modes_n) 
    n_modes = anm._n_modes
    
    anm_Eigvecs_slow = anm[0:slow_num].getEigvecs()
    anm_Eigvecs_all = anm.getEigvecs()
    sqFlucts_all_anm = calcSqFlucts(anm)
    CC_anm = calcCrossCorr(anm[0:slow_num])
    prs_all_anm, effectiveness, sensitivity = calcPerturbResponse(model=anm, atoms=calphas)
    anm_Eigvecs_1_to_3 = anm[0:3].getEigvecs()
    
    #2019/8/27
    stiffness = calcMechStiff(anm, calphas)
    
    f_out_name = out_path + os.sep + pdb_name + '_anm_stiffness.txt'
    savetxt(f_out_name, stiffness, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')

    f_out_name = out_path + os.sep + pdb_name + '_anm_mode_1_to_3.txt'
    savetxt(f_out_name, anm_Eigvecs_1_to_3, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_anm_slow_mode.txt'
    savetxt(f_out_name, anm_Eigvecs_slow, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_anm_all_mode.txt'
    savetxt(f_out_name, anm_Eigvecs_all, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_anm_all_sq.txt'
    savetxt(f_out_name, sqFlucts_all_anm, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')

    f_out_name = out_path + os.sep + pdb_name + '_anm_cc.txt'
    savetxt(f_out_name, CC_anm, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')

    f_out_name = out_path + os.sep + pdb_name + '_anm_all_prs_effector.txt'
    savetxt(f_out_name, effectiveness, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_anm_all_prs_sensor.txt'
    savetxt(f_out_name, sensitivity, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    #f_out_name = out_path + os.sep + pdb_name + '_anm_modes.nmd'
    #writeNMD(f_out_name, anm[:20], calphas)

    #f_out_name = out_path + os.sep + pdb_name + '_gnm_cc.png'
    # fig = showAtomicMatrix(CC, x_array=np.mean(CC, axis=0), y_array=np.mean(CC, axis=0), atoms=calphas)
    #fig = showAtomicMatrix(CC, atoms=calphas)
    #xlabel('Residues')
    #savefig(f_out_name, dpi=200)
    #plt.close()
    close()

f_list.close();




















