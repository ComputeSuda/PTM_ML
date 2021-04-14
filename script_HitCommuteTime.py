# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:58:53 2020

@author: Jixiao He
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 17 15:32:31 2019

@author: jx he
"""

from prody import *
from matplotlib.pylab import *
from IT_HitCommute import *

# ZERO = 1e-6
ioff()


list_file = 'pdb_list.txt'
pdb_path = './'
out_path = './output_6sdc/'
slow_num = 20


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
    calphas = structure.select('calpha')
    print (str(calphas.numAtoms()) + ' nodes are selected.')
    gnm = GNM(pdb_name)
    gnm.buildKirchhoff(calphas)
    K = gnm.getKirchhoff()  
    hc = IT_HitCommute(K)
    H = hc.buildHitTimes(K)  #hitting time matrix
    C = hc.buildCommuteTimes()  #commute time matrix
    
    f_out_name = out_path + os.sep + pdb_name + '_HitTimes.txt'
    savetxt(f_out_name, H, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')
    
    f_out_name = out_path + os.sep + pdb_name + '_CommuteTimes.txt'
    savetxt(f_out_name, C, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='')

    close()

f_list.close();




















