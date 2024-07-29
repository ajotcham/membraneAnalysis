import sys
import numpy as np
# sys.path.append('/Users/ajotcham/Desktop/MembraneAnalysis/membraneAnalysis/lammps')
from mypackage import readlammps

test_path = '../lammps/output/dump_pytest.lammpstrj'

def test_same():
    expected_info = np.array([[1., 1., 1., 1., 1., 1., 1., 1.],
                              [2., 2., 2., 2., 2., 2., 2., 2.]])
    
    test_info = readlammps.test_function(test_path)
    np.testing.assert_array_equal(test_info, expected_info)
    # expectedPos = np.array([[1., 1., 1., 1., 1.],
    #                     [2., 2., 2., 2., 2.]])
    
    # expectedVelo = np.array([[1., 1., 1., 1., 1.],
    #                      [2., 2., 2., 2., 2.]])
    
    # testP, testV = readlammps.test_function(test_path)
    # np.testing.assert_array_equal(testP, expectedPos) 
    # np.testing.assert_array_equal(testV, expectedVelo) 
    

