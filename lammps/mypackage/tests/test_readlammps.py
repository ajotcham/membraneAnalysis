import numpy as np
import mypackage

test_path = 'test_dumps/simple_dump.lammpstrj'

def test_read_atoms():
    """
    Tests if read_atoms function is working properly
    """

    expected_pos = np.array([[1., 1., 1.],
                            [2., 2., 2.]])
    expected_velo = np.array([[1., 1., 1.],
                            [2., 2., 2.]])
    expected_atoms = int(2)

    expected_bounds = np.array([[-20, 20], [-20, 20], [-20, 20]])
    
    test_pos, test_velo, test_map, test_atoms, test_bounds = mypackage.read_atoms(test_path)
    
    np.testing.assert_array_equal(test_pos, expected_pos)
    np.testing.assert_array_equal(test_velo, expected_velo)
    np.testing.assert_array_equal(test_bounds, expected_bounds)
    assert test_atoms == expected_atoms
