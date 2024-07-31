import numpy as np
import mypackage

def test_grid_function():
    """
    tests for find_grid_cell function
    """
    atom_pos = np.array([1.5, 2.5, 3.5])
    boundaries = np.array([[0, 5], [0, 5], [0, 5]])
    discretization = np.array([5, 5, 5])
    lengths = np.array([5, 5, 5])

    expected_result = np.array([1, 2, 3])
    # Call the function
    test_result = mypackage.find_grid_cell(atom_pos, boundaries, discretization, lengths)

    assert np.array_equal(test_result, expected_result), f"Test failed: {test_result} != {expected_result}"
    
