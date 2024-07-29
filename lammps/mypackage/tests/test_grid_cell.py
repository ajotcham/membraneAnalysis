import numpy as np
import mypackage

def test_same_grid():
    expected_result = np.array([1, 2, 3])

    test_result = mypackage.test_grid_function()
    assert np.array_equal(test_result, expected_result), f"Test failed: {test_result} != {expected_result}"
    
