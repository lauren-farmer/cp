import networkx as nx
from submitted_dist_dom_solution import run_ilp
from lecturer_code_sample_dist_dom import distance_dominates

def test_path_graphs():
    """Test path graphs with different distances"""
    print("\n=== Testing Path Graphs ===")
    path = nx.path_graph(6)
    
    test_cases = [
        (1, 2),  # (distance, expected_size)
        (2, 2),
        (5, 1)
    ]
    
    for k, expected_size in test_cases:
        result = run_ilp(path, distance=k)
        if result:
            dom_set = result['dom_set']
            is_valid = distance_dominates(path, dom_set, k)
            print(f"k={k}: dom_set size={len(dom_set)}, expected={expected_size}, valid={is_valid}")
            assert is_valid, f"Solution for k={k} is not valid"
            assert len(dom_set) <= expected_size, f"Solution size {len(dom_set)} exceeds expected size {expected_size}"
        else:
            print(f"k={k}: No solution found")

def test_complete_graphs():
    """Test complete graphs - should always need only one vertex"""
    print("\n=== Testing Complete Graphs ===")
    complete = nx.complete_graph(6)
    
    for k in [1, 2, 3]:
        result = run_ilp(complete, distance=k)
        if result:
            dom_set = result['dom_set']
            is_valid = distance_dominates(complete, dom_set, k)
            print(f"k={k}: dom_set size={len(dom_set)}, valid={is_valid}")
            assert is_valid, f"Solution for k={k} is not valid"
            assert len(dom_set) == 1, f"Complete graph should need only 1 vertex, got {len(dom_set)}"
        else:
            print(f"k={k}: No solution found")

def test_grid_graphs():
    """Test grid graphs with different distances"""
    print("\n=== Testing Grid Graphs ===")
    grid = nx.grid_2d_graph(3, 3)  # 3x3 grid
    
    test_cases = [
        (1, 4),  # (distance, expected_size)
        (2, 1),
        (3, 1)
    ]
    
    for k, expected_size in test_cases:
        result = run_ilp(grid, distance=k)
        if result:
            dom_set = result['dom_set']
            is_valid = distance_dominates(grid, dom_set, k)
            print(f"k={k}: dom_set size={len(dom_set)}, valid={is_valid}")
            assert is_valid, f"Solution for k={k} is not valid"
            assert len(dom_set) <= expected_size, f"Solution size {len(dom_set)} exceeds expected size {expected_size}"
        else:
            print(f"k={k}: No solution found")

def test_timeout():
    """Test timeout functionality"""
    print("\n=== Testing Timeout ===")
    # Create a larger graph that might timeout with very short timeout
    large_grid = nx.grid_2d_graph(10, 10)
    result = run_ilp(large_grid, distance=1, timeout=1)  # 1ms timeout
    print(f"Timeout test: {'Timed out as expected' if result is None else 'Did not timeout'}")

if __name__ == "__main__":
    test_path_graphs()
    test_complete_graphs()
    test_grid_graphs()
    test_timeout()