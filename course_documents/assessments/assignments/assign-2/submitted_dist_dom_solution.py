import random
import networkx as nx
import ortools
from ortools.linear_solver import pywraplp

# THIS FILE IS WHERE STUDENTS SHOULD DO THEIR WORK
def get_distance_matrix(graph, k):
    #Compute distance matrix up to k steps
    nodes = list(graph.nodes())
    dist_matrix = {}
    for source in nodes:
        lengths = nx.single_source_shortest_path_length(graph, source, cutoff=k)
        for target in nodes:
          dist_matrix[(source, target)] = 1 if target in lengths else 0    
    return dist_matrix

#
# This function should run your ILP implementation
# for distance dominating set 
# - instance_graph will be a networkx graph object
# - distance is the distance over which vertices can dominate
# - timeout is the maximum time in ms you should let the search run for
# The function should return a dictionary that has 'dom_set' mapped to 
# a list with the vertex names (as they are in instance_graph)
# that are chosen by your solution
# or None if the model does not halt in the time allowed
# (The dictionary structure is so you can return other things if it's 
# useful for your debugging)
def run_ilp(instance_graph, distance = 1, timeout=1000):
  #  in here you can modify the graph to get whatever format you need, implement your ILP, call your solver
  #  and then translate the result back into a set of nodes from instance_graph   
  
  #  This is obviously not a solution, but just me choosing a single vertex from the graph
  
  
  solver = pywraplp.Solver.CreateSolver('SCIP')
  if not solver:
      return None
  solver.SetTimeLimit(timeout)
    
  # Get nodes and compute distance matrix
  nodes = list(instance_graph.nodes())
  dist_matrix = get_distance_matrix(instance_graph, distance)
    
  x = {v:solver.IntVar(0, 1, f'x_{v}') for v in nodes}  # 1 if v is in dominating set, 0 otherwise
    
  # Constraints: every vertex must be dominated
  for v in nodes:
    # Sum of x[u] for all u that can reach v within distance k must be >= 1
    solver.Add(solver.Sum(x[u] * dist_matrix[(u, v)] for u in nodes) >= 1,)
    
  # Objective: minimize size of dominating set
  solver.Minimize(solver.Sum(x[v] for v in nodes))
    
  # Solve
  status = solver.Solve()
    
  if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
      dom_set = [v for v in nodes if x[v].solution_value() > 0.5]
      return {'dom_set': dom_set}
  else:
      return None
