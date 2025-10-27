import networkx as nx
from ortools.linear_solver import pywraplp

# THIS FILE IS WHERE STUDENTS SHOULD DO THEIR WORK

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
  
  # Ensure undirected simple graph
  G = nx.Graph(instance_graph)
  nodes = list(G.nodes())
  n = len(nodes)

  # Map between node labels and indices
  idx_of = {node: i for i, node in enumerate(nodes)}

  # Handle trivial case k < 0
  if distance < 0:
    raise ValueError("Distance must be non-negative")
    
  solver = pywraplp.Solver.CreateSolver('SCIP')
  if not solver:
      return None
  if timeout is not None and timeout > 0:
    solver.SetTimeLimit(int(timeout))
    
  x = [solver.IntVar(0, 1, f'x_{i}') for i in range(n)]  # 1 if v is in dominating set, 0 otherwise

  neighbourhoods = [] 
  # Constraints: every vertex must be dominated
  for v in nodes:
    reachable = nx.single_source_shortest_path_length(G, v, cutoff=distance)
    neighbourhoods.append(set(reachable.keys()))

  # For each vertex v: must be dominated by at least one chosen node
  for i, neigh in enumerate(neighbourhoods):
    solver.Add(solver.Sum([x[idx_of[u]] for u in neigh]) >= 1)
  
  # Objective: minimize size of dominating set
  solver.Minimize(solver.Sum(x))
    
  # Solve
  status = solver.Solve()
    
  if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
    chosen_nodes = [nodes[i] for i in range(n) if x[i].solution_value() > 0.5]
    return {'dom_set': chosen_nodes}
  else:
      return None
