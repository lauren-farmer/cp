import random
import math
import networkx as nx
from ortools.sat.python import cp_model

BURN = "burn"
OPEN = "open"

def _do_a_spread(graph, state_dict):
  new_burns = []
  for vertex in graph.nodes():
    for other in graph.neighbors(vertex):
      if state_dict[other] == BURN:
        new_burns.append(vertex)
  for vertex in new_burns:
    state_dict[vertex] = BURN

def _is_a_burning_seq(graph, burning_seq):
  state_dict = {}
  for vertex in graph.nodes():
    state_dict[vertex] = OPEN

  for time in range(len(burning_seq)):
    _do_a_spread(graph, state_dict)
    state_dict[burning_seq[time]] = BURN
  # one extra spread after the last ignition
  _do_a_spread(graph, state_dict)

  for vertex in graph.nodes():
    if state_dict[vertex] != BURN:
      return False
  return True


# THIS FILE IS WHERE STUDENTS SHOULD DO THEIR WORK


#
# This function should run your ILP implementation
# for distance dominating set 
# - instance_graph will be a networkx graph object
# - distance is the distance over which vertices can dominate
# - timeout is the maximum time in ms you should let the search run for
# The function should return a dictionary that has 'burn_seq' mapped to 
# a list of vertices in order that your solution sets fires on them. 
# or None if the model does not halt in the time allowed
# That is: if your solution sets fires at vertex "a", then "b", then "c" and that ends the process, then 
# the dicitonary you return should include 'burn_seq': ["a", "b, "c"]
# It is possible you will want your model to have a different decision structure: you are
# welcome to compute in that structure and then translate. 
# (The dictionary structure is so you can return other things if it's 
# useful for your debugging)

def solve_csp1_for_B(G, B, timeout_ms = None, workers=8):
  n = G.number_of_nodes()
  if n == 0:
    return True, []
  
  nodes = list(G.nodes())
  idx_of = {node: i for i, node in enumerate(nodes)}
  node_of = {i: node for i, node in enumerate(nodes)}

  # Precompute adjacency 
  neighbours = {i: [idx_of[u] for u in G.neighbors(node_of[i])] for i in range(n)}

  model = cp_model.CpModel()
  
  # decision[i,j]: vertex i is ignited at round j (1..B)
  decision = {
    (i, j): model.NewBoolVar(f"decision_{i}_{j}")
    for i in range(n) for j in range(1, B + 1)
  }

  # burned[i,j]: vertex i is burned by end of round j (0..B)
  burned = {
    (i, j): model.NewBoolVar(f"burned_{i}_{j}")
    for i in range(n) for j in range(B + 1)
  }

  #Constraints
  #initial state with no burned vertices
  for i in range(n):
    model.Add(burned[i,0] == 0)

  #Constraint 24: If a vertex is burned at turn j-1, it remains burned at turn j
  for i in range(n):
    for j in range(1,B+1):
      model.Add(burned[i,j-1] <= burned[i,j])
  
  #Constraint 25: If a vertex is actively burned, it becomes burned (i.e. vertex is in the burning sequence)
  for i in range(n):
    for j in range(1, B+1):
      model.Add(decision[i,j] <= burned[i,j])

  #Constraint 26: If a vertex is adjacent to a burned vertex at turn j-1, it becomes burned at turn j
  #(spreading)
  for i in range(n):
    for j in range(1, B+1):
      for k in neighbours[i]:
        model.Add(burned[k,j-1] <= burned[i,j])

  #Constraint 27: Vertex can only be burned if it was already burning, was actively chosen or via spread from a burning neighbour
  #(preventing spontaneous combustion)
  #Can actively burn a vertex on turn 0
  for i in range(n):
    for j in range(1, B+1):
      model.Add(burned[i, j] <= burned[i, j-1] + decision[i,j] + sum(burned[k, j-1] for k in neighbours[i]))

  #Constraint 28: Ensure EXACTLY one vertex is actively burned at each turn
  # Must burn exactly one vertex per turn until all are burned
  for j in range(1, B+1):
    model.Add(sum([decision[i,j] for i in range(n)]) == 1)

  #Constraint 29 : All vertices must be burned by turn B
  model.Add(sum([burned[i, B] for i in range(n)]) == n)

  #Constraint 34 : A vertex can only be actively burned once
  for i in range(n):
    model.Add(sum(decision[i,j] for j in range(1, B+1)) <= 1)

  #A vertex can only be actively burned if not already burning
  for i in range(n):
    for j in range(1, B):
      model.Add(decision[i,j] + burned[i,j-1] <= 1)

  solver = cp_model.CpSolver()
  if timeout_ms: 
    solver.parameters.max_time_in_seconds = timeout_ms / 1000.0
  solver.parameters.num_search_workers = workers
  
  status = solver.Solve(model)
  if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    return False, None
  
  burn_seq = []
  for j in range(1, B+1):
    chosen = None
    for i in range(n):
      if solver.Value(decision[i,j]):
        chosen = node_of[i]
        break  
    burn_seq.append(chosen)
  return True, burn_seq

def run_ilp(instance_graph, timeout= 1000):
  G = nx.Graph(instance_graph)
  n = G.number_of_nodes()

  if n == 0:
    return {'burn_seq': []}
  if n == 1:
    return {'burn_seq': list(G.nodes())}
  
  #theoretical bounds on burning number
  upper_bound = math.ceil(math.sqrt(n))
  lower_bound = 1

  best_seq = None

  while lower_bound <= upper_bound:
    B= (lower_bound + upper_bound) // 2
    feasible, seq = solve_csp1_for_B(G, B, timeout_ms=timeout)
    if feasible and seq is not None and _is_a_burning_seq(G, seq):
      best_seq = seq
      upper_bound = B - 1
    else:
      lower_bound = B + 1
  if best_seq is None:
    return None
  return {'burn_seq': best_seq}