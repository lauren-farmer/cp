import random
import math
import networkx as nx
from ortools.sat.python import cp_model

#constants for validating the burning sequence (labelling the vertices)
BURN = "burn"
OPEN = "open"

#perfoms a spread step in the burning process
def _do_a_spread(graph, state_dict):
  new_burns = []
  for vertex in graph.nodes():
    for other in graph.neighbors(vertex):
      if state_dict[other] == BURN:
        new_burns.append(vertex)
  for vertex in new_burns:
    state_dict[vertex] = BURN

#validates whether a given burning sequence burns all vertices in the graph
def _is_a_burning_seq(graph, burning_seq):
  state_dict = {}
  for vertex in graph.nodes():
    state_dict[vertex] = OPEN
  #for every time step in the burning sequence
  for time in range(len(burning_seq)):
    _do_a_spread(graph, state_dict)
    state_dict[burning_seq[time]] = BURN

  #allowing an extra spread step after the last ignition
  _do_a_spread(graph, state_dict)

  #checking if all vertices are burned
  for vertex in graph.nodes():
    if state_dict[vertex] != BURN: #if not burned
      return False #sequence is invalid
  return True

#Build and solve the CSP1 model for a fixed number of rounds B.
#Returns (True, burn_seq) if feasible, else (False, None).
#burn_seq is a list of vertices chosen to ignite at each round 1..B.
def solve_csp1_for_B(G, B, timeout_ms = None, workers=8):
  n = G.number_of_nodes() #number of vertices
  if n == 0: #edge case: empty graph
    return True, [] #trivially feasible with empty burning sequence
  
  nodes = list(G.nodes()) #list of graph nodes
  idx_of = {node: i for i, node in enumerate(nodes)}  #map from node label to index 0..n-1
  node_of = {i: node for i, node in enumerate(nodes)} #map from index 0..n-1 to node label

  # Precompute adjacency 
  #neighbours[i]: list of indices of neighbours of node i
  neighbours = {i: [idx_of[u] for u in G.neighbors(node_of[i])] for i in range(n)}

  #create CP-SAT solver
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

  status = solver.Solve(model) #solve the CSP

  #If neither feasible nor optimal, return infeasible
  if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    return False, None
  
  #Extract burning sequence from the decsion variables
  burn_seq = [] #list of chosen vertices to burn at rounds 1..B
  for j in range(1, B+1): 
    chosen = None
    for i in range(n):
      if solver.Value(decision[i,j]):
        chosen = node_of[i] #translate index back to node label
        break  
    burn_seq.append(chosen)
  return True, burn_seq

#Binary search over B to find the minimum burning number
#Uses binary search over B and CSP1 to find the minimum burning number
#Returns a dictionary with key 'burn_seq' where burn_seq is the optimal burning sequence (list of vertices in ignition order)
def run_ilp(instance_graph, timeout= 1000):
  G = nx.Graph(instance_graph) #ensures simple undirected graph
  n = G.number_of_nodes() #number of vertices

  #handling small graphs directly
  if n == 0: #empty graph
    return {'burn_seq': []}
  if n == 1: #single vertex graph
    return {'burn_seq': list(G.nodes())}
  
  #theoretical bounds on burning number
  upper_bound = math.ceil(math.sqrt(n))
  lower_bound = 1

  best_seq = None

  #binary search over B
  while lower_bound <= upper_bound:
    B= (lower_bound + upper_bound) // 2 #midpoint
    feasible, seq = solve_csp1_for_B(G, B, timeout_ms=timeout)
    #Accept B if solver finds a solution AND the sequence actually burns the entire graph
    if feasible and seq is not None and _is_a_burning_seq(G, seq):
      best_seq = seq #last feasible sequence found
      upper_bound = B - 1 #try smaller B
    else:
      lower_bound = B + 1 #if more rounds needsed, try larger B
  #no valid sequence found
  if best_seq is None:
    return None
  #othererwise return the best sequence found
  return {'burn_seq': best_seq}