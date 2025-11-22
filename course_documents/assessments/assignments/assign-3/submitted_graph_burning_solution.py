import random
import networkx as nx
from ortools.linear_solver import pywraplp

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
def run_ilp(instance_graph, timeout=1000):
  #  in here you can modify the graph to get whatever format you need, implement your ILP, call your solver
  #  and then translate the result back into a set of nodes from instance_graph   

  #ensure undirected simple graph
  G = nx.Graph(instance_graph)
  nodes = list(G.nodes())
  n = len(nodes)

  #handle trivial cases
  if n == 0:
    return {'burn_seq': []}
  if n == 1:
    return {'burn_seq': nodes}
  
  #map between node labels and indices
  idx_of = {node: i for i, node in enumerate(nodes)}
  node_of = {i: node for i, node in enumerate(nodes)}

  #precompute adjcency for efficiency
  neighbours = [[] for _ in range(n)] 
  for i in range(n):
    for j in range(n):
      if i != j and G.has_edge(node_of[i], node_of[j]):
        neighbours[i].append(j)

  #create solver
  solver = pywraplp.Solver.CreateSolver('SCIP')

  #Decision variables
  #decision[i][j] = 1 if vertex i is set on fire at step j, else 0
  #burned[i][j] = 1 if vertex i is burned at the end of step j, else 0
  decision = {}
  burned = {}
    
  for i in range(n):
    for j in range(n):
      decision[i, j] = solver.IntVar(0, 1, f'decision_{i}_{j}')
      burned[i, j] = solver.IntVar(0, 1, f'burned_{i}_{j}')

  #Constraints
  #Constraint 24: If a vertex is burned at turn j-1, it remains burned at turn j
  for i in range(n):
    for j in range(1,n):
      solver.Add(burned[i,j-1] <= burned[i,j])
  
  #Constraint 25: If a vertex is actively burned, it becomes burned (i.e. vertex is in the burning sequence)
  for i in range(n):
    for j in range(n):
      solver.Add(decision[i,j] <= burned[i,j])

  #Constraint 26: If a vertex is adjacent to a burned vertex at turn j-1, it becomes burned at turn j (spreading)
  for i in range(n):
    for j in range(1, n):
      for k in neighbours[i]:
          solver.Add(burned[k,j-1] <= burned[i,j])

  #Constraint 27: Vertex can only be burned if it was already burning, was actively chosen or via spread from a burning neighbour
  #(preventing spontaneous combustion)
  #Can actively burn a vertex on turn 0
  for i in range(n): solver.Add(burned[i,0] <= decision[i,0])
  # At turn j > 0: via staying burned, active burn, or spread
  for i in range(n):
    for j in range(1, n):
    # If burned at turn j, must have been: burned before, actively burned, or has burned neighbor
      valid_ways = [burned[i, j-1], decision[i, j]]
      for k in neighbours[i]:
        valid_ways.append(burned[k, j-1])
            
      # burned[i,j] can only be 1 if at least one valid way is 1
      solver.Add(burned[i, j] <= solver.Sum(valid_ways))

  #Constraint 28: Ensure EXACTLY one vertex is actively burned at each turn
  #Unless it is the last turn and all are already burned
  for j in range(n):
    solver.Add(solver.Sum([decision[i,j] for i in range(n)]) <=1)
  # Must burn exactly one vertex per turn until all are burned
  for j in range(n):
    if j == 0:
      # First turn: must burn exactly 1
      solver.Add(solver.Sum([decision[i, j] for i in range(n)]) == 1)
    else:
      # Subsequent turns: burn 1 if any vertex is unburned, else 0
      # Check if all vertices burned at end of j-1
      sum_burned_prev = solver.Sum([burned[i, j-1] for i in range(n)])
      # If sum_burned_prev < n, then must burn 1 vertex
      # If sum_burned_prev == n, then must burn 0 vertices
      # Use big-M method
      M = n
      is_complete = solver.IntVar(0, 1, f'is_complete_{j}')
      
      # is_complete = 1 iff all vertices are burned
      solver.Add(sum_burned_prev >= n * is_complete)
      solver.Add(sum_burned_prev <= n - 1 + M * is_complete)
      
      # Must burn (1 - is_complete) vertices
      solver.Add(solver.Sum([decision[i, j] for i in range(n)]) == 1 - is_complete)

  #Constraint 29 : All vertices must be burned by turn n 
  for i in range(n):
    solver.Add(solver.Sum([burned[i, j] for j in range(n)]) >= 1)

  #Other Constraints
  #A vertex can only be actively burned once
  for i in range(n):
    solver.Add(solver.Sum([decision[i,j] for j in range(n)]) <= 1)

  #A vertex can only be actively burned if not already burning
  for i in range(n):
    for j in range(1, n):
      solver.Add(decision[i,j] <= 1 - burned[i,j-1])

  #burning_number: the number of turns it takes to burn the graph
  burning_number = solver.IntVar(1, n, 'burning_number')

  for j in range(n):
    turn_used = solver.IntVar(0, 1, f'turn_used_{j}')
    solver.Add(solver.Sum([decision[i,j] for i in range(n)]) >= turn_used)
    solver.Add(solver.Sum([decision[i,j] for i in range(n)]) <= turn_used)
    solver.Add(burning_number >= turn_used * (j + 1))
  
  # Objective: minimize burning_number
  solver.Minimize(burning_number)
    
  # Solve
  status = solver.Solve()
    
  if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
    # Extract the burning sequence
    burn_seq = []
    for j in range(n):
      for i in range(n):
        if decision[i,j].solution_value() > 0.5:
          burn_seq.append(node_of[i])
          break  # Only one vertex per turn
        
    return {'burn_seq': burn_seq}
  else:
    return None
