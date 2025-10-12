%hard code small distances, will merit part marks OR 
%look at the global constraints provided by MiniZinc OR 
%an accessory array that captures at what distance each vertex is reached by each other vertex

int: n;                     % number of vertices
int: m;                     % number of edges
int: k;                     % distance for domination

array[1..m] of int: from;    % edge start vertices
array[1..m] of int: to;      % edge end vertices

% Decision variable: 1 if vertex i is in the dominating set
array[1..n] of var 0..1: decision;

array[1..n, 1..n] of var 0..1: dist_reached; % dist_reached[i,j] = 1 if j is dist_reached from i within distance k
%need to enforce reaching distances with a constraint
%when should dist_reached[vertx, 1] be 1?
%- vertex is adjacent to something in decision
%when should it be dist_reached[vertx, 2] be 1?
%-if it was reached at 1 dist_reached[vertx, 1] = 1
%-if dist_reached[other_vertex, 1] = 1 and other_vertex is adjacent to vertex
%when should dist_reached[vertx, 3] be 1?
%when dist_reached[other_vertex, 2] = 1
%when one of its neighbours was reached at 2

% --- Build adjacency (1 if directly connected, or same vertex) ---
array[1..n, 1..n] of int: adj =
  array2d(1..n, 1..n,
    [ if i = j then 1
      else exists(e in 1..m)( (from[e] = i /\ to[e] = j) \/ (to[e] = i /\ from[e] = j) )
      endif
    | i in 1..n, j in 1..n ]);

% --- Compute reachability within distance ≤ k ---
% Floyd-Warshall-style iterative closure for up to distance k
%
% Initialize dist_reached to 1 for any pair within distance ≤ k
% For small graphs this is fine to express directly.
constraint
  forall(i in 1..n, j in 1..n)(
    dist_reached[i,j] = 
      if adj[i,j] = 1 then 1
      else if k > 1 then
        % There is a path of length ≤ k if there exists a neighbor mid
        % that is adjacent to i and dist_reached[mid,j] = 1 for distance k-1
        max([adj[i,mid] * dist_reached[mid,j] | mid in 1..n]) 
      else 0
      endif
  );

% --- Every vertex must be dominated ---
% Each vertex v must have some selected u within distance ≤ k
constraint
  forall(v in 1..n)(
    sum(u in 1..n)(decision[u] * dist_reached[u,v]) >= 1
  );

% --- Objective: minimize the number of selected vertices ---
solve minimize sum(v in 1..n)(decision[v]);

% --- Output exactly as required ---
output ["decision = ", show(decision)];