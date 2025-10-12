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