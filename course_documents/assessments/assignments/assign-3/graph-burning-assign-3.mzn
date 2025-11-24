%GBP using the CSP1 model
int: n;                     % number of vertices
int: m;                     % number of edges

array[1..m] of int: from;    % edge start vertices
array[1..m] of int: to;      % edge end vertices

set of int: V = 1..n;    % domain constraint (30)
set of int: TURNS = 0..n; % burning turns (0 = initial state, 1..n = burning turns)

% I am allowing two options for decision variables here to allow you to write the model as you want to. 
% Comment out the one that you do not want to use
% As in assignment 1, you may include other decision variables if they help you, but do not rename these
% you may change their range bounds if that is useful to you (that is, restrict their values downwards from all integers

% Decision variable option 1: set to t at position i if vertex i is selected for active burning on turn t, -1 otherwise
%array[0..n] of var int: decision_1;
% Decision variable option 2: a 2D indicator that shows at what time you burn which vertex: 
% decision_2[i, j] should be set to 1 if you choose to burn vertex i at time j, 0 otherwise
array[1..n,0..n] of var 0..1: decision_2;

%Helper variable : burned[i,j] = 1 iff vertex i is burned by the end of turn j
array[1..n,0..n] of var 0..1: burned;

%Helper function : Check if vertices are neighbours
function var bool: adjacent(int: u, int: v) =
    exists(e in 1..m)( (from[e] == u /\ to[e] == v) \/ (from[e] == v /\ to[e] == u) );

%Domain Constraints
%No spreading at turn 0 (33)
constraint forall (i in 1..n) (burned[i,0] = 0);
%Nothing set on fire at turn 0 (34)
constraint forall (i in V) (decision_2[i,0] = 0);

%Constraint 24: If a vertex is burned at turn j-1, it remains burned at turn j
constraint forall(i in V, j in 1..n)(burned[i, j-1] = 1 -> burned[i, j] = 1);

%Constraint 25: If a vertex is actively burned, it becomes burned (i.e. vertex is in the burning sequence)
constraint forall(i in V, j in 1..n)(decision_2[i, j] = 1 -> burned[i, j] = 1);

%Constraint 26: If a vertex is adjacent to a burned vertex at turn j-1, it becomes burned at turn j (spreading)
constraint forall(i in V, j in 1..n)(
    exists(k in 1..n)(adjacent(i,k) /\ burned[k, j-1] = 1) -> burned[i, j] = 1);

%Constraint 27: Vertex can only be burned if it was already burning, was actively chosen or via spread from a burning neighbour
%(preventing spontaneous combustion)
constraint forall(i in V, j in 1..n)(
    burned[i, j] = 1 -> (
        burned[i, j-1] = 1
        \/ decision_2[i, j] = 1  
        \/ exists(k in 1..n)(adjacent(i,k) /\ burned[k, j-1] = 1)));

%Constraint 28: Ensure EXACTLY one vertex is actively burned at each turn
%Unless it is the last turn and all are already burned
constraint forall (j in 1..n) (
    sum(i in V)(decision_2[i, j]) = 1
    \/ forall(i in V)(burned[i, j-1] = 1));

%Constraint 29 : All vertices must be burned by turn n 
constraint forall(i in V)(burned[i, n] = 1);

%Other Constraints
%A vertex can only be actively burned once
constraint forall(i in V)(sum(j in 1..n)(decision_2[i, j]) <= 1);

%A vertex can only be actively burned if not already burning
constraint forall(i in V, j in 1..n)(
    decision_2[i, j] = 1 -> burned[i, j-1] = 0);

%Calculating the burning number 
var 1..n: burning_number = max([j | j in 1..n where sum(i in 1..n)(decision_2[i, j]) > 0]);
%Determining the burning sequence 
array[1..n] of var 0..n: burning_sequence;  
constraint forall(i in V)(burning_sequence[i] = sum(j in 1..n)(j * decision_2[i, j]));

%Objective: Minimize the burning number
solve minimize burning_number;

output [
  "burning_number = ", show(burning_number), "\n",
  "burning_sequence = ", show(burning_sequence), "\n"
];