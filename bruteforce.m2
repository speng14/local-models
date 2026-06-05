loadPackage "TorAlgebra";

-- general case, brute force approach by enumerating all matrix equations. 
n = 5;
r = 3;
I = toList(0,1, 2);

Q = QQ[(n-r)*r*(#I)];
varlist = flatten entries vars Q;

-- enumerate all the F_i's 
F = {}; 
Ir = id_(QQ^r); 

for j from 0 to (#I-1) do (
    i = I#j;
    temp = genericMatrix(Q, varlist#((n-r)*r*j), n-r, r); 
    Fi = Ir || temp; 
    if i != 0 then (
        topRows = submatrix(Fi, {n-i..n-1}, {0..r-1}); 
        botRows = submatrix(Fi, {0..n-i-1}, {0..r-1}); 
        Fi = topRows || botRows; 
    );

    F = append(F, Fi);
);

print F; 

-- enumerate all the M_i's 
M = {}; 
for j from 0 to (#I-1) do (
    i = I#j; 
    test = toList(i..i+r-1); 
    rows = apply(i..i+r-1, x -> x % n); 
    Mi = submatrix(F#((j-1) % #F), rows, {0..r-1}); 
    M = append(M, Mi); 
);

print M; 

-- multiply matrices to get relations 
J = ideal(0); 

for j from 0 to (#I-1) do (
    i = I#j; 
    rows = apply(i..i+r-1, x -> x % n); 
    rows = set(0..n-1) - set(rows); 
    rows = toList(rows);  
    temp = apply(n, k -> if I#((j-1)% #I) <= k and k <= ((i-1)%n) then 0 else 1);
    D = diagonalMatrix temp;

    Fi = submatrix(F#j, rows, {0..r-1}); 
    Gi = submatrix(D*F#((j-1)% #I), rows, {0..r-1}); 
    J = J + ideal(Fi * M#j - Gi);
); 

print J; 
print time isGorenstein(Q/J); 

--- Issue: when there are linear terms, there's no easy way to eliminate
--- so it will error out 