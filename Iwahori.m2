loadPackage "TorAlgebra";

-- Iwahori case 
n = 6;
r = 2;
I = toList (0..n-1);

Q = QQ[n*r];
varlist = vars Q;

-- collect matrices defining equations 
M = {};

for i from 0 to n-1 do (
    -- bottom row is variables
    botRow = flatten entries submatrix(varlist, {0}, {i*r..((i+1)*r-1)});
    -- top parts is identity shifted up by 1 
    topRows = flatten entries submatrix(id_(QQ^r), {1..r-1},{0..r-1});
    temp = matrix{topRows, botRow} ;
    print temp; 
    M = append(M, temp); 
);

-- define the ideal 
I = ideal(0); 
for i from 0 to n-1 do (
    prod = 1; 
    for j from 0 to n-1 do (
        prod = prod* M#((i+j)%n); 
    );
    I = I + ideal(prod); 
);

print I; 

print time isGorenstein(Q/I)


