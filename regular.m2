loadPackage "Depth"; 
loadPackage "MCMApproximations";

-- n>=4, r=2, I = {0, 2} 

Q = QQ[x1, x2, x3, x4, y1, y2, y3, y4]; 

X = genericMatrix(Q, x1, 2, 2); 
Y = genericMatrix(Q, y1, 2, 2); 

J = ideal(X*Y) + ideal(Y*X); 
R = Q/J;

test = {x1 - y1, x2 - y3, x3 - y2, x4 + y4};
time isRegularSequence(test, R)

I = ideal(test); 

R0 = R/I; 
-- print the dimension (as a vector space) of the resulting ring 
degree I


----------------------------------------------------------------
-- n=6, r=3, I = {0, 3}
----------------------------------------------------------------

Q = QQ[x1, x2, x3, x4, x5, x6, x7, x8, x9, y1, y2, y3, y4, y5, y6, y7, y8, y9]; 

X = genericMatrix(Q, x1, 3, 3); 
Y = genericMatrix(Q, y1, 3, 3); 

J = ideal(X*Y) + ideal(Y*X)
R = Q/J;


test = {x1 - y1, x2 - y4, x4 - y2, x5 + y5};

helper = (sequence) -> (
    if #sequence == 9 then (
        return {sequence}; 
    ); 
    output = {}; 
    for a in {x1, x2, x3, x4, x5, x6, x7, x8, x9} do (
        for b in {y1, y2, y3, y4, y5, y6, y7, y8, y9} do (
            newseq = append(sequence, a - b); 
            if isRegularSequence(newseq, R) then (
                print "Success! "; 
                print newseq; 
                print "\n"; 
                output = join(output, helper(newseq)); 
            ); 
            newseq = append(sequence, a+b); 
            if isRegularSequence(newseq, R) then (
                print "Success! "; 
                print newseq; 
                print "\n"; 
                output = join(output, helper(newseq)); 
            ); 
        );
    );
    return output; 
);

helper(test)

time isRegularSequence(test, R)

I = ideal(test); 

R0 = R/I; 
dim R0
-- print the dimension (as a vector space) of the resulting ring 
degree I
