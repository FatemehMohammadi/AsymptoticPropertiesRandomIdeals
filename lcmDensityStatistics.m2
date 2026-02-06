
loadPackage"RandomMonomialIdeals"
loadPackage "Posets"

lcmDensity=method()
-- method to compute lcmDensity of a monomial ideal given the generating set
lcmDensity List := genSet -> (
    m := monomialIdeal genSet; 
    r := numcols  mingens m; 
    lcmI := lcmLattice m;
    return(length(lcmI_*)/2^r)
)


-- example for generating a sample from the Erdos-Renyi random monomial ideal model

n=3; D=4; p={0.0,0.0,0.0,0.5}; N=1000; 
-- n = number of variables 
-- D = max degree of monomial 
-- p = {} probability vector 
-- N = desired sample size 
myModel = ER(n,D,p)

mySample = sample(myModel, N); 
stats =  statistics(mySample, lcmDensity);
print("Sample mean of LCM density is:")
toRR(stats#Mean)
print("Sample standard deviation of LCM density is:")
toRR(stats#StdDev)



-* 
Looking at the RandomMonomialIdeals.m2 package,
there needs to be an update to the basic RMI functions to allow for SQUAREFREE ideals. 
Below, we will copy the original method, and edit it to allow for squarfree option. 
This shold eventually be put to the package itself. 
*- 

squareFree = method(TypicalValue => List)
--returns the list of all squarefree monomials of degree d
squareFree(ZZ,Ring) := List => (d,S) ->
 flatten entries compress(basis(d,S^1) % ideal apply (flatten entries vars S, p -> p^2))

randomMonomialSet = method(TypicalValue => List, Options => {
	Coefficients => QQ,
	VariableName => "x",
	Strategy => "ER"
	})
-- borrowing method from RandomIdeals.m2 but output = list not ideal: 
randomMonomialSet (PolynomialRing,ZZ,List) := List => o -> (R,D,pOrM) -> (
    if #pOrM != D then error "pOrM expected to be a list of length D";
    if not all(pOrM, q->instance(q, ZZ)) and not all(pOrM, q->instance(q,RR))
        then error "pOrM must be a list of all integers or all real numbers";
    B := {};
    if all(pOrM,q->instance(q,ZZ)) then (
        if o.Strategy === "Minimal" then (
            currentRingM := R;
            apply(D, d->(
                chosen := take(random(flatten entries basis(d+1, currentRingM)), pOrM_d);
                B = flatten append(B, chosen/(i->sub(i, R)));
                currentRingM = currentRingM/promote(ideal(chosen), currentRingM)
		)
	    )
	) else B = flatten apply(toList(1..D), d->take(random(flatten entries basis(d,R)), pOrM_(d-1)));
    ) else if all(pOrM,q->instance(q,RR)) then (
        if any(pOrM,q-> q<0.0 or 1.0<q) then error "pOrM expected to be a list of real numbers between 0.0 and 1.0";
        if o.Strategy === "Minimal" then (
            currentRing := R;
            apply(D, d->(
                chosen := select(flatten entries basis(d+1, currentRing), m->random(0.0,1.0)<=pOrM_d);
                B = flatten append(B, chosen/(i->sub(i, R)));
                currentRing = currentRing/promote(ideal(chosen), currentRing)
		)
	    )
	) else if o.Strategy === "Squarefree" then (
	-- squareFree(ZZ,Ring) -- list of all square-free monomials of given degree
	    B = flatten apply(toList(1..D),d-> select(squareFree(d,R),m-> random(0.0,1.0)<=pOrM_(d-1)));
	) else B = flatten apply(toList(1..D),d-> select(flatten entries basis(d,R),m-> random(0.0,1.0)<=pOrM_(d-1)));
    );
    B = apply(B,m->sub(m,R));
    if B==={} then {0_R} else B
)

--- to see how this works, here is an example: 
R = QQ[x_1..x_10];
d=3;
randomMonomialSet(R,d,{0.0,0.0,0.1},Strategy=>"Squarefree")


-* 
below are some simulations and repeated experiments that were used for the paper: 
*- 

studyLCMdensity = method()
studyLCMdensity List := params -> (
    n := params_0; 
    D := params_1;
    p := params_2; 
    N := params_3; 

    myModel = ER(n,D,p);
    mySample = sample(myModel, N); 
    stats =  statistics(mySample, lcmDensity);

    return(toRR(stats#Mean), toRR(stats#StdDev) )
)

-- 1000 ideals each in 4 vars prob param varying 0.1..0.9 by 0.2; D=3"
apply(4,i-> studyLCMdensity({4,3,0.1+0.2*i,1000}))



-*
For those who wish to design other ways to get statistics on monomial ideals, below is some code for generating random monomial ideals
and studying their various properties. 
*- 

loadPackage"RandomMonomialIdeals"

--fix parameter values for the graded model: in the line below, this is degree 2, 5 vars, all quadratic, prob of picking any quadratic monomial = 0.5 (and it's 0 for linear generators): 
n=5; D=2; p = {0,.5}; 

-- define/import the ER model 
myModel = ER(n,D,p)
-- get 10 ideals under the model w/ the above parameters: 
mySample = sample(myModel, 10); 
myMonomialGeneratingSets = getData mySample
statistics(mySample, dim@@monomialIdeal)


-* 
some functions for computing LCM lattices
*- 

loadPackage "Posets"

s = myIdeals_0
lcmI = lcmLattice monomialIdeal s;
-- ok and HOW do I compute the size of lcmI (this is a Poset object) as a set? 
-- is it # lcmI ? or what. 
length(lcmI_*)

displayPoset lcmI -- wow this works! creates a PDF of the LCM lattice! 

r = numcols  mingens monomialIdeal s
2^r
-- Note to self: we can use M2 to get nice tikz pictures of the posets https://arxiv.org/html/2510.22843v1 

dens = length(lcmI_*)/2^r --- this is the LCM lattice density computation! 


-* 
putting it together 
*- 

lcmDensity=method()
lcmDensity List := genSet -> (
    m := monomialIdeal genSet; 
    r := numcols  mingens m; 
    lcmI := lcmLattice m;
    return(length(lcmI_*)/2^r)
)

--what parameters are we using:
myModel#Parameters
-- let's get a sample: 
mySample = sample(myModel, 3); 
-- give me mean, stdev, and tally of all the LCM lattice densities for this sample of ideals: 
statistics(mySample, lcmDensity)

-- hide the tally output: 
mySample = sample(myModel, 100); 
stats =  statistics(mySample, lcmDensity);
toRR(stats#Mean)
toRR(stats#StdDev)

-* 
studying minimal free resolutions of these random objects:
*- 

myModel = ER(n,D,p)
mySample = sample(myModel, 100); 


L = (getData mySample)/ideal; 
bettiStats L 


---Other densities
bettiDensity=method()
bettiDensity List := genSet -> (
    m := monomialIdeal genSet; 
    r := numcols  mingens m; 
    bettis := betti res m;
    return((sum values bettis)/2^r)
)

lengthBettiDensity=method()
lengthBettiDensity List := genSet -> (
    m := monomialIdeal genSet; 
    r := numcols  mingens m; 
    bettis := betti res m;
    return((length values bettis)/2^r)
)
