function [ prob ] = GetProb( MN, HyperPara )

%MN: UT, IT or RT
%HyperPara:the hyperparameters on the Dirichlet priors

[ M, N ] = size( MN );

nume = MN + HyperPara;
sumMN = sum( MN, 2 ) + N * HyperPara;
deno = repmat( sumMN, 1, N );
prob = nume./deno;
    
end