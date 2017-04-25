%% Function GibbsSamplerRLDA
% Runs the Gibbs sampler for the RLDA model
%
% [IT, UT, RT, Z] = GibbsSamplerRLDA(IS, US, RS, T, N, ALPHA, BETA, GAMA, SEED, OUTPUT)
% runs the Gibbs sampler for the RLDA model on an item-user rating data
% provided by the vectors |IS|, |US| and |RS|. |IS(k)|, |US(k)| and |RS(k)|
% contain the item, user and rating indices for the kth token. The maximum 
% of |IS| is |N|, the the number of items. The maximum of |US| is |M|, the 
% number of users. |T| is the number of interests. The first output is the
% sparse matrix |IT|, of size |N| x |T|, where |IT(i,j)| contains the number
% of times item |i| has been assigned to interest |j|. The second output is
% |UT|, a sparse |M| x |T| matrix, where |UT(i,j)| contains the number of 
% times an item in user |u| has been assigned to interest |j|. The third 
% output is |RT|, a sparse |R| x |T| matrix, where |RT(i,j)| contains the 
% number of times an item |i| has been assigned to interest |j|. The fouth 
% output |Z| contains the intereest assignments; |Z(k)| contains the interest 
% assignment for token k.    
%
% NOTES
%
% |IS|, |US| and |RS| should be in double precision
% |N| determines the number of iterations to run the Gibbs sampler.
% |ALPHA|, |BETA| and |GAMA| are the hyperparameters on the Dirichlet priors
% for the interest distributions (|theta|), the interest-item distributions
% (|phi|) and the interest-item distributions (|lambda|) respectively.
% |SEED| sets the seed for the random number generator
%
% |OUTPUT| determines the screen output by the sampler
%   0 = no output provided
%   1 = show the iteration number only
%   2 = show all output
%
% The time to complete the procedure scales linearly with the number of
% interests and the number of iterations. The memory requirements scale
% linearly with the number of interests and users.
%
% A good setting for the number of iterations will depend on the number of
% interests and the complexity of problem. For most problems, 500 to 2000
% iterations will suffice.
%
% Appropriate values for |ALPHA|, |BETA| and |GAMA| depend on the number of
% interests and the number of items, and the number of ratings. For most
% applications, good results can be obtained by setting |ALPHA = 50 / T|,
% |BETA = 200 / W| and |GAMA = 0,01|
%
% The sampler uses its own random number generator and setting the seed for
% this function will not influence the random number seed for Matlab
% functions
%
%%
% REFERENCES
%%
% * Xiuze Zhou and Shunxiang Wu (2016). Rating LDA model for collaborative filtering.
% Knowledge-Based Systems, 110:135-143.
% * D. Blei, A. Ng, and M. Jordan (2003). Latent Dirichlet allocation.
% Journal of Machine Learning Research, 3:993-1022.


