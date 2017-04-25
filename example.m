%% Example of running basic RLDA
%
% This example shows how to run the Gibbs sampler on a small dataset to
% extract a set of interests and shows the most likely items per interest, 
% and ratings per interest. It also writes the results to a file.

%%
load 'train_datas'; % load the nips item-user rating dataset

%%
% Set the number of interests
T=50; 

%%
% Set the hyperparameters
ALPHA=50/T;
BETA=0.01;
GAMA=0.01;

%%
% The number of iterations
N = 1000; 

%%
% The random seed
SEED = 3;

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% This function might need a few minutes to finish
tic
[ IT, UT, RT, Z ] = GibbsSamplerRLDA( IS , US , RS , T , N , ALPHA , BETA , GAMA , SEED , OUTPUT );
toc

%%
% Just in case, save the resulting information from this sample 
temp=['.\output\temp',num2str(T)];
save( temp, 'IT', 'UT' ,'RT', 'Z', 'N' , 'ALPHA' ,'BETA' ,'GAMA' , 'SEED');

pu_t = GetProb( UT , ALPHA ); % user-interest probability
pt_i = GetProb( IT', BETA ); % interest-item probability
pt_r = GetProb( RT', GAMA ); % interest-rating probability
[ pui1, pui2, pui3, pui4, pui5 ] = GetPuir( pu_t, pt_i, pt_r ); % user-interest-rating probability

fprintf( 'OK\n' ); 
