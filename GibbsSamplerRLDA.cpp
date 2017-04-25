#include "mex.h"
#include "cokus.cpp"


// Updated for 64-bit compilers

void GibbsSamplerRLDA( double ALPHA, double BETA, double GAMA, int W, int T, int D, int R, int NN, int OUTPUT, int n, int *z, int *d, int *w, int *r, int *wp, int *dp, int *zr, int *ztot, int *order, double *probs, int startcond )
{
    int wi, di, i, ii, j, interest, rp, temp, iter, wioffset, dioffset, ri, rioffset;
    double totprob, WBETA, RGAMA, rnd, max;

    if (startcond == 1)
    {
        /* start from previously saved state */
        for (i = 0; i < n; i++)
        {
            wi = w[ i ];
            di = d[ i ];
            ri = r[i];
            interest = z[ i ];
            wp[ wi * T + interest ]++; // increment wp count matrix
            dp[ di * T + interest ]++; // increment dp count matrix
            zr[ ri * T + interest ]++;// increment zr count matrix
            ztot[ interest ]++; // increment ztot matrix
        }
    }

    if (startcond == 0)
    {
        /* random initialization */
        if (OUTPUT == 2) mexPrintf( "Starting Random initialization\n" );
        for (i = 0; i < n; i++)
        {
            wi = w[ i ];
            di = d[ i ];
            ri = r[i];
            // pick a random interest 0..T-1
            interest = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) );
            z[ i ] = interest; // assign this item token to this interest
            wp[ wi * T + interest ]++; // increment wp count matrix
            dp[ di * T + interest ]++; // increment dp count matrix
            zr[ ri * T + interest ]++; // increment zr count matrix
            ztot[ interest ]++; // increment ztot matrix
        }
    }

    if (OUTPUT == 2) mexPrintf( "Determining random order update sequence\n" );

    for (i = 0; i < n; i++) order[i] = i; // fill with increasing series
    for (i = 0; i < (n - 1); i++)
    {
        // pick a random integer between i and nw
        rp = i + (int) ((double) (n - i) * (double) randomMT() / (double) (4294967296.0 + 1.0));

        // switch contents on position i and position rp
        temp = order[rp];
        order[rp] = order[i];
        order[i] = temp;
    }

    //for (i=0; i<n; i++) mexPrintf( "i=%3d order[i]=%3d\n" , i , order[ i ] );
    WBETA = (double) (W * BETA);
    RGAMA = (double) (R * GAMA);
    for (iter = 0; iter < NN; iter++)
    {
        if (OUTPUT >= 1)
        {
            if ((iter % 10) == 0) mexPrintf( "\tIteration %d of %d\n" , iter , NN );
            if ((iter % 10) == 0) mexEvalString("drawnow;");
        }
        for (ii = 0; ii < n; ii++)
        {
            i = order[ ii ]; // current item token to assess

            wi  = w[i]; // current item index
            di  = d[i]; // current user index
            ri  = r[i]; // current rating index
            interest = z[i]; // current interest assignment to item token
            ztot[interest]--;  // substract this from counts

            wioffset = wi * T;
            dioffset = di * T;
            rioffset = ri * T;

            wp[wioffset + interest]--;
            dp[dioffset + interest]--;
            zr[rioffset + interest]--;

            //mexPrintf( "(1) Working on ii=%d i=%d wi=%d di=%d interest=%d wp=%d dp=%d\n" , ii , i , wi , di , interest , wp[wi+interest*W] , dp[wi+interest*D] );

            totprob = (double) 0;
            for (j = 0; j < T; j++)
            {
                probs[j] = ((double) wp[ wioffset + j ] + (double) BETA) / ( (double) ztot[j] + (double) WBETA) * ( (double) dp[ dioffset + j ] + (double) ALPHA) * ((double)zr[rioffset + j] + (double)GAMA) / ((double)ztot[j] + (double)RGAMA);
                totprob += probs[j];
            }

            // sample a interest from the distribution
            rnd = (double) totprob * (double) randomMT() / (double) 4294967296.0;
            max = probs[0];
            interest = 0;
            while (rnd > max)
            {
                interest++;
                max += probs[interest];
            }

            z[i] = interest; // assign current item token i to interest j
            wp[wioffset + interest ]++; // and update counts
            dp[dioffset + interest ]++;
            zr[rioffset + interest ]++;
            ztot[interest]++;
            //mexPrintf( "(2) Working on ii=%d i=%d wi=%d di=%d interest=%d wp=%d dp=%d\n" , ii , i , wi , di , interest , wp[wi+interest*W] , dp[wi+interest*D] );
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    double *srwp, *srdp, *probs, *Z, *WS, *DS, *ZIN, *RS, *Ra;
    double ALPHA, BETA, GAMA;
    mwIndex *irwp, *jcwp, *irdp, *jcdp;
    int *z, *d, *w, *r, *order, *wp, *dp, *ztot, *zr;
    int W, T, D, R, NN, SEED, OUTPUT, nzmax, nzmaxwp, nzmaxdp, ntokens;
    int i, j, c, n, nt, wi, di, startcond, ri;

    /* Check for proper number of arguments. */
    if (nrhs < 10)
    {
        mexErrMsgTxt("At least 10 input arguments required");
    }
    else if (nlhs < 4)
    {
        mexErrMsgTxt("4 output arguments required");
    }

    startcond = 0;
    if (nrhs == 11) startcond = 1;

    /* process the input arguments */
    if (mxIsDouble( prhs[ 0 ] ) != 1) mexErrMsgTxt("WS input vector must be a double precision matrix");
    if (mxIsDouble( prhs[ 1 ] ) != 1) mexErrMsgTxt("DS input vector must be a double precision matrix");
    if (mxIsDouble( prhs[ 2 ] ) != 1) mexErrMsgTxt("RS input vector must be a double precision matrix");

    // pointer to item indices
    WS = mxGetPr( prhs[ 0 ] );

    // pointer to user indices
    DS = mxGetPr( prhs[ 1 ] );

    // pointer to rating indices
    RS = mxGetPr( prhs[ 2 ] );

    // get the number of tokens
    ntokens = mxGetM( prhs[ 0 ] ) * mxGetN( prhs[ 0 ] );


    if (ntokens == 0) mexErrMsgTxt("WS vector is empty");
    if (ntokens != ( mxGetM( prhs[ 1 ] ) * mxGetN( prhs[ 1 ] ))) mexErrMsgTxt("WS and DS vectors should have same number of entries");

    T    = (int) mxGetScalar(prhs[3]);
    if (T <= 0) mexErrMsgTxt("Number of interests must be greater than zero");

    NN    = (int) mxGetScalar(prhs[4]);
    if (NN < 0) mexErrMsgTxt("Number of iterations must be positive");

    ALPHA = (double) mxGetScalar(prhs[5]);
    if (ALPHA <= 0) mexErrMsgTxt("ALPHA must be greater than zero");

    BETA = (double) mxGetScalar(prhs[6]);
    if (BETA <= 0) mexErrMsgTxt("BETA must be greater than zero");

    GAMA = (double) mxGetScalar(prhs[7]);
    if (GAMA <= 0) mexErrMsgTxt("GAMA must be greater than zero");

    SEED = (int) mxGetScalar(prhs[8]);

    OUTPUT = (int) mxGetScalar(prhs[9]);

    if (startcond == 1)
    {
        ZIN = mxGetPr( prhs[ 10 ] );
        if (ntokens != ( mxGetM( prhs[ 10 ] ) * mxGetN( prhs[ 10 ] ))) mexErrMsgTxt("WS and ZIN vectors should have same number of entries");
    }

    // seeding
    seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers



    /* allocate memory */
    z  = (int *) mxCalloc( ntokens , sizeof( int ));

    if (startcond == 1)
    {
        for (i = 0; i < ntokens; i++) z[ i ] = (int) ZIN[ i ] - 1;
    }

    d  = (int *) mxCalloc( ntokens , sizeof( int ));
    w  = (int *) mxCalloc( ntokens , sizeof( int ));
    r = (int *)mxCalloc(ntokens , sizeof(int));
    order  = (int *) mxCalloc( ntokens , sizeof( int ));
    ztot  = (int *) mxCalloc( T , sizeof( int ));
    probs  = (double *) mxCalloc( T , sizeof( double ));

    // copy over the item and user indices into internal format
    for (i = 0; i < ntokens; i++)
    {
        w[ i ] = (int) WS[ i ] - 1;
        d[ i ] = (int) DS[ i ] - 1;
        r[i] = (int)RS[i] - 1;
    }

    n = ntokens;

    W = 0;
    D = 0;
    R = 0;
    for (i = 0; i < n; i++)
    {
        if (w[ i ] > W) W = w[ i ];
        if (d[ i ] > D) D = d[ i ];
        if (r[ i ] > R) R = r[ i ];
    }
    W = W + 1;
    D = D + 1;
    R = R + 1;

    wp  = (int *) mxCalloc( T * W , sizeof( int ));
    dp  = (int *) mxCalloc( T * D , sizeof( int ));
    zr  = (int *) mxCalloc( T * R , sizeof( int ));

    //mexPrintf( "N=%d  T=%d W=%d D=%d\n" , ntokens , T , W , D );

    if (OUTPUT == 2)
    {
        mexPrintf( "Running LDA Gibbs Sampler Version 1.0\n" );
        if (startcond == 1) mexPrintf( "Starting from previous state ZIN\n" );
        mexPrintf( "Arguments:\n" );
        mexPrintf( "\tNumber of items      W = %d\n"    , W );
        mexPrintf( "\tNumber of USERs       D = %d\n"    , D );
        mexPrintf( "\tNumber of interests     T = %d\n"    , T );
        mexPrintf( "\tNumber of iterations N = %d\n"    , NN );
        mexPrintf( "\tHyperparameter   ALPHA = %4.4f\n" , ALPHA );
        mexPrintf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
        mexPrintf( "\tHyperparameter    GAMA = %4.4f\n" , GAMA );
        mexPrintf( "\tSeed number            = %d\n"    , SEED );
        mexPrintf( "\tNumber of tokens       = %d\n"    , ntokens );
        mexPrintf( "Internal Memory Allocation\n" );
        mexPrintf( "\tw,d,z,order indices combined = %d bytes\n" , 4 * sizeof( int) * ntokens );
        mexPrintf( "\twp (full) matrix = %d bytes\n" , sizeof( int ) * W * T  );
        mexPrintf( "\tdp (full) matrix = %d bytes\n" , sizeof( int ) * D * T  );
        mexPrintf( "\tzr (full) matrix = %d bytes\n" , sizeof( int ) * T * R  );
        //mexPrintf( "Checking: sizeof(int)=%d sizeof(long)=%d sizeof(double)=%d\n" , sizeof(int) , sizeof(long) , sizeof(double));
    }

    /* run the model */
    GibbsSamplerRLDA( ALPHA, BETA, GAMA, W, T, D, R, NN, OUTPUT, n, z, d, w, r, wp, dp, zr, ztot, order, probs, startcond );

    /* convert the full wp matrix into a sparse matrix */
    nzmaxwp = 0;
    for (i = 0; i < W; i++)
    {
        for (j = 0; j < T; j++)
            nzmaxwp += (int) ( *( wp + j + i * T )) > 0;
    }
    /*if (OUTPUT==2) {
        mexPrintf( "Constructing sparse output matrix wp\n" );
        mexPrintf( "Number of nonzero entries for WP = %d\n" , nzmaxwp );
    }*/

    // MAKE THE WP SPARSE MATRIX
    plhs[0] = mxCreateSparse( W, T, nzmaxwp, mxREAL);
    srwp  = mxGetPr(plhs[0]);
    irwp = mxGetIr(plhs[0]);
    jcwp = mxGetJc(plhs[0]);
    n = 0;
    for (j = 0; j < T; j++)
    {
        *( jcwp + j ) = n;
        for (i = 0; i < W; i++)
        {
            c = (int) * ( wp + i * T + j );
            if (c > 0)
            {
                *( srwp + n ) = c;
                *( irwp + n ) = i;
                n++;
            }
        }
    }
    *( jcwp + T ) = n;

    // MAKE THE DP SPARSE MATRIX
    nzmaxdp = 0;
    for (i = 0; i < D; i++)
    {
        for (j = 0; j < T; j++)
            nzmaxdp += (int) ( *( dp + j + i * T )) > 0;
    }
    /*if (OUTPUT==2) {
        mexPrintf( "Constructing sparse output matrix dp\n" );
        mexPrintf( "Number of nonzero entries for DP = %d\n" , nzmaxdp );
    } */
    plhs[1] = mxCreateSparse( D, T, nzmaxdp, mxREAL);
    srdp  = mxGetPr(plhs[1]);
    irdp = mxGetIr(plhs[1]);
    jcdp = mxGetJc(plhs[1]);
    n = 0;
    for (j = 0; j < T; j++)
    {
        *( jcdp + j ) = n;
        for (i = 0; i < D; i++)
        {
            c = (int) * ( dp + i * T + j );
            if (c > 0)
            {
                *( srdp + n ) = c;
                *( irdp + n ) = i;
                n++;
            }
        }
    }
    *( jcdp + T ) = n;

    plhs[2] = mxCreateDoubleMatrix(R, T, mxREAL);
    Ra = mxGetPr(plhs[2]);
    for (i = 0; i < R; i++)
    {
        for (j = 0; j < T; j++)
            Ra[i * T + j] = (double)zr[i * T + j];
    }

    plhs[ 3 ] = mxCreateDoubleMatrix( 1, ntokens , mxREAL );
    Z = mxGetPr( plhs[ 3 ] );
    for (i = 0; i < ntokens; i++) Z[ i ] = (double) z[ i ] + 1;
}
