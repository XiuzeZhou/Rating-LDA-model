#include "mex.h"

//Syntax
//[Pui1,Pui2,Pui3,Pui4,Pui5]=GetPuir(put,pti,ptr)
//updated for 64-bit compilers

void GetPuir(int U, int I, int R, int T, double *put, double *pti, double *ptr, double *pui1, double *pui2, double *pui3, double *pui4, double *pui5)
{
	int i,j,k,r;
	double p;
	
	r=0;
	for(i=0;i<U;i++)
	{
		for(j=0;j<I;j++)
		{
			p=0;
			for (k=0;k<T;k++)
			{
				p=p+put[k*U+i]*pti[j*T+k]*ptr[r*T+k];
			}
			pui1[j*U+i]=p;
		}
	}

	r=1;
	for(i=0;i<U;i++)
	{
		for(j=0;j<I;j++)
		{
			p=0;
			for (k=0;k<T;k++)
			{
				p=p+put[k*U+i]*pti[j*T+k]*ptr[r*T+k];
			}
			pui2[j*U+i]=p;
		}
	}

	r=2;
	for(i=0;i<U;i++)
	{
		for(j=0;j<I;j++)
		{
			p=0;
			for (k=0;k<T;k++)
			{
				p=p+put[k*U+i]*pti[j*T+k]*ptr[r*T+k];
			}
			pui3[j*U+i]=p;
		}
	}

	r=3;
	for(i=0;i<U;i++)
	{
		for(j=0;j<I;j++)
		{
			p=0;
			for (k=0;k<T;k++)
			{
				p=p+put[k*U+i]*pti[j*T+k]*ptr[r*T+k];
			}
			pui4[j*U+i]=p;
		}
	}

	r=4;
	for(i=0;i<U;i++)
	{
		for(j=0;j<I;j++)
		{
			p=0;
			for (k=0;k<T;k++)
			{
				p=p+put[k*U+i]*pti[j*T+k]*ptr[r*T+k];
			}
			pui5[j*U+i]=p;
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *put,*pti, *ptr, *pui1, *pui2, *pui3,*pui4, *pui5;
	int U,I,R,T,i;
	double *p1, *p2, *p3,*p4 ,*p5;

	put=mxGetPr(prhs[0]);
	U=mxGetM(prhs[0]);
	T=mxGetN(prhs[0]);

	pti=mxGetPr(prhs[1]);
	I=mxGetN(prhs[1]);

	ptr=mxGetPr(prhs[2]);
	R=mxGetN(prhs[2]);

	pui1=(double*)mxCalloc(U*I,sizeof(double));
	pui2=(double*)mxCalloc(U*I,sizeof(double));
	pui3=(double*)mxCalloc(U*I,sizeof(double));
	pui4=(double*)mxCalloc(U*I,sizeof(double));
	pui5=(double*)mxCalloc(U*I,sizeof(double));

	GetPuir(U,I,R,T,put,pti,ptr,pui1,pui2,pui3,pui4,pui5);

	plhs[0]=mxCreateDoubleMatrix(U,I,mxREAL);
	p1=mxGetPr(plhs[0]);
	for(i=0;i<U*I;i++)
		p1[i]=(double)pui1[i];

	plhs[1]=mxCreateDoubleMatrix(U,I,mxREAL);
	p2=mxGetPr(plhs[1]);
	for(i=0;i<U*I;i++)
		p2[i]=(double)pui2[i];

	plhs[2]=mxCreateDoubleMatrix(U,I,mxREAL);
	p3=mxGetPr(plhs[2]);
	for(i=0;i<U*I;i++)
		p3[i]=(double)pui3[i];

	plhs[3]=mxCreateDoubleMatrix(U,I,mxREAL);
	p4=mxGetPr(plhs[3]);
	for(i=0;i<U*I;i++)
		p4[i]=(double)pui4[i];

	plhs[4]=mxCreateDoubleMatrix(U,I,mxREAL);
	p5=mxGetPr(plhs[4]);
	for(i=0;i<U*I;i++)
		p5[i]=(double)pui5[i];
}