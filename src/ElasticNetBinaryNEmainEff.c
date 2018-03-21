//version 4.2 note: initialize the least correlated variable to avoid collinear variables; delete the initial effect at the end of 1st iteration.
//version 4.3 note: initialize the largest correlated variable, but set the largest*0.9 to have correlation = 0;delete initial basis at the end of 1st itaration.
// Jul. 2015
//
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <stdlib.h>
//B: binomial; m: main; Ne: normal + exp. _EN: elasticNet
/*
void printMat(double *a, int M, int N) //MxN
{
	int i,j;
	Rprintf("Printing the matrix\n\n");
	for(i=0;i<M;i++) 
	{
		for(j=0;j<N;j++)
		{
			Rprintf("%f\t", a[j*M +i]); //a[i,j]
		}
		Rprintf("\n");
	}
}

void testFunc(	double *BASIS,
				double *Scales,
				double *PHI2,
				double *beta,
				double *Alpha,
				double *SIGMA2,
				double *H2,
				double *S_in,
				double *Q_in,
				double *S_out,
				double *Q_out,
				double *Mu2,				
				double * Targets,
				double *basisCache,
				int *n, 
				int *kdim,
				int *m,
				int*Used,
				int *nuu,
				double *phi,
				double *newA)
{
	Rprintf("\t Initialized basis %d, Alpha: %f, \n", Used[0],Alpha[0]);
	fEBCatFullStatBmNeEN(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2, 
				Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);

	int nu = *nuu;
	int M_full = *kdim;
	int N = *n;
	int K = M_full;
	int M = *m;
	double newAlpha = *newA;
	int UPDATE_REQUIRED;
	Rprintf("\t add basis %d, Alpha: %f, N:%d\tK:%d\tM:%d\tM_full:%d\n", nu,newAlpha,N,K,M,M_full);
	
	UPDATE_REQUIRED		= ActionAddBmNeEN( BASIS, Scales, PHI2, phi, beta, Alpha,
						newAlpha, SIGMA2, Mu2, S_in, Q_in, nu, M_full, N, K, M);
						
	UPDATE_REQUIRED= ActionResBmNeEN(BASIS,Scales,PHI2, Alpha, SIGMA2, newAlpha,
							Mu2, S_in, Q_in, beta, nu, N, M,M_full);						
	printMat(SIGMA2,(M+1),(M+1));
	UPDATE_REQUIRED = ActionDelBmNeEN(BASIS, Scales, PHI2, Alpha, SIGMA2, 
											Mu2, S_in, Q_in, beta, nu, N, M, M_full);				
	Rprintf("\t AFTER DEL basis %d, Alpha: %f, N:%d\tK:%d\tM:%d\tM_full:%d\n", nu,newAlpha,N,K,M,M_full);
	printMat(SIGMA2,M,M);
	
}

*/

void LinearSolverBmNeEN(double * a, double *logout, int N,int M,double *output);
void fEBInitializationBmNeEN(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, double *Scales, int * initial, int n, int *m, int kdim);

void fEBInitializationBmNeEN2(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, double *Scales, int * initial, int n, int *m, int kdim);

void fEBSigmoidBmNeEN(double * y, double * PHI_Mu,int N);
double fEBDataErrorBmNeEN(double dataError,double *y,double *PHI_Mu,double *Targets,int N);
void MatrixInverseBmNeEN(double * a,int N);
void fEBCatPostModeBmNeEN(double * Mu2, double *beta,double *SIGMA2, double * H2, double *PHI2,
								double *Targets, double *Alpha,int N, int M);
void fEBCatFullStatBmNeEN(double * beta, double * SIGMA2, double * H2, double *S_in, double * Q_in, 
				double * S_out, double * Q_out,  double *BASIS,double * Scales, double *PHI2, 
				double * Targets, int * Used, double *Alpha, double * Mu2, double * BasisCache,
				int *n, int *m, int *kdim);

void fEBDeltaMLBmNeEN(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, double * S_out, double * Q_out, double *Alpha,
				double *a_lmabda,double *b_Alpha, int m, int mBar,double *deltaLogMarginal, int * nu,int N);


void fEBBinaryMexBmNeEN(int *Used, double *Mu2, double *SIGMA2, double *H2, double *Alpha, double *PHI2,
				double *BASIS, double * Targets, double *Scales, double *a_lambda, double *b_Alpha,
				int *iteration, 
				int *n, int *kdim, int *m,double * LOGlikelihood,int basisMax,int verbose);

void ElasticNetBinaryNEmainEff(double *BASIS, double * Targets, double *a_Lambda,double *b_Alpha,
				double * logLIKELIHOOD, double * Beta, double *wald,double *intercept, int *n, int *kdim,int*VB,int*bMax);

int ActionAddBmNeEN( double* BASIS, double*scales, double*PHI, double*Phi,
			double *beta, double* Alpha, double newAlpha, double*SIGMA, double*Mu, double*S_in,
			double*Q_in, int nu, int M_full,int N, int K, int M);
int ActionDelBmNeEN(double* BASIS, double*scales, double*PHI, double*Alpha, double*SIGMA,
		double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full);
int ActionResBmNeEN(double* BASIS, double*scales, double*PHI, double*Alpha, double*SIGMA, double newAlpha,
		double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full);

	
void ElasticNetBinaryNElambdaMax(double *BASIS, double * Targets, double *Lambda_max,
				double * Beta, double *wald,double *intercept, int *n, int *kdim) //,int*VB)
{
	int N					= *n;
	int K					= *kdim;
	int M_full				= K;
	//int verbose = VB[0];
	// set a limit for number of basis

	int basisMax			= 1e6/N;
	if (basisMax>M_full)	basisMax = M_full;

	double temp				= 0;
	int i,j,l;
	double *Scales			= (double * ) R_alloc(M_full, sizeof(double));

	for (i					=0;i<K;i++)
	{
		Beta[i]				= i + 1;
		Beta[M_full + i]	= i + 1;
		Beta[M_full*2 + i]	= 0;
		Beta[M_full*3 + i]	= 0;
		temp				= 0;
		for(l=0;l<N;l++)	temp		= temp + BASIS[i*N + l]*BASIS[i*N + l];
		if(temp ==0) temp	= 1;
		Scales[i]			=sqrt(temp);
	}

	//
	double *Mu2, *SIGMA2, *H2, *Alpha, *PHI2;
	int * Used, *m;
	//Rprintf("max basis: %d\n",basisMax);
	
	Used					= (int* ) R_alloc(basisMax, sizeof(int));
	Mu2						= (double * ) R_alloc(basisMax, sizeof(double));
	SIGMA2					= (double * ) R_alloc(basisMax*basisMax, sizeof(double));
	H2						= (double * ) R_alloc(basisMax*basisMax, sizeof(double));
	Alpha					= (double * ) R_alloc(basisMax, sizeof(double));
	PHI2					= (double * ) R_alloc(N*basisMax, sizeof(double));

	m						= (int* ) R_alloc(1, sizeof(int));
	//Rprintf("outer loop starts");
	m[0]					= 2;
	/***********************/

	int M;

    double *beta			= (double *)R_alloc(N,sizeof(double));
	int *Unused				= (int *) R_alloc(M_full,sizeof(int));

	//Rprintf("Iteration number: %d\n",iter);
    
    //[Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic) 
    int *IniLogic;
	IniLogic				= (int*) R_alloc(1,sizeof(int));

        IniLogic[0]			= 0;
        m[0]				= 2;
		M					= m[0];


	//Rprintf("N_used is: %d; N_unused:%d, M: %d,sample size: %d \n",N_used,N_unused,M,N);

	fEBInitializationBmNeEN(Alpha, PHI2, Used, Unused, Mu2, BASIS, Targets, Scales, IniLogic, N, m, K);

    double *basisCache;
	basisCache				= (double *) R_alloc(N*K,sizeof(double));
    for(i=0;i<K;i++)
    {
        for(j=0;j<N;j++)	basisCache[i*N+j] = BASIS[i*N+j]*BASIS[i*N+j];
    }
	
	double *S_in, *Q_in, *S_out, *Q_out;
	S_in					= (double *) R_alloc(M_full,sizeof(double));
	Q_in					= (double *) R_alloc(M_full,sizeof(double));
	S_out					= (double *) R_alloc(M_full,sizeof(double));
	Q_out					= (double *) R_alloc(M_full,sizeof(double));
    //[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,Intercept] ...
    //                   	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE)
	fEBCatFullStatBmNeEN(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2, 
				Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);

	/****************************/
	int index;
	Lambda_max[0]			= 0;
	temp					= 0;
	for(index=0;index<M_full;index++)
	{
		temp				= (Q_out[index]*Q_out[index] - S_out[index])/2;
		if(temp>Lambda_max[0])			Lambda_max[0] = temp;
	}
	// wald score
	M						= m[0];	
	double *tempW			= (double * ) R_alloc(M,sizeof(double));
	wald[0]					= 0;
	index = 0;
	//Rprintf("number of basis: %d\n",M);
	for(i=0;i<M;i++)
    {
        tempW[i]      		= 0;
        for(j=0;j<M;j++)    tempW[i]     = tempW[i] + Mu2[j]*H2[i*M+j];       
        wald[0]				= wald[0]	 +tempW[i]*Mu2[i];
	}
	for(i=1;i<M;i++)
	{
    // blup collection
		//Rprintf("test Used: %d\n",Used[i-1]);
		index				= Used[i-1] - 1;
		Beta[M_full*2 + index]	= Mu2[i]/Scales[index];
		Beta[M_full*3 + index]  = SIGMA2[i*M + i]/(Scales[index]*Scales[index]);
	}
	//
	intercept[0]			= Mu2[0];
	intercept[1]			= SIGMA2[0];
	//Rprintf("fEB computation compelete!\n");
}


void ElasticNetBinaryNEmainEff(double *BASIS, double * Targets, double *a_Lambda,double *b_Alpha,
				double * logLIKELIHOOD, double * Beta, double *wald,double *intercept, int *n, int *kdim,
				int*VB,int*bMax)
{
	int N					= *n;
	int K					= *kdim;
	//int M_full				= (K+1)*K/2;
	int M_full				= K;
	int verbose = VB[0];
		
	const int iter_max		= 100;
	const double err_max	= 1e-8;
	//------------------------------------------------------------	
	int basisMax			= *bMax;
	//------------------------------------------------------------	
	// set a limit for number of basis

	//int basisMax			= 1e6/N;
	//if (basisMax>M_full)	basisMax = M_full;

	double vk				= 1e-30;
	double vk0				= 1e-30;
	double temp				= 0;
	int i;//,j,l;
	//double scal 			= 1.0;
	//		xixi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
	int inci =1;
	int incj =1;
	int inc0 				= 0;
	double *readPtr1;
	double zero_blas 		= 0;
	
	double *Scales			= (double * ) Calloc(M_full, double);

	for (i					=0;i<K;i++)
	{
		Beta[i]				= i + 1;
		Beta[M_full + i]	= i + 1;
		//Beta[M_full*2 + i]	= 0;
		//Beta[M_full*3 + i]	= 0;
		//temp				= 0;
		//for(l=0;l<N;l++)	temp		= temp + BASIS[i*N + l]*BASIS[i*N + l];
		//if(temp ==0) temp	= 1;
		//Scales[i]			=sqrt(temp);
		
		readPtr1 = &BASIS[i*N];
		temp = F77_CALL(ddot)(&N, readPtr1, &inci,readPtr1, &incj);
		if(temp ==0) temp	= 1;
		Scales[i]			=sqrt(temp);
		
	}
	
	readPtr1 				= &Beta[basisMax*2];
	F77_CALL(dcopy)(&basisMax,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
	readPtr1 				= &Beta[basisMax*3];
	F77_CALL(dcopy)(&basisMax,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	
	//PartII kk
	/*kk						= K;					//index starts at 0;	
	for (i					=0;i<(K-1);i++)
	{
		for (j				=(i+1);j<K;j++)
		{
			Beta[kk]			= i + 1;
			Beta[M_full + kk]	= j + 1;
			Beta[M_full*2 + kk] = 0;
			Beta[M_full*3 + kk] = 0;
			temp				= 0;
			for(l=0;l<N;l++) 	temp		= temp + BASIS[i*N + l]*BASIS[i*N + l]*BASIS[j*N + l]*BASIS[j*N + l];
			if (temp == 0)		temp		= 1;
			Scales[kk]			=sqrt(temp);
			kk++;
		}
	}
	*/
	//
	int iter				= 0;
	double err				= 1000;
	double *Mu2, *SIGMA2, *H2, *Alpha, *PHI2;
	int * Used,*iteration, *m;
	//Rprintf("max basis: %d\n",basisMax);
	
	Used					= (int* ) Calloc(basisMax, int);
	Mu2						= (double * ) Calloc(basisMax, double);
	SIGMA2					= (double * ) Calloc(basisMax*basisMax, double);
	H2						= (double * ) Calloc(basisMax*basisMax, double);
	Alpha					= (double * ) Calloc(basisMax, double);
	PHI2					= (double * ) Calloc(N*basisMax, double);
	iteration				= (int* ) Calloc(1, int);
	m						= (int* ) Calloc(1, int);
	if(verbose >1) Rprintf("Empirical Bayesian Elastic Net outer loop starts\n");
	m[0]			= 2;
	while (iter<iter_max && err>err_max)
	{
		iter				= iter + 1;
		
		vk0					= vk;
		iteration[0]		= iter;
		fEBBinaryMexBmNeEN(Used, Mu2, SIGMA2, H2, Alpha,PHI2,	BASIS, Targets,Scales, a_Lambda,b_Alpha,
						iteration, n, kdim, m,logLIKELIHOOD,basisMax,verbose);

		//vk					= 0;
		//for(i=0;i<m[0]-1;i++)	vk = vk + Alpha[i];
		vk = F77_CALL(dasum)(m, Alpha, &inci);
		err					= fabs(vk - vk0)/m[0];
		if(verbose>2) Rprintf("Iteration number: %d, err: %f\n",iter,err);

	}

	// wald score
	int M					= m[0];	
	double *tempW			= (double * ) Calloc(M,double);

	wald[0]					= 0;
	int index = 0;
	if(verbose >1) Rprintf("EBEN Finished, number of basis: %d\n",M);
	for(i=0;i<M;i++)
    {

        //tempW[i]      		= 0;
        //for(j=0;j<M;j++)    tempW[i]     = tempW[i] + Mu2[j]*H2[i*M+j];       
        //wald[0]				= wald[0]	 +tempW[i]*Mu2[i];
		readPtr1 = &H2[i*M];
		tempW[i] = F77_CALL(ddot)(&M, readPtr1, &inci,Mu2, &incj);
		
	}
		 wald[0] = F77_CALL(ddot)(&M, tempW, &inci,Mu2, &incj);
	for(i=1;i<M;i++)
	{
    // blup collection
		//Rprintf("test Used: %d\n",Used[i-1]);
		index				= Used[i-1] - 1;
		Beta[basisMax*2 + index]	= Mu2[i]/Scales[index];
		Beta[basisMax*3 + index]  = SIGMA2[i*M + i]/(Scales[index]*Scales[index]);
	}
	//

	intercept[0]	= Mu2[0];
	intercept[1]	= SIGMA2[0];
	//Rprintf("fEB computation compelete!\n");
	
	Free(Scales);
	Free(Used);
	Free(Mu2);
	Free(SIGMA2);
	Free(H2);
	Free(Alpha);
	Free(PHI2);
	Free(iteration);	
	Free(m);
	Free(tempW);
	
}




/************** outputs are passed by COPY in R, cann't dynamic realloc memory **************************/
/************** Not a problem in C */
// function [Used,Mu2,SIGMA2,H2,Alpha,PHI2]=fEBBinaryMexBmNe(BASIS,Targets,PHI2,Used,Alpha,Scales,a,b,Mu2,iter)
void fEBBinaryMexBmNeEN(int *Used, double *Mu2, double *SIGMA2, double *H2, double *Alpha, double *PHI2,
				double *BASIS, double * Targets, double *Scales, double *a_Lambda,double *b_Alpha,
				int *iteration, 
				int *n, int *kdim, int *m,double * LOGlikelihood,int basisMax,int verbose)
{
    //basis dimension
   int N,K,M_full,N_unused,M,i,j,iter; //L,
   int N_used = 1;
   N					= *n;			// row number
    K					= *kdim;		// column number
    //M_full				= K*(K+1)/2;
	M_full				= K;
    //kk					= K;

    double *beta	= (double *)Calloc(N,double);
	int *Unused = (int *) Calloc(M_full,int);
    iter				= *iteration;
	//Rprintf("Iteration number: %d\n",iter);
    const int	ACTION_REESTIMATE       = 0;			
	const int	ACTION_ADD          	= 1;
	const int 	ACTION_DELETE        	= -1;
    const int   ACTION_TERMINATE        = 10;    
    	
		
	//Block Coordinate Ascent Algorithm Parameters;
	const double 	nAdd 				= 0.90; //only deltaML >nAdd will be considered as candidate
	double nAdd_real 					= nAdd;
	const double 	MLdelta 			= 0.001; //only deltaML > MLdelta will be selected for action.
	
	//Block Coordinate Ascent Algorithm Parameters;
	double cutoff;
	int *nuAction;
	nuAction = (int*) Calloc(M_full,int);
	int nUpdate, iU;// nTerminate;
	
	
    //[Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic) 
    int *IniLogic;
	IniLogic			= (int*) Calloc(1,int);
    if (iter<=1)    
    {
        IniLogic[0]     = 0;
        m[0]            = 2;
		M				= m[0];
		N_used			= 1;
		N_unused		= M_full -1;

    }else
    {
		IniLogic[0]    = 1;
        M				= *m;          //Used + 1
		N_used			= M -1;
		N_unused		= M_full - N_used;
    }
    //
	//Rprintf("N_used is: %d; N_unused:%d, M: %d,sample size: %d \n",N_used,N_unused,M,N);

	fEBInitializationBmNeEN(Alpha, PHI2, Used, Unused, Mu2, BASIS, Targets, Scales, IniLogic, N, m, K);
	
	//******************************************************
	int initial = Used[0];
	int iniRemoved = 1;
	if (iter<=1)    
    {
		initial = Used[0];
		iniRemoved = 0;
    }
	
	if(verbose >3) Rprintf("\t Initialized basis %d, Alpha: %f, \n", Used[0],Alpha[0]);
	//Rprintf("\t Initialized basis %d, Alpha: %f, Mu: %f \n", Used[0],Alpha[0],Mu2[1]);
	//for(i=0;i<10;i++) Rprintf("PHI2: %f \t  %f; BASIS: %f\n",PHI2[i],PHI2[N+i],BASIS[181*N+i]/Scales[181]);
    double *basisCache;
	basisCache         = (double *) Calloc(N*K,double);
	int inci = 1;
	int incj = 1;
	double *readPtr1;//,*readPtr2
	double scal;
	
    for(i=0;i<K;i++)
    {
        for(j=0;j<N;j++)            basisCache[i*N+j] = BASIS[i*N+j]*BASIS[i*N+j];
    }
	
	double *S_in, *Q_in, *S_out, *Q_out;
	S_in				= (double *) Calloc(M_full,double);
	Q_in				= (double *) Calloc(M_full,double);
	S_out				= (double *) Calloc(M_full,double);
	Q_out				= (double *) Calloc(M_full,double);
    //[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,Intercept] ...
    //                   	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE)
	fEBCatFullStatBmNeEN(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2, 
				Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);
//Rprintf("\t SIGMA2: %f, %f \t H2: %f, %f\n\t %f, %f, \t %f, %f\n",SIGMA2[0],SIGMA2[1], H2[0], H2[1],SIGMA2[2],SIGMA2[3],H2[2],H2[3]);
//Rprintf("\t 182th S_out: %f, Q_out: %f\n",S_out[181],Q_out[181]);
//for(i=0;i<10;i++) Rprintf("S_out: %f \t Q_out: %f\n",S_out[i],Q_out[i]);
    //              For: [DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLBmNeEN(Used,Unused,S_out,Q_out,Alpha,a,b);
    double *DeltaML, *AlphaRoot,deltaLogMarginal,*phi,newAlpha;//oldAlpha,
    double *tmp,*tmpp;//s_ii,,mu_iMujj,deltaInv,kappa,
	
	//double oldAlpha,s_ii,*readPtr2,mu_i,Mujj,deltaInv,kappa;
    //
	int *Action, *anyToDelete;
		int selectedAction = -10;
	anyToDelete			= (int*) Calloc(1,int);
	DeltaML				=	(double *) Calloc(M_full,double);
	AlphaRoot			=	(double *) Calloc(M_full,double);
	Action				= (int *) Calloc(M_full,int);
  	phi					= (double *) Calloc(N,double);
    
  	tmp					= (double *) Calloc(basisMax,double);
  	tmpp				= (double *) Calloc(basisMax,double);

    int nu,jj,index;
    jj					= -1;
    int anyWorthwhileAction,UPDATE_REQUIRED;
    // mexPrintf("cols:%d \n",M);
  	//
    int i_iter;
    i_iter              = 0;
    int LAST_ITERATION  = 0;
	double logLikelihood,dL,logL0;
	logLikelihood		= 1e-30;
	dL					= 1e-30;
	double *PHI_Mu;
	PHI_Mu				= (double*) Calloc(N,double);
	if(verbose>3 && iter ==1)
	{
		Rprintf("check point 3: before loop,initial number of basis:%d\t Basis initialized: ",M);
		for(i=0;i<M;i++) Rprintf("%dth:%d\t",i,Used[i]);
		Rprintf("\n");
		
	}
int itMax = 100;
if(iter==1) itMax = 10;	
    while(LAST_ITERATION!=1)
    {
        i_iter						= i_iter + 1;

		//if(verbose >4) Rprintf("\t inner loop %d; number of basis: %d ; Alpha:\n",i_iter, M);
		
		if(verbose >4) Rprintf("\t inner loop %d; number of basis: %d \t actionStatus: %d \tiniRemoved: %d\n",i_iter, M,selectedAction,iniRemoved);
        //M							= N_used + 1;
     	//N_unused					= M_full - N_used;
		logL0						= logLikelihood;
        //[DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLBmNeEN(Used,Unused,S_out,Q_out,Alpha,a,b);
		fEBDeltaMLBmNeEN(DeltaML, Action, AlphaRoot,anyToDelete,Used, Unused, S_out, Q_out, Alpha,
				a_Lambda,b_Alpha, N_used, N_unused,
				&deltaLogMarginal,&nu,N);
				
		//printMat(Alpha,1,N_used);
		//
        //deltaLogMarginal			= MLdelta;
        //nu							= -1;
        //for(i=0;i<M_full;i++)
        //{
        //    if(DeltaML[i]>deltaLogMarginal)
         //   {
         //       deltaLogMarginal    = DeltaML[i];
         //       nu                  = i;
         //   }
        //}
        //
		if(selectedAction          == ACTION_TERMINATE && iniRemoved ==0 &&M>2)
		{
			nu = -1;
			if(verbose>4)  Rprintf("************************** set to Remove the initial basis **************** \n");
		}
			
		
        if(nu==-1 && iniRemoved ==1)
		{
			anyWorthwhileAction     = 0;
			selectedAction          = ACTION_TERMINATE;
		}else if(nu==-1 &&iniRemoved ==0 && M>2) //remove initialized basis
		{
			if(verbose>4)  Rprintf("************************** Removing the initial basis **************** \n");
			anyWorthwhileAction	= 1;
			nu = initial-1;
			Action[nu] = ACTION_DELETE;		
			nUpdate = 1;	
			nuAction[0]= initial - 1;
			iniRemoved = 1;
			selectedAction          = ACTION_DELETE;
		}else
		{
			anyWorthwhileAction	= 1;
			if(Action[nu] == ACTION_ADD)
			{
				nAdd_real			= nAdd;				
			}else
			{
				nAdd_real 			= 1;

			}
			cutoff 					= deltaLogMarginal * nAdd_real;
			if(cutoff< MLdelta)
			{
				cutoff 				= MLdelta;
			}
			nUpdate 				= 0;
			for(i=0;i<M_full;i++)
			{

				if(DeltaML[i]>=cutoff)
				{
					nuAction[nUpdate]= i;
					nUpdate++;
				}
				
			}
			if(Action[nu] == ACTION_DELETE && nUpdate>1) nUpdate = 1; // when deleting co-linear basis: delete one at a time
if(nUpdate==0) anyWorthwhileAction=0;	
		}
if(anyWorthwhileAction==0)  selectedAction = ACTION_TERMINATE;	
		if(anyWorthwhileAction	==1)
		{
			//Block Coordinate Ascent Algorithm Parameters;
			// find a block to update;

			//if(verbose >4) Rprintf("\t deltaML cutoff: %f ; number to be updated: %d \n",cutoff, nUpdate);
			

				
			for(iU=0;iU<nUpdate;iU++) //update each basis
			{
				nu 					= nuAction[iU];
				deltaLogMarginal 	= DeltaML[nu];
				
				selectedAction              = Action[nu];
				newAlpha                    = AlphaRoot[nu];

				//Rprintf("\t ActionOn nu= %d, deltaML: %f, selectedAction: %d\n",nu+1, DeltaML[nu],selectedAction);
				if(selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE)
				{
					index                   = nu + 1; 
					for(i=0;i<N_used;i++)
					{
						if (Used[i]==index)	jj  = i;
					}
				}
				//kk                          = K;                          
				for(i=0;i<K;i++)
				{
					if (i==nu)
					{
						//for(L=0;L<N;L++)    phi[L]  = BASIS[i*N+L]/Scales[i];
						readPtr1 = &BASIS[i*N];
						scal = 1.0/Scales[i];
						F77_CALL(dcopy)(&N,readPtr1,&inci,phi,&incj);
						F77_CALL(dscal)(&N,&scal,phi,&inci);
						
					}
					/*else if(i<(K-1))
					{
						for(j=i+1;j<K;j++)
						{
							if(kk==nu)
							{
								for(L=0;L<N;L++)    phi[L] =BASIS[i*N+L]*BASIS[j*N+L]/Scales[kk];
							}
							kk              = kk + 1;
						}
					}*/
				}
				//newAlpha                    = AlphaRoot[nu];
				//if(anyWorthwhileAction==0)  selectedAction = ACTION_TERMINATE;
				if(selectedAction==ACTION_REESTIMATE)
				{
					if (fabs(log(newAlpha)-log(Alpha[jj]))<=1e-3 && anyToDelete[0] ==0)
					{	
						//cout<<"Terminated by small deltaAlpha on basis:"<<nu<<endl;
						selectedAction		= ACTION_TERMINATE;
					}
				}
				//Rprintf("\t selectedAction: %d\n",selectedAction);
				//
				UPDATE_REQUIRED				= 0;
				if(selectedAction==ACTION_REESTIMATE)
				{
					if(verbose>4) Rprintf("\t\t Action: Reestimate : %d \t deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
										nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					UPDATE_REQUIRED= ActionResBmNeEN(BASIS,Scales,PHI2, Alpha, SIGMA2, newAlpha,
							Mu2, S_in, Q_in, beta, jj, N, M,M_full);			
							

					//UPDATE_REQUIRED			= 1;
				}
				/////////////////////////////////////////////////////////////////////////////////
				else if(selectedAction==ACTION_ADD)
				{
					if(verbose>4) Rprintf("\t\t Action:add : %d \t deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
										nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					//Rprintf("\t\t newAlpha: %f\n",newAlpha);
					// B_phi*PHI2*SIGMA2        tmp = B_phi*PHI2 
					index					= N_used + 1;
					if(index >basisMax && iter>1) {
						Rprintf("bases: %d, warning: out of Memory, alloc more to Neffect!\n",index);
						return;
					}
					
					//BASIS_B_Phi	=BASIS_Phi*beta;
					UPDATE_REQUIRED		= ActionAddBmNeEN( BASIS, Scales, PHI2, phi, beta, Alpha,
						newAlpha, SIGMA2, Mu2, S_in, Q_in, nu, M_full, N, K, M);
					
					Used[N_used]			= nu + 1;						//new element
					N_used					= N_used + 1;

					//
					N_unused				= N_unused - 1;
					for(i=0;i<N_unused;i++)
					{                
						if(Unused[i]== (nu + 1))		Unused[i] =Unused[N_unused];
					}
					m[0]					= N_used + 1;
					M						= m[0];
					
					//for(i=1;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i-1],Mu2[i]);

				}
				//////////////////////////////////////////////////////////////////////////////////
				else if(selectedAction==ACTION_DELETE)
				{
					if(verbose>4) Rprintf("\t\t Action: delete : %d deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
										nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					
//Rprintf("jj: %d, Alpha: %f;M: %d; \t Nused: %d\n",jj, Alpha[jj], M, N_used);
					M = N_used + 1;
					//printMat(SIGMA2,M,M);
					UPDATE_REQUIRED = ActionDelBmNeEN(BASIS, Scales, PHI2, Alpha, SIGMA2, 
											Mu2, S_in, Q_in, beta, jj, N, M, M_full);
					//printMat(SIGMA2,M,M);
					index = N_used  -1;	
					//
				
					//Used; Unused; jj max = 480;
					Used[jj]				= Used[index];
					N_used					= index;

					//
					N_unused				= N_unused + 1;
					Unused[N_unused -1]		= nu + 1;

					m[0]					= N_used + 1;
					M						= m[0];
					//UPDATE_REQUIRED			= 1;
				//for(i=1;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i-1],Mu2[i]);
				//*********************************************************
				if((nu+1)==initial) iniRemoved =1;
				}
				//Rprintf("\t\t Update_required: %d\n",UPDATE_REQUIRED);
				//

				if(iU==(nUpdate-1) )
				{
					UPDATE_REQUIRED			= 1;
				}else
				{
					UPDATE_REQUIRED			=0;
				}
				
				if(UPDATE_REQUIRED==1)
				{
					fEBCatFullStatBmNeEN(beta, SIGMA2, H2, S_in, Q_in, S_out, Q_out,  BASIS,Scales, PHI2, 
					Targets, Used, Alpha, Mu2, basisCache,n, m, kdim);

				}
				//printMat(SIGMA2,M,M);
			} // //end of for(iU)
		}
	
		//Rprintf("\t\t selected Action: %d\n",selectedAction);
		//

		
        if(selectedAction==ACTION_TERMINATE && iniRemoved == 1) LAST_ITERATION =1;
		if((i_iter==itMax && M==2) || i_iter>itMax) LAST_ITERATION =1;
//*******************************************************************


        if(i_iter==itMax )
		{
			selectedAction=ACTION_TERMINATE;
			if(verbose>4) Rprintf("\t\t Reaching the maximum iteration\n");
			
		}
		logLikelihood = 0;
		for(i = 0;i<N;i++)
		{
			PHI_Mu[i]				= 0;
			//for(j = 0;j<M;j++)		PHI_Mu[i]	= PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
			readPtr1 = &PHI2[i];
			PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj); 
			
			
			logLikelihood			= logLikelihood + Targets[i]*log(exp(PHI_Mu[i])/(1+exp(PHI_Mu[i]))) + 
										(1-Targets[i])*log(1/(1+exp(PHI_Mu[i])));
		}
		dL							= fabs((logLikelihood - logL0)/logL0);
		//Rprintf("\t\t likelihoodChange: %f\n",dL);
		if(dL <1e-3 )
		{
			selectedAction=ACTION_TERMINATE;
			//LAST_ITERATION = 1;
			if(verbose>4) Rprintf("\t\t small global likelihood change.\n");		
		}		//Rprintf("\t\t Last_iteration value: %d\n",LAST_ITERATION);

    }

	LOGlikelihood[0] = logLikelihood;
    	
		
	Free(beta);
	Free(Unused);
	Free(IniLogic);
	Free(basisCache);
	Free(S_in);
	Free(Q_in);
	Free(S_out);
	Free(Q_out);
	Free(anyToDelete);
	Free(DeltaML);	
	Free(AlphaRoot);	
	Free(Action);	
	Free(phi);	
	Free(tmp);	
	Free(tmpp);	
	Free(PHI_Mu);	
	//Block Coordinate Ascent Algorithm Parameters;
	Free(nuAction);
	
}

//888888888888888888888888888888888888888888888888888888888888888888888888888888888888
int ActionAddBmNeEN( double* BASIS, double*scales, double*PHI, double*Phi,
			double *beta, double* Alpha, double newAlpha, double*SIGMA, double*Mu, double*S_in,
			double*Q_in, int nu, int M_full,int N, int K, int M)
{
	int N_used = M - 1;
	int index = M + 1;
	int M1M1 = index*index;
	double *SIGMANEW = (double *) Calloc(M1M1,double);
	
	//Rprintf("\t add basis %d, Alpha: %f, N:%d\tK:%d\tM:%d\tM_full:%d \tN_used:%d\tM1:%d\n",
	//				nu,newAlpha,N,K,M,M_full,N_used,index);
	
	double *B_Phi		= (double *) Calloc(N,double);//beta.*Phi;
	double *BASIS_B_PHI		= (double *) Calloc(M,double); //BASIS'*B*PHI: --> 1xM for one column of X in a K loop; total K x M
	double *BASIS_B_Phi		= (double *) Calloc(M_full,double);
	double *mCi				= (double *) Calloc(M_full,double);
	double *z				= (double *) Calloc(N,double);
	//int kk					= K;
	int i,j,h;

	double*   	tmp			= (double *) Calloc(M,double);
  	double*		tmpp		= (double *) Calloc(M,double);
	double * s_i			= (double *) Calloc(M,double);
		
	double s_ii,mu_i,TAU;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1;
	double b_blas = 1.0;

	//lapack end
	//B_Phi_lowcase       = Phi_lowcase.*beta;
	//BASIS_B_Phi         = BASIS'*B_Phi;
	//Rprintf("\t\t Inside ActionAddGmNeEN: index: %d\n",index);

	for(j=0;j<N;j++)		B_Phi[j]			= beta[j]*Phi[j]; //N   x 1
	
	
	for (i					= 0;i<K;i++)
	{
		BASIS_B_Phi[i]		= 0;

		for(h=0;h<N;h++)
		{
			z[h]			= BASIS[i*N+h]*B_Phi[h];
			BASIS_B_Phi[i]	= BASIS_B_Phi[i] + z[h];
		}
		BASIS_B_Phi[i] 		= BASIS_B_Phi[i]/scales[i];
		
		//if(i<K-1)
		//{
		//	for (j					=(i+1);j<K;j++)
		//	{
		//		BASIS_Phi[kk]		= 0; 
		//		for(h=0;h<N;h++)	BASIS_Phi[kk]	= BASIS_Phi[kk] + BASIS[j*N + h] *z[h];
		//		BASIS_Phi[kk]		= BASIS_Phi[kk]/scales[kk];
		//		BASIS_B_Phi[kk]		= beta[0]*BASIS_Phi[kk];
		//		kk					= kk +1;
		//	}
		//}
	}

	
	
	//tmp						= ((B_Phi_lowcase'*PHI)*SIGMA)';				
	for(i=0;i<M;i++)
    {
        tmp[i]					= 0;
        //for(j=0;j<N;j++)		tmp[i]			= tmp[i] + beta[0]*Phi[j]*PHI[i*N+j]; //M+1   x 1
		readPtr1 	= &PHI[i*N];
		tmp[i]  = F77_CALL(ddot)(&N, readPtr1, &inci,B_Phi,&incj);
    }	
    for(i=0;i<M;i++)
    {
        tmpp[i]					= 0;
        //for(j=0;j<M;j++)	tmpp[i]				= tmpp[i] + tmp[j]*SIGMA[i*M+j];
		readPtr1 	= &SIGMA[i*M];
		tmpp[i]  = F77_CALL(ddot)(&M, readPtr1, &inci,tmp,&incj);	
    }//tmpp is tmp in matlab.
	//Alpha: M -> M + 1
	Alpha[N_used]				= newAlpha;                 //new element
	//PHI2	M	-> M + 1
	//for(i=0;i<N;i++)		PHI[M*N + i]		= Phi[i];		//new column
	readPtr1 	= &PHI[M*N];
	F77_CALL(dcopy)(&N,Phi,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x

	s_ii					= 1.0/(newAlpha + S_in[nu]);   
	
	mu_i					= s_ii*Q_in[nu];
    //for(i=0;i<M;i++)		Mu[i]				= Mu[i] - mu_i*tmpp[i];
	b_blas = -mu_i;
	F77_CALL(daxpy)(&M, &b_blas,tmpp, &inci,Mu, &incj); 
	//daxpy(n, a, x, incx, y, incy) y := a*x + y
	
    Mu[M]					= mu_i;							//new element
	

	//for(i=0;i<M;i++)		s_i[i]				= - tmpp[i]	*s_ii;
	F77_CALL(dcopy)(&M,tmpp,&inci,s_i,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	b_blas = -s_ii;
	F77_CALL(dscal)(&M,&b_blas,s_i,&inci); 		//dscal(n, a, x, incx) x = a*x
	

	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			TAU				= -s_i[i]*tmpp[j];
			SIGMANEW[j*index+i]					= SIGMA[j*M+i] + TAU;
		}
	}

	//for(i=0;i<M;i++)	
	//{
	//	SIGMANEW[M*index+i] = s_i[i];
	//	SIGMANEW[i*index+M] = s_i[i];
	//}
	readPtr1 = &SIGMANEW[M*index];
	F77_CALL(dcopy)(&M,s_i,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	readPtr1 = &SIGMANEW[M];
	F77_CALL(dcopy)(&M,s_i,&inci,readPtr1,&index);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	SIGMANEW[M*index+M]		= s_ii;

	//SIGMA = SIGMANEW;
	F77_CALL(dcopy)(&M1M1,SIGMANEW,&inci,SIGMA,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	//S_in                = S_in - s_ii * mCi.^2;
	//Q_in                = Q_in - mu_i * mCi;
	//mCi                 = BASIS_B_Phi_lowcase - beta*BASIS_PHI*tmp;
	//mCi 	= X'Bx - X'BZSZ'Bx;
	//BASIS_B_PHI = BASIS'*(PHI.*beta*ones(1,M));
	
	//BASIS_B_PHI;
	
	double temp;
	for(i=0;i<M_full;i++)
	{
		//xi.*beta;
		for(j=0;j<N;j++) z[j] = BASIS[i*N+j]*beta[j];
		for(j=0;j<M;j++)
		{
			readPtr1 = &PHI[j*N];
			BASIS_B_PHI[j] = F77_CALL(ddot)(&N, z, &inci,readPtr1,&incj);
			BASIS_B_PHI[j] = BASIS_B_PHI[j]/scales[i];
			
		}

		temp				= F77_CALL(ddot)(&M, BASIS_B_PHI, &inci,tmpp,&incj); 
		//readPtr1 			= &BASIS_PHI[i];
		//temp 				= F77_CALL(ddot)(&M, readPtr1, &M_full,tmpp,&incj); 
		
		mCi[i]				= BASIS_B_Phi[i] -  temp;
		//BASIS_PHI[M*M_full + i]					= BASIS_Phi[i];
		S_in[i]				= S_in[i] - mCi[i]*mCi[i]*s_ii;
		Q_in[i]				= Q_in[i] - mu_i*mCi[i];
	}


	
	int UPDATE_REQUIRED		= 1;
	Free(B_Phi);
	Free(BASIS_B_Phi);	
	Free(BASIS_B_PHI);	
	Free(mCi);
	Free(z);	
	Free(tmp);
	Free(tmpp);	
	Free(s_i);
	Free(SIGMANEW);
	
	return  UPDATE_REQUIRED; 
 }

 //PHI: M x 1
 //Alpha N_used x 1
 //SIGMA: M x M
 //Mu: M x 1
 
int ActionDelBmNeEN(double* BASIS, double*scales,double*PHI, double*Alpha, double*SIGMA,
		double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full)
 {
	 int i,j;
	int index				= M - 1;
	//int N_used 			= M - 1;
	int M1M1 				= (M-1)*(M-1);
	double *tempSIGMA 		= (double *) Calloc((M*M),double);
	double *SIGMANEW 		= (double *) Calloc(M1M1,double);
	double *z				= (double *) Calloc(N,double);
	double *BASIS_B_PHI		= (double *) Calloc(M,double); //BASIS'*B*PHI: --> 1xM for one column of X in a K loop; total K x M
	
//Rprintf("\t INSIDE DELETION FUNCTION\n");
	//lapack
	int inci 				=1;
	int incj 				=1;
	double *readPtr1, *readPtr2;
	//lapack end	

    int j1  				= jj + 1;
	//PART 1. Mu;
    double Mujj;        
	Mujj					= Mu[j1];
	for(i=0;i<M;i++)  Mu[i] = Mu[i] - Mujj*SIGMA[j1*M +i]/SIGMA[j1*M + j1];
	//PART 2. S, Q
	//jPm         = (BASIS_B_PHI * s_j); 		BASIS'*B*PHI--> 1xM for one column of X in a K loop;
    //S_in        = S_in + jPm.^2 / s_jj;
    //Q_in        = Q_in + jPm * mu_j / s_jj;
	// By definition, PHI, SIGMA in calculating S_in,Q_in should be before removing jj.
	double temp;
	for(i=0;i<M_full;i++)
	{
		//BASIS_B_PHI: xi.*beta;
		for(j=0;j<N;j++) z[j] = BASIS[i*N+j]*beta[j];
		for(j=0;j<M;j++)
		{
			readPtr1 = &PHI[j*N];
			BASIS_B_PHI[j] = F77_CALL(ddot)(&N, z, &inci,readPtr1,&incj);
			BASIS_B_PHI[j] = BASIS_B_PHI[j]/scales[i];
			
		}
		//jPm = (BASIS_B_PHI * s_j); 	
		temp				= 0;
		//for(j=0;j<M;j++)	temp				= temp + BASIS_PHI[j*M_full + i]*SIGMA[jj*M + j];
		for(j=0;j<M;j++)	temp				= temp + BASIS_B_PHI[j]*SIGMA[j1*M + j];
		S_in[i]				= S_in[i] +  pow(temp,2)/SIGMA[j1*M + j1];
		Q_in[i]				= Q_in[i] +  temp *Mujj /SIGMA[j1*M + j1];
		//jPm = beta*temp
	}
    //PART 3. SIGMA.
    //tmp				= s_j/s_jj;
    //SIGMANEW		= SIGMA - tmp*s_j'; = SIGMA - SIGMA[,jj]*SIGMA[,jj]'/SIGMA[jj,jj];
    //SIGMANEW(j,:)	= [];
    //SIGMANEW(:,j)	= [];
	
	//------------------------------------------------------------------------------
	//BLOCK MIGRATION OF SIGMANEW ------------------JUN142013 -------IN EBELASTICNET
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)	tempSIGMA[j*M + i]	= SIGMA[j*M + i] - SIGMA[j1*M+i]*SIGMA[j1*M+j]/SIGMA[j1*M+j1];
	}
	//printMat(tempSIGMA,M,M);
	for(i=0;i<index;i++)
	{
		for(j=0;j<index;j++) SIGMANEW[j*index + i]	= tempSIGMA[j*M + i];
	}
	if(j1 != index)// incase the one to be deleted is the last column.
	{
		//
		Alpha[jj]			= Alpha[index-1];           //Alpha: M -> M - 1
		//
		Mu[j1]				= Mu[index];    
		//for(i=0;i<N;i++)		PHI[jj*N + i]		= PHI[index*N+i];		//M -> M -1
		readPtr1 			= &PHI[j1*N];
		readPtr2 			= &PHI[index*N];
		F77_CALL(dcopy)(&N,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
		
		
		//step 1: last column
		readPtr1 			= &SIGMANEW[j1*index];
		readPtr2 			= &tempSIGMA[index*M];
		F77_CALL(dcopy)(&index,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
		//step 2: prepare for last row
		tempSIGMA[j1*M + M-1] = tempSIGMA[M*M-1];

		//last step 	
		readPtr1 = &SIGMANEW[j1];
		readPtr2 = &tempSIGMA[M-1];
		F77_CALL(dcopy)(&index,readPtr2,&M,readPtr1,&index);  //dcopy(n, x, incx, y, incy) ---> y = x
	}

	//BLOCK MIGRATION OF SIGMANEW ------------------JUN142013 -------IN EBELASTICNET	

	//After removing jj, PHI: N x (M-1)
		//SIGMA: (M-1)*(M-1)
	
	//SIGMA = SIGMANEW;	
	F77_CALL(dcopy)(&M1M1,SIGMANEW,&inci,SIGMA,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
    
	//BASIS_PHI(:,jj) = [];
	//for(i=0;i<M_full;i++) BASIS_PHI[jj*M_full + i]	= BASIS_PHI[index*M_full+i];
	//readPtr1 = &BASIS_PHI[jj*M_full];
	//readPtr2 = &BASIS_PHI[index*M_full];
	//F77_CALL(dcopy)(&M_full,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	int UPDATE_REQUIRED		=1;
	Free(tempSIGMA);
	Free(SIGMANEW);
	Free(z);
	Free(BASIS_B_PHI);
	return  UPDATE_REQUIRED;
}//end of ACTION_DELETE


//					UPDATE_REQUIRED= ActionResBmNeEN(BASIS,Scales,PHI2, Alpha, SIGMA2, newAlpha,
//								Mu2, S_in, Q_in, beta, jj, N, M,M_full);
//jj: index from 0 - N_used;
int ActionResBmNeEN(double* BASIS, double*scales, double*PHI, double*Alpha, double*SIGMA, double newAlpha,
		double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full)
 {
	int i,j;
	//int index				= M - 1;
	//int N_used 				= M - 1;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1;//, *readPtr2
	//lapack end
	//Rprintf("\t\t INSIDE RES, jj: %d \t M: %d\n",jj,M);
	int MM = M *M;
	double *SIGMANEW = (double *) Calloc(MM,double);
	double *BASIS_B_PHI		= (double *) Calloc(M,double); //BASIS'*B*PHI: --> 1xM for one column of X in a K loop; total K x M
	double *z				= (double *) Calloc(N,double);
	int j1 = jj +1;
	//printMat(Alpha,N_used,1);
	//Rprintf("\t\t INSIDE RES, jj: %d \t M: %d \tx: %f \t y: %f\n",jj,M,Alpha[jj], newAlpha);
	double oldAlpha			= Alpha[jj];
	Alpha[jj]				= newAlpha;

	double deltaInv			= 1.0/(newAlpha-oldAlpha);
	double kappa			= 1.0/(SIGMA[j1*M+j1] + deltaInv);
	double Mujj				= Mu[j1];
	
	//for(i=0;i<M;i++)		Mu[i]    = Mu[i] - Mujj *kappa * SIGMA[jj*M+i];
	readPtr1 				= &SIGMA[j1*M];
	double b_blas 			= -Mujj * kappa;

	F77_CALL(daxpy)(&M, &b_blas,readPtr1, &inci,Mu, &incj); //daxpy(n, a, x, incx, y, incy) y := a*x + y

	//	SIGMANEW	=SIGMA-sjj;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)	SIGMANEW[j*M + i] = SIGMA[j*M + i] - kappa * SIGMA[j1*M+i]*SIGMA[j1*M+j];
	}
	
	//SIGMA = SIGMANEW;	
	F77_CALL(dcopy)(&MM,SIGMANEW,&inci,SIGMA,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	//S_in		= S_in + kappa*(beta*BASIS_PHI * s_j).^2;
	//Q_in		= Q_in - beta*BASIS_PHI*(deltaMu);           
	double temp;
	for(i=0;i<M_full;i++)
	{
		//BASIS_B_PHI: xi.*beta;
		for(j=0;j<N;j++) z[j] = BASIS[i*N+j]*beta[j];
		for(j=0;j<M;j++)
		{
			readPtr1 = &PHI[j*N];
			BASIS_B_PHI[j] = F77_CALL(ddot)(&N, z, &inci,readPtr1,&incj);
			BASIS_B_PHI[j] = BASIS_B_PHI[j]/scales[i];
			
		}
		
		temp	= 0;
		//tmp                     = kappa*s_j;
		//deltaMu                 = -Mu(jj)*tmp;
		
		//for(j=0;j<M;j++) temp = temp + BASIS_PHI[j*M_full + i]*SIGMA[jj*M + j];
		//readPtr1 			= &BASIS_PHI[i];
		//readPtr2 			= &SIGMA[jj*M];
		//temp = F77_CALL(ddot)(&M, readPtr1, &M_full,readPtr2, &inci);
		for(j=0;j<M;j++) temp = temp + BASIS_B_PHI[j]*SIGMA[j1*M + j];
		S_in[i]				= S_in[i] +  pow(temp,2)*kappa;
		Q_in[i]				= Q_in[i] + Mujj *kappa*temp;
	}
	

	int UPDATE_REQUIRED		=1;

	Free(SIGMANEW);
	Free(z);
	Free(BASIS_B_PHI);
	return UPDATE_REQUIRED;
}//end of ACTION_res




 
 
 
 /****************************************************************************/

// [Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic)  //IniLogic: whether input is empty or not

void fEBInitializationBmNeEN(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, double *Scales, int * initial, int N, int *m, int K)
{
    //basis dimension
    int M,M_full,N_used,i,j,kk;
  	
    //M_full				= K*(K+1)/2;
	M_full					= K;
	int IniLogic			= *initial;
    //INPUT
    if(IniLogic==0)							// is empty
    {
		m[0]				= 2;
		M					= m[0];
        N_used				= 1;
    }else									// not empty
    {
       	N_used              = m[0]-1;
		M					= m[0];
    }
    //output
    const double init_alpha_max     = 1e3;
    const double init_alpha_min     = 1e-3;    
    
	if(IniLogic==0)            // is empty
    {
		//Rprintf("\t Inside Initialization, M: %d, K: %d\n",M, K);
        double *TargetPseudo,proj_ini,proj;
        int loc1 = 0;
		//int loc2 = 0;
		int inci =1;
		int incj = 1;
		double *readPtr1;
		double scal;
		//*********************************************************************
		//double *project   	= (double *) Calloc(K,double);
		//const double 	nAdd 				= 0.9; 
		
        TargetPseudo		= (double *) Calloc(N,double);
        //for(i=0;i<N;i++)            TargetPseudo[i]     = 2*Targets[i] -1;
		scal = -1;
		F77_CALL(dcopy)(&N,&scal,&loc1,TargetPseudo,&inci);
		scal = 2;
		F77_CALL(daxpy)(&N, &scal,Targets, &inci,TargetPseudo, &incj);
		
        //kk					= K;
        proj_ini			= 0;
		Used[0]				= 1;
        for(i=0;i<K;i++)
        {

         	//proj			= 0;
            //for(j=0;j<N;j++)                proj    = proj + BASIS[i*N+j]*TargetPseudo[j];
			readPtr1 = &BASIS[i*N];
			proj = F77_CALL(ddot)(&N, readPtr1, &inci,TargetPseudo, &incj);
			
            proj			= proj/Scales[i];
			//*********************************************************************
			//project[i] 		= fabs(proj);
            if(fabs(proj) < fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                loc1        = i;
                //loc2        = i;
                Used[0]		= i + 1;
            }
        }
		//*********************************************************************		
		//printMat(project,1,K);
		/*
		double cutoff = proj_ini * nAdd;
		for(i=0;i<K;i++)
		{
			if(project[i] >=cutoff) 	project[i] = 0;
			
		}
		proj_ini = 0;
		for(i=0;i<K;i++)
        {
		
            proj			= project[i];
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                loc1        = i;
                //loc2        = i;
                Used[0]		= i + 1;
            }
        }
		*/
		//printMat(project,1,K);
        /*for(i=0;i< (K - 1);i++)
        {
            for(j= (i+ 1); j<K; j++)
            {
                proj			= 0;
                for(k=0;k<N;k++)            proj    = proj + BASIS[i*N+k]*BASIS[j*N+k]*TargetPseudo[k];
                proj            = proj/Scales[kk];
                if(fabs(proj) > fabs(proj_ini))
                {
                    proj_ini    = proj;
                    loc1        = i;
                    loc2        = j;
                    Used[0]		= kk  + 1;
                }
                kk              = kk + 1;
            }
        }*/
        //PHI2, duplicate for linear solver
		double *PHIqr;
		PHIqr					= (double *) Calloc(N*M,double);
		scal = 1;
		i = 0;
		F77_CALL(dcopy)(&N,&scal,&i,PHI2,&inci);
		F77_CALL(dcopy)(&N,&scal,&i,PHIqr,&inci);
     	//for(i=0;i<N;i++)		
		//{
		//	PHI2[i]         = 1;
		//	PHIqr[i]		= 1;
		//}
        
        //
        double * PHI;
        PHI						= (double *) Calloc(N,double);
        /*if(loc1==loc2)
        {
            for(i=0;i<N;i++)
            {
                PHI[i]			= BASIS[loc1*N+i]/Scales[loc1];
			      PHI2[N+i]		= PHI[i];
					PHIqr[N+i]		= PHI[i];
            }*/
		readPtr1 = &BASIS[loc1*N];
		F77_CALL(dcopy)(&N,readPtr1,&inci,PHI,&incj);
		scal = 1/(double)Scales[loc1];
		F77_CALL(dscal)(&N,&scal, PHI,&inci);
		readPtr1 = &PHI2[N];
		F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
		readPtr1 = &PHIqr[N];
		F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
			
        /*}else
        {
            index				= Used[0] -1;
           	for(i=0;i<N;i++)
            {
                PHI[i]			= BASIS[loc1*N+i]*BASIS[loc1*N+i]/Scales[index];
                PHI2[N+i]		= PHI[i];
				PHIqr[N+i]		= PHI[i];
            }
        }*/
		double *logout;
        logout                      = (double *) Calloc(N,double);
        for (i=0;i<N;i++)            logout[i]               = log(((TargetPseudo[i] * 0.9 + 1)/2)/(1-(TargetPseudo[i] * 0.9 + 1)/2));
		// Call function
		LinearSolverBmNeEN(PHIqr, logout, N, M, Mu2);

        if(Mu2[1] == 0) Alpha[0] 	= 1;
        else            Alpha[0]    = 1/(Mu2[1]*Mu2[1]);
		//Alpha[0] 		= 8.6102;
        if(Alpha[0]< init_alpha_min) Alpha[0]				= init_alpha_min;
        if(Alpha[0]> init_alpha_max) Alpha[0]				= init_alpha_max;
		
		//Free(project);
		Free(TargetPseudo);
		Free(PHIqr);
		Free(PHI);
		Free(logout);
		
		
	}

	int IsUsed			= 0;
    kk                  = 0;
    for(i=0;i<M_full;i++)
    {
        IsUsed          = 0;
        for(j=0;j<N_used;j++)
        {
            //index   = Used[j];
            if ((i+1)==Used[j])     IsUsed  = 1;
        }
        if(IsUsed==0)     
        {
            Unused[kk]  = (i+ 1);
            kk          = kk + 1;
        }
    }
	
	
	
}






void fEBInitializationBmNeEN2(double *Alpha, double * PHI2, int *Used, int *Unused, double *Mu2,
				double *BASIS, double *Targets, double *Scales, int * initial, int N, int *m, int K)
{
    //basis dimension
    int M,M_full,N_used,i,j,kk;
  	
    //M_full				= K*(K+1)/2;
	M_full					= K;
	int IniLogic			= *initial;
    //INPUT
    if(IniLogic==0)							// is empty
    {
		m[0]				= 2;
		M					= m[0];
        N_used				= 1;
    }else									// not empty
    {
       	N_used              = m[0]-1;
		M					= m[0];
    }
    //output
    const double init_alpha_max     = 1e3;
    const double init_alpha_min     = 1e-3;    
    
	if(IniLogic==0)            // is empty
    {
		//Rprintf("\t Inside Initialization, M: %d, K: %d\n",M, K);
        double *TargetPseudo,proj_ini,proj;
        int loc1 = 0;
		//int loc2 = 0;
		int inci =1;
		int incj = 1;
		double *readPtr1;
		double scal;
		//*********************************************************************
		double *project   	= (double *) Calloc(K,double);
		const double 	nAdd 				= 0.9; 
		
        TargetPseudo		= (double *) Calloc(N,double);
        //for(i=0;i<N;i++)            TargetPseudo[i]     = 2*Targets[i] -1;
		scal = -1;
		F77_CALL(dcopy)(&N,&scal,&loc1,TargetPseudo,&inci);
		scal = 2;
		F77_CALL(daxpy)(&N, &scal,Targets, &inci,TargetPseudo, &incj);
		
        //kk					= K;
        proj_ini			= 0;
		Used[0]				= 1;
        for(i=0;i<K;i++)
        {

         	//proj			= 0;
            //for(j=0;j<N;j++)                proj    = proj + BASIS[i*N+j]*TargetPseudo[j];
			readPtr1 = &BASIS[i*N];
			proj = F77_CALL(ddot)(&N, readPtr1, &inci,TargetPseudo, &incj);
			
            proj			= proj/Scales[i];
			//*********************************************************************
			project[i] 		= fabs(proj);
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                loc1        = i;
                //loc2        = i;
                Used[0]		= i + 1;
            }
        }
		//*********************************************************************		
		//printMat(project,1,K);
		
		double cutoff = proj_ini * nAdd;
		for(i=0;i<K;i++)
		{
			if(project[i] >=cutoff) 	project[i] = 0;
			
		}
		proj_ini = 0;
		for(i=0;i<K;i++)
        {
		
            proj			= project[i];
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                loc1        = i;
                //loc2        = i;
                Used[0]		= i + 1;
            }
        }
		
		//printMat(project,1,K);
        /*for(i=0;i< (K - 1);i++)
        {
            for(j= (i+ 1); j<K; j++)
            {
                proj			= 0;
                for(k=0;k<N;k++)            proj    = proj + BASIS[i*N+k]*BASIS[j*N+k]*TargetPseudo[k];
                proj            = proj/Scales[kk];
                if(fabs(proj) > fabs(proj_ini))
                {
                    proj_ini    = proj;
                    loc1        = i;
                    loc2        = j;
                    Used[0]		= kk  + 1;
                }
                kk              = kk + 1;
            }
        }*/
        //PHI2, duplicate for linear solver
		double *PHIqr;
		PHIqr					= (double *) Calloc(N*M,double);
		scal = 1;
		i = 0;
		F77_CALL(dcopy)(&N,&scal,&i,PHI2,&inci);
		F77_CALL(dcopy)(&N,&scal,&i,PHIqr,&inci);
     	//for(i=0;i<N;i++)		
		//{
		//	PHI2[i]         = 1;
		//	PHIqr[i]		= 1;
		//}
        
        //
        double * PHI;
        PHI						= (double *) Calloc(N,double);

		
			readPtr1 = &BASIS[loc1*N];
			F77_CALL(dcopy)(&N,readPtr1,&inci,PHI,&incj);
			scal = 1/(double)Scales[loc1];
			F77_CALL(dscal)(&N,&scal, PHI,&inci);
			readPtr1 = &PHI2[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
			readPtr1 = &PHIqr[N];
			F77_CALL(dcopy)(&N,PHI,&inci,readPtr1,&incj);
			

		double *logout;
        logout                      = (double *) Calloc(N,double);
        for (i=0;i<N;i++)            logout[i]               = log(((TargetPseudo[i] * 0.9 + 1)/2)/(1-(TargetPseudo[i] * 0.9 + 1)/2));
		// Call function
		LinearSolverBmNeEN(PHIqr, logout, N, M, Mu2);

        if(Mu2[1] == 0) Alpha[0] 	= 1;
        else            Alpha[0]    = 1/(Mu2[1]*Mu2[1]);
		//Alpha[0] 		= 8.6102;
        if(Alpha[0]< init_alpha_min) Alpha[0]				= init_alpha_min;
        if(Alpha[0]> init_alpha_max) Alpha[0]				= init_alpha_max;
		
		Free(project);
		Free(TargetPseudo);
		Free(PHIqr);
		Free(PHI);
		Free(logout);
		
		
	}

	int IsUsed			= 0;
    kk                  = 0;
    for(i=0;i<M_full;i++)
    {
        IsUsed          = 0;
        for(j=0;j<N_used;j++)
        {
            //index   = Used[j];
            if ((i+1)==Used[j])     IsUsed  = 1;
        }
        if(IsUsed==0)     
        {
            Unused[kk]  = (i+ 1);
            kk          = kk + 1;
        }
    }
}

void LinearSolverBmNeEN(double * a, double *logout, int N, int M,double *output)
{
	const int nrhs		= 1;
	const double Rcond	= 1e-5;
	int rank			= M;
	int *jpvt;
	jpvt				= (int * ) Calloc(M,int);
	const int lwork	= M*N + 4*N;
	double * work;
	work				= (double *) Calloc(lwork,double);

	int info			= 0;
	// *************************Call LAPACK library ************************
	F77_CALL(dgelsy)(&N, &M, &nrhs, a, &N, logout, &N,jpvt,&Rcond, &rank, work, &lwork,&info);
	if(info!=0) 	
	{
		Rprintf("Call linear solver degls error!\n");
		return;
	}
	
	//Rprintf("Matrix inversed!\n");	
	//output				= (double *) realloc(NULL,M*sizeof(double));
	int inci = 1;
	int incj = 1;
	F77_CALL(dcopy)(&M,logout,&inci,output,&incj);
//for (i=0;i<M;i++) output[i] = logout[i];
	
	
	Free(jpvt);
	Free(work);	
}

 /// ***********************************************************************************************   


//[Mu2 beta SIGMA2] 	= fEBCatPostModeBmNeEN(PHI2,Targets,Alpha,Mu2);


/// *********[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,BASIS_B_PHI,Intercept] ...
///                       	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE) ************
    //Mu2 is the same size of M in PHI2; one dimension more than Alpha
    //Function fEBCatPostModeBmNeEN:
    //    [Mu2 U beta]     	= fEBCatPostModeBmNeEN(PHI2,Targets,Alpha,Mu2);
    //Targets: nx1

void fEBCatFullStatBmNeEN(double * beta, double * SIGMA2, double * H2, double *S_in, double * Q_in, 
				double * S_out, double * Q_out,  double *BASIS,double * Scales, double *PHI2, 
				double * Targets, int * Used, double *Alpha, double * Mu2, double * BasisCache,
				int *n, int *m, int* kdim)
{
    //basis dimension
    int N,K,M,i,j,p;
   	N					= *n;			// row number
    K					= *kdim;		// column number
    //M_full				= K*(K+1)/2;
    M					= *m;
    //kk					= K;
    
    //output mxArrays; BASIS_B_PHI is not counte
  
    //[Mu2 Ui beta SIGMA2]     	= fEBCatPostModeBmNeEN(PHI2,Targets,Alpha,Mu2);
	fEBCatPostModeBmNeEN(Mu2, beta,SIGMA2, H2, PHI2,Targets, Alpha,N, M);
        
    //	y				= fEBSigmoidBmNeEN(PHI2 * Mu2);
    double *PHI_Mu,*y;  
    PHI_Mu        		= (double *) Calloc(N,double);
	y					= (double *) Calloc(N,double);
	int inci =1;
	int incj =1;
	double *readPtr1; // *readPtr2;
	double scal;
	
	
    for(i=0;i<N;i++)
    {
        //PHI_Mu[i]		= 0;
        //for(j=0;j<M;j++)            PHI_Mu[i]           = PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
		readPtr1 = &PHI2[i];
		PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj); 
		
    }
	fEBSigmoidBmNeEN(y, PHI_Mu,N);
        
    //e=Targets-y
    double *e;
    e					= (double *) Calloc(N,double);
    //for(i=0;i<N;i++)    e[i]        = Targets[i] - y[i];
    F77_CALL(dcopy)(&N,Targets,&inci,e,&incj);
	scal = -1;
	F77_CALL(daxpy)(&N, &scal,y, &inci,e, &incj);    
	
    //Main loop
        //temp parameters: BPvector
    double *BPvector,*temp,*temp1,*temp2;
    BPvector			= (double *) Calloc(M,double);
    temp				= (double *) Calloc(M,double);
    temp1				= (double *) Calloc(M*N,double);
    temp2				= (double *) Calloc(N,double);
    double tempSum,BBsquare,tempZE;
    //Cache BASIS.^2	outside the Inner loop of the program: save computation
	//double Print;
	//Rprintf("K is: %d\n",K);
    for(i=0;i<K;i++)
    {
		scal = Scales[i];
        for(p=0;p<M;p++)
        {
            BPvector[p]     = 0;
            for(j=0;j<N;j++)
            {
                temp1[p*N+j]= BASIS[i*N+j]*PHI2[p*N+j]*beta[j];
                BPvector[p] = BPvector[p] + temp1[p*N+j];
            }
            BPvector[p]     = BPvector[p]/Scales[i];
        }
        //temp             	= (BPvector*Ui).^2; 
		//temp				= BPvector*SIGMA*BPvector'

        tempSum             = 0;
		for(p=0;p<M;p++)
        {
            //temp[p]      	= 0;
            //for(j=0;j<M;j++)                temp[p]     = temp[p] + BPvector[j]*SIGMA2[p*M+j];
            //temp[p]         = temp[p]*BPvector[p];
            //tempSum         = tempSum + temp[p];
			readPtr1  = &SIGMA2[p*M];
			temp[p]      	= F77_CALL(ddot)(&M, readPtr1, &inci,BPvector, &incj);
			
        }
		tempSum = F77_CALL(ddot)(&M, temp, &inci,BPvector, &incj);
				
        //S_in(i)           = beta'*BASIS2(:,i)/Scales(i)^2 -sum(temp);
        BBsquare            = 0;
        tempZE              = 0;
        //for(p=0;p<N;p++)
        //{
         //   BBsquare        = BBsquare + beta[p]*BasisCache[i*N+p];
        //    temp2[p]           = BASIS[i*N+p]*e[p];
        //    tempZE        	= tempZE + temp2[p];
        //}
		readPtr1 = &BasisCache[i*N];
		BBsquare = F77_CALL(ddot)(&N, beta, &inci,readPtr1, &incj);
		
		readPtr1 = &BASIS[i*N];
		tempZE = F77_CALL(ddot)(&N, readPtr1, &inci,e, &incj);
		
        S_in[i]             = BBsquare/(Scales[i]*Scales[i])-tempSum;
        Q_in[i]             = tempZE/Scales[i];
		//Print				= S_in[i];
		//Rprintf("S_in: %f\n",Print);
        S_out[i]            = S_in[i];
        Q_out[i]            = Q_in[i];
        //Interactions
        /*if(i!=(K-1))
        {
            for(L=(i+1);L<K;L++)
            {
                for(p=0;p<M;p++)
                {
                    BPvector[p]     = 0;
                    for(j=0;j<N;j++)	BPvector[p] = BPvector[p] + temp1[p*N+j]*BASIS[L*N+j];
                    BPvector[p]     = BPvector[p]/Scales[kk];
                }
                //temp             	= (BPvector*Ui).^2; 
                tempSum             = 0;
                for(p=0;p<M;p++)
                {
                    temp[p]      	= 0;
                    for(j=0;j<M;j++)	temp[p]     = temp[p] + BPvector[j]*SIGMA2[p*M+j];
                    temp[p]         = temp[p]*BPvector[p];
                    tempSum         = tempSum + temp[p];
                }
                //S_in(kk)       	= beta'*(BASIS2(:,i).*BASIS2(:,L))/Scales(kk)^2 -sum(temp);      
                BBsquare            = 0;
                tempZE              = 0;
                for(p=0;p<N;p++)
                {
                    BBsquare        = BBsquare + beta[p]*BasisCache[i*N+p]*BasisCache[L*N+p];
                    tempZE        	= tempZE + temp2[p]*BASIS[L*N+p];
                }
                S_in[kk]            = BBsquare/(Scales[kk]*Scales[kk])-tempSum;
                Q_in[kk]            = tempZE/Scales[kk];
                    
                S_out[kk]         	= S_in[kk];
                Q_out[kk]          	= Q_in[kk];
                kk                  = kk +1;
            }
        }*/
    }// Main loop ends
                    
    //S_out(Used),Q_out(Used)  
    //S_out(Used)				= (Alpha .* S_in(Used)) ./ (Alpha - S_in(Used));
    //_out(Used)       			= (Alpha .* Q_in(Used)) ./ (Alpha - S_in(Used));
    int N_used,index;
    N_used						= M-1;
    for(i=0;i<N_used;i++)
    {
        index					= Used[i] -1;
        // mexPrintf("index: \t %d\n",index);
        S_out[index]			= Alpha[i]*S_in[index]/(Alpha[i]-S_in[index]);
        Q_out[index]			= Alpha[i]*Q_in[index]/(Alpha[i]-S_in[index]);
    	/*Print				= S_out[index];
		Rprintf("S_out: %f\n",Print);
			Print				= Alpha[i];
			Rprintf("Alpha: %f\n",Print);*/
	}
	
	Free(PHI_Mu);
	Free(y);
	Free(e);
	Free(BPvector);
	Free(temp);
	Free(temp1);
	Free(temp2);

}


// #*******************************************************************
//[Mu2 beta SIGMA2] 			= fEBCatPostModeBmNeEN(PHI2,Targets,Alpha,Mu2);
void fEBCatPostModeBmNeEN(double * Mu2, double *beta,double *SIGMA2, double * H2, double *PHI2,
								double *Targets, double *Alpha,int N, int M)
{
	int i,j,k,L;
    //control parameters
    const double GRADIENT_MIN  	= 1e-6;
    double temp                 = pow(2.0,8);
    const double STEP_MIN       = 1/temp;
    const int   itsMax    		= 25;

    //Function fEBDataErrorBmNeEN [error,y]   =fEBDataErrorBmNeEN(PHI_Mu,Targets)
    //Calculate PHI_Mu: size: Nx1
    double *PHI_Mu, *y;
	PHI_Mu						= (double *) Calloc(N,double);
	int inci =1;
	int incj =1;
	int inc0 = 0;
	double *readPtr1;
	double scal;	
	
    for(i = 0;i<N;i++)
    {
        //PHI_Mu[i]				= 0;
        //for(j = 0;j<M;j++)		PHI_Mu[i]	= PHI_Mu[i] + PHI2[j*N+i]*Mu2[j];
		readPtr1 = &PHI2[i];
		PHI_Mu[i] = F77_CALL(ddot)(&M, readPtr1, &N,Mu2, &incj); 
    }
	double dataError			= 0;
	y							= (double *) Calloc(N,double);
	dataError					= fEBDataErrorBmNeEN(dataError,y,PHI_Mu,Targets,N);    
    //
    double regulariser			= 0;
    for(i = 1;i<M;i++)			regulariser		= regulariser + Alpha[i-1]*Mu2[i]*Mu2[i]/2;
    double newTotalError,g0,h0;
    newTotalError				= regulariser + dataError;
    double * errorLog,*e,*g2;
    errorLog					= (double *) Calloc(itsMax,double);
	e							= (double *) Calloc(N,double);
	g2							= (double *) Calloc(M,double);
    //
    g0							= 0;
    h0							= 0;
    //setup H
	int countGrad;
    double *DeltaMu, *Mu_new2;
    DeltaMu						= (double *) Calloc(M,double);
    Mu_new2						= (double *) Calloc(M,double);
    double step					= 1.0;

    //*****************************************************************************
    // Loops start here 
    for(i = 0;i<itsMax;i++)
    {
        errorLog[i]				= newTotalError;    //log error value
        //Gradient
        //g0						= 0;
        //h0						= 0;
        //for(j = 0;j<N;j++)
        //{
        //    e[j]				= Targets[j]-y[j];
        //    g0					= g0 + e[j];
        //    beta[j]				= y[j]*(1-y[j]);
        //    h0					= h0 + beta[j];
        //}
		
		F77_CALL(dcopy)(&N,Targets,&inci,e,&incj);
		scal = -1;
		F77_CALL(daxpy)(&N, &scal,y, &inci,e, &incj);
		scal = 1;
		g0 = F77_CALL(ddot)(&N, &scal, &inc0,e, &incj); 
		for(j = 0;j<N;j++)
        {
            beta[j]				= (y[j])*(1-y[j]); //exact 0 means separation: result in infinite values;
			if(beta[j]<1e-10) beta[j] = 1e-5;
			if(beta[j]>1e10) beta[j] = 1e5;
        }
		h0 = F77_CALL(ddot)(&N, &scal, &inc0,beta, &incj);

//		if(h0==0) h0 = 100;
        g2[0]					= g0;
        H2[0]					= h0;
        // first row and first column of H
        for(j=1;j<M;j++)
        {
            //g2[j]				= 0;
            //H2[j]				= 0;
            //for(k=0;k<N;k++)
            //{
            //    g2[j]			= g2[j] + PHI2[j*N+k]*e[k]; //Beta2 =Mu2
             //   H2[j]			= H2[j] +beta[k]*PHI2[j*N+k];
            //}
            //g2[j]				= g2[j]-Alpha[j-1]*Mu2[j];
            //H2[j*M]				= H2[j];
			readPtr1 = &PHI2[j*N];
            g2[j]				= F77_CALL(ddot)(&N, e, &inci,readPtr1, &incj);
            H2[j]				= F77_CALL(ddot)(&N, beta, &inci,readPtr1, &incj);

            g2[j]				= g2[j]-Alpha[j-1]*Mu2[j];
            H2[j*M]				= H2[j];
        }
        //Hessian
        //h01
        //H;
        for(j=1;j<M;j++)
        {
            for(k = 1;k<M;k++)
            {
                H2[k*M+j]       = 0;
                for(L = 0;L<N;L++)		H2[k*M+j]   = H2[k*M+j] + PHI2[j*N+L]*beta[L]*PHI2[k*N+L];
                if(j==k)                H2[k*M+j]   = H2[k*M+j] + Alpha[k-1];
            }
        }
		//save a copy of H
		//for(j=0;j<M;j++)
		//{
		//	for (k=0;k<M;k++)	SIGMA2[k*M+j] = H2[k*M+j];
		//}
		j = M*M;
		F77_CALL(dcopy)(&j,H2,&inci,SIGMA2,&incj);
		
		MatrixInverseBmNeEN(SIGMA2,M);				//inverse of H2 is needed for wald score
//		if(SIGMA2[0]>1e2) SIGMA2[0] = 0.01;
		//if(i==itsMax-1) 
		//{
		//	printMat(SIGMA2,M,M);		
		//	printMat(H2,M,M);

		//	printMat(beta,M,1);
		//}
		
		//Rprintf("INSIDE fEBPostMode: \n");
		//Rprintf("\t SIGMA2: %f, %f \t H2: %f, %f\n\t %f, %f, \t %f, %f\n",SIGMA2[0],SIGMA2[1], H2[0], H2[1],SIGMA2[2],SIGMA2[3],H2[2],H2[3]);
        countGrad               = 0;
        for(j = 1;j<M;j++) 
		{
            if(fabs(g2[j])<GRADIENT_MIN)	countGrad++;
		}
        if(countGrad==M-1)						//end loop
        {
            for (k = 0;k<M;k++)
            {
               // DeltaMu[k]      = 0;
               // for(L = 0;L<M;L++)	DeltaMu[k]  = DeltaMu[k] + SIGMA2[L*M+k]*g2[L];
               readPtr1 =&SIGMA2[k];
				DeltaMu[k]      = F77_CALL(ddot)(&M, g2, &inci,readPtr1, &M);
			   
            }
         
            break;
        }
        //Comput Full Newton step H^(-1)*g;            
         //sigma*g
		for (k = 0;k<M;k++)
		{
                //DeltaMu[k]      = 0;
                //for(L=0;L<M;L++)	DeltaMu[k]  = DeltaMu[k] + SIGMA2[L*M+k]*g2[L];
				readPtr1 =&SIGMA2[k];
				DeltaMu[k]      = F77_CALL(ddot)(&M, g2, &inci,readPtr1, &M);	
				
		}
        //
        step					= 1;
        while (step>STEP_MIN)
        {
            //for(j=0;j<M;j++)	Mu_new2[j]      = Mu2[j] + step*DeltaMu[j];
            F77_CALL(dcopy)(&M,Mu2,&inci,Mu_new2,&incj);
			F77_CALL(daxpy)(&M,&step,DeltaMu,&inci,Mu_new2,&inci);
			
			//
            for(j=0;j<N;j++)
            {
                //PHI_Mu[j]       = 0;
                //for(k=0;k<M;k++)	PHI_Mu[j]   = PHI_Mu[j] + PHI2[k*N+j]*Mu_new2[k];
				readPtr1 = &PHI2[j];
				PHI_Mu[j] = F77_CALL(ddot)(&M, readPtr1, &N,Mu_new2, &incj); 
				
            }
            //Function 2:   fEBDataErrorBmNeEN
            dataError			= fEBDataErrorBmNeEN(dataError,y,PHI_Mu,Targets,N);
            //regulariser update            
            regulariser			= 0;
            for(j=1;j<M;j++)	regulariser  	= regulariser + Alpha[j-1]*Mu_new2[j]*Mu_new2[j]/2;
            newTotalError       = dataError + regulariser;
            if(newTotalError>=errorLog[i])      step		= step/2;
            else
            {
                //for(j=0;j<M;j++)				Mu2[j]		= Mu_new2[j];
				F77_CALL(dcopy)(&M,Mu_new2,&inci,Mu2,&incj);
                step            =0;
            }
        }//end of while
        if(step==1) break;
    }
	
	Free(PHI_Mu);
	Free(y);	
	Free(errorLog);
	Free(e);	
	Free(g2);
	Free(DeltaMu);	
	Free(Mu_new2);
	
}

 //************************************************************************************************
double fEBDataErrorBmNeEN(double dataError,double *y,double *PHI_Mu,double *Targets,int N)
{
	int i;
    dataError					= 0;		// scalar
    // Call function fEBSigmoidBmNeEN;    
	fEBSigmoidBmNeEN(y,PHI_Mu,N);
    for (i = 0;i<N;i++)
    {
        if(y[i]!=0)		dataError   = dataError - Targets[i]*log(y[i]);
        if (y[i]!=1)	dataError   = dataError - (1-Targets[i])*log(1-y[i]);
    }   
	return dataError;
}

void fEBSigmoidBmNeEN(double * y, double * PHI_Mu,int N) // from left to right 
{
	int i;
    for (i=0;i<N;i++)	y [i]	= 1/(1+exp(-PHI_Mu[i]));

}

//****************************************************************
//** Matrix Inverse by call lapack library of cholesky decomposition and linear solver *
void MatrixInverseBmNeEN(double * a,int N)
{
	const char uplo			= 'U';
	int info				= 0;
	
	//*************************Call LAPACK library ************************
	F77_CALL(dpotrf)(&uplo, &N,a, &N,&info);
	if(info!=0) 	
	{
		Rprintf("Call 1st function. dpotrf error, Ill-conditioned Hessian!\n");
		return;
	}
	F77_CALL(dpotri)(&uplo,&N,a,&N,&info);
	if(info!=0) 	
	{
		Rprintf("Call 2nd function dpotri error!\n");
		return;
	}
	//a is upper triangular
	int i,j;
	for (i=1;i<N;i++)
	{
		for(j=0;j<i;j++)	a[j*N+i]	= a[i*N+j];
	}
}

//***************************************************************************************************
void fEBDeltaMLBmNeEN(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, double * S_out, double * Q_out, double *Alpha,
				double *a_lambda, double *b_Alpha, int N_used, int N_unused,
				double *deltaLogMarginal, int *nu,int N)
{
    //basis dimension
    int M_full,i,index;
anyToDelete[0] = 0;
    M_full								= N_used + N_unused;
    //
    //Parameters setup      //Action is return as int
    const int	ACTION_REESTIMATE		= 0;			
	const int 	ACTION_ADD				= 1;
	const int	ACTION_DELETE			= -1;
    //const double    ACTION_OUT		= 0;      //for Unused index and stay unused. not needed: mxArray initialized as zeros
    const double    logL_inf			= 0.0;
    double alpha, beta, gamma,delta,root2,logL,oldAlpha,lambda1,lambda2; //lambda
	lambda1 							= *a_lambda * (b_Alpha[0]);
	lambda2 							= *a_lambda * (1 - b_Alpha[0]);
    //int   anyToDelete					= 0;
    int   anyToAdd						= 0;
    //const int CNPriorityAddition		= 0;
    //const int CNPriorityDeletion		= 1;
	//Rprintf("Inside fEBDeltaMLBmNeEN: a: %f, b %f \n",a, b);
	int CNPriorityAddition		= 0;
	int CNPriorityDeletion		= 0;		
	if(N_used<10)
	{
		CNPriorityAddition		= 1;
		CNPriorityDeletion		= 0;
		//Rprintf("Inside fEBDeltaMLBmNeEN: a: %f, b %f \n",a, b);
	}
	if(N_used>100 || N_used >=N)
	{
		CNPriorityAddition		= 0;
		CNPriorityDeletion		= 1;
	}	
	//Block Coordinate Ascent Algorithm Parameters;
	for(i=0;i<M_full;i++) Action[i] = -10;
	double deltaMLmax = 0;
	int 	max = 0;
	
    for(i=0;i<N_used;i++)
    {
        //anyToDelete				= false;
        index						= Used[i] -1;
		DeltaML[index]				= 0;
        //    mexPrintf("N_used: \t %d,\t%d\n",N_used,index);
        alpha						= S_out[index] - Q_out[index]*Q_out[index] + 2*lambda1 + lambda2;
        beta						= (S_out[index]+lambda2)*(S_out[index] + 4 *lambda1 + lambda2);
        gamma						= 2*lambda1*(S_out[index] + lambda2)*(S_out[index]+lambda2);
        delta						= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha<0 && delta>0)
        {
            root2					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root2/(root2+S_out[index] + lambda2))+pow(Q_out[index],2)/(root2+S_out[index] + lambda2))*0.5 -lambda1/root2;
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root2 + lambda2;
                Action[index]       = ACTION_REESTIMATE;
                //
                oldAlpha            = Alpha[i]-lambda2;				
                DeltaML[index]  	= 0.5*(log(root2*(oldAlpha + S_out[index]+lambda2)/(oldAlpha*(root2 + S_out[index] + lambda2))) +
                                    Q_out[index]*Q_out[index]*(1/(root2 + S_out[index] +lambda2) - 1/(oldAlpha+S_out[index]+ lambda2))) -
                                    lambda1*(1/root2 - 1/oldAlpha);
            }
        }
        //case 2
        //case 3

        //DELETE
        else if (N_used>1)
        {
            anyToDelete[0]      = 1;
            Action[index]       = ACTION_DELETE;
            oldAlpha            = Alpha[i] - lambda2;
            logL                = (log(oldAlpha/(oldAlpha+S_out[index] + lambda2)) + pow(Q_out[index],2)/(oldAlpha + S_out[index] + lambda2))*0.5 - lambda1/oldAlpha;
            DeltaML[index]      = - logL;
        }
		
		//Block Coordinate Ascent Algorithm Parameters;
		if(DeltaML[index]>deltaMLmax)
		{
				max			= index;
				deltaMLmax 	= DeltaML[index];
		}
		
    }
    //ADDITION
    for(i=0;i<N_unused;i++)
    {
        index					= Unused[i] -1;
		DeltaML[index]			= 0;
		alpha						= S_out[index] - Q_out[index]*Q_out[index] + 2*lambda1 + lambda2;
        beta						= (S_out[index]+ lambda2)*(S_out[index] +lambda2 + 4 *lambda1);
        gamma						= 2*lambda1*(S_out[index]+lambda2)*(S_out[index] + lambda2);
        delta						= beta*beta - 4*alpha*gamma;
        //case1
        if(alpha<0 && delta>0)
        {
            root2					= (- beta-sqrt(delta))/(2*alpha);
            logL					= (log(root2/(root2 + S_out[index] + lambda2)) + pow(Q_out[index],2)/(root2 + S_out[index] + lambda2))*0.5 - lambda1/root2;
            if (logL > logL_inf)
            {
                AlphaRoot[index]    = root2 + lambda2;
                Action[index]       = ACTION_ADD;
                //
                DeltaML[index]  	= (log(root2/(root2 + S_out[index] + lambda2)) + pow(Q_out[index],2)/(root2 + S_out[index] +lambda2))*0.5 - lambda1/root2;
            }
        }
        //case 2
        //case 3  
		//Block Coordinate Ascent Algorithm Parameters;
		if(DeltaML[index]>deltaMLmax)
		{
				max			= index;
				deltaMLmax 	= DeltaML[index];
		} 		
    }
    
	//Block Coordinate Ascent Algorithm Parameters; 
	// only one action for a block
	//int selectedAction = Action[max];
	
	/*
	if(i_iter <=(1*N_used)) // re-estimate
	{
		
		for(i=0;i<M_full;i++)
		{
			if(Action[i] !=ACTION_REESTIMATE) DeltaML[i] = 0;
		}
	}else
	{
		for(i=0;i<M_full;i++)
		{
			if(Action[i] !=selectedAction) DeltaML[i] = 0;
		}
	}
	*/
	//	for(i=0;i<M_full;i++)
	//	{
	//		if(Action[i] !=selectedAction) DeltaML[i] = 0;
	//	}
	
   if((anyToAdd==1 && CNPriorityAddition==1) || (anyToDelete[0]==1 && CNPriorityDeletion==1))
    {
        for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_REESTIMATE)											DeltaML[i]     = 0;
			else if (Action[i] == ACTION_DELETE)
            {
                    if(anyToAdd==1 && CNPriorityAddition==1 && CNPriorityDeletion!=1)    DeltaML[i]     = 0;
            }else if (Action[i] == ACTION_ADD)
            {
                    if(anyToDelete[0] ==1 && CNPriorityDeletion==1 && CNPriorityAddition!=1) DeltaML[i] = 0;
            }
        }
		
		deltaMLmax = 0;
		max = 0;		
		for(i=0;i<M_full;i++)
		{	
			if(DeltaML[i]>deltaMLmax)
			{
					max			= i;
					deltaMLmax 	= DeltaML[i];
			} 	
		}	
		
    }  

deltaLogMarginal[0] = deltaMLmax;
nu[0] 	= max; 
}




