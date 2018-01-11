//version 4.2 note: initialize the least correlated variable to avoid collinear variables; delete the initial effect at the end of 1st iteration.
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
//G: gauss; m:main; Ne: normal exp
void LinearSolverGmNeEN(double * a, double *logout, int N,int M,double *output);
void fEBInitializationGmNeEN(double *Alpha, double * PHI, int *Used, int *Unused, 
				double *BASIS, double *Targets, double *Scales, int * initial, 
				int n, int *m, int kdim, double *beta);
				
void fEBInitializationGmNeEN2(double *Alpha, double *PHI, int *Used, int *Unused, double *BASIS, 
			double *Targets, double *Scales, int * initial, int N, int *m, int K, double *beta);
			
void MatrixInverseGmNeEN(double * a,int N);
void CacheBPGmNeEN(double **BASIS_PHI, double *BASIS_Targets, double *BASIS, double *PHI,
				double *Targets,double *scales,int N,int K,int M, int M_full);

void fEBLinearFullStatGmNeEN(double *beta,double * SIGMA,double *H, double *S_in, double * Q_in, double * S_out, 
				double * Q_out,  double *BASIS, double * Scales, double *PHI, double **BASIS_PHI,
				double *BASIS_Targets, double * Targets, int * Used, double *Alpha, double * Mu, 
				 double *Gamma,int *n, int *m, int* kdim, int *iteration,int *i_iter);

void fEBDeltaMLGmNeEN(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, double * S_out, double * Q_out, double *Alpha,
				double *a_lmabda, double *b_Alpha,
				int m, int mBar, int N,
				double *deltaLogMarginal,int *nu,double residual,double varY, int iter, int i_iter);

void LinearFastEmpBayesGmNeEN(int *Used, double *Mu, double *SIGMA, double *H, double *Alpha, double *PHI,
				double *BASIS, double * Targets, double *Scales, double *a_lambda,double *b_Alpha,
				int *iteration, int *n, int *kdim, int *m, int basisMax, double *b, double *beta,double * C_inv,
				int verbose,double residual,double varY);

int ActionAddGmNeEN(double **BASIS_PHI, double* BASIS, double*scales, double*PHI, double*Phi,
			double *beta, double* Alpha, double newAlpha, double*SIGMA, double*Mu, double*S_in,
			double*Q_in, int nu, double*SIGMANEW, int M_full, int N, int K, int M);


int ActionDelGmNeEN(double*PHI, double*Alpha, double*SIGMA, double*SIGMANEW, double**BASIS_PHI,
				double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full);

double varTargetsGmNeEN(double* Target,int N);

void FinalUpdateGmNeEN(double *PHI,double * H,double*SIGMA,double *Targets,double *Mu,double *Alpha,
				double *beta,int N, int M);

//
void elasticNetLinearNeMainEff(double *BASIS, double *y, double *a_lambda, double *b_Alpha,
				double *Beta, 
				double *wald, double *intercept, int *n, int *kdim, int *verb,double *residual)
{
	int N					= *n;
	int K					= *kdim;
	int verbose 				= *verb;
	//int M_full				= (K+1)*K/2;
	int M_full = K;
	const int iter_max		= 100;
	const double err_max	= 1e-8;
	// set a limit for number of basis
	
	int basisMax			= 1e7/M_full;
	if (basisMax>M_full)	basisMax = M_full;
	if(verbose>0) Rprintf("basisMax: %d",basisMax);
	if(verbose >1) Rprintf("start EB-elasticNet with alpha: %f, lambda: %f\n",*b_Alpha, *a_lambda);
	double vk				= 1e-30;
	double vk0				= 1e-30;
	double temp				= 0;
	int i;
	double *Scales			= (double * ) Calloc(M_full, double);
	//lapack
	int inci 				=1;
	int incj 				=1;
	double *readPtr1, *readPtr2;
	int inc0 				= 0;
	double a_blas 			= 1;
	double b_blas 			= 1;
	double zero_blas 		= 0;
	//lapack end
	
	for (i					=0;i<K;i++)
	{
		Beta[i]				= i + 1;
		Beta[M_full + i]	= i + 1;
		//Beta[M_full*2 + i]	= 0;
		//Beta[M_full*3 + i]	= 0;
		temp				= 0;
		//for(l=0;l<N;l++)	temp		= temp + BASIS[i*N + l]*BASIS[i*N + l];
		readPtr1  			= &BASIS[i*N];
		temp  				= F77_CALL(ddot)(&N, readPtr1, &inci,readPtr1, &incj);
		if(temp ==0) temp	= 1;
		Scales[i]			=sqrt(temp);
	}
	readPtr1 				= &Beta[K*2];
	F77_CALL(dcopy)(&K,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
	readPtr1 				= &Beta[K*3];
	F77_CALL(dcopy)(&K,&zero_blas,&inc0,readPtr1,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
		
	//PartII kk
	//kk						= K;					//index starts at 0;	
	/*for (i					=0;i<(K-1);i++)
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
			kk					= kk + 1;
		}
	}*/

	//
	int iter				= 0;
	double err				= 1000;
	double *Mu, *SIGMA, *H, *Alpha, *PHI,*Targets,*C_inv;
	int * Used,*iteration, *m;
	
	Used					= (int* ) Calloc(basisMax, int);
	Mu						= (double * ) Calloc(basisMax, double);							  
	SIGMA					= (double * ) Calloc(basisMax*basisMax, double);
	H						= (double * ) Calloc(basisMax*basisMax, double);
	Alpha					= (double * ) Calloc(basisMax, double);
	PHI						= (double * ) Calloc(N*basisMax, double);
	Targets					= (double * ) Calloc(N, double);
	iteration				= (int* ) Calloc(1, int);
	m						= (int* ) Calloc(1, int);
	C_inv					= (double * ) Calloc(N*N, double);
	if(verbose>0) Rprintf("outer loop starts");
	m[0]			= 1;
	int M					= m[0];	
	//Fixed Effect
	double b				= 0;
	//for(i=0;i<N;i++)	b	= b + y[i];
	F77_CALL(daxpy)(&N, &a_blas,y, &inci,&b, &inc0); //daxpy(n, a, x, incx, y, incy) y := a*x + y

	b						= b/N;
	double beta;
	double *Csum			= (double *) Calloc(N,double);
	double Cinv,Cinvy;
	
	double varT = varTargetsGmNeEN(y,N);
	double epsilon = varT*0.01;
	double residvar = 1e10;
	while (iter<iter_max && err>err_max && residvar>=epsilon)
	{
		iter				= iter + 1;
		
		vk0					= vk;
		iteration[0]		= iter;
		//for(i = 0;i<N;i++)	Targets[i]	= y[i] - b;
		b_blas = -b;
		F77_CALL(dcopy)(&N,&b_blas,&inc0,Targets,&inci);  //dcopy(n, x, incx, y, incy) ---> y = x
		F77_CALL(daxpy)(&N, &a_blas,y, &inci,Targets, &incj); //daxpy(n, a, x, incx, y, incy) y := a*x + y
		
//
		LinearFastEmpBayesGmNeEN(Used, Mu, SIGMA, H, Alpha,PHI,	BASIS, Targets,Scales, a_lambda,b_Alpha,
						iteration, n, kdim, m,basisMax,&b,&beta,C_inv,verbose,residvar,varT); //residvar is current residual variance.

		//b				=1/(C_inv.Transpose()*x*x)*(C_inv.Transpose()*x*phenotype);

		 for(i=0;i<N;i++)
		 {
			 Csum[i]	= 0;
			 //for(j=0;j<N;j++)	Csum[i]		= Csum[i] + C_inv[i*N+j];
			 readPtr1 		= &Csum[i];
			 readPtr2 		= &C_inv[i*N];
			F77_CALL(daxpy)(&N, &a_blas,readPtr2, &inci,readPtr1, &inc0); //daxpy(n, a, x, incx, y, incy) y := a*x + y

		 }
		 Cinv = 0;
		 //for(i=0;i<N;i++)	Cinv = Cinv + Csum[i];
		F77_CALL(daxpy)(&N, &a_blas,Csum, &inci,&Cinv, &inc0); //daxpy(n, a, x, incx, y, incy) y := a*x + y
	
		 Cinvy = 0;
		 //for(i=0;i<N;i++)	Cinvy = Cinvy + Csum[i]*y[i];
		 Cinvy 				= F77_CALL(ddot)(&N, Csum, &inci,y, &incj);
		 b		= Cinvy/(Cinv+ 1e-10);
		vk					= 0;
		//for(i=0;i<m[0];i++)	vk = vk + Alpha[i];
		M 					= m[0];
		F77_CALL(daxpy)(&M, &a_blas,Alpha, &inci,&vk, &inc0); //daxpy(n, a, x, incx, y, incy) y := a*x + y
	
		err					= fabs(vk - vk0)/m[0];
		residvar = 1/(beta + 1e-10);
		if(verbose >2) Rprintf("Iteration number: %d, err: %f;\t mu: %f\tsigma0:%f.\n",iter,err,b,residvar);
	}//outerloop

	// wald score
	M					= m[0];	
	double *tempW			= (double * ) Calloc(M,double);

	wald[0]					= 0;
	int index = 0;
	if(verbose >1) Rprintf("EBEN Finished, number of basis: %d\n",M);
	for(i=0;i<M;i++)
    {

        tempW[i]      		= 0;
        //for(j=0;j<M;j++)    tempW[i]     = tempW[i] + Mu[j]*H[i*M+j];       
		readPtr1 			= &H[i*M];		
        tempW[i] 			= F77_CALL(ddot)(&M, Mu, &inci,readPtr1, &incj);     
        //wald[0]				= wald[0]	 +tempW[i]*Mu[i];
	}
	wald[0] 				= F77_CALL(ddot)(&M, tempW, &inci,Mu, &incj);	
	for(i=0;i<M;i++)
	{
    // blup collection
		//Rprintf("test Used: %d\n",Used[i-1]);
		index				= Used[i] - 1;
		Beta[M_full*2 + index]	= Mu[i]/Scales[index];
		Beta[M_full*3 + index]  = SIGMA[i*M + i]/(Scales[index]*Scales[index]);
	}
	//

	intercept[0]	= b;
	residual[0] 	= 1/(beta + 1e-10);
	//Rprintf("fEB computation compelete!\n");
	Free(Scales);
	Free(Used);	
	Free(Mu);
	Free(SIGMA);	
	Free(H);
	Free(Alpha);	
	Free(PHI);
	Free(Targets);	
	Free(iteration);	
	Free(m);
	Free(C_inv);	
	Free(tempW);
	Free(Csum);
}


/************** outputs are passed by COPY in R, cann't dynamic realloc memory **************************/
/************** Not a problem in C */
// function [Used,Mu2,SIGMA2,H2,Alpha,PHI2]=fEBBinaryMex(BASIS,Targets,PHI2,Used,Alpha,Scales,a,b,Mu2,iter)
void LinearFastEmpBayesGmNeEN(int *Used, double *Mu, double *SIGMA, double *H, double *Alpha, double *PHI,
				double *BASIS, double * Targets, double *Scales, double *a_lambda, double *b_Alpha,
				int *iteration, int *n, int *kdim, int *m,int basisMax,double *b,double *beta,double * C_inv,
				int verbose,double residual,double varY)
{
    //basis dimension
   int N,K,M_full,N_unused,M,i,j,iter;
   	N					= *n;			// row number
    K					= *kdim;		// column number
    //M_full				= K*(K+1)/2;
	M_full = K;
	//kk					= K;

	int *Unused				= (int *) Calloc(M_full,int);
    iter				= *iteration;
	//Rprintf("Iteration number: %d\n",iter);
    const int	ACTION_REESTIMATE       = 0;			
	const int	ACTION_ADD          	= 1;
	const int 	ACTION_DELETE        	= -1;
    const int   ACTION_TERMINATE        = 10;    
    
	//
	const int		CNBetaUpdateStart	=10;
	const double	BetaMaxFactor		=1e6;
	const double	MinDeltaLogBeta		=1e-6;
	
	//Block Coordinate Ascent Algorithm Parameters;
	const double 	nAdd 				= 0.9; //only deltaML >nAdd will be considered as candidate
	double nAdd_real 					= nAdd;
	const double 	MLdelta 			= 0.001; //only deltaML > MLdelta will be selected for action.
	
    //[Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic) 
    int *IniLogic;
	IniLogic				= (int*) Calloc(1,int);
    if (iter<=1)    
    {
        IniLogic[0]     = 0;
        m[0]            = 1;
		M				= m[0];
		N_unused		= M_full -1;

    }else
    {
		IniLogic[0]    = 1;
        M				= *m;          //Used + 1
		N_unused		= M_full - M;
    }
    //
	//lapack
	int inci 				=1;
	int incj 				=1;
	double *readPtr1;//, *readPtr2;
	int inc0 				= 0;
	double a_blas 			= 1;
	double b_blas 			= 1;
	double c_blas 			= 1;
	int MM;
	char transa 			= 'N';
	char transb 			= 'N';
	int lda,ldb,ldc,ldk;
	//lapack end
	
	//Rprintf("N_used is: %d; N_unused:%d, M: %d,sample size: %d \n",N_used,N_unused,M,N);

	fEBInitializationGmNeEN(Alpha, PHI, Used, Unused, BASIS, Targets, Scales, IniLogic, N, m, K,beta);
	M = m[0];
	N_unused = M_full - M;
	
	//******************************************************	
	int initial = Used[0];
	int iniRemoved = 1;
	if (iter<=1)    
    {
		initial = Used[0];
		iniRemoved = 0;
    }
	//Rprintf("\t Initialized basis %d, Alpha: %f, \n", Used[0],Alpha[0]);
	//Rprintf("\t beta: %f\n",beta[0]);
	//Rprintf("\t m: %d\n",M);
	//for(i=0;i<N;i++) Rprintf("\tPHI: \t %f\n",PHI[i]);
	//for(i=0;i<10;i++) Rprintf("PHI2: %f \t  %f; BASIS: %f\n",PHI2[i],PHI2[N+i],BASIS[181*N+i]/Scales[181]);
    //CACHE MATRIX
	double *BASIS_Targets,**BASIS_PHI;
	BASIS_Targets		= (double *) Calloc(M_full,double);
	BASIS_PHI			= (double **) Calloc(basisMax,double);
	for(i=0;i<M;i++)
	{
		BASIS_PHI[i] 	= (double *) Calloc(M_full,double);
	}
	CacheBPGmNeEN(BASIS_PHI, BASIS_Targets, BASIS, PHI,	Targets,Scales,N,K,M,M_full);
	
	double *S_in, *Q_in, *S_out, *Q_out,*gamma;
	S_in				= (double *) Calloc(M_full,double);
	Q_in				= (double *) Calloc(M_full,double);
	S_out				= (double *) Calloc(M_full,double);
	Q_out				= (double *) Calloc(M_full,double);
	gamma				= (double *) Calloc(basisMax,double);
    //[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,Intercept] ...
    //                   	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE)
	int i_iter = 0;
	fEBLinearFullStatGmNeEN(beta,SIGMA, H, S_in, Q_in, S_out,Q_out,  BASIS, Scales, 
			PHI, BASIS_PHI,BASIS_Targets, Targets, Used, Alpha, Mu, 
				 gamma, n, m, kdim, iteration,&i_iter);
//Rprintf("\t gamma: %f\n",gamma[0]);
//Rprintf("\t 182th S_in: %f, Q_in: %f, S_out: %f, Q_out: %f\n",S_in[181],Q_in[181],S_out[181],Q_out[181]);
//for(i=0;i<100;i++) Rprintf("BASIS_T: %f\n",BASIS_Targets[i]);
    //              For: [DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLGmNeEN(Used,Unused,S_out,Q_out,Alpha,a,b);
   double *DeltaML, *AlphaRoot,deltaLogMarginal,*phi,newAlpha,oldAlpha;
    double deltaInv,kappa,Mujj;
    //
	int *Action, *anyToDelete;
	int selectedAction = -10;
	anyToDelete			= (int*) Calloc(1,int);
	DeltaML				=	(double *) Calloc(M_full,double);
	AlphaRoot			=	(double *) Calloc(M_full,double);
	Action				= (int *) Calloc(M_full,int);
  	phi					= (double *) Calloc(N,double);

    int nu,jj,index;
    jj					= -1;
    int anyWorthwhileAction,UPDATE_REQUIRED;
    // mexPrintf("cols:%d \n",M);
  	//

    int LAST_ITERATION  = 0;
	//Gauss update
	double *PHI_Mu,*e;
	PHI_Mu				= (double*) Calloc(N,double);
	e					= (double*) Calloc(N,double);
	double betaZ1;
	double deltaLogBeta;
	double ee;
	double varT;
	
	
	double temp;			// for action_reestimate
	double * SIGMANEW	= (double * ) Calloc(basisMax*basisMax, double);

//Block Coordinate Ascent Algorithm Parameters;
	double cutoff = 0;
	int *nuAction;
	nuAction = (int*) Calloc(M_full,int);
	int nUpdate, iU;// nTerminate;
	
	if(verbose>3 && iter ==1)
	{
		Rprintf("check point 3: before loop,initial number of basis:%d\t Basis initialized: ",M);
		for(i=0;i<M;i++) Rprintf("%dth:%d\t Alpha: %f \t w: %f\t",i,Used[i],Alpha[i],Mu[i]);
		Rprintf("beta:%f\n",beta[0]);

	}

int itMax = 100;
if(iter==1) itMax = 10;
   while(LAST_ITERATION!=1)
    {
        i_iter						= i_iter + 1;
		if(verbose >4) Rprintf("\t inner loop %d; number of basis: %d \t actionStatus: %d \tiniRemoved: %d\n",i_iter, M,selectedAction,iniRemoved);
        //M							= N_used + 1;
     	//N_unused					= M_full - N_used;

        //[DeltaML,Action,AlphaRoot,anyToDelete]     = fEBDeltaMLGmNeEN(Used,Unused,S_out,Q_out,Alpha,a,b);
		fEBDeltaMLGmNeEN(DeltaML, Action, AlphaRoot,anyToDelete,Used, Unused, S_out, Q_out, Alpha,
				a_lambda, b_Alpha, M, N_unused, N,
				&deltaLogMarginal,&nu,residual,varY,iter,i_iter);
		//
        //deltaLogMarginal			= MLdelta;
        //nu							= -1;
        //for(i=0;i<M_full;i++)
        //{
        //    if(DeltaML[i]>deltaLogMarginal)
         //   {
        //        deltaLogMarginal    = DeltaML[i];
         //       nu                  = i;
         //   }
        //}
        //selectedAction          = Action[nu];
		//*******************************************************************************************
		if(selectedAction          == ACTION_TERMINATE && iniRemoved ==0 &&M>1)
		{
			nu = -1;
			if(verbose >4) Rprintf("\t\t************************ set to Remove the initial basis **************** \n");
		}
		
		if(nu==-1 && iniRemoved ==1)
		{
			anyWorthwhileAction     = 0;
			selectedAction          = ACTION_TERMINATE;
			if(verbose >4) Rprintf("\t\t********Terminate inner loop due to no worthwhile action basis\n");
		}else if(nu==-1 &&iniRemoved ==0 && M>1) //remove initialized basis
		{
			if(verbose >4) Rprintf("\t\t************************** Removing the initial basis **************** \n");
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
if(nUpdate==0)
{
	anyWorthwhileAction=0;	
if(verbose >4) Rprintf("\t\t********Terminate inner loop due to no worthwhile action basis: nu: %d \tdeltaLogMarginal: %f\tcutoff:%f\tAction:%d\tdeltaML:%f\n",nu,deltaLogMarginal,cutoff,Action[nu],DeltaML[nu]);
}
		}
if(anyWorthwhileAction==0)  selectedAction = ACTION_TERMINATE;		
		if(anyWorthwhileAction	==1)
		{
		
			if(verbose >5) Rprintf("\t deltaML cutoff: %f ; number to be updated: %d \n",cutoff, nUpdate);
			

				
			for(iU=0;iU<nUpdate;iU++) //update each basis
			{
				nu 					= nuAction[iU];
								
				selectedAction          = Action[nu];
				newAlpha                = AlphaRoot[nu];
				deltaLogMarginal 	= DeltaML[nu];

				if(verbose>5) Rprintf("\t ActionOn nu= %d, deltaML: %f, selectedAction: %d, Alpha: %f\n",nu+1, DeltaML[nu],selectedAction,AlphaRoot[nu]);
				if(selectedAction==ACTION_REESTIMATE || selectedAction==ACTION_DELETE)
				{
					index                   = nu + 1; 
					for(i=0;i<M;i++)
					{
						if (Used[i]==index)	
						{
								jj  = i;
								break;
						}
					}
				}

				//kk                          = K;                          
				for(i=0;i<K;i++)
				{
					if (i==nu)
					{
						//for(L=0;L<N;L++)    phi[L]  = BASIS[i*N+L]/Scales[i];
						readPtr1 		= &BASIS[i*N];
						F77_CALL(dcopy)(&N,readPtr1,&inci,phi,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
						b_blas 			= 1/Scales[i];
						F77_CALL(dscal)(&N,&b_blas,phi,&inci); 		//dscal(n, a, x, incx) x = a*x
						//break;
					}
					/*else if(i<(K-1))
					{
						for(j=(i+1);j<K;j++)
						{
							if(kk==nu)
							{
								for(L=0;L<N;L++)    phi[L] =BASIS[i*N+L]*BASIS[j*N+L]/Scales[kk];
								//break;
							}
							kk              = kk + 1;
						}
					}*/
				}
				//Rprintf("N: %d, kk: %d,i: %d, j: %d, K: %d\n",N,kk,i,j,K);
				//for(i=0;i<N;i++) Rprintf("phi: %f\n",phi[i]);
				//newAlpha                    = AlphaRoot[nu];
				//Rprintf("\tNew alpha: %f\n",newAlpha);
				
				if(selectedAction==ACTION_REESTIMATE)
				{
					if (fabs(log(newAlpha)-log(Alpha[jj]))<=1e-3 && anyToDelete[0] ==0)
					{	
						//cout<<"Terminated by small deltaAlpha on basis:"<<nu<<endl;
						selectedAction		= ACTION_TERMINATE;
						if(verbose >4) Rprintf("\t\t********Terminate inner loop due to too small changes in re-estimation while no other actions\n");
					}
				}
				//Rprintf("\t selectedAction: %d\n",selectedAction);
				//
				UPDATE_REQUIRED				= 0;
				if(selectedAction==ACTION_REESTIMATE)
				{
					if(verbose>4) Rprintf("\t\t Action: Reestimate : %d \t deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
								nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					oldAlpha				= Alpha[jj];
					Alpha[jj]				= newAlpha;

					deltaInv				= 1.0/(newAlpha-oldAlpha);
					kappa					= 1.0/(SIGMA[jj*M+jj] + deltaInv);
					Mujj					= Mu[jj];
					//for(i=0;i<M;i++)		Mu[i]    = Mu[i] - Mujj *kappa * SIGMA[jj*M+i];
					readPtr1 				= &SIGMA[jj*M];
					b_blas 					= -Mujj * kappa;
					F77_CALL(daxpy)(&M, &b_blas,readPtr1, &inci,Mu, &incj); //daxpy(n, a, x, incx, y, incy) y := a*x + y

					//	SIGMANEW	=SIGMA-sjj;
					for(i=0;i<M;i++)
					{
						for(j=0;j<M;j++)	SIGMANEW[j*M + i] = SIGMA[j*M + i] - kappa * SIGMA[jj*M+i]*SIGMA[jj*M+j];
					}
					
					//S_in		= S_in + kappa*(beta*BASIS_PHI * s_j).^2;
					//Q_in		= Q_in - beta*BASIS_PHI*(deltaMu);           
					
					for(i=0;i<M_full;i++)
					{
						temp	= 0;
						//for(j=0;j<M;j++) temp = temp + BASIS_PHI[j*M_full + i]*SIGMA[jj*M + j];
						//readPtr1 			= &BASIS_PHI[i];
						//readPtr2 			= &SIGMA[jj*M];
						//temp = F77_CALL(ddot)(&M, readPtr1, &M_full,readPtr2, &inci);
						for(j=0;j<M;j++) temp = temp + BASIS_PHI[j][i]*SIGMA[jj*M + j];
						S_in[i]				= S_in[i] +  pow(beta[0]*temp,2)*kappa;
						Q_in[i]				= Q_in[i] +  beta[0]*Mujj *kappa*temp;
					}
					//
					//for(i=0;i<M_full;i++)
					//{
					//	temp	= 0;
					//	for(j=0;j<M;j++) temp = temp + BASIS_PHI[j*M_full + i]* SIGMA[jj*M+j];
					//	Q_in[i]				= Q_in[i] +  beta[0]*Mujj *kappa*temp;
					//}
					UPDATE_REQUIRED			= 1;
				}
				/////////////////////////////////////////////////////////////////////////////////
				else if(selectedAction==ACTION_ADD)
				{
					if(verbose>4) Rprintf("\t\t Action:add : %d \t deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
								nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					//Rprintf("\t\t newAlpha: %f\n",newAlpha);
					// B_phi*PHI2*SIGMA2        tmp = B_phi*PHI2 
					index					= M + 1;
					if(index > (basisMax -10) && iter>1) {
						Rprintf("bases: %d, warning: out of Memory!\n",index);
					}//return;
					if(index > (basisMax -1) && iter>1) {
						Rprintf("bases: %d, out of Memory,exiting program!\n",index);
						//exit(1);			
					}
					//BASIS_B_Phi	=BASIS_Phi*beta;
					UPDATE_REQUIRED		= ActionAddGmNeEN(BASIS_PHI, BASIS, Scales, PHI, phi, beta, Alpha,
						newAlpha, SIGMA, Mu, S_in, Q_in, nu, SIGMANEW, M_full, N, K, M);

					//			
					Used[M]			= nu + 1;						//new element

					//

					N_unused				= N_unused - 1;
					for(i=0;i<N_unused;i++)
					{                
						if(Unused[i]== (nu + 1))		Unused[i] =Unused[N_unused];
					}
					m[0]					= M + 1;
					M						= m[0];
					//for(i=0;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i],Mu[i]);
					//Rprintf(" \t\t S_in[73]: %f, Q_in[73]: %f \n",S_in[nu],Q_in[nu]);
					//for(i=0;i<4;i++) Rprintf(" \t\t SIGMANEW: %f \n",SIGMANEW[i]);
					//for(i=0;i<N;i++) Rprintf("phi: %f\n",phi[i]);
				}
				//
				else if(selectedAction==ACTION_DELETE)
				{
					if(verbose>4) Rprintf("\t\t Action: delete : %d deltaML: %f Alpha: %f s_i: %f q_i: %f S_i: %f Q_i: %f\n",
								nu + 1, deltaLogMarginal,newAlpha,S_out[nu],Q_out[nu],S_in[nu],Q_in[nu]);
					UPDATE_REQUIRED = ActionDelGmNeEN(PHI, Alpha, SIGMA, SIGMANEW, BASIS_PHI,
						Mu, S_in, Q_in, beta, jj, N, M, M_full);
					index					= M -1;
					//free deleted row of BASIS_PHI
					Free(BASIS_PHI[index]);
					//Used; Unused;
					Used[jj]				= Used[index];

					//
					N_unused				= N_unused + 1;
					Unused[N_unused -1]		= nu + 1;

					m[0]					= M -1;
					M						= m[0];
				//for(i=1;i<M;i++) Rprintf(" \t\t basis: %d :new weight: %f \n",Used[i-1],Mu2[i]);
				}

				//Rprintf("\t\t Update_required: %d\n",UPDATE_REQUIRED);
				//
				if(UPDATE_REQUIRED==1)
				{
					//for(i=0;i<M_full;i++)
					//{
					//	S_out[i]	= S_in[i];
					//	Q_out[i]	= Q_in[i];
					//}
					F77_CALL(dcopy)(&M_full,S_in,&inci,S_out,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
					F77_CALL(dcopy)(&M_full,Q_in,&inci,Q_out,&incj); 
					for(i=0;i<M;i++)
					{
						index					= Used[i] -1;
						S_out[index]			= Alpha[i]*S_in[index]/(Alpha[i]-S_in[index]);
						Q_out[index]			= Alpha[i]*Q_in[index]/(Alpha[i]-S_in[index]);
					}
					//for(i=0;i<M;i++)
					//{
					//	for(j=0;j<M;j++) SIGMA[i*M+j]	= SIGMANEW[i*M+j];
					//}
					MM = M*M;			
					F77_CALL(dcopy)(&MM,SIGMANEW,&inci,SIGMA,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x

					
					for(i=0;i<M;i++) gamma[i]		= 1- Alpha[i]*SIGMA[i*M+i];
				}
			} // //end of for(iU)
		}	
		
		if ((selectedAction==ACTION_TERMINATE)||(i_iter<=CNBetaUpdateStart)||(i_iter%5==0) || (nUpdate>=2))
		{	
			ee					= 0;
			//for(i=0;i<N;i++)
            //{
            //    PHI_Mu[i]       = 0;
            //    for(j=0;j<M;j++)	PHI_Mu[i]   = PHI_Mu[i] + PHI[j*N+i]*Mu[j];
			//	e[i]			= Targets[i] - PHI_Mu[i];
			//	ee				= ee + e[i]*e[i];
			//}
			lda 				= N;
			b_blas 				= 0;
			F77_CALL(dgemv)(&transa, &N, &M,&a_blas, PHI, &lda, Mu, &inci, &b_blas,PHI_Mu, &incj); 
			//y := alpha*A*x + beta*y, dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
			F77_CALL(dcopy)(&N,Targets,&inci,e,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
			b_blas 				= -1;
			F77_CALL(daxpy)(&N, &b_blas,PHI_Mu, &inci,e, &incj);//daxpy(n, a, x, incx, y, incy) y := a*x + y
	
			 ee = F77_CALL(ddot)(&N, e, &inci,e, &incj);
			 //Rprintf("Gaussian updat before: beta: %f\n",beta[0]);
			 betaZ1				= beta[0];
			 temp				= 0;
			//for(i=0;i<M;i++) temp = temp  + gamma[i];
			F77_CALL(daxpy)(&M, &a_blas,gamma, &inci,&temp, &inc0);//daxpy(n, a, x, incx, y, incy) y := a*x + y
			
			 beta[0]			= (N-temp)/ee;
			 varT				= varTargetsGmNeEN(Targets,N);
			 //Rprintf("Gaussian updat after: beta: %f\n",beta[0]);
			 if(beta[0]>(BetaMaxFactor/varT))	beta[0] = BetaMaxFactor/varT;
			 deltaLogBeta			= log(beta[0]) - log(betaZ1);
			 //
			 if (fabs(deltaLogBeta)>MinDeltaLogBeta)
			 {
				//Rprintf("final update\n");
				FinalUpdateGmNeEN(PHI,H,SIGMA,Targets,Mu,Alpha,beta,N, M);
				//for(i=0;i<M;i++) Rprintf("Mu: %f \n",Mu[i]);
				if (selectedAction!=ACTION_TERMINATE)
				{
					fEBLinearFullStatGmNeEN(beta,SIGMA,H, S_in, Q_in, S_out,Q_out,  BASIS, Scales, 
							PHI, BASIS_PHI, BASIS_Targets, Targets, Used, Alpha, Mu, 
							gamma, n, m, kdim, iteration,&i_iter);
				}
			 }
			 //Rprintf("\tend of loop values N_used: %d, N_unused %d\n", M,N_unused);
        }

		//Rprintf("\t\t selected Action: %d\n",selectedAction);
		//
        if(selectedAction==ACTION_TERMINATE && iniRemoved == 1) LAST_ITERATION =1;
		if((i_iter==itMax && M==1) || i_iter>itMax) LAST_ITERATION =1;
 //       if(i_iter==100)   LAST_ITERATION = 1;
 if(i_iter==itMax)   selectedAction=ACTION_TERMINATE;

		//Rprintf("\t\t Last_iteration value: %d\n",LAST_ITERATION);
    }
	
	//C_inv                       = beta*eye(N)-beta^2*PHI*SIGMA*PHI';
	double*PHIsig	= (double *) Calloc(N*M,double); // PHI *SIGMA
	transb = 'N';
	lda = N;
	ldb = M;
	ldc = N;
	ldk = M; //b copy
	b_blas = 1;
	c_blas = 0;
	F77_CALL(dgemm)(&transa, &transb,&N, &M, &ldk,&b_blas, PHI, &lda, SIGMA, &ldb, &c_blas, PHIsig, &ldc);
	//dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	//	C := alpha*op(A)*op(B) + beta*C,
	//for(i=0;i<N;i++)
	//{
	//	 for(j=0;j<M;j++)
	//	 { 
	//		PHIsig[j*N+i]	= 0;
	//		for(h=0;h<M;h++) PHIsig[j*N+i]		= PHIsig[j*N+i] + PHI[h*N+i]*SIGMA[j*M+h];
	//	 }
	 //}
	  //C_inv = PHIsig * PHI';
	transb = 'T';
	ldk = N;
	lda = N;
	ldb = N;
	ldc = N;
	b_blas = -beta[0]*beta[0];
	c_blas = 0;
	F77_CALL(dgemm)(&transa, &transb,&N, &ldk, &M,&b_blas, PHIsig, &lda, PHI, &ldb, &c_blas, C_inv, &ldc);

	
	//for(i=0;i<N;i++)
	//{
	//	 for(j=0;j<N;j++)
	//	 { 
	//		C_inv[j*N+i]	= 0;
	//		for(h=0;h<M;h++) C_inv[j*N+i]		= C_inv[j*N+i] + PHIsig[h*N+i]*PHI[h*N+j];
	//		C_inv[j*N+i]	= 0-C_inv[j*N+i]*beta[0]*beta[0];
	//	 }
	 //}
	 for(i=0;i<N;i++) C_inv[i*N+i]	=C_inv[i*N+i]	+ beta[0];
	Free(Unused);	
	Free(IniLogic);
	Free(BASIS_Targets);
//Free(BASIS_PHI);	
	for(i=0;i<M;i++)
	{
		Free(BASIS_PHI[i]);
	}	
	Free(BASIS_PHI);	
	Free(S_in);
	Free(Q_in);
	Free(S_out);	
	Free(Q_out);
	Free(anyToDelete);	
	Free(DeltaML);	
	Free(AlphaRoot);
	Free(Action);	
	Free(phi);	
	Free(PHI_Mu);
	Free(e);
	Free(SIGMANEW);
	Free(PHIsig);
	Free(gamma);
		//
	//Block Coordinate Ascent Algorithm Parameters;
	Free(nuAction);
	
}

/****************************************************************************/

// [Alpha,PHI2,Used,Unused,Mu2]=InitialCategory(BASIS,Targets,Scales,PHI2,Used,Alpha,Mu2,IniLogic)  //IniLogic: whether input is empty or not

void fEBInitializationGmNeEN2(double *Alpha, double *PHI, int *Used, int *Unused, double *BASIS, 
			double *Targets, double *Scales, int * initial, int N, int *m, int K, double *beta)
{
    //basis dimension
    int M,M_full,i,j,kk,index;
  	//kk					= K;
    //M_full				= K*(K+1)/2;
	M_full					= K;
	int IniLogic			= *initial;
    //INPUT
    if(IniLogic==0)							// is empty
    {
		m[0]				= 1;
		M					= m[0];
    }else									// not empty
    {
       	M					= m[0];
    }
    //output
    const double init_alpha_max     = 1e3;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1;
	double b_blas = 1;
	//lapack end
    
	if(IniLogic==0)            // is empty
    {
		//Rprintf("\t Inside Initialization, M: %d, K: %d\n",M, K);
        double proj_ini,proj;
        int loc1			= 0;
		int loc2			= 0;
        proj_ini			= 0;
		//Used[0]				= 1;
		
		//Block Coordinate Ascent Algorithm Parameters;
		double *projection = (double *) Calloc(M_full,double);
		double 	nAdd 		= 0.9; //only deltaML >nAdd will be considered as candidate
		
        for(i=0;i<K;i++)
        {
         	proj			= 0;
            //for(j=0;j<N;j++)                proj    = proj + BASIS[i*N+j]*Targets[j];
			readPtr1 = &BASIS[i*N];
			proj = F77_CALL(ddot)(&N, readPtr1, &inci,Targets, &incj);
            proj			= proj/Scales[i];
			projection[i] 	= fabs(proj);
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                //loc1        = i;
                // loc2        = i;
                //Used[0]		= i + 1;
            }
        }
		//Block Coordinate Ascent Algorithm Parameters;

		double cutoff = proj_ini * nAdd;
		for(i=0;i<K;i++)
		{
			if(projection[i] >=cutoff) 	projection[i] = 0;
			
		}
		proj_ini = 0;
		for(i=0;i<K;i++)
        {
		
            proj			= projection[i];
            if(fabs(proj) > fabs(proj_ini))
            {
                proj_ini    = fabs(proj);
                loc1        = i;
                //loc2        = i;
                Used[0]		= i + 1;
            }
        }
		
		
        /*for(i=0;i< (K - 1);i++)
        {
            for(j= (i+ 1); j<K; j++)
            {
                proj			= 0;
                for(k=0;k<N;k++)            proj    = proj + BASIS[i*N+k]*BASIS[j*N+k]*Targets[k];
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
        
        //
        if(loc1==loc2)
        {
            //for(i=0;i<N;i++)                PHI[i]			= BASIS[loc1*N+i]/Scales[loc1];
			readPtr1 		= &BASIS[loc1*N];
			F77_CALL(dcopy)(&N,readPtr1,&inci,PHI,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
			b_blas 			= 1/Scales[loc1];
			F77_CALL(dscal)(&N,&b_blas,PHI,&inci); 		//dscal(n, a, x, incx) x = a*x	
        }else
        {
            index				= Used[0] -1;
           	for(i=0;i<N;i++)                PHI[i]			= BASIS[loc1*N+i]*BASIS[loc2*N+i]/Scales[index];
        }
		
		// beta
		double stdT;
		stdT					= varTargetsGmNeEN(Targets,N);
		stdT					= sqrt(stdT);
		if (stdT<1e-6)	stdT	= 1e-6;
		beta[0]					= 1/pow((stdT*0.1),2);

		//Alpha
		//p                           = diag(PHI'*PHI)*beta;
        //q                           = (PHI'*Targets)*beta;
        //Alpha	= p^2/(q^2-p);
		double p,q;
		p						= 0;
		q						= 0;
		//for(i=0;i<N;i++) 
		//{
		//	p					= p + PHI[i]*PHI[i];
		//	q					= q + PHI[i]*Targets[i];
		//}
		p = F77_CALL(ddot)(&N, PHI, &inci,PHI, &incj);
		q = F77_CALL(ddot)(&N, PHI, &inci,Targets, &incj);
		p						= p*beta[0];
		q						= q*beta[0];
		Alpha[0]				= p*p/(q*q-p);
		
        if(Alpha[0]< 0) Alpha[0]				= init_alpha_max;
        if(Alpha[0]> init_alpha_max) Alpha[0]				= init_alpha_max;
Free(projection);
	}

	int IsUsed					= 0;
    kk							= 0;
    for(i=0;i<M_full;i++)
    {
        IsUsed					= 0;
        for(j=0;j<M;j++)
        {
            //index   = Used[j];
            if ((i+1)==Used[j])     IsUsed  = 1;
        }
        if(IsUsed==0)     
        {
            Unused[kk]		= (i+ 1);
            kk				= kk + 1;
        }
    }
}


void fEBInitializationGmNeEN(double *Alpha, double *PHI, int *Used, int *Unused, double *BASIS, 
			double *Targets, double *Scales, int * initial, int N, int *m, int K, double *beta)
{
    //basis dimension
    int M,M_full,i,j,kk,index;
  	//kk					= K;
    //M_full				= K*(K+1)/2;
	M_full					= K;
	int IniLogic			= *initial;
    //INPUT
    if(IniLogic==0)							// is empty
    {
		m[0]				= 1;
		M					= m[0];
    }else									// not empty
    {
       	M					= m[0];
    }
    //output
    const double init_alpha_max     = 1e2;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1;
	double b_blas = 1;
	//lapack end
    
	if(IniLogic==0)            // is empty
    {
		//Rprintf("\t Inside Initialization, M: %d, K: %d\n",M, K);
        double proj_ini,proj;
        int loc1			= 0;
		int loc2			= 0;
        proj_ini			= 0;		
		
		
		Used[0]				= 1;
        for(i=0;i<K;i++)
        {
         	proj			= 0;
            //for(j=0;j<N;j++)                proj    = proj + BASIS[i*N+j]*Targets[j];
			readPtr1 = &BASIS[i*N];
			proj = F77_CALL(ddot)(&N, readPtr1, &inci,Targets, &incj);
            proj			= proj/Scales[i];
			
            if(fabs(proj) < fabs(proj_ini))
            {
                proj_ini    = proj;
                loc1        = i;
                loc2        = i;
                Used[0]		= i + 1;
            }
        }
        /*for(i=0;i< (K - 1);i++)
        {
            for(j= (i+ 1); j<K; j++)
            {
                proj			= 0;
                for(k=0;k<N;k++)            proj    = proj + BASIS[i*N+k]*BASIS[j*N+k]*Targets[k];
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
        
        //
        if(loc1==loc2)
        {
            //for(i=0;i<N;i++)                PHI[i]			= BASIS[loc1*N+i]/Scales[loc1];
			readPtr1 		= &BASIS[loc1*N];
			F77_CALL(dcopy)(&N,readPtr1,&inci,PHI,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
			b_blas 			= 1/Scales[loc1];
			F77_CALL(dscal)(&N,&b_blas,PHI,&inci); 		//dscal(n, a, x, incx) x = a*x	
        }else
        {
            index				= Used[0] -1;
           	for(i=0;i<N;i++)                PHI[i]			= BASIS[loc1*N+i]*BASIS[loc2*N+i]/Scales[index];
        }
		
		// beta
		double stdT;
		stdT					= varTargetsGmNeEN(Targets,N);
		//stdT					= sqrt(stdT);
		//if (stdT<1e-6)	stdT	= 1e-6;
		//beta[0]					= 1/pow((stdT*0.1),2);
		beta[0] 				= 1/(stdT*0.01+1e-10);
		//Alpha
		//p                           = diag(PHI'*PHI)*beta;
        //q                           = (PHI'*Targets)*beta;
        //Alpha	= p^2/(q^2-p);
		double p,q;
		p						= 0;
		q						= 0;
		//for(i=0;i<N;i++) 
		//{
		//	p					= p + PHI[i]*PHI[i];
		//	q					= q + PHI[i]*Targets[i];
		//}
		p = F77_CALL(ddot)(&N, PHI, &inci,PHI, &incj);
		q = F77_CALL(ddot)(&N, PHI, &inci,Targets, &incj);
		p						= p*beta[0];
		q						= q*beta[0];
		Alpha[0]				= p*p/(q*q-p);
		
        if(Alpha[0]< 0) Alpha[0]				= init_alpha_max;
        if(Alpha[0]> init_alpha_max) Alpha[0]				= init_alpha_max;

	}

	int IsUsed					= 0;
    kk							= 0;
    for(i=0;i<M_full;i++)
    {
        IsUsed					= 0;
        for(j=0;j<M;j++)
        {
            //index   = Used[j];
            if ((i+1)==Used[j])     IsUsed  = 1;
        }
        if(IsUsed==0)     
        {
            Unused[kk]		= (i+ 1);
            kk				= kk + 1;
        }
    }
}


void LinearSolverGmNeEN(double * a, double *logout, int N, int M,double *output)
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
	//int i;
	//output				= (double *) realloc(NULL,M*sizeof(double));
	//for (i=0;i<M;i++) output[i] = logout[i];	
	int inci = 1;
	int incj = 1;
	F77_CALL(dcopy)(&M,logout,&inci,output,&incj);
	Free(jpvt);
	Free(work);
}

 /// ***********************************************************************************************   


void CacheBPGmNeEN(double **BASIS_PHI, double *BASIS_Targets, double *BASIS, double *PHI,
				double *Targets, double *scales,int N,int K,int M,int M_full)
{
	double	zTargets;
	double *z2					= (double *) Calloc(M,double);
	double *cache1				= (double *) Calloc(N,double);
	double *cache2				= (double *) Calloc(N*M,double);


	int i,h,l;
	//int kk						= K;

	//part 1 1-k
	for (i						=0;i<K;i++)
	{
		for(l=0;l<M;l++)
		{
			z2[l]					= 0;
			for(h=0;h<N;h++)
			{
				cache2[h*M+l]		= PHI[l*N+h] * BASIS[i*N+ h];
				z2[l]				= z2[l] + cache2[h*M+l];
			}
			//BASIS_PHI[l*M_full+i]	= z2[l]/scales[i];
			BASIS_PHI[l][i]	= (z2[l]/scales[i]);
		}

		zTargets					= 0;
		for(l=0;l<N;l++) 
		{
			cache1[l]				= BASIS[i*N+ l]*Targets[l];
			zTargets				= zTargets + cache1[l];	
		}
		BASIS_Targets[i]			= zTargets/scales[i];
	
		//part 2 k+1 to kk
		/*if(i<(K-1))
		{
			for (j					= (i+1);j<K;j++)
			{
				for(h=0;h<M;h++)
				{
					z2[h]	= 0;
					for(l =0;l<N;l++) z2[h] = z2[h] + cache2[l*M+h]*BASIS[j*N+l];
					BASIS_PHI[h*M_full+kk]	= z2[h]/scales[kk];
				}

				zTargets			= 0;
				for(h=0;h<N;h++)	zTargets		= zTargets + BASIS[j*N+h]*cache1[h];
				BASIS_Targets[kk]	= zTargets/scales[kk];
				kk					= kk + 1;
			}
		}*/
	}
	Free(z2);
	Free(cache1);
	Free(cache2);
}


/// *********[beta,SIGMA2,Mu2,S_in,Q_in,S_out,Q_out,BASIS_B_PHI,Intercept] ...
///                       	= FullstatCategory(BASIS,Scales,PHI2,Targets,Used,Alpha,Mu2,BASIS_CACHE) ************
    //Mu2 is the same size of M in PHI2; one dimension more than Alpha
    //Targets: nx1

void fEBLinearFullStatGmNeEN(double *beta, double * SIGMA, double *H, double *S_in, double * Q_in, double * S_out, 
				double * Q_out,   double *BASIS, double * Scales, double *PHI, double **BASIS_PHI,
				double *BASIS_Targets, double * Targets, int * Used, double *Alpha, double * Mu, 
				 double *gamma,int *n, int *m, int* kdim, int *iteration,int *i_iter)
{
    //basis dimension
    int N,K,M,i,j,p;
   	N					= *n;			// row number
    K					= *kdim;		// column number
	int M_full;
	//M_full				= (K+1)*K/2;
	M_full 				= K;
    M					= *m;
	//lapack
	int inci 			=1;
	int incj 			=1;
	//double *readPtr1;
	double a_blas 		= 1;
	double b_blas 		= 1;

	char transa 		= 'N';

	//int lda;
	//lapack end
	
	//Rprintf("Inside test: beta= %f\n",beta[0]);
//H2 = t(PHI2)*PHI2*beta + diag(A)
	if(iteration[0]		== 1 && i_iter[0]==0)			// initialize SIGMA
	{
		H[0]			= 0;
		//for(p=0;p<N;p++)	H[0]	= H[0] + PHI[p]*PHI[p];
		H[0] = F77_CALL(ddot)(&N, PHI, &inci,PHI, &incj);
		H[0]			=	H[0]* beta[0] + Alpha[0];
		
		//save a copy of H
		SIGMA[0]		= 1/H[0];
		//MatrixInverseGmNeEN(SIGMA,M);				//inverse of H2 is needed for wald score	
	}

		/*for(i=0;i<M;i++)
		{
			for (j=0;j<M;j++)	Rprintf("sigma:  %f\t", SIGMA[j*M+i]);
			Rprintf("\n");
		}*/


	//Muu				=SIGMA*(PHI.Transpose()*Targets)*beta;
	double * PHIt		= (double *) Calloc(M,double);
	//for(i=0;i<M;i++)
	//{	
	//	PHIt[i]		= 0;
	//	for(j=0;j<N;j++)	PHIt[i] =	PHIt[i] + PHI[i*N +j]*Targets[j];
	//}
	transa = 'T';
	a_blas = 1;
	b_blas = 0;
	//lda 	= N;
	F77_CALL(dgemv)(&transa, &N, &M,&a_blas, PHI, &N, Targets, &inci, &b_blas,PHIt, &incj);
	//y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
	
	
	//for(i=0;i<M;i++)
	//{
	//	Mu[i] = 0;
	//	for(j = 0;j<M;j++)		Mu[i] = Mu[i]  + SIGMA[j*M+i]*PHIt[j];
	//	Mu[i] = Mu[i]*beta[0];
	//}
	transa = 'N';
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &M,&a_blas, SIGMA, &M, PHIt, &inci, &b_blas,Mu, &incj);
	b_blas = beta[0];
	F77_CALL(dscal)(&M,&b_blas,Mu,&inci); //dscal(n, a, x, incx)
	
	//gamma
	for(i=1;i<M;i++)	gamma[i]	= 1- SIGMA[i*M+i] *Alpha[i];

    //Main loop
        //temp parameters: BPvector
    double *BPvector;
    BPvector			= (double *) Calloc(M,double);
    double tempSum,tempBPMu;

    for(i=0; i<M_full; i++)
    {
		for(j=0;j<M;j++)
		{
			BPvector[j]			= 0;
		//	for(p=0;p<M;p++)	BPvector[j]		= BPvector[j] + BASIS_PHI[p*M_full + i]*SIGMA[j*M+p]; 
		for(p=0;p<M;p++)	BPvector[j]		= BPvector[j] + BASIS_PHI[p][i]*SIGMA[j*M+p]; 
		}
		//readPtr1 	= &BASIS_PHI[i];
		//b_blas 		= 0;
		//F77_CALL(dgemv)(&transa, &M, &M,&a_blas, SIGMA, &M, readPtr1, &M_full, &b_blas,BPvector, &incj);
		
        tempSum					= 0;
		//for(j=0;j<M;j++)		tempSum			= tempSum + BPvector[j]*BASIS_PHI[j*M_full+i];
		for(j=0;j<M;j++)		tempSum			= tempSum + BPvector[j]*BASIS_PHI[j][i];
		//readPtr1 	= &BASIS_PHI[i];
		//tempSum		= F77_CALL(ddot)(&M, BPvector, &inci,readPtr1, &M_full);
		
		S_in[i]					= beta[0] - beta[0]*tempSum*beta[0];
		tempBPMu				= 0;
		//for(p=0;p<M;p++) tempBPMu = tempBPMu + BASIS_PHI[p*M_full + i] *Mu[p];
		for(p=0;p<M;p++) tempBPMu = tempBPMu + BASIS_PHI[p][i] *Mu[p];
		//readPtr1 	= &BASIS_PHI[i];
		//tempBPMu 	= F77_CALL(ddot)(&M, Mu, &inci,readPtr1, &M_full);	
		
		Q_in[i]					= beta[0]*(BASIS_Targets[i] - tempBPMu);
        //S_out[i]         		= S_in[i];
        //Q_out[i]          		= Q_in[i];
    }// Main loop ends
	F77_CALL(dcopy)(&M_full,S_in,&inci,S_out,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	F77_CALL(dcopy)(&M_full,Q_in,&inci,Q_out,&incj);  
	
    //S_out(Used),Q_out(Used)  
    //S_out(Used)				= (Alpha .* S_in(Used)) ./ (Alpha - S_in(Used));
    //_out(Used)       			= (Alpha .* Q_in(Used)) ./ (Alpha - S_in(Used));
    int index;

    for(i=0;i<M;i++)
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
	Free(PHIt);	
	Free(BPvector);
}


//****************************************************************
//** Matrix Inverse by call lapack library of cholesky decomposition and linear solver *
void MatrixInverseGmNeEN(double * a,int N)
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
void fEBDeltaMLGmNeEN(double *DeltaML, int *Action, double *AlphaRoot, int *anyToDelete,
				int *Used, int * Unused, double * S_out, double * Q_out, double *Alpha,
				double *a_lambda, double *b_Alpha, int N_used, int N_unused, int N,
				double *deltaLogMarginal, int *nu,double residual,double varY,int iter,int i_iter)//Block Coordinate Ascent Algorithm Parameters;
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
	int CNPriorityAddition		= 0;
	int CNPriorityDeletion		= 0;		
//****************version4.5************************* if residual variance <=10%of overall variance, switch priority to deletion.

	double epsilon = varY*0.1;// when >90% of residual variance explained, priority is deletion.
	double epsilonUp = varY*0.95;
	

	if(N_used<10)
	{
		CNPriorityAddition		= 1;
		CNPriorityDeletion		= 0;
		//Rprintf("Inside fEBDeltaMLBmNeEN: a: %f, b %f \n",a, b);
	}
	if(N_used>100 || N_used >=N || residual<=epsilon)
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
				anyToAdd = 1;
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

		//if(DeltaML[index]>0) Rprintf("there are additions:%d\tdeltaML:%f!\n",index,DeltaML[index]);
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
		//for(i=0;i<M_full;i++)
		//{
		//	if(Action[i] !=selectedAction) DeltaML[i] = 0;//
		//}


    

    if((anyToAdd==1 && CNPriorityAddition==1) || (anyToDelete[0]==1 && CNPriorityDeletion==1))
    {
        for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_REESTIMATE)											DeltaML[i]     = 0;
			else if (Action[i] == ACTION_DELETE)
            {
                    if(anyToAdd==1 && CNPriorityAddition==1 && CNPriorityDeletion!=1)
					{
						DeltaML[i]     = 0;
						//Rprintf("THERE ARE ADDTIONS, SET DELETION TO ZEROs!\n");
					}
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
	if((anyToAdd==0 &&iter==1 && i_iter<10) || (anyToAdd==0 && residual>=epsilonUp))
	{
		for(i=0;i<M_full;i++)
        {
            if (Action[i] == ACTION_DELETE)	DeltaML[i]     = 0;
		}
		//Rprintf("Force EBEN not to remove basis !\n");
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

//888888888888888888888888888888888888888888888888888888888888888888888888888888888888
int ActionAddGmNeEN(double **BASIS_PHI, double* BASIS, double*scales, double*PHI, double*Phi,
			double *beta, double* Alpha, double newAlpha, double*SIGMA, double*Mu, double*S_in,
			double*Q_in, int nu, double*SIGMANEW, int M_full,int N, int K, int M)
{
	double *BASIS_Phi		= (double *) Calloc(M_full,double);
	double *BASIS_B_Phi		= (double *) Calloc(M_full,double);
	double *mCi				= (double *) Calloc(M_full,double);
	double *z				= (double *) Calloc(N,double);
	//int kk					= K;
	int i,j,h;
	int index				= M + 1;
	double*   	tmp			= (double *) Calloc(M,double);
  	double*		tmpp		= (double *) Calloc(M,double);
	double s_ii,mu_i,TAU;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1;
	double b_blas = 1.0;

	//lapack end
	
	//Rprintf("\t\t Inside ActionAddGmNeEN: index: %d\n",index);
	for (i					= 0;i<K;i++)
	{
		BASIS_Phi[i]		= 0;
		for(h=0;h<N;h++)
		{
			z[h]			= BASIS[i*N+h]*Phi[h];
			BASIS_Phi[i]	= BASIS_Phi[i] + z[h];
		}
		BASIS_Phi[i]		= BASIS_Phi[i]/scales[i];
		BASIS_B_Phi[i]		= beta[0]*BASIS_Phi[i];
		//		
		/*if(i<K-1)
		{
			for (j					=(i+1);j<K;j++)
			{
				BASIS_Phi[kk]		= 0; 
				for(h=0;h<N;h++)	BASIS_Phi[kk]	= BASIS_Phi[kk] + BASIS[j*N + h] *z[h];
				BASIS_Phi[kk]		= BASIS_Phi[kk]/scales[kk];
				BASIS_B_Phi[kk]		= beta[0]*BASIS_Phi[kk];
				kk					= kk +1;
			}
		}*/
	}
	//tmp						= ((B_Phi_lowcase'*PHI)*SIGMA)';				
	for(i=0;i<M;i++)
    {
        tmp[i]					= 0;
        //for(j=0;j<N;j++)		tmp[i]			= tmp[i] + beta[0]*Phi[j]*PHI[i*N+j]; //M+1   x 1
		readPtr1 	= &PHI[i*N];
		tmp[i]  = F77_CALL(ddot)(&N, readPtr1, &inci,Phi,&incj);
    }
	b_blas = beta[0];
	F77_CALL(dscal)(&M,&b_blas,tmp,&inci); 		//dscal(n, a, x, incx) x = a*x
	
	
    for(i=0;i<M;i++)
    {
        tmpp[i]					= 0;
        //for(j=0;j<M;j++)	tmpp[i]				= tmpp[i] + tmp[j]*SIGMA[i*M+j];
		readPtr1 	= &SIGMA[i*M];
		tmpp[i]  = F77_CALL(ddot)(&M, readPtr1, &inci,tmp,&incj);	
    }//tmpp is tmp in matlab.
	//Alpha: M -> M + 1
	Alpha[M]				= newAlpha;                 //new element
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
	
	double * s_i			= (double *) Calloc(M,double);
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

	//S_in                = S_in - s_ii * mCi.^2;
	//Q_in                = Q_in - mu_i * mCi;
	//mCi                 = BASIS_B_Phi_lowcase - beta*BASIS_PHI*tmp;
	double temp;
	for(i=0;i<M_full;i++)
	{
		temp				= 0;
		//for(j=0;j<M;j++)	temp				= temp + BASIS_PHI[j*M_full + i]*tmpp[j];
		for(j=0;j<M;j++)	temp				= temp + BASIS_PHI[j][i]*tmpp[j];
		//readPtr1 			= &BASIS_PHI[i];
		//temp 				= F77_CALL(ddot)(&M, readPtr1, &M_full,tmpp,&incj); 
		
		mCi[i]				= BASIS_B_Phi[i] -  beta[0]*temp;
		//BASIS_PHI[M*M_full + i]					= BASIS_Phi[i];
		S_in[i]				= S_in[i] - mCi[i]*mCi[i]*s_ii;
		Q_in[i]				= Q_in[i] - mu_i*mCi[i];
	}
	BASIS_PHI[M]=BASIS_Phi;
	int UPDATE_REQUIRED		= 1;
	//Free(BASIS_Phi);
	Free(BASIS_B_Phi);	
	Free(mCi);
	Free(z);	
	Free(tmp);
	Free(tmpp);	
	Free(s_i);
	
	return  UPDATE_REQUIRED; 
 }

int ActionDelGmNeEN(double*PHI, double*Alpha, double*SIGMA, double*SIGMANEW, double**BASIS_PHI,
		double*Mu, double*S_in, double*Q_in, double *beta, int jj, int N, int M, int M_full)
 {
	 int i,j;
	int index				= M - 1;
	//lapack
	int inci =1;
	int incj =1;
	double *readPtr1, *readPtr2;
	//lapack end
	
	
	//
	Alpha[jj]				= Alpha[index];           //Alpha: M -> M - 1
            
	//
	//for(i=0;i<N;i++)		PHI[jj*N + i]		= PHI[index*N+i];		//M -> M -1
	readPtr1 = &PHI[jj*N];
	readPtr2 = &PHI[index*N];
	F77_CALL(dcopy)(&N,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x

    int Mujj;        
	Mujj					= Mu[jj];
	for(i=0;i<M;i++)  Mu[i] = Mu[i] - Mujj*SIGMA[jj*M +i]/SIGMA[jj*M + jj];

	//
	Mu[jj]					= Mu[index];        
	//------------------------------------------------------------------------------
	//BLOCK MIGRATION OF SIGMANEW ------------------JUN142013 -------IN EBELASTICNET
	double *tempSIGMA = (double *) Calloc((M*M),double);
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)	tempSIGMA[j*M + i]	= SIGMA[j*M + i] - SIGMA[jj*M+i]/SIGMA[jj*M+jj]*SIGMA[jj*M+j];
	}
	for(i=0;i<index;i++)
	{
		for(j=0;j<index;j++) SIGMANEW[j*index + i]	= tempSIGMA[j*M + i];
	}
	if(jj != index)// incase the one to be deleted is the last column.
	{
	//step 1: last column
	readPtr1 = &SIGMANEW[jj*index];
	readPtr2 = &tempSIGMA[index*M];
	F77_CALL(dcopy)(&index,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	//step 2: prepare for last row
	tempSIGMA[jj*M + M-1] = tempSIGMA[M*M-1];

	//last step 	
	readPtr1 = &SIGMANEW[jj];
	readPtr2 = &tempSIGMA[M-1];
	F77_CALL(dcopy)(&index,readPtr2,&M,readPtr1,&index);  //dcopy(n, x, incx, y, incy) ---> y = x
	}


	
	
	//for(i=0;i<M;i++)
	//{
	//	for(j=0;j<M;j++)	SIGMANEW[j*M + i]	= SIGMA[j*M + i] - SIGMA[jj*M+i]/SIGMA[jj*M+jj]*SIGMA[jj*M+j];
	//}
	//
	//for(i=0;i<M;i++)		SIGMANEW[jj*M+i]	= SIGMANEW[index*M + i];			//M isnot index here
	//for(i=0;i<M;i++)		SIGMANEW[i*M+jj]	= SIGMANEW[i*M + index];	
	//readPtr1 = &SIGMANEW[jj*M];
	//readPtr2 = &SIGMANEW[index*M];
	//F77_CALL(dcopy)(&M,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
		
	//readPtr1 = &SIGMANEW[jj];
	//readPtr2 = &SIGMANEW[index];
	//F77_CALL(dcopy)(&M,readPtr2,&M,readPtr1,&M);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	//BLOCK MIGRATION OF SIGMANEW ------------------JUN142013 -------IN EBELASTICNET	
	
	
	double temp;
	for(i=0;i<M_full;i++)
	{
		temp				= 0;
		//for(j=0;j<M;j++)	temp				= temp + BASIS_PHI[j*M_full + i]*SIGMA[jj*M + j];
		for(j=0;j<M;j++)	temp				= temp + BASIS_PHI[j][i]*SIGMA[jj*M + j];
		S_in[i]				= S_in[i] +  pow(beta[0]*temp,2)/SIGMA[jj*M + jj];
		Q_in[i]				= Q_in[i] +  beta[0]*temp *Mujj /SIGMA[jj*M + jj];
		//jPm = beta*temp
	}
	double *ptr;
	ptr=BASIS_PHI[jj];
	BASIS_PHI[jj]=BASIS_PHI[index];
	BASIS_PHI[index]=ptr;
	//BASIS_PHI(:,jj) = [];
	//for(i=0;i<M_full;i++) BASIS_PHI[jj*M_full + i]	= BASIS_PHI[index*M_full+i];
	//readPtr1 = &BASIS_PHI[jj*M_full];
	//readPtr2 = &BASIS_PHI[index*M_full];
	//F77_CALL(dcopy)(&M_full,readPtr2,&inci,readPtr1,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	
	int UPDATE_REQUIRED		=1;
	Free(tempSIGMA);
	return  UPDATE_REQUIRED;
}//end of ACTION_DELETE



double varTargetsGmNeEN(double* Targets,int N)
{
	int i;
	double meanT			= 0;
	double stdT				= 0;
	double varT;
	for(i=0;i<N;i++) meanT	= meanT + Targets[i];
	meanT					= meanT/N;
	for(i=0;i<N;i++) stdT	= stdT + pow((Targets[i] - meanT),2);
	varT					= stdT/(N-1);
	stdT					= sqrt(varT);
	return varT;
}


void FinalUpdateGmNeEN(double *PHI, double *H,double *SIGMA, double *Targets,
			double *Mu, double *Alpha, double *beta, int N, int M)
{
	int i;
	//PHI'*PHI*beta + diag(Alpha)
	//lapack
	int inci =1;
	int incj =1;

	double a_blas = 1;
	double b_blas = 1;

	int MM;
	char transa = 'T';
	char transb = 'N';
	int lda,ldb,ldc,ldk;
	//lapack end
	
	//PHI'*PHI*beta + diag(Alpha)
	//for(i=0;i<M;i++)
	//{
	//	for(j=0;j<M;j++)
	//	{
	//	H[j*M+i]				= 0;
	//	for(p=0;p<N;p++)		H[j*M+i]		= H[j*M+i] + PHI[i*N + p]*PHI[j*N+p];
	//	H[j*M+i]				= H[j*M+i]* beta[0];
	//	}
	//}
	ldk  	= N;
	lda 	= N;
	ldb 	= N;
	ldc 	= M;
	b_blas 	= 0;	
	F77_CALL(dgemm)(&transa, &transb,&M, &M, &ldk,&a_blas, PHI, &lda, PHI, &ldb, &b_blas, H, &ldc);
	//dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	//	C := alpha*op(A)*op(B) + beta*C,
	b_blas 	= beta[0];
	MM	 	= M*M;
	F77_CALL(dscal)(&MM,&b_blas,H,&inci); 		//dscal(n, a, x, incx) x = a*x
	
	for(i=0;i<M;i++)			H[i*M+i]		= H[i*M+i] + Alpha[i];
		
	//save a copy of H
	//for(i=0;i<M;i++)
	//{
	//	for (j=0;j<M;j++)		SIGMA[j*M+i]	= H[j*M+i];
	//}
	F77_CALL(dcopy)(&MM,H,&inci,SIGMA,&incj);  //dcopy(n, x, incx, y, incy) ---> y = x
	

	MatrixInverseGmNeEN(SIGMA,M);				//inverse of H2 is needed for wald score	

	//Muu						=SIGMA*(PHI.Transpose()*Targets)*beta;
double * PHIt				= (double *) Calloc(M,double);
	//for(i=0;i<M;i++)
	//{	
	//	PHIt[i]					= 0;
	//	for(j=0;j<N;j++)		PHIt[i]			=	PHIt[i] + PHI[i*N +j]*Targets[j];
	//}
	transa = 'T';
	a_blas = 1;
	b_blas = 0;
	//lda 	= N;
	F77_CALL(dgemv)(&transa, &N, &M,&a_blas, PHI, &N, Targets, &inci, &b_blas,PHIt, &incj);
	
	//for(i=0;i<M;i++)
	//{
	//	Mu[i]					= 0;
	//	for(j = 0;j<M;j++)		Mu[i]			= Mu[i]  + SIGMA[j*M+i]*PHIt[j];
	//	Mu[i]					= Mu[i]*beta[0];
	//}
	transa = 'N';
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &M,&a_blas, SIGMA, &M, PHIt, &inci, &b_blas,Mu, &incj);
	b_blas = beta[0];
	F77_CALL(dscal)(&M,&b_blas,Mu,&inci); //dscal(n, a, x, incx)
	
	
	Free(PHIt);	

}

