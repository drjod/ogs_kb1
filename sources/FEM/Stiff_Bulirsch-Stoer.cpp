#include "Stiff_Bulirsch-Stoer.h"
#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <fstream>



std::ofstream odeint_derivs;
std::ofstream odeint_yscal;
std::ofstream odeint_ystart;
std::ofstream rkqs_ytemp;
std::ofstream rkqs_yerr;
std::ofstream rkqs_errmax;
std::ofstream rkck_derivs; 
std::ofstream rkck_ytemp;
std::ofstream tstepinfo;
int idcs[30];

extern size_t aktueller_zeitschritt;


/* interne Deklarationen */
void nrerror(char *error_text);
double *dvector(long nl, long nh);
int *ivector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
bool ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
bool simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double [], int, long, double), long);
void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double [], int, long, double));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
double hmin, double *hnext, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], int, long, double),
	bool (*stifbs)(double [], double [], int, double *, double, double, double [],
	  double *, double *, void (*)(double, double [], double [], int, long, double)), 
  bool (*rkqs)(double [], double [], int, double *, double , double , double [],
	  double *, double *, void (*)(double , double [], double [], int, long, double), long),
  long, int);

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double , double [], double [], int, long, double), long node);
/*************************************************************************************/
/* Numerical Recipes Signum                                                          */
/*************************************************************************************/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*************************************************************************************/
/* Numerical Recipes Square                                                          */
/*************************************************************************************/

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

/*************************************************************************************/
/* Numerical Recipes double-max                                                      */
/*************************************************************************************/

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))

/*************************************************************************************/
/* Numerical Recipes double-min                                                      */
/*************************************************************************************/

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
(dminarg1) : (dminarg2))

/*************************************************************************************/
/* Numerical Recipes standard error handler                                           */
/*************************************************************************************/

void nrerror(const char *error_text)
//char error_text[];
{
   //void exit();

   fprintf(stderr,"\n \n Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   std::cout.flush();
   exit(1);
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate a double vector with subscript range v[nl..nh]                           */
/*************************************************************************************/

double *dvector(long nl, long nh)
{
   double *v;

   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v)                                        //CB
   {
      std::cout <<  "allocation failure in dvector()" << "\n";
      nrerror("allocation failure in dvector()");
   }
   return v-nl+NR_END;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate an int vector with subscript range v[nl..nh]                             */
/*************************************************************************************/

int *ivector(long nl, long nh)
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v)                                        //CB
   {
      std::cout <<  "allocation failure in ivector()" << "\n";
      nrerror("allocation failure in ivector()");
   }
   return v-nl+NR_END;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]               */
/*************************************************************************************/

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m)                                        //CB
   {
      std::cout <<  "allocation failure 1 in matrix()" << "\n";
      nrerror("allocation failure 1 in matrix()");
   }
   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl])                                   //CB
   {
      std::cout <<  "allocation failure 2 in matrix()" << "\n";
      nrerror("allocation failure 2 in matrix()");
   }
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free a double vector allocated with dvector()                                     */
/*************************************************************************************/

void free_dvector(double *v, long nl, long nh)
{
   long i = nh; i++;                              // to avoid warning on nh
   free((FREE_ARG) (v+nl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free an int vector allocated with ivector()                                       */
/*************************************************************************************/

void free_ivector(int *v, long nl, long nh)
{
   long i = nh; i++;                              // to avoid warning on nh
   free((FREE_ARG) (v+nl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free a double matrix allocated by dmatrix()                                       */
/*************************************************************************************/

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   long i = nch; i++;                             // to avoid warning on nh
   i=nrh;                                         // to avoid warning on nrh
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes Chapter 2.3                                                     */
/* LU Decomposition of matrix a[1..n][1..n]                                          */

/* a[1..n][1..n] input-Matrix und output-Matrix										 */
/* n input: Dimension des arrays a													 */
/* indx output: vector zeichnet row permutation auf (kann weg)						 */
/* d output:																		 */

/* Diese Funktion entspricht der Rockflow-Routine LU-Decomposition (solver.c),       */
/* indiziert aber ALLE arrays von 1-nn (Rockflow 0-(nn-1))                           */

/*************************************************************************************/

extern double **d,*x;

bool ludcmp(double **a, int n, int *indx, double *d)
{
   int i,imax=0,j,k;                              //imax=0 SB warning avoided
   double big,dum,sum,temp;
   double *vv;
  bool success = true;

   vv=dvector(1,n);
   *d=1.0;

   for (i=1;i<=n;i++)
   {
      big=0.0;
      for (j=1;j<=n;j++)
         if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0)                             //CB
      {
         std::cout << "Singular matrix in routine ludcmp" << "\n";
    //nrerror("Singular matrix in routine ludcmp");
    success = false;
   	free_dvector(vv,1,n);
    return success;
      }
      vv[i]=1.0/big;
   }

   for (j=1;j<=n;j++)
   {
      for (i=1;i<j;i++)
      {
         sum=a[i][j];
         for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      big=0.0;

      for (i=j;i<=n;i++)
      {
         sum=a[i][j];
         for (k=1;k<j;k++)
            sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (dum=vv[i]*fabs(sum)) >= big)
         {
            big=dum;
            imax=i;
         }
      }

      if (j != imax)
      {
         for (k=1;k<=n;k++)
         {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = -(*d);
         vv[imax]=vv[j];
      }

      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n)
      {
         dum=1.0/(a[j][j]);
         for (i=j+1;i<=n;i++) a[i][j] *= dum;
      }
   }

   free_dvector(vv,1,n);
  return success;
}


/*************************************************************************************/
/* Numerical Recipes Chapter 2.3                                                     */
/* Solves Set of n linear equations                                                  */

/* Diese Funktion entspricht der Rockflow-Routine lubksb_3 (solver.c),               */
/* indiziert aber ALLE arrays von 1-nn (Rockflow 0-(nn-1))                           */

/*************************************************************************************/

void lubksb(double **a, int n, int *indx, double b[])
{
   int i,ii=0,ip,j;
   double sum;

   for (i=1;i<=n;i++)
   {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
   }

   for (i=n;i>=1;i--)
   {
      sum=b[i];
      for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.4                                                    */
/* Use polynominal extrapolation to evaluate nv functions at x=0 by fitting a        */
/* polynominal to a sequence of estimates with progressively smaller values x=xest   */
/*************************************************************************************/

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
   int k1,j;
   double q,f2,f1,delta,*c;

   c=dvector(1,nv);
   x[iest]=xest;
   for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];

   if (iest == 1)
   {
      for (j=1;j<=nv;j++) d[j][1]=yest[j];

   } else
   {
      for (j=1;j<=nv;j++) c[j]=yest[j];
      for (k1=1;k1<iest;k1++)
      {
         delta=1.0/(x[iest-k1]-xest);
         f1=xest*delta;
         f2=x[iest-k1]*delta;

         for (j=1;j<=nv;j++)
         {
            q=d[j][k1];
            d[j][k1]=dy[j];
            delta=c[j]-q;
            dy[j]=f1*delta;
            c[j]=f2*delta;
            yz[j] += dy[j];
         }
      }

      for (j=1;j<=nv;j++) d[j][iest]=dy[j];
   }

   free_dvector(c,1,nv);
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.6                                                    */
/* Performs one step of semi-implicit midpoint rule                                  */
/*************************************************************************************/

bool simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double [], double [], int, long, double), long node)
{
   int i,j,nn,*indx;
   double d,h,x,**a,*del,*ytemp;
  bool success = true;

   indx=ivector(1,n);
   a=dmatrix(1,n,1,n);
   del=dvector(1,n);
   ytemp=dvector(1,n);
	h=htot/nstep;        // current stepsize

   for (i=1;i<=n;i++)
   {
      for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
      ++a[i][i];
   }

  // LU decomp of matrix a
	success = ludcmp(a,n,indx,&d);   
  if(success == 0){
  	free_dvector(ytemp,1,n);
	  free_dvector(del,1,n);
	  free_dmatrix(a,1,n,1,n);
	  free_ivector(indx,1,n);  
    return success;
  }
	// set up rhs for 1st step, yout for temp storage 
   for (i=1;i<=n;i++)
      yout[i]=h*(dydx[i]+h*dfdx[i]);
   lubksb(a,n,indx,yout);

  // 1st step
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+(del[i]=yout[i]);
   x=xs+h;

  // calculate derivatives with external user-provided function 
  // yout for temp storage of derivatives
	(*derivs)(x,ytemp,yout,n,node, h);

  // General step	
   for (nn=2;nn<=nstep;nn++)
   {
     	// set up rhs for general step
      for (i=1;i<=n;i++)
         yout[i]=h*yout[i]-del[i];

      /* Rockflow: void lubksb_3(double *a, long n, long *indx, double *b)*/
      /* achtung, 1D/2D array*/
      lubksb(a,n,indx,yout);
      for (i=1;i<=n;i++)
         ytemp[i] += (del[i] += 2.0*yout[i]);
      x += h;
		(*derivs)(x,ytemp,yout,n,node, h);
   }
	// set up rhs for last step
   for (i=1;i<=n;i++)
      yout[i]=h*yout[i]-del[i];
   lubksb(a,n,indx,yout);
  //take last step
   for (i=1;i<=n;i++)
      yout[i] += ytemp[i];

   free_dvector(ytemp,1,n);
   free_dvector(del,1,n);
   free_dmatrix(a,1,n,1,n);
   free_ivector(indx,1,n);
  return success;
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.6                                                    */
/* Semi-implicit extrapolation step for integrating stiff ODEs                       */
/*************************************************************************************/

#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define SCALMX 0.1

double **d,*x;

/* Input: y=current_conc, dydx=their_derivs, xx=current_time, htry=suggested_stepsize     */
/* Ouput: y=updated_conc,  xx=end_time, hdid=achieved_stepsize , hnext= estimated_next_ss */
bool stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double [], int, long, double), long node)
{
   int i,iq,k,kk,km;
   static int first=1,kmax,kopt,nvold = -1;
   static double epsold = -1.0,xnew;
   double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
   double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
   static double a[IMAXX+1];
   static double alf[KMAXX+1][KMAXX+1];
   static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
   //	static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70,98,138,194,274,386};
   // Num Recip S. 744
   // Differenz zwischen zwei Werten muss ein Vielfaches von 4 sein
   // So wählen, dass Verhältnis der Werte <= 5/7 ist
   // z.B. 10, 14 -> 10/14 <= 5/7
   // nächster wäre 18, aber 14/18 > 5/7, deshalb 22 mit 14/22 <= 5/7
   int reduct,exitflag=0;
   km=0;                                          //SB avoid warnings
   red = 0.0;                                     //SB avoid warning
   errmax = 0.0;                                  // SB avoid warning
   scale = 0.0;                                   // SB avoid warning
  bool success = true;

   d=dmatrix(1,nv,1,KMAXX);
   dfdx=dvector(1,nv);
   dfdy=dmatrix(1,nv,1,nv);
   err=dvector(1,KMAXX);
   x=dvector(1,KMAXX);
   yerr=dvector(1,nv);
   ysav=dvector(1,nv);
   yseq=dvector(1,nv);

  // reinitialize as eps is a new tolerance
  // or if nv have changed
  	if(eps != epsold || nv != nvold) 
   {
		*hnext = xnew = -1.0e29;  // impossible values
		eps1=SAFE1*eps;
    // compute work coefficients Ak
		a[1]=nseq[1]+1;           
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		// compute alpha(k,q)
      for (iq=2;iq<=KMAXX;iq++)
      {
         for (k=1;k<iq;k++)
            alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
               ((a[iq+1]-a[1]+1.0)*(2*k+1))));
      }
      epsold=eps;
    // save nv
      nvold=nv;
    // add cost of Jacobian evals to work coeffs a[]
      a[1] += nv;
      for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
    // Determine opt. row no. for convergence
      for (kopt=2;kopt<KMAXX;kopt++)
         if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
      kmax=kopt;
   }
   h=htry;
  // save starting values
   for (i=1;i<=nv;i++) ysav[i]=y[i];
  // evaluate jacobian matrix, update dfdx, dfdy
   jacobn(*xx,y,dfdx,dfdy,nv,node);
	// new stepsize or new integration --> re-establish order window
   if (*xx != xnew || h != (*hnext))
   {
      first=1;
      kopt=kmax;
   }
   reduct=0;

  // start stepping
   for (;;)
   {
    // evaluate the sequence of modified midpoint integrations
      for (k=1;k<=kmax;k++)
      {
         xnew=(*xx)+h;
         if (xnew == (*xx))                       //CB
         {
            std::cout << "step size underflow in stifbs" << "\n";
            //nrerror("step size underflow in stifbs");
				        success = false;
					      	free_dvector(yseq,1,nv);
						      free_dvector(ysav,1,nv);
						      free_dvector(yerr,1,nv);
						      free_dvector(x,1,KMAXX);
						      free_dvector(err,1,KMAXX);
						      free_dmatrix(dfdy,1,nv,1,nv);
						      free_dvector(dfdx,1,nv);
						      free_dmatrix(d,1,KMAXX,1,KMAXX);
					        return success; 
         }
      // semi-implicite midpoint algorithm
			success = simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs,node);
      if(success==0){
      	free_dvector(yseq,1,nv);
	      free_dvector(ysav,1,nv);
	      free_dvector(yerr,1,nv);
	      free_dvector(x,1,KMAXX);
	      free_dvector(err,1,KMAXX);
	      free_dmatrix(dfdy,1,nv,1,nv);
	      free_dvector(dfdx,1,nv);
	      free_dmatrix(d,1,KMAXX,1,KMAXX);
        return success;       
      }
      // squared since error is even
         xest=DSQR(h/nseq[k]);

         pzextr(k,xest,yseq,y,yerr,nv);
      // compute normalized error estimate eps(k)

         if (k != 1)
         {
            errmax=TINY;
            for (i=1;i<=nv;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
				// scale error relative to tolerance
            errmax /= eps;
            km=k-1;
            err[km]=pow(errmax/SAFE1,1.0/(double)(2*km+1));
         }

         if (k != 1 && (k >= kopt-1 || first))
         {
            if (errmax < 1.0)
            {
               exitflag=1;
               break;
            }
            if (k == kmax || k == kopt+1)
            {
               red=SAFE2/err[km];
               break;
            }
            else if (k == kopt && alf[kopt-1][kopt] < err[km])
            {
               red=1.0/err[km];
               break;
            }
            else if (kopt == kmax && alf[km][kmax-1] < err[km])
            {
               red=alf[km][kmax-1]*SAFE2/err[km];
               break;
            }
            else if (alf[km][kopt] < err[km])
            {
               red=alf[km][kopt-1]/err[km];
               break;
            }
         }
      }
      //		if (exitflag) std::cout << " Exitflag > 0 in stifbs of biodegradation" << "\n";
      if (exitflag) break;
		// reduce stepsize by at least REDMIN and at most by REDMAX
      red=DMIN(red,REDMIN);
      red=DMAX(red,REDMAX);
      h *= red;
      reduct=1;
	} // try again

  // successfull step was taken
   *xx=xnew;
   *hdid=h;
   first=0;
   wrkmin=1.0e35;
  // compute optimal row for convergence and corresponding stepsize
   for (kk=1;kk<=km;kk++)
   {
      fact=DMAX(err[kk],SCALMX);
      work=fact*a[kk+1];
      if (work < wrkmin)
      {
         scale=fact;
         wrkmin=work;
         kopt=kk+1;
      }
   }
   *hnext=h/scale;
   if (kopt >= k && kopt != kmax && !reduct)
   {
    // check for possible order increse but not if stepsize was just reduced
      fact=DMAX(scale/alf[kopt-1][kopt],SCALMX);
      if (a[kopt+1]*fact <= wrkmin)
      {
         *hnext=h/fact;
         kopt++;
      }
   }

   free_dvector(yseq,1,nv);
   free_dvector(ysav,1,nv);
   free_dvector(yerr,1,nv);
   free_dvector(x,1,KMAXX);
   free_dvector(err,1,KMAXX);
   free_dmatrix(dfdy,1,nv,1,nv);
   free_dvector(dfdx,1,nv);
   free_dmatrix(d,1,KMAXX,1,KMAXX);

  return success;
}


#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef SCALMX
#undef NRANSI

#define NRANSI
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define DEBUGRK 0
#define TSTEPINFO 0


/*************************************************************************************/
/* Numerical Recipes Chapter 16.2                                                    */
/* stepper routine                                                                   */
/* Adaptive step size controlled 5th order Runge Kutta method                        */
/*************************************************************************************/

bool rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double , double [], double [], int, long, double), long node)
{

	void rkck(double y[], double dydx[], int n, double x, double h,
		double yout[], double yerr[], void (*derivs)(double , double [], double [], int, long, double), long node);
	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;
 bool success = false;
 double *epsi;

 // y[]    - c(step)
 // dydx[] - dcdt(step)
 // yscal  - scaling of errors

	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry; //set stepsize to initial trial value
	epsi=dvector(1,n);

 for (i=1;i<=n;i++) 
   epsi[i]=eps;
 //epsi[2]= eps*1000;
 //epsi[20]= eps*1000;
 //epsi[21]= eps*1000;
	
 for (;;) {
	  rkck(y,dydx,n,*x,h,ytemp,yerr,derivs, node); // take a step --> update yerr, ytemp
	  //evaluate accuracy
   errmax=0.0;

   //for (i=1;i<=n;i++) 
     //yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; // for Runge Kutta, this is done in driver routine odeint
     //yscal[i]= eps + eps*DMAX(fabs(y[i]),fabs(ytemp[i]));  
     //yscal[i]= 0 + eps*fabs(y[i]);  
   for (i=1;i<=n;i++) 
     //if(i==2 || i==20 || i==21 ) 
     //  errmax=DMAX(errmax,errmax/10); // find max rel. error |d0/d1|
     //  errmax=DMAX(errmax, fabs(yerr[i]/yscal[i])); // find max rel. error |d0/d1|
       errmax=DMAX(errmax, fabs(yerr[i]/yscal[i])/(epsi[i])); // find max rel. error |d0/d1|
       //errmax += pow(yerr[i]/yscal[i],2);
       //errmax = DMAX(errmax, pow(yerr[i]/yscal[i],2));
   //errmax /= eps;              // SCALE relative to required tolerance
   //errmax = sqrt(errmax/n);
   //errmax = sqrt(errmax);

   if(DEBUGRK>0){
  
   rkqs_yerr << "node" << node <<" " << *x << " " << h << " ";
   for (i=1;i<=25;i++) rkqs_yerr << " " << fabs(yerr[ idcs[i] ]) ;
   rkqs_yerr << "\n";

   rkqs_ytemp << "node" << node <<" " << *x << " " << h << " ";
   for (i=1;i<=25;i++) rkqs_ytemp << " " << ytemp[ idcs[i] ] ;
   rkqs_ytemp << "\n";
   
   rkqs_errmax << "node" << node <<" " << *x << " " << h << " ";
   for (i=1;i<=25;i++) {
       rkqs_errmax << " " << fabs((yerr[ idcs[i] ] / yscal[ idcs[i] ]))/epsi[i];
   }
   rkqs_errmax << " " << errmax; 
   rkqs_errmax << "\n";  
  }

   if (errmax <= 1.0) break;   // success, break and compute next step size

   // failed, truncation error too large
   htemp=SAFETY*h*pow(errmax,PSHRNK);                 // reduce step size
   h=(h >= 0.0 ?  (htemp,0.1*h) : DMIN(htemp,0.1*h)); // no more than factor of 1/10
	  xnew=(*x)+h;
    if (xnew == *x){ 
      std::cout << "step size underflow in rkqs" << "\n";
      //nrerror("stepsize underflow in rkqs");
   	  free_dvector(ytemp,1,n);
      free_dvector(yerr,1,n);
      success = false;
      return success;
   }
 }

  if(DEBUGRK>0){
   rkqs_yerr << "\n";
   rkqs_ytemp << "\n";
   rkqs_errmax << "\n";  
  }

 if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;

//if(*hnext > 3500) *hnext=3500;

	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];  // update y[i] concentration vector
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
	free_dvector(epsi,1,n);
 success = true;
 return success;

}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI


#define NRANSI

/*************************************************************************************/
/* Numerical Recipes Chapter 16.2                                                    */
/* 5th order Runge Kutta Cash Karp step                                              */
/*************************************************************************************/
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double , double [], double [], int, long, double), long node)
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	ak2=dvector(1,n);
	ak3=dvector(1,n);
	ak4=dvector(1,n);
	ak5=dvector(1,n);
	ak6=dvector(1,n);
	ytemp=dvector(1,n);

	for (i=1;i<=n;i++) ytemp[i]=y[i]+b21*h*dydx[i]; // step 1 

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << dydx[ idcs[i] ] ;
   rkck_derivs << "\n";

   rkck_ytemp << node << " " << x << " " << h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << ytemp[ idcs[i] ] ;
   rkck_ytemp << "\n";
 }

	(*derivs)(x+a2*h,ytemp,ak2, n, node, a2*h);     // step 2
	for (i=1;i<=n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]); 

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << a2*h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << ak2[ idcs[i] ] ;
   rkck_derivs << "\n";

   rkck_ytemp << node << " " << x << " " << a2*h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << ytemp[ idcs[i] ] ;
   rkck_ytemp << "\n";
 }

	(*derivs)(x+a3*h,ytemp,ak3, n, node, a3*h);     // step 3
	for (i=1;i<=n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << a3*h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << ak3[ idcs[i] ] ;
   rkck_derivs << "\n";

   rkck_ytemp << node << " " << x << " " << a3*h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << ytemp[ idcs[i] ] ;
   rkck_ytemp << "\n";
 }

	(*derivs)(x+a4*h,ytemp,ak4, n, node, a4*h);     // step 4
	for (i=1;i<=n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << a4*h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << ak4[ idcs[i] ] ;
   rkck_derivs << "\n";

   rkck_ytemp << node << " " << x << " " << a4*h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << ytemp[ idcs[i] ] ;
   rkck_ytemp << "\n";
 } 
 
 (*derivs)(x+a5*h,ytemp,ak5, n, node, a5*h);     // step 5
	for (i=1;i<=n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << a5*h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << ak5[ idcs[i] ] ;
   rkck_derivs << "\n";

   rkck_ytemp << node << " " << x << " " << a5*h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << ytemp[ idcs[i] ] ;
   rkck_ytemp << "\n";
 } 
 
 
 (*derivs)(x+a6*h,ytemp,ak6, n, node, a6*h);     // step 6
 // accumulate increments with proper weights
 for (i=1;i<=n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);

 if(DEBUGRK>0){
   rkck_derivs << node << " " << x << " " << a6*h << " " ;
   for (i=1;i<=25;i++) rkck_derivs << " " << ak6[ idcs[i] ] ;
   rkck_derivs << "\n" << "\n";

   rkck_ytemp << node << " " << x << " " << a6*h << " " ;
   for (i=1;i<=25;i++) rkck_ytemp << " " << yout[ idcs[i] ] ;
   rkck_ytemp << "\n" << "\n";
 }


 // error estimate delta=y2-y1 = difference between 4th and 5th order method
 for (i=1;i<=n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

	free_dvector(ytemp,1,n);
	free_dvector(ak6,1,n);
	free_dvector(ak5,1,n);
	free_dvector(ak4,1,n);
	free_dvector(ak3,1,n);
	free_dvector(ak2,1,n);
}
#undef NRANSI


/*************************************************************************************/
/* Numerical Recipes Chapter 16.2                                                    */
/* Driver for Runge-Kutta Algorithm fitted to stiff Burlish-Stoer Method (stifbs)    */
/* Or 5th order Runge-Kutta Cash-Karp method (rkqs)                                  */
/* Integrates starting values ystart[1..nvar] from x1 to x2 using an adaptive step   */
/* size control                                                                      */
/*************************************************************************************/

/*extern int kmax,kount;
extern double *xp,**yp,dxsav;
*/
bool odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, double *nexth, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], int, long, double),
	bool (*stifbs)(double [], double [], int, double *, double, double, double [],
	  double *, double *, void (*)(double, double [], double [], int, long, double), long), 
  bool (*rkqs)(double [], double [], int, double *, double , double , double [],
	  double *, double *, void (*)(double , double [], double [], int, long, double),long),
  long node, int SolverType)
{
	int nstp,i;
//	double xsav,x,hnext,hdid,h;
	double x,hdid,hnext,h;
	double *yscal,*y,*dydx;
  bool success=false;


  //if(TSTEPINFO>0){
  //  if(aktueller_zeitschritt = 1)
  //    tstepinfo.open("tstepinfo.txt");  
  //  else
  //    tstepinfo.open("tstepinfo.txt", std::ios::app);  
  //  tstepinfo << "Timestep: " << aktueller_zeitschritt << " Time: " << x1 << " Node: " << node << "\n";
  //}


//  if(DEBUGRK>0){
//
//  odeint_derivs.precision(12);
//  odeint_derivs.open("odeint_derivs.txt");
//  odeint_yscal.precision(12);
//  odeint_yscal.open("odeint_yscal.txt");
//  odeint_ystart.precision(12);
//  odeint_ystart.open("odeint_ystart.txt");
//  rkqs_yerr.precision(12);
//  rkqs_yerr.open("rkqs_yerr.txt");
//  rkqs_ytemp.precision(12);
//  rkqs_ytemp.open("rkqs_ytemp.txt");
//  rkqs_errmax.precision(12);
//  rkqs_errmax.open("rkqs_errmax.txt");
//  rkck_derivs.precision(12);
//  rkck_derivs.open("rkck_derivs.txt"); 
//  rkck_ytemp.precision(12);
//  rkck_ytemp.open("rkck_ytemp.txt");
//}

  // ystart[1..nvar] starting values 
  // nvar = no. species
  // x1 ... x2 = time interval
  // eps required accuracy
  // h1 first trial step size 
  // hmin min step size
  // nexth next estimated step size
  // nok, no of good steps
  // nbad no of bad steps

	yscal=dvector(1,nvar); // scaling of error
	y=dvector(1,nvar);     // to store start concentrations
	dydx=dvector(1,nvar);  // to store derivatives dcdt
	x=x1;                  // current time during time stepping, start with x1 = t0
	h=SIGN(h1,x2-x1);

 //	*nok = (*nbad) = kount = 0;
	*nok = (*nbad) = 0;
  // to store start concentrations
   for (i=1;i<=nvar;i++) y[i]=ystart[i];

  //if(DEBUGRK>0){

  // odeint_derivs << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  // odeint_ystart << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  // odeint_yscal << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  //rkqs_yerr << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  //rkqs_ytemp << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  //rkqs_errmax << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";
  //rkck_derivs  << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";  
  //rkck_ytemp << "T_start: " << x1 << ", h0: " << h  << ", T_end: " << x2 << "\n";

  // odeint_ystart << "node" << node <<" " << x << " " << h << " " ;
  // for (i=1;i<=25;i++) odeint_ystart << " " << y[ idcs[i] ] ;
  // odeint_ystart << "\n";
  //}
  //no storage of intermediate results
  //if (kmax > 0) xsav=x-dxsav*2.0;

   for (nstp=1;nstp<=MAXSTEP;nstp++)
   {

    // calculate derivatives with external user-provided function */
		 (*derivs)(x,y,dydx,nvar,node, h);

  //if(DEBUGRK>0){
  // odeint_derivs << "node" << node << " " << x << " "; 
  // odeint_derivs << h <<" "; 
  // for (i=1;i<=25;i++) odeint_derivs << " " << dydx[ idcs[i] ] ;
  // odeint_derivs << "\n";
  //}
  
  // Preconditioning 
    // Scaling used to monitor accuracy
  if(SolverType==1) { 
      for (i=1;i<=nvar;i++) yscal[i]=DMAX(1,fabs(y[i])); // recommended for stiff problems
      //for (i=1;i<=nvar;i++) yscal[i]=DMAX(1.0,fabs(y[i])); // recommended for stiff problems
  }
  else if(SolverType==2){ 
      //for (i=1;i<=nvar;i++) yscal[i]=DMAX(1,fabs(y[i])); // recommended for stiff problems
      for (i=1;i<=nvar;i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; // for Runge Kutta
  }

  //if(DEBUGRK>0){
  //  
  //  odeint_yscal << "node" << node <<" " << x << " "; 
  // odeint_yscal << h <<" "; 
  // for (i=1;i<=25;i++) odeint_yscal << " " << yscal[ idcs[i] ] ;
  // odeint_yscal << "\n";
  //}
    
    //no storage of intermediate results
    //if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
	  	//	xp[++kount]=x;
		  //	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
		  //	xsav=x;
		  //}

		  // if stepsize can overshoot, limit to x2-x
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

    // call stepper 
    if(SolverType==1){ // stiff burlisch-stoer  
      success = (*stifbs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs, node);
      if(success==0) {
			     free_dvector(dydx,1,nvar);
			     free_dvector(y,1,nvar);
			     free_dvector(yscal,1,nvar);
        return success;
      }
    }
    else if(SolverType==2){ // Runge Kutta
      success = (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,node);
      if(success==0) {
			     free_dvector(dydx,1,nvar);
			     free_dvector(y,1,nvar);
			     free_dvector(yscal,1,nvar);
        return success;
      }
    }

    //if(DEBUGRK>0){
    //  odeint_ystart << "node" << node <<" " << x << " " << hdid << " " ;
    //  for (i=1;i<=25;i++) odeint_ystart << " " << y[ idcs[i] ] ;
    //  odeint_ystart << "\n";
    //}

    // save info on success or failure
      if (hdid == h) ++(*nok); else ++(*nbad);
      if ((x-x2)*(x2-x1) >= 0.0)
      {
         for (i=1;i<=nvar;i++) ystart[i]=y[i];
			      //if (kmax) {
			      //	xp[++kount]=x;
			      //	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			      //}
         free_dvector(dydx,1,nvar);
         free_dvector(y,1,nvar);
         free_dvector(yscal,1,nvar);
         *nexth=hnext;
			      success = true;

      //if(DEBUGRK>0){

      //   odeint_ystart.close();
      //   odeint_derivs.close();
      //   odeint_yscal.close();
      //   rkqs_yerr.close();
      //   rkqs_ytemp.close();
      //   rkqs_errmax.close();
      //   rkck_derivs.close();  
      //   rkck_ytemp.close();
      //}

         return success;
      }
      if (fabs(hnext) <= hmin)                    //CB
      {
        std::cout << "Step size too small in odeint" << "\n";
        //nrerror("Step size too small in odeint");
        free_dvector(dydx,1,nvar);
			     free_dvector(y,1,nvar);
			     free_dvector(yscal,1,nvar);
			     success = false;
        return success;
      }
      h=hnext;
   }
                                                  //CB
   std::cout << "Too many steps in routine odeint" << "\n";
	  //nrerror("Too many steps in routine odeint");
   free_dvector(dydx,1,nvar);
	  free_dvector(y,1,nvar);
	  free_dvector(yscal,1,nvar);
	  success = false;
   return success;
}


#undef MAXSTEP
#undef TINY
#undef DEBUGRK
