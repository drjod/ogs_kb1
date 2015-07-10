#ifndef _STEPPERBULISCHSTOER_H 
#define _STEPPERBULISCHSTOER_H

#include <limits>
#include "Eigen/Eigen"

#ifndef SQR
template<class T>
inline T SQR(const T a) {return a*a;}
#endif

#ifndef MAX
template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}
#endif

#ifndef MIN
template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}
#endif

class StepperBase
{
public:
	StepperBase(Eigen::VectorXd &yy, 
		        Eigen::VectorXd &dydxx, 
				double          &xx,
				const double    atoll, 
				const double    rtoll, 
				bool            dens): x(xx), y(yy), dydx(dydxx), atol(atoll), rtol(rtoll), dense(dens), n(y.size()), neqn(n), yout(n), yerr(n) {};

	double          &x; 
	double          xold; 
	Eigen::VectorXd &y; 
	Eigen::VectorXd &dydx; 
	double          atol; 
	double          rtol;
	bool            dense;
	double          hdid;
	double          hnext;
	double          EPS;
	size_t          n, neqn;
	Eigen::VectorXd yout, yerr; 
};



template <class D>
class StepperBulischStoer : StepperBase
{
public:
	StepperBulischStoer(Eigen::VectorXd &yy, Eigen::VectorXd &dydxx, double &xx, double atoll, double rtoll, bool dens)
	: StepperBase(yy, dydxx, xx, atoll, rtoll, dens), KMAXX(8), IMAXX(KMAXX+1), 
	  nseq(IMAXX), cost(IMAXX), table(KMAXX, n), dydxnew(n), coeff(IMAXX, IMAXX), errfac(2*IMAXX+2),
	  ysave(IMAXX,n), fsave(IMAXX*(2*IMAXX+1),n), ipoint(IMAXX+1), dens((2*IMAXX+5)*n)
	{
		size_t i,j,k,l; 
		EPS=std::numeric_limits<double>::epsilon();

		if (dense) // Choose the sequence (17.3.23) ...
			for (i=0; i<IMAXX; i++ )
				nseq[i]=4*i+2;
 		else // ... or (17.3.6).
			for (i=0; i<IMAXX; i++)
				nseq[i]=2*(i+1);

		cost[0]=nseq[0]+1; // Equation (17.3.12).
	
		for (k=0;k<KMAXX;k++) 
			cost[k+1]=cost[k]+nseq[k+1];
		hnext=-1.0e99; // Impossible value.
	
		double logfact=-log10( MAX(1.0e-12,rtol) )*0.6+0.5;

		k_targ= MAX((size_t)1, MIN( KMAXX-1,(size_t)logfact) ); // Initial estimate of optimal k.

		for ( k=0; k<IMAXX;  k++) // Coecients in equation (17.3.8).
		{
			for  ( l=0; l<k;  l++)
			{
				double ratio=double(nseq[k])/nseq[l];
				coeff(k,l)=1.0/(ratio*ratio-1.0);
			}
		}

		size_t ip5;
		double e;
		for ( i=0; i<2*IMAXX+1; i++)
		{
			ip5=i+5;
			errfac[i]=1.0/(ip5*ip5);
			e = 0.5*sqrt(double(i+1)/ip5);
			for  ( j=0; j<=i; j++)
			{
				errfac[i] *= e/(j+1);
			}
		}

		ipoint[0]=0;

		size_t njadd;
		for ( i=1; i<=IMAXX; i++)
		{
			njadd=4*i-2;
			if  (nseq[i-1] > (int)njadd)
				njadd++;		
			ipoint[i]=ipoint[i-1]+njadd;
		}
	};
	
	~StepperBulischStoer(void){};

	void dy (Eigen::VectorXd &y, 
		     const double    htot, 
			 const size_t    k, 
			 Eigen::VectorXd &yend, 
			 int             &ipt, 
			 D               *derivs_class)
	{
		size_t i,nn; 
		Eigen::VectorXd ym(n),yn(n);
		size_t  nstep=nseq[k];
		double  h=htot/nstep; // Stepsize this trip.
		for (i=0; i<n; i++) // First step.
		{ 
			ym[i]=y[i];
			yn[i]=y[i]+h*dydx[i];
		}
	
		double xnew=x+h;
	
		derivs_class->eval(xnew,yn,yend); // Use yend for temporary storage of derivatives. 
	
		double h2=2.0*h;
	
		for ( nn=1; nn<nstep; nn++)
		{
			if (dense && nn == nstep/2)
			{
				for ( i=0; i<n; i++)
					ysave(k,i)=yn[i];
			}
		
			if (dense && abs( (int)(nn-nstep/2)) <= (int)(2*k+1) )
			{
				ipt++;
				for ( i=0; i<n; i++)
					fsave(ipt,i)=yend[i];
			}
		
			for ( i=0; i<n; i++ )
			{   //General step.
				double swap=ym[i]+h2*yend[i];
				ym[i]=yn[i];
				yn[i]=swap;
			}
		
			xnew  +=  h;
			derivs_class->eval(xnew,yn,yend);
		}
	
		if ( dense && nstep/2 <= 2*k+1)
		{
			ipt++;
			for  ( i=0;i<n;i++)
				fsave(ipt,i)=yend[i];
		}
	
		for  ( i=0; i<n; i++) // Last step.
			yend[i]=0.5*(ym[i]+yn[i]+h*yend[i]);
	};

	void step(const double htry, D *derivs_class)
	{
		const double STEPFAC1=0.65, STEPFAC2=0.94, STEPFAC3=0.02, STEPFAC4=4.0, KFAC1=0.8, KFAC2=0.9;
		bool  first_step=true; 
		bool  last_step=false;
		bool  forward;
		bool  reject=false;
		bool  prev_reject=false;
	
		size_t  i,k;
		double  fac, h, hnew, hopt_int, err;
		bool  firstk;
	
		Eigen::VectorXd  hopt(IMAXX),work(IMAXX);
		Eigen::VectorXd  ysav(n),yseq(n);
		Eigen::VectorXd  ymid(n),scale(n);
	
		work[0]=0;
		h=htry;
	
		if ( h>0 )
			forward = true;
		else
			forward = false;

		for (i=0;i<n;i++) 
			ysav[i]=y[i]; // Save the starting values.

		if (h != hnext && !first_step )
		{
			// h gets reset in Odeint for the last step.
			last_step = true;
		}
	
		if (reject) //Previous step was rejected.
		{
			prev_reject=true;
			last_step=false;
		}
	
		reject=false;
		firstk=true;
		hnew=fabs(h);
    
	interp_error: // Restart here if interpolation error too big.

		while  (firstk || reject) // Loop until step accepted
		{
			if ( forward )
				h = hnew;
			else
				h = -1.0*hnew;

			firstk=false;
			reject=false;


			// if ( fabs(h) <= fabs(x)*EPS )
			// 	throw( "step size underflow  in  StepperBS");
		

			int ipt=-1; // Initialize counter for saving stu.

			for (k=0; k<=k_targ+1; k++) // Evaluate the sequence of modied midpoint integrations.
			{
				dy(ysav,h,k,yseq,ipt, derivs_class);
				if (k == 0)
					y=yseq;
				else // Store result in tableau.
					for (i=0;i<n;i++)
						table(k-1,i)=yseq[i];
				    // debug
				    // std::cout << table; 
					if (k  !=  0)  
					{
						polyextr(k,table,y); // Perform extrapolation.
						err=0.0; // Compute normalized error estimate errk.
					
						for (i=0;i<n;i++)
						{
							scale[i]=atol+rtol*MAX(fabs(ysav[i]),fabs(y[i]));
							err+=SQR((y[i]-table(0,i))/scale[i]);
						}
					
						err=sqrt(err/n);
						double  expo=1.0/(2*k+1); // Compute optimal step size for this order.
						double  facmin=pow(STEPFAC3,expo);
					
						if (err  ==  0.0)
							fac=1.0/facmin;
						else
						{
							fac=STEPFAC2/pow(err/STEPFAC1,expo);
							fac=MAX(facmin/STEPFAC4, MIN(1.0/facmin,fac));
						}
					
						hopt[k]=fabs(h*fac);
						work[k]=cost[k]/hopt[k];  // Work per unit step (17.3.13).
					
						if ((first_step ||  last_step)  &&  err  <=  1.0)
							break;
						if  (k  ==  k_targ-1  &&  !prev_reject &&  !first_step  &&  !last_step)  
						{
							if  (err  <=  1.0) // Converged within order window.
								break;
							else if (err>SQR(nseq[k_targ]*nseq[k_targ+1]/(nseq[0]*nseq[0])))
							{
								reject=true; // Criterion (17.3.17) predicts step will fail.
								k_targ=k;
								if  (k_targ>1  &&  work[k-1]<KFAC1*work[k])
									k_targ--;
							
								hnew=hopt[k_targ];
								break;
							}
						}
					
						if  (k == k_targ)  {
							if  (err  <=  1.0)
								break; // Converged within order window.
							else if (err>SQR(nseq[k+1]/nseq[0])) 
							{
								reject=true; // Criterion (17.3.20) predicts step will fail.
								if  (k_targ>1  &&  work[k-1]<KFAC1*work[k])
									k_targ--;
								hnew=hopt[k_targ];
								break;
							}
						}
					
						if  ( k == k_targ+1 )
						{
							if  ( err > 1.0 )
							{
								reject=true;
								if  (k_targ>1  &&  work[k_targ-1]<KFAC1*work[k_targ])
									k_targ--;
								hnew=hopt[k_targ];
							}
							break;
						}
					}
			} // Go back and try next k.
		
			if (reject) // Arrive here from any break in for loop.
				prev_reject=true;
		} // Go back if step was rejected.
	
		derivs_class->eval(x+h,y,dydxnew); // Used for start of next step and in dense out-put.
		if (dense)
		{
			prepare_dense(h,dydxnew,ysav,scale,k,err);
			hopt_int=h/MAX(pow(err,1.0/(2*k+3)),0.01);
			// Stepsize based on interpolation error. 
			if (err > 10.0)  // Interpolation error too big, reject step.
			{
				hnew=fabs(hopt_int);
				reject=true;
				prev_reject=true;
				goto  interp_error;
			}
		}

		dydx=dydxnew; // Update for start of next step.
		xold=x; // For dense output.
		x+=h;
		hdid=h;
		first_step=false;
	
		size_t kopt; // Determine optimal order for next step.
	
		if (k == 1)
			kopt=2;
		else if (k  <=  k_targ)
		{
			kopt=k;
			if  (work[k-1] <  KFAC1*work[k])
				kopt=k-1;
			else if (work[k]  <  KFAC2*work[k-1])
				kopt=MIN(k+1,KMAXX-1);
		}
		else
		{
			kopt=k-1;
			if ( k > 2  &&  work[k-2] < KFAC1*work[k-1])
				kopt=k-2;
			if ( work[k] < KFAC2*work[kopt])
				kopt=MIN(k,KMAXX-1);
		}
	
		if  (prev_reject) // After a rejected step neither order nor step-
		{ 
			// size should increase. 
			k_targ=MIN(kopt,k);
			hnew=MIN((double)fabs(h),(double)hopt[k_targ]);
			prev_reject=false;
		}
		else // Stepsize control for next step.
		{
			if  (kopt  <=  k)
				hnew=hopt[kopt];
			else
			{
				if  (k<k_targ  &&  work[k]<KFAC2*work[k-1])
					hnew=hopt[k]*cost[kopt+1]/cost[k];
				else
					hnew=hopt[k]*cost[kopt]/cost[k];
			}
			k_targ=kopt;
		}
	
		if (dense) // Keep interpolation error small enough.
		  hnew=MIN((double)hnew,(double)fabs(hopt_int));

		if (forward)
			hnext=hnew;
		else
			hnext=-hnew;
	}; 

	void polyextr(const size_t k, 
		          Eigen::MatrixXd &table, 
				  Eigen::VectorXd &last)
	{
		size_t l=last.size();
		size_t i,j;
	
		for(j=k-1; j>0; j--) // Update the current row using the Neville recursive formula. 
			for  (i=0; i<l; i++)
				table(j-1,i)=table(j,i)+coeff(k,j)*(table(j,i)-table(j-1,i));
		for (i=0; i<l; i++) // Update the last element.
			last[i]=table(0,i)+coeff(k,0)*(table(0,i)-last[i]);
	};

	void prepare_dense(const double h, 
		               Eigen::VectorXd &dydxnew, 
					   Eigen::VectorXd &ysav, 
					   Eigen::VectorXd &scale, 
					   const size_t k, 
					   double &error)
	{
		size_t i,j,l,kmi,kk;

		mu=2*k-1; //Degree of interpolating polynomial is mu C 4.
	
		for (i=0; i<n; i++)  //Store y and y 0 at both ends of interval.
		{
			dens[i]=ysav[i];
			dens[n+i]=h*dydx[i];
			dens[2*n+i]=y[i];
			dens[3*n+i]=dydxnew[i]*h;
		}

		for ( j=1; j<=k; j++)
		{   //Compute solution at midpoint.
			double dblenj=nseq[j];
			for (l=j; l>=1; l--)
			{
				double factor=SQR(dblenj/nseq[l-1])-1.0;
				for  (size_t  i=0;  i<n;  i++)
					ysave(l-1,i)=ysave(l,i)+(ysave(l,i)-ysave(l-1,i))/factor;
			}
		}

		for ( i=0; i<n; i++)
			dens[4*n+i]=ysave(0,i);

		double facnj;
		double dblenj;

		for ( kmi=1;  kmi<=mu;  kmi++)
		{   //Compute kmi-th derivative at midpoint.
			size_t kbeg=(kmi-1)/2;
			for ( kk=kbeg;  kk<=k;  kk++)
			{
				facnj=pow((double)(nseq[kk]/2.0),(int)(kmi-1));
				size_t ipt=ipoint[kk+1]-2*kk+kmi-3;
				for  ( i=0;  i<n;  i++)
					ysave(kk,i)=fsave(ipt,i)*facnj;
			}

			for ( j=kbeg+1;  j<=k;  j++)  {
				dblenj=nseq[j];
				for  ( l=j; l>=kbeg+1; l--)
				{
					double  factor=SQR(dblenj/nseq[l-1])-1.0;
					for  ( i=0;  i<n;  i++)
						ysave(l-1,i)=ysave(l,i)+(ysave(l,i)-ysave(l-1,i))/factor;
				}
			}

			for ( i=0; i<n; i++)
				dens[(kmi+4)*n+i]=ysave(kbeg,i)*h;
    
			if (kmi  ==  mu)
				continue;

			size_t lbeg;
			size_t lend;

			for ( kk=kmi/2;  kk<=k;  kk++)
			{   
				//Compute dierences.
				lbeg=ipoint[kk+1]-1;
				lend=ipoint[kk]+kmi;
				if (kmi  ==  1)
					lend  +=  2;
				for ( l=lbeg;  l>=lend;  l-=2)
					for  ( i=0;  i<n;  i++)
						fsave(l,i)=fsave(l,i)-fsave(l-2,i);
		       
				if (kmi == 1)
				{
					l=lend-2;
					for  ( i=0;  i<n;  i++)
						fsave(l,i)=fsave(l,i)-dydx[i];
				}
			}

			for  ( kk=kmi/2;  kk<=k;  kk++) 
			{
				lbeg=ipoint[kk+1]-2;
				lend=ipoint[kk]+kmi+1;
				for  ( l=lbeg;  l>=lend;  l-=2)
					for  ( i=0;  i<n;  i++)
						fsave(l,i)=fsave(l,i)-fsave(l-2,i);
			}

		}

		dense_interp(n,dens,mu); // Compute the interpolation coecients in dens.
		error=0.0; // Estimate the interpolation error.
	
		if ( mu >= 1)
		{
			for  ( i=0;  i<n;  i++)
				error += SQR(dens[(mu+4)*n+i]/scale[i]);
			error=sqrt(error/n)*errfac[mu-1];
		}

	}; 

	double dense_out(const size_t i, 
		             const double x, 
					 const double h)
	{
		double theta=(x-xold)/h;
		double theta1=1.0-theta;
		double yinterp=dens[i]+theta*(dens[n+i]+theta1*(dens[2*n+i]*theta + dens[3*n+i]*theta1));

		// std::cout << dens;
		if (mu<0)
			return yinterp;

		double  theta05=theta-0.5;
		double  t4=SQR(theta*theta1);
		double  c=dens[n*(mu+4)+i];
		for (size_t j=mu; j>0; j--)
			c=dens[n*(j+3)+i]+c*theta05/j;
		yinterp += t4*c;
	
	return  yinterp;
	}; 

	void dense_interp(const size_t n, 
		              Eigen::VectorXd &y, 
					  const size_t imit)
	{
		double  y0,y1,yp0,yp1,ydiff,aspl,bspl,ph0,ph1,ph2,ph3,fac1,fac2;
		Eigen::VectorXd a(31);
		for (size_t i=0; i<n; i++)
		{
			y0=y[i];
			y1=y[2*n+i];
			yp0=y[n+i];
			yp1=y[3*n+i];
			ydiff=y1-y0;
			aspl=-yp1+ydiff;
			bspl=yp0-ydiff;
			y[n+i]=ydiff;
			y[2*n+i]=aspl;
			y[3*n+i]=bspl;
		
			if  (imit  <  0) 
				continue;
			ph0=(y0+y1)*0.5+0.125*(aspl+bspl);
			ph1=ydiff+(aspl-bspl)*0.25;
			ph2=-(yp0-yp1);
			ph3=6.0*(bspl-aspl);
			if  ( imit >= 1 )
			{
				a[1]=16.0*(y[5*n+i]-ph1);
				if  (imit  >=  3)
				{
					a[3]=16.0*(y[7*n+i]-ph3+3*a[1]);
					for  (size_t  im=5;  im  <=imit;  im+=2) 
					{
						fac1=im*(im-1)/2.0;
						fac2=fac1*(im-2)*(im-3)*2.0;
						a[im]=16.0*(y[(im+4)*n+i]+fac1*a[im-2]-fac2*a[im-4]);
					}
				}
			}
		
			a[0]=(y[4*n+i]-ph0)*16.0;
		
			if  (imit  >=  2) 
			{
				a[2]=(y[n*6+i]-ph2+a[0])*16.0;
				for  (size_t  im=4;  im  <=imit;  im+=2) 
				{
					fac1=im*(im-1)/2.0;
					fac2=im*(im-1)*(im-2)*(im-3);
					a[im]=(y[n*(im+4)+i]+a[im-2]*fac1-a[im-4]*fac2)*16.0;
				}
			}
		
			for  (size_t  im=0;  im<=imit;  im++)
				y[n*(im+4)+i]=a[im];
		}

	}; 

	Eigen::VectorXd& get_y(void){return y;};
	Eigen::VectorXd& get_dydx(void){return dydx;};
	void set_y(Eigen::VectorXd& new_y){y = new_y;};
	void set_dydx(Eigen::VectorXd& new_dydx){dydx = new_dydx;};

	const size_t     KMAXX, IMAXX;   // KMAXX is the maximum number of rows used in the extrapolation
	size_t           k_targ;         // optimal row number for convergence
	Eigen::VectorXi  nseq;           // stepsize sequence
	Eigen::VectorXi  cost;           // A_k
	Eigen::MatrixXd  table;          // extrapolation tableau
	Eigen::VectorXd  dydxnew;        
	size_t           mu;             // used for dense output
	Eigen::MatrixXd  coeff;          // coefficients used in extrapolation tableau 
	Eigen::VectorXd  errfac;         // used to compute dense interpolation error
	Eigen::MatrixXd  ysave;          // ysave and fsave store values and derivatives
	Eigen::MatrixXd  fsave;          // to be used for dense output
	Eigen::VectorXi  ipoint;         // keep track of where values are storeed in fsave
	Eigen::VectorXd  dens;           // store quantities for dense interpolating polynomial

};

#endif
