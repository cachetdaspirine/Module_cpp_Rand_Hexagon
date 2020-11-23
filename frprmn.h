#ifndef FRPRMN_H_
#define FRPRMN_H_

#include "dlinmin.h"

template <class T>
struct Frprmn : Dlinemethod<T> {
//	Multidimensional minimization by the Fletcher-Reeves-Polak-Ribiere method.
	Int iter;
	Doub fret; //Value of the function at the minimum.
	using Dlinemethod<T>::func; //Variables from a templated base class
	using Dlinemethod<T>::linmin; //are not automatically inherited.
	using Dlinemethod<T>::p;
	using Dlinemethod<T>::xi;
	const Doub ftol;
 Frprmn(T &funcd, const Doub ftoll=1e-10) : Dlinemethod<T>(funcd),
	  ftol(ftoll) {}
//	Constructor arguments are funcd, the function or functor to be minimized, and an optional
//	argument ftoll, the fractional tolerance in the function value such that failure to decrease
//	by more than this amount on one iteration signals doneness.
	VecDoub minimize(VecDoub_I &pp)
//	Given a starting point pp[0..n-1], performs the minimization on a function whose value
//	and gradient are provided by a functor funcd (see text).
	{
		const Int ITMAX=1000000000;
		const Doub EPS=1.0e-10;//18;
		const Doub GTOL=1.0e-8;
				
//		Here ITMAX is the maximum allowed number of iterations; EPS is a small number to
//		rectify the special case of converging to exactly zero function value; and GTOL is the
//		convergence criterion for the zero gradient test.
		Doub gg,dgg;
		Int n=pp.size(); //Initializations.
		p=pp;
		
		VecDoub g(n),h(n);
		xi.resize(n);
		Doub fp=func(p);
		func.df(p,xi);
		for (Int j=0;j<n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j];
		}
		for (Int its=0;its<ITMAX;its++) { //Loop over iterations.

			iter=its;
			
			fret=linmin(); //Next statement is one possible return:
			/////
	//		cout << "the function value at iteration "<<its<<" = " <<fret<<"\n";
			////
						
			if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) {

			return p;}
				
			fp=fret;
						
			func.df(p,xi);
			Doub test=0.0; //Test for convergence on zero gradient.
				Doub den=MAX(fp,1.0);
				for (Int j=0;j<n;j++) {
					Doub temp=abs(xi[j])*MAX(abs(p[j]),1.0)/den;
					if (temp > test) test=temp;
						}
			if (test < GTOL){ 				

			return p;} //The other possible return.
			dgg=gg=0.0;
			for (Int j=0;j<n;j++) {
				gg += g[j]*g[j];
				
				// dgg += xi[j]*xi[j]; //This statement for Fletcher-Reeves.
				dgg += (xi[j]+g[j])*xi[j]; //This statement for Polak-Ribiere.
					}
			if (gg == 0.0){// Unlikely. If gradient is exactly zero, then

			return p;} //we are already done.
			Doub gam=dgg/gg;

			for (Int j=0;j<n;j++) {
				g[j] = -xi[j];
				xi[j]=h[j]=g[j]+gam*h[j];
			}

		}
		throw("Too many iterations in frprmn");
		//cout<<"CG not converged!!!\n";
	}
};

#endif /*FRPRMN_H_*/
