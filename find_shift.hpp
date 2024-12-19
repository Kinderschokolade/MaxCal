#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <string> 
#include <math.h>
#include <float.h>
//#include <boost/math/special_functions/expint.hpp>

using namespace std;




double u(double & x,double & beta, double & U , double & k)
{
	return (-beta * U / (1.+exp(-2.*k*x)));
} 

//double g(double x,double & beta, double & U , double & k)
//{

	//cout<<x <<"\t"<< beta << "\t" << U << "\t" << k << endl;
	//cout<<beta << "\t" <<  u(x,beta,U,k)<< "\t" << beta*U << endl; 

//	return(1./(2.*k) *boost::math::expint(u(x,beta,U,k)) - exp(-beta*U) / (2.*k)* boost::math::expint(u(x,beta,U,k) + beta*U ));
//}


double func1(double & x, double & beta1, double & beta2, double & U1, double & U2, double & k1 ,double & k2, double & c)
{
	double out;
	if (x > 0)
	{
		out = - g(x,beta2,U2,k2) + g(c,beta2,U2,k2) - exp(-beta2*U2) * (-x+c);
	}
	else
	{
		out = -g(0.,beta2,U2,k2) + g(c,beta2,U2,k2) - exp(-beta2*U2) * c + g(-x,beta1,U1,k1) - g(0.,beta1,U1,k1) + exp(-beta1*U1) *x;	
	}
	return out;
}


double func2(double & x, double & beta1, double & beta2, double & U1, double & U2, double & k1 ,double & k2, double & c)
{
	double out;
	if (x < 0)
	{
		out = - g(-x,beta1,U1,k1) + g(c,beta1,U1,k1) - exp(-beta1*U1) * (x+c);
	}
	else
	{
		out = -g(0.,beta1,U1,k1) + g(c,beta1,U1,k1) - exp(-beta1*U1) * c + g(x,beta2,U2,k2) - g(0.,beta2,U2,k2) - exp(-beta2*U2) *x;	
	}
	return out;
}

double f(double & x, double & beta1, double & beta2, double & U1, double & U2, double & k1 ,double & k2, double & c)
{
	return (func1(x,beta1,beta2,U1,U2,k1,k2,c) - func2(x,beta1,beta2,U1,U2,k1,k2,c));
}





double findroot(double & beta1, double & beta2, double & U1, double & U2, double & k1 ,double & k2, double & c)
{
	double x1=-c,x2=c,x3;
	double precision = 0.00001;
	int count = 0; 
	int maxIt = pow(10,5);

	if (fabs(beta1 * U1 - beta2*U2) < precision) {x1 =0; x2 =0; x3 =0;}
		//algorithm fails for same exponents, solve right here. 
	
    	while (fabs(x1 - x2) > precision )
    	{
 
	       	 x3 = (x1 + x2)/2;
			//stuck for x3 = 0 
	
        	if( f(x1,beta1,beta2,U1,U2,k1,k2,c) * f(x3,beta1,beta2,U1,U2,k1,k2,c) < 0 )
        	{
            		x2 = x3;
        	}
        	else
        	{
            		x1 = x3;
       		}
        	count++;
		if (count > maxIt) {break; cout<< "root finding failed"<< endl;}
    	}
	return x3;
}

