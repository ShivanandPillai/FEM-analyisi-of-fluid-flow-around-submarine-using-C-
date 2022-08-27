 //OE21M027//
//Term Paper//

#include <iostream>
#include <fstream>
#include <cmath>

#define N 40000

using namespace std ;

double co[N][2] , u[N] , v[N] , P[N] ;
double m , a1 , a2 , b1 , b2 , c1 , c2 , x1 , x2 , y3 , y2 , u1 , u2 , v1 , v2 , K11 , K12 , K13 , K14 , K22 , K23 , K24 , K33 , K34 , K44 , X1 , X2 , X3 , X4 ;
double A1 , A2 , B1 , B2 , C1 , C2 , D1 , D2 , E1 , E2 , F1 , F2 , e1 , e2 , f1 , f2 , g1 , g2 ;
double rho = 1025 ; 
int p[N] , r , s , i , j , k , l , n ;

struct mesh
{
	double x1 , y1 , x2 , y2 , x3 , y3 ;
};

double find_coefficients(double x1 , double y1 , double x2 , double y2 ,double x3 , double y3 , double u1 , double v1 , double u2 , double v2)
{
	if(p[s]%2 == 1)
	{
		m = (x3-x1)/(y3-y1) ;
	    K11 = 0.5*pow(x1,2)*(y1-y3)+pow(x3,3)/(6*m) ;
	    K12 = 0.25*pow(x1,2)*(pow(y1,2)-pow(y3,2))+(y3*pow(x3,3))/(6*m)-pow(x3,4)/(24*pow(m,2)) ;
	    K13 = K11 ;
	    K14 = (1/3)*pow(x1,3)*(y1-y3)+pow(x3,4)/(12*m) ;
	    K22 = x1*y1-x1*(K12/K11) ;
	    K23 = -x1*(K13/K11) ;
	    K24 = -x1*(K14/K11) ;
	    K33 = K13*((x2*(y2-y1))/(y1*K11-K12)) ;
	    K34 = K14*((x2*(y2-y1))/(y1*K11-K12)) ;
	    K44 = x1*y1-y1*(K14/K13) ;
	    X1 = 0 ;
	    X2 = u1 ;
	    X3 = u2-u1*((x2*y2*K11-x2*K12)/(x1*y1*K11-x1*K12)) ;
	    if(x1==0)
	    {
		    X3 = u2-u1*((x2*y2*K11-x2*K12)/(x1*y1*K11-(x1-0.00001)*K12)) ;
	    }
	    X4 = v1-y1*((u2*(x1*y1*K11-x1*K12)-u1*(x2*y2*K11-x2*K12))/(x1*x2*K13*(y2-y1))) ;
	    if((x1==0)||(x2==0)||((y2-y1)==0))
	    {
		    X4 = v1-y1*((u2*(x1*y1*K11-x1*K12)-u1*(x2*y2*K11-x2*K12))/((x1+0.00001)*(x2+0.00001)*K13*(y2-y1+0.00001))) ;
	    }
	    b2 = X4/K44 ;
	    if(K44==0)
	    {
		    b2 = X4/(K44+0.00001) ;
	    }
	    b1 = (X3-K34*b2)/K33 ;
	    if(K33==0)
	    {
	 	    b1 = (X3-K34*b2)/(K33+0.00001) ;
	    }
	    a2 = (X2-K24*b2-K23*b1)/K22 ;
	    if(K22==0)
	    {
		    a2 = (X2-K24*b2-K23*b1)/(K22+0.00001) ;
	    }
	    a1 = (X1-K14*b2-K13*b1-K12*a2)/K11 ;
	    
		A1 = (pow(a1,2)/3)*(pow(x1,2)*(y1-y3)+pow(x3,4)/(4*m)) ;
	    B1 = (pow(a2,2)/3)*((pow(x1,3)/3)*(pow(y1,3)-pow(y3,3))+(pow(y3,2)*pow(x3,4))/(4*m)-(y3*pow(x3,5))/(10*pow(m,2))+pow(x3,6)/(60*pow(m,3))) ;
	    C1 = ((2*a1*a2)/3)*((pow(x1,3)/2)*(pow(y1,2)-pow(y3,2))+(y3*pow(x3,4))/(4*m)-pow(x3,5)/(20*pow(m,2))) ;
	    D1 = (b1/(2*a1))*C1 ;
	    if(a1==0)
	    {
	    	D1 = (b1/(2*(a1+0.00001)))*C1 ;
		}
	    E1 = 0.25*a2*b2*((pow(x1,4)/2)*(pow(y1,2)-pow(y3,2))+(y3*pow(x3,5))/(5*m)-pow(x3,6)/(30*pow(m,2))) ;
	    F1 = (2/(rho*pow(a1,2)))*A1 ;
	    if(a1==0)
	    {
	    	F1 = (2/(rho*pow(a1+0.00001,2)))*A1 ;
		}
	    A2 = m*pow(b1,2)*((y1/3)*(pow(y1,3)-pow(y3,3))+0.25*(pow(y1,4)-pow(y3,4))) ;
	    B2 = (pow(b2,2)/pow(a2,2))*B1 ;
	    if(a2==0)
	    {
	    	B2 = (pow(b2,2)/pow(a2+0.00001,2))*B1 ;
		}
	    C2 = b1*b2*((pow(x1,2)/3)*(pow(y1,3)-pow(y3,3))+(pow(y3,2)*pow(x3,3))/(3*m)-(y3*pow(x3,4))/(6*pow(m,2))+pow(x3,5)/(30*pow(m,3))) ;
	    D2 = 0.5*a1*b2*((pow(x1,2)/2)*(pow(y1,3)-pow(y3,3))+(pow(y3,2)*pow(x3,3))/(3*m)-(y3*pow(x3,4))/(6*pow(m,2))+pow(x3,5)/(30*pow(m,3))) ;
	    E2 = 0.5*a2*b2*((pow(x1,2)/4)*(pow(y1,4)-pow(y3,4))+(pow(y3,4)*pow(x3,2))/4-m*x3*(pow(y3,5)/10)+(pow(m,2)/60)*(pow(y3,6)-pow(y1,6))) ;
	    F2 = (2/(rho*pow(b1,2)))*A2 ;
	    if(b1==0)
	    {
	    	F2 = (2/(rho*pow(b1+0.00001,2)))*A2 ;
		}
	    c1 = -(A1+B1+C1+D1+E1)/F1 ;
	    if(F1==0)
	    {
	    	c1 = -(A1+B1+C1+D1+E1)/(F1+0.00001) ;
		}
	    c2 = -(A2+B2+C2+D2+E2)/F2 ;
	    if(F2==0)
	    {
	    	c2 = -(A2+B2+C2+D2+E2)/(F2+0.00001) ;
		}
	}
	else
	{
		u2 = u[k+3] ;
		v2 = v[k+3] ;
		m = (x3-x1)/(y3-y1) ;
	    K11 = 0.5*pow(x3,2)*(y3-y1)+pow(x1,3)/(6*m) ;
	    K12 = 0.25*pow(x3,2)*(pow(y3,2)-pow(y1,2))+(y1*pow(x1,3))/(6*m)-pow(x1,4)/(24*pow(m,2)) ;
	    K13 = K11 ;
	    K14 = (1/3)*pow(x3,3)*(y3-y1)+pow(x1,4)/(12*m) ;
	    K22 = x1*y1-x1*(K12/K11) ;
	    K23 = -x1*(K13/K11) ;
	    K24 = -x1*(K14/K11) ;
	    K33 = K13*((x2*(y2-y1))/(y1*K11-K12)) ;
	    K34 = K14*((x2*(y2-y1))/(y1*K11-K12)) ;
	    K44 = x1*y1-y1*(K14/K13) ;
	    X1 = 0 ;
	    X2 = u1 ;
	    X3 = u2-u1*((x2*y2*K11-x2*K12)/(x1*y1*K11-x1*K12)) ;
	    if(x1==0)
	    {
		    X3 = u2-u1*((x2*y2*K11-x2*K12)/(x1*y1*K11-(x1-0.00001)*K12)) ;
	    }
	    X4 = v1-y1*((u2*(x1*y1*K11-x1*K12)-u1*(x2*y2*K11-x2*K12))/(x1*x2*K13*(y2-y1))) ;
	    if((x1==0)||(x2==0)||((y2-y1)==0))
	    {
		    X4 = v1-y1*((u2*(x1*y1*K11-x1*K12)-u1*(x2*y2*K11-x2*K12))/((x1+0.00001)*(x2+0.00001)*K13*(y2-y1+0.00001))) ;
	    }
	    b2 = X4/K44 ;
	    if(K44==0)
	    {
		    b2 = X4/(K44+0.00001) ;
	    }
	    b1 = (X3-K34*b2)/K33 ;
	    if(K33==0)
	    {
	 	    b1 = (X3-K34*b2)/(K33+0.00001) ;
	    }
	    a2 = (X2-K24*b2-K23*b1)/K22 ;
	    if(K22==0)
	    {
		    a2 = (X2-K24*b2-K23*b1)/(K22+0.00001) ;
	    }
	    a1 = (X1-K14*b2-K13*b1-K12*a2)/K11 ;
	    
	    A1 = (pow(a1,2)/3)*(pow(x3,2)*(y3-y1)+pow(x1,4)/(4*m)) ;
	    B1 = (pow(a2,2)/3)*((pow(x3,3)/3)*(pow(y3,3)-pow(y1,3))+(pow(y1,2)*pow(x1,4))/(4*m)-(y1*pow(x1,5))/(10*pow(m,2))+pow(x1,6)/(60*pow(m,3))) ;
	    C1 = ((2*a1*a2)/3)*((pow(x3,3)/2)*(pow(y3,2)-pow(y1,2))+(y1*pow(x1,4))/(4*m)-pow(x1,5)/(20*pow(m,2))) ;
	    D1 = (b1/(2*a1))*C1 ;
	    if(a1==0)
	    {
	    	D1 = (b1/(2*(a1+0.00001)))*C1 ;
		}
	    E1 = 0.25*a2*b2*((pow(x3,4)/2)*(pow(y3,2)-pow(y1,2))+(y1*pow(x1,5))/(5*m)-pow(x1,6)/(30*pow(m,2))) ;
	    F1 = (2/(rho*pow(a1,2)))*A1 ;
	    if(a1==0)
	    {
	    	F1 = (2/(rho*pow(a1+0.00001,2)))*A1 ;
		}
	    A2 = m*pow(b1,2)*((y3/3)*(pow(y3,3)-pow(y1,3))+0.25*(pow(y3,4)-pow(y1,4))) ;
	    B2 = (pow(b2,2)/pow(a2,2))*B1 ;
	    if(a2==0)
	    {
	    	B2 = (pow(b2,2)/pow(a2+0.00001,2))*B1 ;
		}
	    C2 = b1*b2*((pow(x3,2)/3)*(pow(y3,3)-pow(y1,3))+(pow(y1,2)*pow(x1,3))/(3*m)-(y1*pow(x1,4))/(6*pow(m,2))+pow(x1,5)/(30*pow(m,3))) ;
	    D2 = 0.5*a1*b2*((pow(x3,2)/2)*(pow(y3,3)-pow(y1,3))+(pow(y1,2)*pow(x1,3))/(3*m)-(y1*pow(x1,4))/(6*pow(m,2))+pow(x1,5)/(30*pow(m,3))) ;
	    E2 = 0.5*a2*b2*((pow(x3,2)/4)*(pow(y3,4)-pow(y1,4))+(pow(y1,4)*pow(x1,2))/4-m*x1*(pow(y1,5)/10)+(pow(m,2)/60)*(pow(y1,6)-pow(y3,6))) ;
	    F2 = (2/(rho*pow(b1,2)))*A2 ;
	    if(b1==0)
	    {
	    	F2 = (2/(rho*pow(b1+0.00001,2)))*A2 ;
		}
	    c1 = -(A1+B1+C1+D1+E1)/F1 ;
	    if(F1==0)
	    {
	    	c1 = -(A1+B1+C1+D1+E1)/(F1+0.00001) ;
		}
	    c2 = -(A2+B2+C2+D2+E2)/F2 ;
	    if(F2==0)
	    {
	    	c2 = -(A2+B2+C2+D2+E2)/(F2+0.00001) ;
		}
	}
}

    
int main()
{
	struct mesh element[N] ;
	double dx = 25  ;
	double dy = 25 ;
	double rel_speed = 10 ;
	double x = 0 ;
	double y = 0 ;
	double x_nose = 0 ;
	double x_tail = 0 ;
	i = 0;
	l = 0;
	while(y < 222)
	{
		if(y < 100)
		{
			while(x <= 500)
			{
					
			    element[i].x1 = x ;
	            element[i].y1 = y ;
	            element[i].x2 = element[i].x1 + ((106-y)/105)*dx ;
	            element[i].y2 = element[i].y1 ;
	            element[i].x3 = element[i].x2 ;
	            element[i].y3 = element[i].y2 + ((106-y)/105)*dy ;
	            element[i+1].x1 = x ;
	            element[i+1].y1 = y ;
	            element[i+1].x2 = element[i+1].x1 ;
	            element[i+1].y2 = element[i+1].y1 + ((106-y)/105)*dy ;
	            element[i+1].x3 = element[i+1].x2 + ((106-y)/105)*dx ;
	            element[i+1].y3 = element[i+1].y2 ;
				x = element[i].x2 ;
		        y = element[i].y2 ;
		        i += 2 ; 
			}
		}
		else if((y>=100)&&(y<110))
		{
			if(y<=105)
			{
				x_nose = -sqrt(100-4*pow((y-105),2)) + 210 ;
			    while(x <= x_nose)
			    {
			        element[i].x1 = x ;
	                element[i].y1 = y ;
	                element[i].x2 = element[i].x1 + ((106-y)/105)*dx ;
	                element[i].y2 = element[i].y1 ;
	                element[i].x3 = element[i].x2 ;
	                element[i].y3 = element[i].y2 + ((106-y)/105)*dy ;
	                element[i+1].x1 = x ;
	                element[i+1].y1 = y ;
	                element[i+1].x2 = element[i+1].x1 ;
	                element[i+1].y2 = element[i+1].y1 + ((106-y)/105)*dy ;
	                element[i+1].x3 = element[i+1].x2 + ((106-y)/105)*dx ;
	                element[i+1].y3 = element[i+1].y2 ;
				    x = element[i].x2 ;
		            y = element[i].y2 ;
		            i += 2 ; 
			    }
			
				x_tail = 4*y - 140 ;
			    x = x_tail ; 
			    while(x <= 270)
			    {
			        element[i].x1 = x ;
	                element[i].y1 = y ;
	                element[i].x2 = element[i].x1 + ((106-y)/105)*dx ;
	                element[i].y2 = element[i].y1 ;
	                element[i].x3 = element[i].x2 ;
	                element[i].y3 = element[i].y2 + ((106-y)/105)*dy ;
	                element[i+1].x1 = x ;
	                element[i+1].y1 = y ;
	                element[i+1].x2 = element[i+1].x1 ;
	                element[i+1].y2 = element[i+1].y1 + ((106-y)/105)*dy ;
	                element[i+1].x3 = element[i+1].x2 + ((106-y)/105)*dx ;
	                element[i+1].y3 = element[i+1].y2 ;
				    x = element[i].x2 ;
		            y = element[i].y2 ;
		            i += 2 ;  
			    }
			    
			    if(y<=103.75)
			    {
			    	x = 275 ; 
			        while(x <= 500)
			        {
			            element[i].x1 = x ;
	                    element[i].y1 = y ;
	                    element[i].x2 = element[i].x1 + ((106-y)/105)*dx ;
	                    element[i].y2 = element[i].y1 ;
	                    element[i].x3 = element[i].x2 ;
	                    element[i].y3 = element[i].y2 + ((106-y)/105)*dy ;
	                    element[i+1].x1 = x ;
	                    element[i+1].y1 = y ;
	                    element[i+1].x2 = element[i+1].x1 ;
	                    element[i+1].y2 = element[i+1].y1 + ((106-y)/105)*dy ;
	                    element[i+1].x3 = element[i+1].x2 + ((106-y)/105)*dx ;
	                    element[i+1].y3 = element[i+1].y2 ;
				        x = element[i].x2 ;
		                y = element[i].y2 ;
		                i += 2 ; 
			        }
				}
				else if(y>103.75)
				{
					x_tail = 4*y - 140 ;
			        x = x_tail ; 
			        while(x <= 500)
			        {
			            element[i].x1 = x ;
	                    element[i].y1 = y ;
	                    element[i].x2 = element[i].x1 + ((106-y)/105)*dx ;
	                    element[i].y2 = element[i].y1 ;
	                    element[i].x3 = element[i].x2 ;
	                    element[i].y3 = element[i].y2 + ((106-y)/105)*dy ;
	                    element[i+1].x1 = x ;
	                    element[i+1].y1 = y ;
	                    element[i+1].x2 = element[i+1].x1 ;
	                    element[i+1].y2 = element[i+1].y1 + ((106-y)/105)*dy ;
	                    element[i+1].x3 = element[i+1].x2 + ((106-y)/105)*dx ;
	                    element[i+1].y3 = element[i+1].y2 ;
				        x = element[i].x2 ;
		                y = element[i].y2 ;
		                i += 2 ; 
			        }
				}
			}
			else if(y>105)
			{
				x_nose = -sqrt(100-4*pow((y-105),2)) + 210 ;
			    while(x <= x_nose)
			    {
			        element[i].x1 = x ;
	                element[i].y1 = y ;
	                element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	                element[i].y2 = element[i].y1 ;
	                element[i].x3 = element[i].x2 ;
	                element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	                element[i+1].x1 = x ;
	                element[i+1].y1 = y ;
	                element[i+1].x2 = element[i+1].x1 ;
	                element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	                element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	                element[i+1].y3 = element[i+1].y2 ;
				    x = element[i].x2 ;
		            y = element[i].y2 ;
		            i += 2 ; 
			    }
			
				x_tail = 700 - 4*y ;
			    x = x_tail ; 
			    while(x <= 270)
			    {
			        element[i].x1 = x ;
	                element[i].y1 = y ;
	                element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	                element[i].y2 = element[i].y1 ;
	                element[i].x3 = element[i].x2 ;
	                element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	                element[i+1].x1 = x ;
	                element[i+1].y1 = y ;
	                element[i+1].x2 = element[i+1].x1 ;
	                element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	                element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	                element[i+1].y3 = element[i+1].y2 ;
				    x = element[i].x2 ;
		            y = element[i].y2 ;
		            i += 2 ;
			    }
			    
			    if(y>106.25)
			    {
			    	x = 275 ; 
			        while(x <= 500)
			        {
			            element[i].x1 = x ;
	                    element[i].y1 = y ;
	                    element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	                    element[i].y2 = element[i].y1 ;
	                    element[i].x3 = element[i].x2 ;
	                    element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	                    element[i+1].x1 = x ;
	                    element[i+1].y1 = y ;
	                    element[i+1].x2 = element[i+1].x1 ;
	                    element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	                    element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	                    element[i+1].y3 = element[i+1].y2 ;
				        x = element[i].x2 ;
		                y = element[i].y2 ;
		                i += 2 ;
			        }
				}
				else if(y<=106.25)
				{
					x_tail = 700 - 4*y ;
			        x = x_tail ;
			        while(x <= 500)
			        {
			            element[i].x1 = x ;
	                    element[i].y1 = y ;
	                    element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	                    element[i].y2 = element[i].y1 ;
	                    element[i].x3 = element[i].x2 ;
	                    element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	                    element[i+1].x1 = x ;
	                    element[i+1].y1 = y ;
	                    element[i+1].x2 = element[i+1].x1 ;
	                    element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	                    element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	                    element[i+1].y3 = element[i+1].y2 ;
				        x = element[i].x2 ;
		                y = element[i].y2 ;
		                i += 2 ;
			        }
				}
			}
		}
		else if((y>=110)&&(y<122))
		{
			x_nose = 0.41667*y + 169.1667 ;
			while(x <= x_nose)
			{
			    element[i].x1 = x ;
	            element[i].y1 = y ;
	            element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	            element[i].y2 = element[i].y1 ;
	            element[i].x3 = element[i].x2 ;
	            element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	            element[i+1].x1 = x ;
	            element[i+1].y1 = y ;
	            element[i+1].x2 = element[i+1].x1 ;
	            element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	            element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	            element[i+1].y3 = element[i+1].y2 ;
				x = element[i].x2 ;
		        y = element[i].y2 ;
		        i += 2 ;
			}
			
			x = 227 ; 
			while(x <= 500)
			{
			    element[i].x1 = x ;
	            element[i].y1 = y ;
	            element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	            element[i].y2 = element[i].y1 ;
	            element[i].x3 = element[i].x2 ;
	            element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	            element[i+1].x1 = x ;
	            element[i+1].y1 = y ;
	            element[i+1].x2 = element[i+1].x1 ;
	            element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	            element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	            element[i+1].y3 = element[i+1].y2 ;
				x = element[i].x2 ;
		        y = element[i].y2 ;
		        i += 2 ;
			}
		}
		else if(y>=122)
		{
			while(x<=500)
			{
				element[i].x1 = x ;
	            element[i].y1 = y ;
	            element[i].x2 = element[i].x1 + ((y-104)/105)*dx ;
	            element[i].y2 = element[i].y1 ;
	            element[i].x3 = element[i].x2 ;
	            element[i].y3 = element[i].y2 + ((y-104)/105)*dy ;
	            element[i+1].x1 = x ;
	            element[i+1].y1 = y ;
	            element[i+1].x2 = element[i+1].x1 ;
	            element[i+1].y2 = element[i+1].y1 + ((y-104)/105)*dy ;
	            element[i+1].x3 = element[i+1].x2 + ((y-104)/105)*dx ;
	            element[i+1].y3 = element[i+1].y2 ;
				x = element[i].x2 ;
		        y = element[i].y2 ;
		        i += 2 ;
			}
		}
		x = 0 ;
		y = element[i-1].y3 ;
		l+=1 ;
    }
    
    cout << " No. of elements = " << i-1 << "\n" ;
	
	ofstream aStream ;
    aStream.open("Mesh.txt") ;
    
	for(j=0 ; j<i ; j++)
	{
		aStream << element[j].x1 << " " << element[j].y1 << " \n" << element[j].x2 << " " << element[j].y2 << " \n" << element[j].x3 << " " << element[j].y3 << " \n" ;
	}
	
	aStream.close();
    cout<<"\n Mesh details are stored in Mesh.txt file"<<endl;
    
    ofstream bStream ;
    bStream.open("Triangle_connectivity.txt") ;
    n = (i-1)*3 ;
    
	for(k=0 ; k<=n ; k++)
	{
		bStream << k+1 << " " << k+2 << " " << k+3 << " \n" ;
		k+=2 ;
	}
	
	bStream.close();
    cout<<"\n Triangle connectivity details are stored in Triangle_connectivity.txt file"<<endl;
    
    j = 0 ;
    k = 0 ;
    while(j<i)
    {
        co[k][0] = element[j+1].x1 ;
        co[k][1] = element[j+1].y1 ;
        co[k+1][0] = element[j+1].x2 ;
        co[k+1][1] = element[j+1].y2 ;
        co[k+2][0] = element[j+1].x3 ;
        co[k+2][1] = element[j+1].y3 ;
        co[k+3][0] = element[j].x2 ;
        co[k+3][1] = element[j].y2 ;
        j+=2 ;
        k+=4 ;
	}
	
	j = 0 ;
	r = 0 ;
	u[0] = rel_speed ;
	v[0] = 0 ;
	u[1] = rel_speed ;
	v[1] = 0 ;
	u[3] = rel_speed ;
	v[3] = 0 ;
	e1 = 0 ;
	e2 = 0 ;
	f1 = 0 ;
	f2 = 0 ;
	g1 = 0 ;
	g2 = 0 ;
    while(k <= i)
	{
		while(j<i)
		{
			if((element[j].x1==co[k][0])&&(element[j].y1==co[k][1]))
			{
				p[r] = j ;
				r+=1;
			}
			else if((element[j].x2==co[k][0])&&(element[j].y2==co[k][1]))
			{
				p[r] = j ;
				r+=1 ;
			}
			else if((element[j].x3==co[k][0])&&(element[j].y3==co[k][1]))
			{
				p[r] = j ;
				r+=1 ;
			}
			j+=1;
		}
		s = 0 ;
		while(s<r)
		{
			find_coefficients(element[p[s]].x1,element[p[s]].y1,element[p[s]].x2,element[p[s]].y2,element[p[s]].x3,element[p[s]].y3,10,0,10,0);
			e1 = e1 + a1 ;
			e2 = e2 + a2 ;
			f1 = f1 + b1 ;
			f2 = f2 + b2 ;
			g1 = g1 + c1 ;
			g2 = g2 + c2 ;
			
        }
        u[k] = e1*co[k][0]+e2*co[k][1]*co[k][0] ;
		v[k] = f1*co[k][1]+f2*co[k][1]*co[k][0] ;
		P[k] = u[k] = g1*pow(co[k][0],2)+g2*pow(co[k][1],2) ; 
	}
    
    ofstream cStream ;
    cStream.open("u-velocity.txt") ;
    
	for(j=0 ; j<i ; j++)
	{
		cStream << u[j] << "\n" ;
	}
	
	cStream.close();
    cout<<"\n Horizontal velocities are stored in u-velocity.txt file"<<endl;
    
    ofstream dStream ;
    dStream.open("v-velocity.txt") ;
    
	for(j=0 ; j<i ; j++)
	{
		dStream << v[j] << "\n" ;
	}
	
	dStream.close();
    cout<<"\n Vertical velocities are stored in v-velocity.txt file"<<endl;
    
    ofstream eStream ;
    eStream.open("Pressure.txt") ;
    
	for(j=0 ; j<i ; j++)
	{
		eStream << u[j] << "\n" ;
	}
	
	eStream.close();
    cout<<"\n Pressure values are stored in Pressure.txt file"<<endl;
    cout<<"\n -Submitted by Shivanand Paramashivam Pillai ( Roll.No.OE21M027 )" ;
    
	return 0 ;
}
