 /** create vector of weights for Finite difference - second derivative 
  * 4th order difference scheme - does not include boundary **/
 
void d2ydx2( vector<vector<double> > &vec2d, // 2d- vector of coefficients
             vector<vector<int> > &vec2dh){  // 2d - vector of grid positions - factors for grid step size h
  
  int n=vec2d.size();// n = total number of equations NOT number of grid points
  
  vector<double> coef1{ (5./6), (-5./4), (-1./3), (7./6), (-1./2), (1./12) }; // 4th order coefs. next step from boundary
  
  vector<double> coef2{ (-1./12), (4./3), (-5./2), (4./3), (-1./12), };// 4th order coefs. middle range grid points
  
    int r=0,c=0;
 
    r=1;// second from left edge in x - next step from lower boundary. r="row" 0 to n-1
    for(c=0; c<6; c++){
          vec2d[r][c]=coef1[c];// set coefs in a row which we will sum up later in the main ODE system
          vec2dh[r][c]=c; // this is the h, e.g. f[0], f[1],...f[5] is part of the finite different equation
    }
    
    //IMPORTANT NOTE: n here is number of total equations NOT the number of grid steps, e.g. 10 grid points has n=11
    // 0,...,10 = 11 slots, so [n-2] = [9] and so the second last equation. 
    for ( r = 2; r <(n-2); r++) { 
        for(int c=0; c<5; c++){
          vec2d[r][c]=coef2[c];
          vec2dh[r][c]=c+r-2;
        }
    }
 
    r= n-2;// second from right edge in x - REVERSE ORDER, note r=9 (on a n=10 grid)
    for(c=0; c<6; c++){
          vec2d[r][c]=coef1[5-c];
         vec2dh[r][5-c]=(n-1)-c;
        }
  
}