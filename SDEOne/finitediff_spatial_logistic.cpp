// to run to sourceCpp("finitediff_spatial_logistic.cpp")
// then call m<-ode(inits,0.2,50.0,150.0,100,(100- -0)/100)
// note rep(0,11) as 11 equations, grid size = 10
// needs libraries Rcpp and BH
// THIS HAS BEEN CHECKED WITH MATHEMATICA - it's correct
// finitediff_ex1.nb

// We can now use the BH package
// [[Rcpp::depends(BH)]]

#ifdef IGNORE
inits<-rep(0,101);
inits[34:48]<-c(0.0000888416, 0.00071113, 0.0041334, 0.0174458, 0.0534689,  0.118997,0.192308, 
 0.225676,  0.192308, 
 0.118997,0.0534689, 
 0.0174458, 0.0041334, 
 0.00071113,  0.0000888416)


library(plotly)
fig <- plot_ly(
  type = 'surface',#heatmap',colors = topo.colors(100),#colorRamp(c("blue", "green")),
#contours = list(
#    z = list(show = TRUE, start = 0.0, end = 0.07, size = 0.01)),
  z = ~m[,-1],showscale = TRUE, colorbar=list(title = "Allele Density"))
fig <- fig %>% layout(
    scene = list(
      yaxis = list(showticklabels=FALSE, title='Time'),
      xaxis = list(showticklabels=FALSE, title="Distance"),
      zaxis = list(showticklabels=FALSE, title="Density")))
      #zaxis = list(range = list(0, 0.07),nticks = 4)))#,#zaxis = list(nticks = 4)
      #camera = list(eye = list(x = 0, y = -1, z = 0.05)),
#aspectratio = list(x = .5, y = .5, z = 0.02)))

fig

    
#endif  


#include <Rcpp.h>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

using namespace Rcpp;
using namespace std; 
using namespace boost::numeric::odeint;

// if defined then print out extra things
#define DEBUGyes

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

// Don't change - boost boilerplate
//[ integrate_observer - this saves the state at each timestep used by the solver
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
//]

/** FINITE DIFFERENT SPECIFICATION
    second derivatie wrt space, partial f^2/partial x^2  e.g. u_t = (1/8) u_xx
    In Mathematica
    NDSolve`FiniteDifferenceDerivative[Derivative[2], h Range[0, 10], Map[f, h Range[0, 10]]]
   
    Basic idea is that finite difference has a consistent pattern, boundaries are special, [0], [n]
    and dealt with elsewhere, and then the next steps from each bounary - [1] and [n-1] - are mirror version of each other and 
    not symmetrical about the mid point - these are included below, and then the remaining grid is symmetrical from [2],...[n-2]
 
    The coefficients are hard coded, one set for [1],[n-1] and one set for [2],...,[n-2]. 
 
*/ 
 
  /** create vector of weights for Finite difference - second derivative and including endpoints which should be overwritten **/
void d2ydx2( vector<vector<double> > &vec2d, 
             vector<vector<int> > &vec2dh){ //n=20 steps 0,...20
  
  int n=vec2d.size();// n = total number of equations NOT number of grid points
  
  vector<double> coef1{ (5./6), (-5./4), (-1./3), (7./6), (-1./2), (1./12) }; // next step from boundary
  
  vector<double> coef2{ (-1./12), (4./3), (-5./2), (4./3), (-1./12), };// middle range grid points
  
    //
    int r=0,c=0;
 
    r=1;// second from left edge in x - next step from lower boundary
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
 
    r= n-2;// second from right edge in x - REVERSE ORDER, note r=9 (on a n=10 grid
    for(c=0; c<6; c++){
          vec2d[r][c]=coef1[5-c];
         vec2dh[r][5-c]=(n-1)-c;
        }
  
}


  /** create vector of weights for Finite difference - first derivative and including endpoints which should be overwritten **/
void dydx( vector<vector<double> > &vec2d, 
             vector<vector<int> > &vec2dh){ //n=20 steps 0,...20
  
  int n=vec2d.size();

  vector<double> coef1{ (-1./4), (-5./6), (3./2), (-1./2), (1./12) };
  
  vector<double> coef2{ (1./12), (-2./3), (2./3), (-1./12) };
  
    //
    int r=0,c=0;
 
    r=1;// second from left edge in x
    for(c=0; c<5; c++){
          vec2d[r][c]=coef1[c];
          vec2dh[r][c]=c;
    }
    
 
    for ( r = 2; r <(n-2); r++) { 
        for(int c=0; c<4; c++){
          vec2d[r][c]=coef2[c];
        }
        for(int c=0; c<2; c++){
          vec2dh[r][c]=c+r-2;
        }
        for(int c=3; c<5; c++){
          vec2dh[r][c-1]=c+r-2;
        }
        
    }
 
    r= n-2;// second from right edge in x - REVERSE ORDER
    for(c=0; c<5; c++){
          vec2d[r][c]= -coef1[4-c];
         vec2dh[r][4-c]=(n-1)-c;
        }
  
}



//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class f_ode {

    //members
    double m_alpha,m_beta;                     //misc model parameters
    vector<vector<double> > &m_vec2d_xx; //second derivative FD 
    vector<vector<int> > &m_vec2dh_xx;   //second derivative FD
    vector<vector<double> > &m_vec2d_x;  //second derivative FD 
    vector<vector<int> > &m_vec2dh_x;    //second derivative FD
    double m_h;                          //spatial grid step size, e.g. h=1./10;
    int m_n;                             //number of equations, e.g. if h=1/10 then m_n=11
    vector<double> &m_vec_xgrid;         //grid of x values - boundary to boundary 
    
    // other useful variables
    double coef1_0=0.0;
    double coef2_0=0.0;
    double coef3_0=0.0;
    double coef4_0=0.0;
    double coef1_1=0.0;
    double coef2_1=0.0;
    double coef3_1=0.0;
    double coef4_1=0.0;

    
public:
    f_ode( double alpha,
           double beta,
           double h,
           int n,
           vector<vector<double> > &vec2d_xx, 
           vector<vector<int> > &vec2dh_xx,
           vector<vector<double> > &vec2d_x, 
           vector<vector<int> > &vec2dh_x,
           vector<double> &vec_xgrid): m_alpha(alpha), 
                                       m_beta(beta),
                                             m_h(h),
                                             m_n(n),
                                             m_vec2d_xx(vec2d_xx),
                                             m_vec2dh_xx(vec2dh_xx),
                                             m_vec2d_x(vec2d_x),
                                             m_vec2dh_x(vec2dh_x),
                                             m_vec_xgrid(vec_xgrid){ }

    void operator() ( const state_type &x , state_type &dxdt , const double  t  )
    {
      
    /******************* START OF KEY SECTION = THE MODELLING PART WHICH SHOULD BE CUSTOMIZED **************************/
        // CUSTOM BOUNDARY - LOWER in space dimension

        //Rcout<<m_mu<<"  "<<m_sigma2<<" "<<m_theta<<" "<<m_vec_xgrid[0]<<" "<<m_vec_xgrid[100]<<endl;
        //Rcpp::stop("stopping\n");
      
      dxdt[0] = -x[0];// f+f'=0
      
        // FINITE DIFFERENCE inside boundaries
        // build up equation as a sum of terms
        for(int r=1;r<(m_n-1);r++){//e.g. m_n=11 for n=10 steps
          
            dxdt[r]=0;// initialize to zero
          
          for(int c=0;c<m_vec2d_xx[r].size();c++){ //add up all the relevant terms using the cofficients stored in m_vec2d[] and m_vec2dh
            
            // this is PDE specific - tailor to the equation, this line is for f_xx 
              dxdt[r] += m_vec2d_xx[r][c]*x[m_vec2dh_xx[r][c]]/(m_h*m_h); // f_t = f_xx ...
            }
          
          // no finite diff here
          dxdt[r] += m_alpha*x[r]*(1.-x[r]/m_beta); 
          
          
        } // end of row loop, i.e. ODEs in the body
        
       // CUSTOM BOUNDARY - UPPER in space dimension

                  
        
      dxdt[m_n-1] = -x[m_n-1];// f+f'=0
                       
        
      
    /******************* END OF MODELLING SECTION *********************************************************************/  
    }
};
//]


/** This last part is the calling function from R */

// [[Rcpp::export]]
NumericMatrix ode(NumericVector x0,double alpha, double beta, double tend,int nn, double h) {
                 // initial values                      time end point, number of grid steps, grid step size
  
  int n=nn+1;//nn is the grid steps, but n is the number of equations, i.e n=10 then nn=11 as 0,...10 inclusive
  int j=0;
  
  //Allocate storage for the finite difference vectors of coefficients
  //These are u_xx hence _xx subscript, myvec is the coefficients, e.g. myvec[i]*f(.), and myvech is the grid point, e.g. f(myvech[i])
  vector<vector<double> > myvec_xx(n); //coeffs must be 11 because need 0,...,10 inclusive
  vector<vector<int> > myvech_xx(n);   // h factors
  // allocate number of cols - note ragged as LHS and RHS near boundary are different patterns from in the middle of grid
  for(int i=2;i<(n-2);i++){ // these 5 and 6 are hardcoded based on the FD order of accuracy - fourth order as per Mathematica default scheme
    myvec_xx[i].resize(5);
    myvech_xx[i].resize(5);
    }
  myvec_xx[1].resize(6);
  myvec_xx[n-2].resize(6);
  myvech_xx[1].resize(6);
  myvech_xx[n-2].resize(6);
  
  //As above but vectors for the finite difference coefficients for u_x - first order partial derivative
  vector<vector<double> > myvec_x(n);// must be 11 because need 0,...,10 inclusive
  vector<vector<int> > myvech_x(n);
  // as above - ragged and hardcoded
  for(int i=2;i<(n-2);i++){
    myvec_x[i].resize(4);
    myvech_x[i].resize(4);
    }
  myvec_x[1].resize(5);
  myvec_x[n-2].resize(5);
  myvech_x[1].resize(5);
  myvech_x[n-2].resize(5);
  
  
  // Call functions to fill in the vectors created above with the necessary FD coefficients.
   d2ydx2(myvec_xx,myvech_xx); // second derivative FD
   dydx(myvec_x,myvech_x);     // first derivative FD
   
  //create a vector for storing the values of x across the grid, e.g. x[0]=Lower+0*h,x[1]=Lower+1*h,x[2]=Lower+2*h etc   
  vector<double> myvec_xgrid(n);
  double lowerX=0;//lower end of x
  int i=0;
  // so this will be e.g. from -4,...,+4
  for (auto it = myvec_xgrid.begin();it != myvec_xgrid.end(); it++){*it=lowerX+(i++)*h;}
  
   
#ifdef DEBUGyes
   
   Rcout<<"MYVEC_XX"<<endl;
    for (auto it = myvec_xx.begin();it != myvec_xx.end(); it++){
       for (auto it2 =(*it).begin();it2 != (*it).end(); it2++){
            cout << *it2 << " ";
       }
        Rcout << endl;
    }
 Rcout<<endl;
  
    Rcout<<"MYVECH_XX"<<endl;
 for (auto it = myvech_xx.begin();it != myvech_xx.end(); it++){
       for (auto it2 =(*it).begin();it2 != (*it).end(); it2++){
            Rcout << *it2 << " ";
       }
        Rcout << endl;
    }
  
  Rcout<<"MYVEC_X"<<endl;
    for (auto it = myvec_x.begin();it != myvec_x.end(); it++){
       for (auto it2 =(*it).begin();it2 != (*it).end(); it2++){
            Rcout << *it2 << " ";
       }
        Rcout << endl;
    }
 Rcout<<endl;
    Rcout<<"MYVECH_X"<<endl;
 for (auto it = myvech_x.begin();it != myvech_x.end(); it++){
       for (auto it2 =(*it).begin();it2 != (*it).end(); it2++){
            Rcout << *it2 << " ";
       }
        Rcout << endl;
    }
 Rcout<<"xgrid=="<<endl;
  for (auto it = myvec_xgrid.begin();it != myvec_xgrid.end(); it++){
   Rcout<<"==>"<<*it<<endl;}
  
#endif
  
 //Rcpp::stop("stopping\n"); 
  
  /** NOW FOR THE BOOST ODEINT PART **/
    //[ state_initialization
    state_type x(n);//n is the dimension of the problem
    // set initial state at t=tstart for the system
    for(j=0;j<n;j++){ 
    x[j] = x0[j];} 


    //[ integration_class
   // f_ode myf(k,r);
  //  size_t steps = integrate( myf ,
  //          x , 0.0 , 4.0 , 0.1 );
    //]

     vector<state_type> x_vec;
    vector<double> times;
f_ode myf(alpha,beta,h,n,myvec_xx,myvech_xx,myvec_x,myvech_x,myvec_xgrid);

    size_t steps = integrate( myf,
            x , 0.0 , tend , 0.01 ,
            push_back_state_and_time( x_vec , times ) );

//  bulirsch_stoer_dense_out< state_type > stepper( 1E-5 , 1E-5 , 1E-5 , 1E-5 );
  
  //make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ); 
  
 /*size_t steps =    integrate_adaptive(stepper, myf,
            x , 0.0 , tend , 0.01 ,
            push_back_state_and_time( x_vec , times ) );
   */ 
    /* output */
    /*for( size_t i=0; i<=steps; i++ )
    {
        Rcout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << endl;
    }
    */
    
    NumericMatrix results (steps,x.size()+1);
    
    for (auto r=0;r<steps;r++){
      for(auto c=0;c<(x.size()+1);c++){
        if(c==0){results(r,c)=times[r];
        }else{
    /*Rcout<<"result="<<i<<endl;*/
    results(r,c)=x_vec[r][c-1];
        }}}
    
    return(results);
}

