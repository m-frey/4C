/*!-----------------------------------------------------------------------------------------------------------
\file beam2_evaluate.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*-----------------------------------------------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "beam2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam2::ActionType act = Beam2::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Beam2::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Beam2::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam2::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Beam2::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Beam2::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam2::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = Beam2::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Beam2::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Beam2::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam2::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam2::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = Beam2::calc_struct_reset_istep;
  else if (action=="calc_brownian")       act = Beam2::calc_brownian;
  else if (action=="calc_struct_ptcstiff")        act = Beam2::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam2");

  switch(act)
  {
    case Beam2::calc_struct_ptcstiff:
    {
      //Beam2 element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size
    }
    break;
    //action type for evaluating statistical forces
    case Beam2::calc_brownian:
    {   
    	/*
      // get element displacements (for use in shear flow fields)
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      */
      
       
      /*in case of parallel computing the random forces have to be handled in a special way: normally
       * in the frame of evaluation each processor evaluates forces for its column elements (including
       * ghost elements); later on in the assembly each processor adds the thereby gained forces only
       * for those DOF of whose owner it is. In case of random forces such an assembly would render it
       * impossible to establish certain correlations between forces related to nodes with different ownerss
       * (for each nodes the random forces would be evaluated in an identical process, but due to 
       * independent random numbers); as correlation between forces is restricted to the support of at
       * the maximum one element a solution to this problem is, to evaluate all the forces of one element
       * only by means of one specific processor (here we employ the elemnet owner processor); these
       * forces are assembled in a column map vector and later exported to a row map force vector; this 
       * export is carried out additively so that it is important not to evaluate any forces at all if
       * this processor is not owner of the element;
       * note: the crucial difference between this assembly and the common one is that for certain nodal
       * forces not the owner of the node is responsible, but the owner of the element*/
      
      //test whether this processor is row map owner of the element (otherwise no forces added)
      //if(this->Owner() != discretization.Comm().MyPID()) return 0;

      
      //compute stochastic forces in local frame
      ComputeLocalBrownianForces(params);
      

    }
    break;
    /*in case that only linear stiffness matrix is required b2_nlstiffmass is called with zero dispalcement and
     residual values*/
    case Beam2::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam2::calc_struct_nlnstiffmass:
    case Beam2::calc_struct_nlnstifflmass:
    case Beam2::calc_struct_nlnstiff:
    case Beam2::calc_struct_internalforce:
    {
      int lumpedmass = 0;  // 0=consistent, 1=lumped
      if (act==Beam2::calc_struct_nlnstifflmass) lumpedmass = 1;

      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // get element velocities
      RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==null) dserror("Cannot get state vectors 'velocity'");
      vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      // determine element matrices and forces
      if (act == Beam2::calc_struct_nlnstiffmass)
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      else if (act == Beam2::calc_struct_nlnstifflmass)
      {
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      }
      else if (act == Beam2::calc_struct_nlnstiff)
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,NULL,&elevec1,lumpedmass);
      else if  (act ==  calc_struct_internalforce)
        nlnstiffmass(params,lm,myvel,mydisp,NULL,NULL,&elevec1,lumpedmass);
      
      /*at the end of an iteration step the geometric ocnfiguration has to be updated: the starting point for the
       * next iteration step is the configuration at the end of the current step */
      numperiodsold_ = numperiodsnew_;
      alphaold_ = alphanew_;


      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      /*if(Id() == 2) //limiting the following tests to certain element numbers
      {     
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(6,6);

        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(6,6);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-8;

        //flag indicating whether approximation lead to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<3; i++)
        {
          for(int k=0; k<2; k++)
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(6);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            vector<double> vel_aux(6);
            vector<double> disp_aux(6);
            for(int id = 0;id<6;id++)
            {
                DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
                DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);
            }
            
            //modifying displacment artificially (for numerical derivative of internal forces):
            disp_aux[i + 3*k] = disp_aux[i + 3*k] + h_rel;
             vel_aux[i + 3*k] =  vel_aux[i + 3*k] + h_rel * params.get<double>("gamma",1.0) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",1.0) );
            
            //std::cout << " approximation force calc " <<  i << k;
             //calc forces depend on vel_aux and disp_aux
            nlnstiffmass(params,lm,vel_aux,disp_aux,NULL,NULL,&force_aux,lumpedmass);
            

            
            //calc approx stiffmatrix
            for(int u = 0;u<6;u++)
            {
              stiff_approx(u,i+k*3)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
            }
            
          }
        }

        
       for(int line=0; line<6; line++)
       {
         for(int col=0; col<6; col++)
         {
           stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

           //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
           if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0)
             stiff_relerr(line,col) = 0;

           if ( stiff_relerr(line,col) > 0)
             outputflag = 1;
         }
       }

       if(outputflag ==1)
       {
         std::cout<<"\n\n actually calculated stiffness matrix"<< elemat1;
         std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;
         std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
       }
      }
      //end of section in which numerical approximation for stiffness matrix is computed
      */

    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geomatric
       * status of the beam have to be updated;*/
      numperiodsconv_ = numperiodsnew_;
      alphaconv_ = alphanew_;
    }
    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      numperiodsold_ = numperiodsconv_;
      alphaold_ = alphaconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam2 elements");
    default:
      dserror("Unknown type of action for Beam2 %d", act);
  }
  return 0;

}//DRT::ELEMENTS::Beam2::Evaluate


/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  // element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  //jacobian determinant
  double det = lrefe_/2;


  // no. of nodes on this element
  const int iel = NumNode();
  const int numdf = 3;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D  intpoints(gaussrule_);


  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector      funct(iel);

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val   = condition.Get<vector<double> >("val");

  //integration loops
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[ip][0];
    const double wgt = intpoints.qwgt[ip];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * det;

    // load vector ar
    double ar[numdf];
    // loop the dofs of a node

    for (int i=0; i<numdf; ++i)
    {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    }

    //sum up load components
    for (int node=0; node<iel; ++node)
      for (int dof=0; dof<numdf; ++dof)
         elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | Compute forces for Brownian motion (public)                                                       06/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2::ComputeLocalBrownianForces(ParameterList& params)
{
  /*creating a random generator object which creates random numbers with mean = 0 and standard deviation
   * (2kT/dt)^0,5 with thermal energy kT, time step size dt; using Blitz namespace "ranlib" for random number generation*/
  ranlib::Normal<double> normalGen(0,pow(2.0 * params.get<double>("KT",0.0) / params.get<double>("delta time",0.0),0.5));
  
  //fstoch consists of 4 values, two forces at each node
  LINALG::Matrix<4,1> aux;
  aux(0) = normalGen.random();
  aux(1) = normalGen.random();
  aux(2) = normalGen.random();
  aux(3) = normalGen.random();
  

  

  /*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
  * length is approximated by the one of an infinitely long staff for friciton orthogonal to staff axis*/
  double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);
  

  int stochasticorder = params.get<int>("STOCH_ORDER",-1);

  
  switch(stochasticorder)
  {
  case -1:
  {
    //multiply S_loc (cholesky decomposition of C_loc) for C_loc only diagonal
    floc_(0,0) = pow(zeta/2.0,0.5)*aux(0);
    floc_(1,0) = pow(zeta/2.0,0.5)*aux(1);
    floc_(2,0) = pow(zeta/2.0,0.5)*aux(2);
    floc_(3,0) = pow(zeta/2.0,0.5)*aux(3); 
    
    
  }
  break;
  case 0:
  {
    //multiply S_loc(cholesky decomposition of C_loc) for gamma_parallel=gamma_perp
  	floc_(0,0) = pow(zeta/3.0,0.5)*aux(0);
  	floc_(1,0) = pow(zeta/3.0,0.5)*aux(1);
  	floc_(2,0) = pow(zeta/12.0,0.5)*aux(0)+pow(zeta/4.0,0.5)*aux(2);
  	floc_(3,0) = pow(zeta/12.0,0.5)*aux(1)+pow(zeta/4.0,0.5)*aux(3);  
    
    
  }
  break;
  case 1:
  {
    //multiply S_loc(cholesky decomposition of C_loc) for gamma_parallel=gamma_perp/2
  	floc_(0,0) = pow(zeta/6.0,0.5)*aux(0);
  	floc_(1,0) = pow(zeta/3.0,0.5)*aux(1);
  	floc_(2,0) = pow(zeta/24.0,0.5)*aux(0)+pow(zeta*3.0/24.0,0.5)*aux(2);
  	floc_(3,0) = pow(zeta/12.0,0.5)*aux(1)+pow(zeta/4.0,0.5)*aux(3);     

  }
  break;
  default:
      dserror("Unknown type of stochasticorder for Beam2 %d", stochasticorder);
  }

return 0; 
}//DRT::ELEMENTS::Beam2::ComputeLocalBrownianForces
/*-----------------------------------------------------------------------------------------------------------*
 | Assemble statistical forces and damping matrix according to fluctuation dissipation theorem (public) 06/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam2::CalcBrownian(ParameterList& params,
                              vector<int>&              lm,
                              vector<double>&           vel,  //!< element velocity vector
                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                              Epetra_SerialDenseVector* force)  //!< element internal force vector
{
	//define parameters
	 double dt = params.get<double>("delta time",0.0);
   int stochasticorder = params.get<int>("STOCH_ORDER",0);// polynomial order for interpolation of stochastic line load
   double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);

   switch(stochasticorder)
   {
   //simple isotropic model of Brownian motion with uncorrelated nodal forces
   case -1:
   {
   	 
  	 //calc brownian damp matrix 
  	 if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
  	 {
  		 (*stiffmatrix)(0,0) += zeta/(2.0*dt);
  		 (*stiffmatrix)(1,1) += zeta/(2.0*dt);
  		 (*stiffmatrix)(3,3) += zeta/(2.0*dt);
  		 (*stiffmatrix)(4,4) += zeta/(2.0*dt);
  	 }
  	    
  	 if (force != NULL)
  	 {
  		 //calc internal brownian forces  
  		 (*force)(0) +=zeta/(2.0)*vel[0];
  		 (*force)(1) +=zeta/(2.0)*vel[1];
  		 (*force)(3) +=zeta/(2.0)*vel[3];
  		 (*force)(4) +=zeta/(2.0)*vel[4];
      
  		 //calc external brownian forces 
  		 (*force)(0) -=floc_(0,0);
  		 (*force)(1) -=floc_(1,0);
  		 (*force)(3) -=floc_(2,0);
  		 (*force)(4) -=floc_(3,0); 
  	 }
     

    }
   break;
   
   //isotropic model of Brownian motion with correlated forces
   case 0:
   {   
	 
     //calc brownian damp matrix 
  	 if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
  	 {
  		 (*stiffmatrix)(0,0) += zeta/(3.0*dt);
  		 (*stiffmatrix)(1,1) += zeta/(3.0*dt);
  		 (*stiffmatrix)(3,3) += zeta/(3.0*dt);
  		 (*stiffmatrix)(4,4) += zeta/(3.0*dt);
  		 (*stiffmatrix)(0,3) += zeta/(6.0*dt);
  		 (*stiffmatrix)(1,4) += zeta/(6.0*dt);
  		 (*stiffmatrix)(3,0) += zeta/(6.0*dt);
  		 (*stiffmatrix)(4,1) += zeta/(6.0*dt);
  	 } 
     
 	 
  	 //calc int brownian forces
  	 if (force !=  NULL)
  	 {
  		 (*force)(0) +=zeta/(3.0)*vel[0]+zeta/(6.0)*vel[3];
  		 (*force)(1) +=zeta/(3.0)*vel[1]+zeta/(6.0)*vel[4];
  		 (*force)(3) +=zeta/(6.0)*vel[0]+zeta/(3.0)*vel[3];
  		 (*force)(4) +=zeta/(6.0)*vel[1]+zeta/(3.0)*vel[4];

  		 //calc ext brownian forces

  		 (*force)(0) -=floc_(0,0);
  		 (*force)(1) -=floc_(1,0);
  		 (*force)(3) -=floc_(2,0);
  		 (*force)(4) -=floc_(3,0); 
  	 }
     
     
     }
   break;
   //anisotropic model of Brownian motion with correlated nodal forces
   case 1:
   { 
		//triad to rotate local configuration into global
		LINALG::Matrix<2,2> T(true);
		T(0,0) =  cos(alphanew_);
		T(0,1) = -sin(alphanew_);
		T(1,0) =  sin(alphanew_);
		T(1,1) =  cos(alphanew_);

    //local damping matrix
    LINALG::Matrix<2,2> dampbasis(true);
    dampbasis(0,0) = zeta/2.0;
    dampbasis(1,1) = zeta;


    //turning local into global damping matrix (storing intermediate result in variable "aux1")
    LINALG::Matrix<2,2> aux1;
		aux1.Multiply(T,dampbasis);
		dampbasis.MultiplyNT(aux1,T);
  
    	
    //calc brownian damp
    if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
    {
    	
    	//complete first term due to variation of velocity     
    	for(int i=0; i<2; i++)
    	{
    		for(int j=0; j<2; j++)
    		{
    			(*stiffmatrix)(i,j) += dampbasis(i,j)/(3.0*dt);
    			(*stiffmatrix)(i+3,j+3) += dampbasis(i,j)/(3.0*dt);
    			(*stiffmatrix)(i,j+3) += dampbasis(i,j)/(6.0*dt);
    			(*stiffmatrix)(i+3,j) += dampbasis(i,j)/(6.0*dt);
    		}
    	}
    	
  
    	//calc second term due to variation of Triad
    	LINALG::Matrix<2,2> Spin(true);//Spinmatrix
    	Spin(0,1)=-1.0;
    	Spin(1,0)=1.0;
     
    	//multiply Spin from right side to C_global
    	aux1.Multiply(Spin,dampbasis);
     
    	//multiply Spin from left side to C_global
    	LINALG::Matrix<2,2> aux2(true);
    	aux2.Multiply(dampbasis,Spin);

    	for (int i=0; i<2; i++)//difference
    	{
    		for (int j=0; j<2; j++)
    		{
    			aux1(i,j)=aux1(i,j)-aux2(i,j);
    		}
    	} 
    	
    	LINALG::Matrix<4,4> aux3(true);
    	for(int i=0; i<2; i++)//multiply velcoefficients
    	{
    		for(int j=0; j<2; j++)
    		{
    			aux3(i,j) += aux1(i,j)/(3.0);
    			aux3(i+2,j+2) += aux1(i,j)/(3.0);
    			aux3(i,j+2) += aux1(i,j)/(6.0);
    			aux3(i+2,j) += aux1(i,j)/(6.0);
    		}
    	}   
	   
     LINALG::Matrix<4,1> aux4(true);
     for (int i=0; i<4; i++)//multiply velocity
     {
    	 for (int j=0; j<2; j++)
    	 {
    		 aux4(i,0)+=aux3(i,j)*vel[j];
    		 aux4(i,0)+=aux3(i,j+2)*vel[j+3];
    	 }
     }
     

     double cos_alpha = cos(alphanew_);
     double sin_alpha = sin(alphanew_);
         	
     //vector z according to Crisfield, Vol. 1, (7.66)(reduced to transl. dofs)
     LINALG::Matrix<4,1> z;
     z(0) = sin_alpha/lrefe_;
     z(1) = -cos_alpha/lrefe_;
     z(2) = -sin_alpha/lrefe_;
     z(3) = cos_alpha/lrefe_;
     
     
     //add calculated contribution to stiffness
     for(int i=0; i<2; i++)
     {
    	 for(int j=0; j<2; j++)
    	 {
    		 (*stiffmatrix)(i,j) +=aux4(i,0)*z(j,0);
    		 (*stiffmatrix)(i+3,j) +=aux4(i+2,0)*z(j,0);
    		 (*stiffmatrix)(i+3,j+3) +=aux4(i+2,0)*z(j+2,0);
    		 (*stiffmatrix)(i,j+3) +=aux4(i,0)*z(j+2,0);			
    	 }
     }
     //end internal stiffness
    
  
     //calc ext. stiffness			
     aux1.Multiply(Spin,T);			
     
     LINALG::Matrix<4,1>aux5(true);//multiply fstoch with ST, aux5 vector now for two nodes
     for (int i=0; i<2; i++)
     {
    	 for (int j=0; j<2; j++)
    	 {
    		 aux5(i,0)+=aux1(i,j)*floc_(j,0);
    		 aux5(i+2,0)+=aux1(i,j)*floc_(j+2,0);
    	 }
     }
     		
     //subtract ext stiffness
     for(int i=0; i<2; i++)
     {
    	 for(int j=0; j<2; j++)
    	 {
    		 (*stiffmatrix)(i,j) -=aux5(i,0)*z(j,0);
    		 (*stiffmatrix)(i+3,j) -=aux5(i+2,0)*z(j,0);
    		 (*stiffmatrix)(i+3,j+3) -=aux5(i+2,0)*z(j+2,0);
    		 (*stiffmatrix)(i,j+3) -=aux5(i,0)*z(j+2,0);
    	 }
     }//end ext stiffness 	      
    }//end stiffmatrix calc

    
    //calc forces   
    if (force != NULL)
    {
    	//calc internal brownian forces
    	LINALG::Matrix<4,4> aux6(true);
	   	for(int i=0; i<2; i++)
	   	{
	   		for(int j=0; j<2; j++)
	   		{
	   			aux6(i,j) = dampbasis(i,j)/(3.0);
	   			aux6(i+2,j+2) = dampbasis(i,j)/(3.0);
	   			aux6(i,j+2) = dampbasis(i,j)/(6.0);
	   			aux6(i+2,j) = dampbasis(i,j)/(6.0);
	   		}
	   	}
    	
    	
    	for (int i=0; i<2; i++)		  
    	{
    		for (int j=0; j<2; j++)
    		{
    			(*force)(i)+=aux6(i,j)*vel[j];
    			(*force)(i)+=aux6(i,j+2)*vel[j+3];
    			(*force)(i+3)+=aux6(i+2,j)*vel[j];
    			(*force)(i+3)+=aux6(i+2,j+2)*vel[j+3];
    		}
    	}

    
    	//calc ext brownian forces
    	for (int i=0; i<2; i++)
    	{
    		for (int j=0; j<2; j++)
    		{
    			(*force)(i)-=T(i,j)*floc_(j,0);
    			(*force)(i+3)-=T(i,j)*floc_(j+2,0);
    		}
    	} 

    }//end calc forces
   
  }//end case 1 
   break;
   
  }//end calc brownian

  
  return;
  
}//DRT::ELEMENTS::Beam2::CalcBrownian

/*-----------------------------------------------------------------------------------------------------------*
 | compute current rotation absolute rotation angle of element frame out of x-axisrr              cyron 03/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam2::updatealpha(const LINALG::Matrix<3,2>& xcurr,const double& lcurr)
{
  /*befor computing the absolute rotation angle of the element frame we first compute an angle beta \in [-PI;PI[
   * from the current nodal positions; beta denotes a rotation out of the x-axis in a x-y-plane; note that 
   * this angle may differ from the absolute rotation angle alpha by a multiple of 2*PI;*/
  double beta;
  
  // beta is the rotation angle out of x-axis in a x-y-plane
  double cos_beta = (xcurr(0,1)-xcurr(0,0))/lcurr;
  double sin_beta = (xcurr(1,1)-xcurr(1,0))/lcurr;
  
  //computation of beta according to Crisfield, Vol. 1, (7.60)
  
  //if coc_beta >= 0 we know -PI/2 <= beta <= PI/2
  if (cos_beta >= 0)
    beta = asin(sin_beta);
  //else we know  beta > PI/2 or beta < -PI/2
  else
  { 
    //if sin_beta >=0 we know beta > PI/2
    if(sin_beta >= 0)
      beta =  acos(cos_beta);
    //elss we know beta > -PI/2
    else
      beta = -acos(cos_beta);
  }
  
  /* by default we assume that the difference between beta and the absolute rotation angle alpha is the same
   * multiple of 2*PI as in the iteration step before; then beta + numperiodsnew_*2*PI would be the new absolute
   * rotation angle alpha; if the difference between this angle and the absolute angle in the last converged step
   * is smaller than minus PI we assume that beta, which is evaluated in [-PI; PI[ has exceeded the upper limit of 
   * this interval in positive direciton from the last to this iteration step; then alpha can be computed from 
   * beta by adding (numperiodsnew_ + 1)*2*PI; analogously with a difference greater than +PI we assume that beta 
   * has exceeded the lower limit of the interval [-PI; PI[ in negative direction so that alpha can be computed
   * adding (numperiodsnew_ - 1)*2*PI  */
  numperiodsnew_ = numperiodsold_;
  
  if(beta + numperiodsnew_*2*PI - alphaold_ < -PI)
    numperiodsnew_ ++;
  else if(beta + numperiodsnew_*2*PI - alphaold_ > PI)
    numperiodsnew_ --;
    
 
  alphanew_ = beta + 2*PI*numperiodsnew_;
  
}

/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for corotational formulation                           cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1;
inline void DRT::ELEMENTS::Beam2::local_aux(LINALG::Matrix<3,6>& Bcurr,
                    			              LINALG::Matrix<6,1>& rcurr,
                    			              LINALG::Matrix<6,1>& zcurr,
                                        const double& lcurr,
                                        const double& lrefe_)
{
  double cos_alpha = cos(alphanew_);
  double sin_alpha = sin(alphanew_);
  
  //vector r according to Crisfield, Vol. 1, (7.62)
  rcurr(0) = -cos_alpha;
  rcurr(1) = -sin_alpha;
  rcurr(2) = 0;
  rcurr(3) = cos_alpha;
  rcurr(4) = sin_alpha;
  rcurr(5) = 0;

  //vector z according to Crisfield, Vol. 1, (7.66)
  zcurr(0) = sin_alpha;
  zcurr(1) = -cos_alpha;
  zcurr(2) = 0;
  zcurr(3) = -sin_alpha;
  zcurr(4) = cos_alpha;
  zcurr(5) = 0;

  //assigning values to each element of the Bcurr matrix, Crisfield, Vol. 1, (7.99)
  for(int id_col=0; id_col<6; id_col++)
  	{
  	  Bcurr(0,id_col) = rcurr(id_col);
  	  Bcurr(1,id_col) = 0;
  	  Bcurr(2,id_col) = (lrefe_ / lcurr) * zcurr(id_col);
  	}
    Bcurr(2,2) -= (lrefe_ / 2);
    Bcurr(2,5) -= (lrefe_ / 2);
    Bcurr(1,2) += 1;
    Bcurr(1,5) -= 1;

  return;
} /* DRT::ELEMENTS::Beam2::local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::nlnstiffmass(ParameterList& params,
                                        vector<int>&              lm,
                                        vector<double>&           vel,
                                        vector<double>&           disp,
                                        Epetra_SerialDenseMatrix* stiffmatrix,
                                        Epetra_SerialDenseMatrix* massmatrix,
                                        Epetra_SerialDenseVector* force,
                                        int lumpedmass)
{
  const int numdf = 3;
  const int iel = NumNode();
  //coordinates in current configuration of all the nodes in two dimensions stored in 3 x iel matrices
  LINALG::Matrix<3,2> xcurr;

  //current length of beam in physical space
  double lcurr = 0;

  //some geometric auxiliary variables according to Crisfield, Vol. 1
  LINALG::Matrix<6,1> zcurr;
  LINALG::Matrix<6,1> rcurr;
  LINALG::Matrix<3,6> Bcurr;
  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::Matrix<3,6> aux_CB;
  //declaration of local internal forces
  LINALG::Matrix<3,1> force_loc;
  //declaration of material parameters
  double ym; //Young's modulus
  double sm; //shear modulus
  double density; //density

  //calculating refenrence configuration xrefe and current configuration xcurr
  for (int k=0; k<iel; ++k)
  {
    xcurr(0,k) = Nodes()[k]->X()[0] + disp[k*numdf+0];
    xcurr(1,k) = Nodes()[k]->X()[1] + disp[k*numdf+1];

    /*note that xcurr(2,0),xcurr(2,1) are local angles; in Crisfield, Vol. 1, (7.98) they 
     *are denoted by theta_{l1},theta_{l2} in contrast to the local  angles theta_{1},theta_{2};
     *the global director angle is not used at all in the present element formulation*/
    xcurr(2,k) = disp[k*numdf+2];
  }

  //current length
  lcurr = pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );
  
  //update absolute rotation angle alpha of element frame
  updatealpha(xcurr,lcurr);

  //calculation of local geometrically important matrices and vectors
  local_aux(Bcurr,rcurr,zcurr,lcurr,lrefe_);


  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
	{
        case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = ym / (2*(1 + actmat->PoissonRatio()));
      density = actmat->Density();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
 }

  //Crisfield, Vol. 1, (7.52 - 7.55)
  force_loc(0) = ym*crosssec_*(lcurr*lcurr - lrefe_*lrefe_)/(lrefe_*(lcurr + lrefe_));

  //local internal bending moment, Crisfield, Vol. 1, (7.97)
  force_loc(1) = -ym*mominer_*(xcurr(2,1)-xcurr(2,0))/lrefe_;

  //local internal shear force, Crisfield, Vol. 1, (7.98)
  /*in the following internal forces are computed; note that xcurr(2,0),xcurr(2,1) are global angles;
   * in Crisfield, Vol. 1, (7.98) they are denoted by theta_{1},theta_{2} in contrast to the local
   * angles theta_{l1},theta_{l2}; note also that these variables do not represent the absolute director
   * angle (since they are zero in the beginning even for an initially rotated beam), but only
   * the absolute director angle minus the reference angle alpha0_; as a consequence the shear force
   * has to be computed by substraction of (alphanew_ - alpha0_)*/
  force_loc(2) = -sm*crosssecshear_*( ( xcurr(2,1) + xcurr(2,0) )/2 - (alphanew_ - alpha0_) );

  if (force != NULL)
  {
  //declaration of global internal force
  LINALG::Matrix<6,1> force_glob;
  //calculation of global internal forces from Crisfield, Vol. 1, (7.102): q_i = B^T q_{li}
  force_glob.MultiplyTN(Bcurr,force_loc);

    for(int k = 0; k<6; k++)
      (*force)(k) = force_glob(k);
  }


  //calculating tangential stiffness matrix in global coordinates, Crisfield, Vol. 1, (7.107)
  if (stiffmatrix != NULL)
  {
    //declaration of fixed size matrix for global tiffness
    LINALG::Matrix<6,6> stiff_glob;
    
    //linear elastic part including rotation: B^T C_t B / l_0
    for(int id_col=0; id_col<6; id_col++)
    {
      aux_CB(0,id_col) = Bcurr(0,id_col) * (ym*crosssec_/lrefe_);
      aux_CB(1,id_col) = Bcurr(1,id_col) * (ym*mominer_/lrefe_);
      aux_CB(2,id_col) = Bcurr(2,id_col) * (sm*crosssecshear_/lrefe_);
    }
    
    stiff_glob.MultiplyTN(aux_CB,Bcurr);

    //adding geometric stiffness by shear force: N z z^T / l_n
    double aux_Q_fac = force_loc(2)*lrefe_ / pow(lcurr,2);
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
        {
          stiff_glob(id_lin,id_col) -= aux_Q_fac * rcurr(id_lin) * zcurr(id_col);
          stiff_glob(id_lin,id_col) -= aux_Q_fac * rcurr(id_col) * zcurr(id_lin);
        }

    //adding geometric stiffness by axial force: Q l_0 (r z^T + z r^T) / (l_n)^2
    double aux_N_fac = force_loc(0)/lcurr;
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
          stiff_glob(id_lin,id_col) += aux_N_fac * zcurr(id_lin) * zcurr(id_col);
    
    //shfting values from fixed size matrix to epetra matrix *stiffmatrix
    for(int i = 0; i < 6; i++)
      for(int j = 0; j < 6; j++)
        (*stiffmatrix)(i,j) = stiff_glob(i,j);
    
  }
  


  //calculating mass matrix (local version = global version)
  if (massmatrix != NULL)
  {
      //if lumped_flag == 0 a consistent mass Timoshenko beam mass matrix is applied
      if (lumpedmass == 0)
      {
        //assignment of massmatrix by means of auxiliary diagonal matrix aux_E stored as an array
        double aux_E[3]={density*lrefe_*crosssec_/6.0, density*lrefe_*crosssec_/6.0, density*lrefe_*mominer_/6.0};
        for(int id=0; id<3; id++)
        {
        	    (*massmatrix)(id,id)     = 2.0*aux_E[id];
              (*massmatrix)(id+3,id+3) = 2.0*aux_E[id];
              (*massmatrix)(id,id+3)   = aux_E[id];
              (*massmatrix)(id+3,id)   = aux_E[id];
        }
      }
      /*if lumped_flag == 1 a lumped mass matrix is applied where the cross sectional moment of inertia is
       * assumed to be approximately zero so that the 3,3 and 5,5 element are both zero */

      else if (lumpedmass == 1)
      {
        //note: this is not an exact lumped mass matrix, but it is modified in such a way that it leads
        //to a diagonal mass matrix with constant diagonal entries
        (*massmatrix)(0,0) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(1,1) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(2,2) = density*lrefe_*mominer_/2.0;
        (*massmatrix)(3,3) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(4,4) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(5,5) = density*lrefe_*mominer_/2.0;
       }
      else
        dserror("improper value of variable lumpedmass");
  }
  
  
  
  /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian(params,lm,vel,stiffmatrix,force);
   

  return;
} // DRT::ELEMENTS::Beam2::nlnstiffmass


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2

