/*!----------------------------------------------------------------------
\file randomfield_fft.cpp
Created on: 15 November, 2011
\brief Class for generating samples of gaussian and non gaussian random fields based on spectral representation
using FFT algorithms

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "gen_randomfield.H"
#include "mlmc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <cmath>
// boost currently not in use, use blitz instead
//#include <boost/random.hpp>
// include fftw++ tsuff for multidimensional FFT
#include"fftw3.h"
//using namespace DRT;
#include <fstream>


#include <boost/math/distributions/beta.hpp> // for beta_distribution.
#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/math/distributions/lognormal.hpp>
using boost::math::beta_distribution;
using boost::math::lognormal_distribution;
using boost::math::normal_distribution;



/*----------------------------------------------------------------------*/
/* standard constructor */
GenRandomField::GenRandomField(unsigned int  seed,double sigma, double corr_length,Teuchos::RCP<DRT::Discretization> discret)
{

  // Init the necessesary stuff
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  // Dimension
  dim_ = mlmcp.get<int>("RANDOM_FIELD_DIMENSION");
  if(dim_!=3&&dim_!=2)
      dserror("Dimension of random field must be 2 or 3, fix your input file");
  N_= mlmcp.get<int>("NUM_COS_TERMS");
  seed_ = seed;
  d_ = corr_length;
  sigma_0_ = sigma;
  pi_=M_PI;
  periodicity_=mlmcp.get<double>("PERIODICITY");
  M_=N_*8;
  dx_=periodicity_/M_;
  marginal_pdf_=lognormal;
  //marginal_pdf_=normal;

  // create Multidimesional array to store the values
  values_ = new double[M_*M_];
  // The StoPro will have a period of 2*pi / Deltakappa == 2*pi*N*d / 6.
  // We want this to be >= 200.
  // ceil:= return next largest integer
  //N_ = (int)ceil( 600.0 / ( pi_ * d_ ) );
// Heuristic: PSD is of `insignificant magnitude' for
   //   abs(kappa) <= 6/d

  dkappa_ = 2*pi_/(mlmcp.get<double>("PERIODICITY"));
  cout << "dkappa " << dkappa_ << endl;
  //dserror("stop herer");
  //dkappa_ = 2*pi_/N_;

   switch(dim_){
   case 3:
     Phi_0_.reserve( N_ * N_ * N_ );
     Phi_1_.reserve( N_ * N_ * N_ );
     Phi_2_.reserve( N_ * N_ * N_ );
     Phi_3_.reserve( N_ * N_ * N_ );
     break;
   case 2:
        Phi_0_.reserve( N_ * N_ );
        Phi_1_.reserve( N_ * N_ );
        break;
   default:
     dserror("Dimension of random field must be 2 or 3, fix your input file");
     break;
   }
   ComputeBoundingBox(discret);
   CreateNewPhaseAngles(seed_);
   CalcDiscretePSD();
   SimGaussRandomFieldFFT();
   TranslateToNonGaussian();
}
void GenRandomField::CreateNewSample(unsigned int seed)
{
  CreateNewPhaseAngles(seed);
  SimGaussRandomFieldFFT();
  TranslateToNonGaussian();

}


void GenRandomField::CreateNewPhaseAngles(unsigned int seed)
{

  // This defines a random number genrator
  //boost::mt19937 mt;
  // Defines random number generator for numbers between 0 and 2pi
  //boost::uniform_real<double> random( 0, 2*pi_ );

  // set seed of random number generator
  // same seed produces same string of random numbers
  //mt.seed(seed);
  //try out time in seconds since janaury
  uniformclosedgen_.seed(seed);
  switch (dim_){

  case 3:
    Phi_0_.clear();
    Phi_1_.clear();
    Phi_2_.clear();
    Phi_3_.clear();

    for ( int k5 = 0; k5 < N_ * N_ * N_; ++k5 )
      {
        //Phi_0_.push_back( random( mt ) );
        //Phi_1_.push_back( random( mt ) );
        //Phi_2_.push_back( random( mt ) );
        //.push_back( random( mt ) );
        // blitz
        Phi_0_.push_back( uniformclosedgen_.random()*2*pi_ );
        Phi_1_.push_back( uniformclosedgen_.random()*2*pi_ );
        Phi_2_.push_back( uniformclosedgen_.random()*2*pi_ );
        Phi_3_.push_back( uniformclosedgen_.random()*2*pi_ );
      }
    break;
  case 2:
    Phi_0_.clear();
    Phi_1_.clear();
    for ( int k5 = 0; k5 < N_ * N_ ; ++k5 )
         {
           //Phi_0_.push_back( random( mt ) );
           //Phi_1_.push_back( random( mt ) );
           Phi_0_.push_back( uniformclosedgen_.random()*2*pi_ );
           Phi_1_.push_back( uniformclosedgen_.random()*2*pi_ );

         }
    break;
  default:
    dserror("Dimension of random field must be 2 or 3, fix your input file");
    break;
  }

}

// compute power spectral density
void GenRandomField::CalcDiscretePSD()
{
  // check wether pdf is gaussian
  //if(marginal_pdf_==normal)
  //{
    // just compute PSD
    for (int j=0;j<N_;j++)
      {
        for (int k=0;k<N_;k++)
        {
          discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2)))));
        }
      }

  //}
  if(marginal_pdf_==lognormal)
  {
    // compute underlying gaussian distribution based on shields2011
    SpectralMatching();
    //dserror("Beta Distribution not supported yet");
  }
  else if (marginal_pdf_==normal)
  {
    cout << " Nothing to to marginal pdf gaussian " << endl;
  }
  else
  {
    dserror("Only normal and beta distribution supported fix your input file");
  }
}

// HERE comes the experimental FFT Stuff
void GenRandomField::SimGaussRandomFieldFFT()
{
// Lets see if we can speed up things with FFTW
  //int M =N_*8; // Define number of points
  // double for loops to compute coefficients
  double A; // store some stuff
  // store coefficients


  complex<double>* b1;
  complex<double>* b2;
  b1 = new complex<double>[M_*M_];
  b2 = new complex<double>[M_*M_];
  // define complex i
  complex<double> i_comp (0,1);


  for (int j=0;j<M_;j++)
  {
    for (int k=0;k<M_;k++)
    {
      //A=sqrt(2*(pow(sigma_0_,2)*d_/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
      // sort entries row major style
      // set first elements to zero
      if(k==0||j==0||j>(N_-2)||k>(N_-2))
      {
        b1[k+M_*j]=0.0;
        b2[k+M_*j]=0.0;
      }
      else
      {
        //A=sqrt(2*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
        A=sqrt(2*(discrete_PSD_[k+j*N_]*(pow(dkappa_,2))));
        real(b1[k+M_*j])=A*sqrt(2)*cos(Phi_1_[k+M_*j]);
        imag(b1[k+M_*j])= A*sqrt(2)*sin(Phi_1_[k+M_*j]);
        real(b2[k+M_*j])= A*sqrt(2)*cos(Phi_2_[k+M_*j]);
        imag(b2[k+M_*j])= A*sqrt(2)*sin(Phi_2_[k+M_*j]);
      }
     }
  }



  int rank = 1; /* not 2: we are computing 1d transforms */
  /* 1d transforms of length M_ */
  int N_fftw = M_;
  int howmany = M_; // same here
  int idist = M_;
  int    odist = M_;
  int istride =1;
  int ostride = 1; /* distance between two elements in the same column */
  //int *inembed = n;
  //int *onembed = n;

  // allocate output arrays
  complex<double>* d1;
  complex<double>* d2;
  d1 = new complex<double>[M_*M_];
  d2 = new complex<double>[M_*M_];
  fftw_plan ifft_of_rows;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_collums;

  ifft_of_rows = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(b1)),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d1)),
                  NULL,
                  ostride,
                  odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);
  ifft_of_rows2 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(b2)),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d2)),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
  //ifft_of_rows = fftw_plan_dft_1d(1024, (reinterpret_cast<fftw_complex*>(b1)), (reinterpret_cast<fftw_complex*>(b1)), FFTW_FORWARD, FFTW_ESTIMATE);

  istride =M_;
  ostride=M_;
  idist=1;
  odist=1;
  ifft_of_collums = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(d1)),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d2)),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);

  fftw_execute(ifft_of_rows);
  fftw_execute(ifft_of_rows2);

  complex<double> scaling (M_,M_);
  // transpose d1
   for (int k=0;k<M_*M_;k++)
   {
     d2[k]=conj(d2[k]);
     d1[k]=d1[k]+d2[k];
   }


   fftw_execute(ifft_of_collums);

   // move values into class variable
   for(int i=0;i<M_*M_;i++)
   {
     values_[i]=real(d2[i]);
   }

   // free memory
   delete b1;
   delete b2;
   delete d1;
   delete d2;
   fftw_destroy_plan(ifft_of_rows);
   fftw_destroy_plan(ifft_of_collums);

}
void GenRandomField::ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret)
{
  // root bounding Box
  vector<double> maxrbb;
  maxrbb.push_back(10.0e-19);
  maxrbb.push_back(10.0e-19);
  maxrbb.push_back(10.0e-19);
 // maxrbb.push_back(10.0e-19);
  vector<double> minrbb;
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);


  bb_max_.push_back(10.0e-19);
  bb_max_.push_back(10.0e-19);
  bb_max_.push_back(10.0e-19);

  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  {
      for (int lid = 0; lid <discret->NumMyColNodes(); ++lid)
      {
        const DRT::Node* node = discret->lColNode(lid);
        // check if greater than maxrbb
        if (maxrbb[0]<node->X()[0])
          maxrbb[0]=node->X()[0];
        if (maxrbb[1]<node->X()[1])
          maxrbb[1]=node->X()[1];
        if (maxrbb[2]<node->X()[2])
           maxrbb[2]=node->X()[2];
        // check if smaller than minrbb
        if (minrbb[0]>node->X()[0])
          minrbb[0]=node->X()[0];
        if (minrbb[1]>node->X()[1])
          minrbb[1]=node->X()[1];
        if (minrbb[2]>node->X()[2])
          minrbb[2]=node->X()[2];
      }


  }

  discret->Comm().MaxAll(&maxrbb[0],&bb_max_[0],3);
  discret->Comm().MinAll(&minrbb[0],&bb_min_[0],3);

  discret->Comm().Barrier();
  const int myrank = discret->Comm().MyPID();
  if (myrank == 0)
  {
    cout << "min " << bb_min_[0] << " "<< bb_min_[1]  << " "<< bb_min_[2] << endl;
    cout << "max " << bb_max_[0] << " "<< bb_max_[1]  << " "<< bb_max_[2] << endl;
  }

}
double GenRandomField::EvalFieldAtLocation(vector<double> location, bool writetofile)
{
  int index_x;
  int index_y;
  int index_z;

  // Compute indices
  index_x=int(floor((location[0]-bb_min_[0])/dx_));
  // HACH SET z to y
  //cout << "hack in use" << endl;
  index_y=int(floor((location[1]-bb_min_[1])/dx_));
  index_z=int(floor((location[2]-bb_min_[2])/dx_));
  if (writetofile)
  {
    ofstream File;
    File.open("RFatPoint.txt",ios::app);
    File << setprecision (9) << values_[index_x+M_*index_y]<< endl;
    File.close();
  }

  return values_[index_x+M_*index_y];


}
// Translate Gaussian to nonGaussian process based on Mircea Grigoriu's translation process
// theory
void GenRandomField::TranslateToNonGaussian()
{
  normal_distribution<>  my_norm(0,1);
  // check wether pdf is gaussian
  switch(marginal_pdf_)
  {

    case normal:
      //cout << RED_LIGHT << "WARNING: Target marginal PDF is gaussian so nothing to do here"<< END_COLOR << endl;
    break;

    case beta:
    {
      // init a prefixed beta distribution
      double a_beta = 4.0;
      double b_beta = 2.0;
      beta_distribution<> my_beta(a_beta,b_beta);
      //translate
      for(int i=0;i<M_*M_;i++)
      {
        values_[i]=quantile(my_beta,cdf(my_norm, values_[i]))*5.61-3.74;;
      }
    }
    break;

    case lognormal:
    {
      // init params for logn distribution based on shields paper
      double mu_bar=1.8;
      double sigma_N = sqrt(log(1+1/(pow(1.8,2))));
      double mu_N = log(mu_bar)-pow(sigma_N,2)/2;
      lognormal_distribution<>  my_lognorm(mu_N,sigma_N);
      //Translate the data
      for(int i=0;i<M_*M_;i++)
      {
        values_[i]=quantile(my_lognorm,cdf(my_norm, values_[i]))-1.8;
      }
    }
    break;
    default:
      dserror("Only lognormal and beta distribution supported fix your input file");
  }
}

// Transform PSD of underlying gauusian process
void GenRandomField::SpectralMatching()
{
  // PSD of underlying gaussian field
  //vector<double> PSD_ul_g;
  // Target PSD of non gassian field
  vector<double> PSD_ng_target;

  vector<double> PSD_ng(N_*N_,0.0);
  vector<double> PSD_ul_g(N_*N_,0.0);

  //vector<double> autocorr_ng(2*N_*2*N_,0.0);
  vector<double> rho(2*N_*2*N_,0.0);

  PSD_ng_target=discrete_PSD_;
  PSD_ng_target[0]=0.0;
  complex<double>* autocorr;
  autocorr = new complex<double>[N_*2*N_*2];
  complex<double>* almost_autocorr;
  almost_autocorr = new complex<double>[N_*2*N_*2];
  complex<double>* autocorr_ng;
  autocorr_ng= new complex<double>[N_*2*N_*2];

  //temp complex storage
  complex<double>* PSD_ul_g_complex;
  PSD_ul_g_complex = new complex<double>[2*N_*2*N_];
  complex<double>* almost_PSD_ng_complex;
  almost_PSD_ng_complex= new complex<double>[2*N_*2*N_];
  complex<double>* PSD_ng_complex;
  PSD_ng_complex= new complex<double>[2*N_*2*N_];

  for (int j=0;j<N_*2;j++)
    {
      for (int k=0;k<N_*2;k++)
      {
        //A=sqrt(2*(pow(sigma_0_,2)*d_/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
        // sort entries ro w major style
        // set first elements to zero
        //if(k==0||j==0||j>(N_-1)||k>(N_-1))
        if(j>(N_-1)||k>(N_-1))
        {
          real(PSD_ul_g_complex[k+N_*2*j])=0.0;
          imag(PSD_ul_g_complex[k+N_*2*j])=0.0;
        }
        else
        {
          real(PSD_ul_g_complex[k+N_*2*j])=discrete_PSD_[k+j*N_];
          imag(PSD_ul_g_complex[k+N_*2*j])=0.0;
          if(k==0||j==0)
          {
           real(PSD_ul_g_complex[k+N_*2*j])=real(PSD_ul_g_complex[k+N_*2*j])*0.5;
           imag(PSD_ul_g_complex[k+N_*2*j])=0.0;
          }
          else if (k==0&&j==0)
          {
            real(PSD_ul_g_complex[k+N_*2*j])=0;
            imag(PSD_ul_g_complex[k+N_*2*j])=0;
          }
        }
      }
    }

  for(int i=0;i<(2*N_*N_*2);i++)
  {
    real(autocorr[i])=0.0;
    imag(autocorr[i])=0.0;
    real(almost_autocorr[i])=0.0;
    imag(almost_autocorr[i])=0.0;
  }


   // TWO DIM FFTS for spectral matching
   int rank = 1; /* not 2: we are computing 1d transforms */
   //int n[] = {1024}; /* 1d transforms of length 10 */
   int N_fftw = 2*N_;
   int howmany = 2*N_; // same here
   int idist = 2*N_;
   int    odist = 2*N_;
   int istride =1;
   int ostride = 1; /* distance between two elements in the same column */


   fftw_plan ifft_of_rows_of_psd;
   fftw_plan ifft_of_columns_of_psd;
  // ifft for autocorr
  fftw_plan ifft_of_rows_of_autocorr_ng;
  fftw_plan ifft_of_columns_of_autocorr_ng;

   ifft_of_rows_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(PSD_ul_g_complex)),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(almost_autocorr)),
                   NULL,
                   ostride,
                   odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);

   ifft_of_rows_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(autocorr_ng)),
          NULL,istride, idist,(reinterpret_cast<fftw_complex*>(almost_PSD_ng_complex)),
                      NULL,
                      ostride,
                      odist,
                     FFTW_BACKWARD,
                     FFTW_ESTIMATE);


   istride =2*N_;
   ostride=2*N_;
   idist=1;
   odist=1;
   ifft_of_columns_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(almost_autocorr)),
        NULL,istride, idist,(reinterpret_cast<fftw_complex*>(autocorr)),
                    NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
   ifft_of_columns_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(almost_PSD_ng_complex)),
           NULL,istride, idist,(reinterpret_cast<fftw_complex*>(PSD_ng_complex)),
                       NULL,
                      ostride,
                      odist,
                      FFTW_BACKWARD,
                      FFTW_ESTIMATE);
 // end of two dim ffts for spectral matching


  // iteration without errorcheck for now

  for(int i =0;i<10;i++)
  {
    if(i!=0)// set new psd_ul_g
    {
      for (int j=0;j<N_*2;j++)
         {
           for (int k=0;k<N_*2;k++)
           {
             if(j>(N_-1)||k>(N_-1))
             {
               real(PSD_ul_g_complex[k+N_*2*j])=0.0;
               imag(PSD_ul_g_complex[k+N_*2*j])=0.0;
             }
             else
             {
               real(PSD_ul_g_complex[k+N_*2*j])=PSD_ul_g[k+j*N_];
               imag(PSD_ul_g_complex[k+N_*2*j])=0.0;
               if(k==0||j==0)
               {
               // real(PSD_ul_g_complex[k+N_*2*j])=real(PSD_ul_g_complex[k+N_*2*j])*0.5;
               }
             }
           }
         }
    }

    fftw_execute(ifft_of_rows_of_psd);
   // delete real part
    for(int k=0;k<2*N_*2*N_;k++)
    {
      //imag(almost_autocorr[k])=0.0;
      //real(almost_autocorr[k])= real(almost_autocorr[k])*2;
    }

    fftw_execute(ifft_of_columns_of_psd);
    double scaling_fac= 2*pi_;
     scaling_fac=dkappa_*N_;
         //scaling_fac=1.5707963;
    // loop over vectorlength
    for(int k=0;k<2*N_*2*N_;k++)
    {
      //cout << "N_ " << N_ << endl;

     rho[k] = real(autocorr[k])*2*2*pow((scaling_fac),2)/(2*N_*2*N_*(pow(sigma_0_,2)));
     real(autocorr_ng[k])=Integrate(-18.0,7.0,-18.0,7.0,rho[k]);
    }
    // test
    //autocorr_ng[0]=0;

    fftw_execute(ifft_of_rows_of_autocorr_ng);
    fftw_execute(ifft_of_columns_of_autocorr_ng);
    for (int j=0;j<N_*2;j++)
        {
          for (int k=0;k<N_*2;k++)
          {
            // sort entries ro w major style
            // set first elements to zero
            //if(k==0||j==0||j>(N_-1)||k>(N_-1))
            if(j>(N_-1)||k>(N_-1))
            {}
            //else if (j==0&& k== 0)
             // PSD_ng[k+j*N_]=0.0;
            else
            {
              PSD_ng[k+j*N_]=real(PSD_ng_complex[k+j*2*N_])/(pow(scaling_fac,2));//*(pow((N_*2.0),2));
              //PSD_ng[k+j*N_]=real(PSD_ng_complex[k+j*2*N_])/(4*pow(pi_,2));//*(pow((N_*2.0),2));
              PSD_ul_g[k+j*N_]=real(PSD_ul_g_complex[k+j*2*N_]);
            }
          }
        }

    PSD_ng[0]=0.0;
    for(int k=0;k<N_*N_;k++)
    {
      if(PSD_ng[k]>10e-10)
        PSD_ul_g[k]=pow(PSD_ng_target[k]/PSD_ng[k],1.4)*PSD_ul_g[k];
       // do not set to zero because if once zero you'll never get it non-zero again
      else
        PSD_ul_g[k]=10e-10;
    }


         ofstream File2;
            File2.open("Psd_coplex.txt",ios::app);
            for (int j=0;j<2*N_;j++)
                    {
                      for (int k=0;k<2*N_;k++)
                      {
                        File2 << k << " " << j << " " << real(PSD_ng_complex[k+j*2*N_])/(4*pow(pi_,2)) << endl;

                      }

                    }
             File2.close();

  }



  // free memory
    delete autocorr;
    delete almost_autocorr;
    delete PSD_ul_g_complex;
    delete PSD_ng_complex;
    fftw_destroy_plan(ifft_of_columns_of_psd);
    fftw_destroy_plan(ifft_of_rows_of_psd);
    fftw_destroy_plan(ifft_of_rows_of_autocorr_ng);
    fftw_destroy_plan(ifft_of_columns_of_autocorr_ng);


    // Write PSD_ul_g PSD_ng_target and PSD_ng to a file
    // Dimensioj is 128* 128

     for(int h=0;h<N_*N_;h++)
     {
       if(PSD_ng[h]>10e-10)// change that
         discrete_PSD_[h]=PSD_ul_g[h];
       else
         //remove all the very small entries to get rid of the wiggles
         discrete_PSD_[h]=0.0;
     }
     ofstream File;
        File.open("Psd.txt",ios::app);
        for (int j=0;j<N_;j++)
                {
                  for (int k=0;k<N_;k++)
                  {
                    File << k << " " << j << " " << PSD_ng_target[k+j*N_] << " " << PSD_ng[k+j*N_]<< " "  << PSD_ul_g[k+j*N_] << endl;

                  }

                }
         File.close();


}

// Routine to calculate
double GenRandomField::Integrate(double xmin, double xmax, double ymin, double ymax, double rho)
{
  // get trillios gausspoints with hig order
  Teuchos::RCP<DRT::UTILS::GaussPoints> gp = DRT::UTILS::GaussPointCache::Instance().Create( DRT::Element::quad4, 30 );
  // needed for transformation in [-1,1];[-1,1] space
  double hx= abs(xmax-xmin);
  double hy= abs(ymax-ymin);
  double jdet= hx*hy/4;
  double integral_value =0.0;

  for (int i=0; i<gp->NumPoints(); i++)
  {
    integral_value+=gp->Weight(i)*jdet*Testfunction(xmin+hx/2*(1+gp->Point(i)[0]),ymin+hy/2*(1+gp->Point(i)[1]),rho);
  }
  return integral_value;
}

double GenRandomField::Testfunction(double argument_x ,double argument_y, double rho)
{
  double result;
  double mu_bar = 1.8;
  double sigma_N = sqrt(log(1+1/(pow(1.8,2))));
  double mu_N = log(mu_bar)-pow(sigma_N,2)/2;
  lognormal_distribution<>  my_lognorm(mu_N,sigma_N);
  normal_distribution<> my_normal(0,1);
  // calc function value
  result=quantile(my_lognorm,(cdf(my_normal,argument_x)))*quantile(my_lognorm,(cdf(my_normal,argument_y)))*
      (1/(2*pi_*pow(sigma_0_,2)*sqrt(1-pow(rho,2))))*exp(-(pow(argument_x,2)+pow(argument_y,2)-2*rho*argument_x*argument_y)
          /(2*pow(sigma_0_,2)*(1-pow(rho,2))));
  return result;
}

// Write Random Field to file
void GenRandomField::WriteRandomFieldToFile()
{
  ofstream File;
  File.open("RandomField.txt",ios::out);
  for(int i=0;i<N_*N_;i++)
  {
    File << values_[i]<< endl;
  }
  File.close();
}

#endif
