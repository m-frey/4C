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


#include <boost/random.hpp>
// include fftw++ stuff for multidimensional FFT
#include"fftw3.h"
//using namespace DRT;
#include <fstream>


#include <boost/math/distributions/beta.hpp> // for beta_distribution.
#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/math/distributions/lognormal.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/accumulators.hpp>


#include "../drt_inpar/inpar_mlmc.H"
//include <boost/math/distributions.hpp>
using boost::math::beta_distribution;
using boost::math::lognormal_distribution;
using boost::math::normal_distribution;
using  boost::accumulators::variance;
using  boost::accumulators::mean;
using  boost::accumulators::stats;

/*----------------------------------------------------------------------*/
/* standard constructor */
GenRandomField::GenRandomField(unsigned int  seed,Teuchos::RCP<DRT::Discretization> discret)
{
  init_psd_=false;
   myrank_ = discret->Comm().MyPID();
  // Init the necessesary stuff
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  // Dimension
  dim_ = mlmcp.get<int>("RANDOM_FIELD_DIMENSION");
  if(dim_!=3&&dim_!=2)
      dserror("Dimension of random field must be 2 or 3, fix your input file");
  N_= mlmcp.get<int>("NUM_COS_TERMS");
  kappa_u_=mlmcp.get<double>("KAPPA_U");
  seed_ = seed;
  d_ = mlmcp.get<double>("CORRLENGTH");
  sigma_0_= mlmcp.get<double>("SIGMA");
  // Sigma of target nongaussian pdf
  sigma_ul_g_cur_it_ = 0.0;
  pi_=M_PI;
  M_=mlmcp.get<int>("SIZE_PER_DIM");
  dkappa_=kappa_u_/N_;
  periodicity_=2.*pi_/dkappa_;
  dx_=periodicity_/M_;

  // distribution parameters of non gauss pdf
  distribution_params_.push_back(mlmcp.get<double>("NONGAUSSPARAM1"));
  distribution_params_.push_back(mlmcp.get<double>("NONGAUSSPARAM2"));
  // Get correlation structure
  INPAR::MLMC::CorrStruct cstruct = DRT::INPUT::IntegralValue<INPAR::MLMC::CorrStruct>(mlmcp,"CORRSTRUCT");
  switch(cstruct){
    case INPAR::MLMC::corr_gaussian:
      //blebla
      break;
    default:
      dserror("Unknown Correlation structure");
  }
  double upper_bound;
  INPAR::MLMC::MarginalPdf mpdf = DRT::INPUT::IntegralValue<INPAR::MLMC::MarginalPdf>(mlmcp,"MARGINALPDF");
  switch(mpdf){
    case INPAR::MLMC::pdf_gaussian:
      marginal_pdf_=normal;
      // Hack
      cout << "remove a line here" << endl;
      sigma_0_ =sqrt((exp(pow(distribution_params_[1],2))-1)*exp(2*distribution_params_[0]+pow(distribution_params_[1],2)));
      break;
    case INPAR::MLMC::pdf_beta:
      marginal_pdf_=beta;
      // compute bounds of distribution
      // using mu_b = 0 and sigma_b = 1 the lower and upper bounds can be computed according to Yamazaki1988
      //lower bound
      distribution_params_.push_back(-1.0*sqrt(distribution_params_[0]*(distribution_params_[0]+distribution_params_[1]+1)/distribution_params_[1]));
      // upper bound
      upper_bound =(-1.0*sqrt(distribution_params_[1]*(distribution_params_[0]+distribution_params_[1]+1)/distribution_params_[0]));
      // I need abs(lB)+abs(uB) hence
      distribution_params_.push_back(abs(upper_bound)+abs(distribution_params_[2]));
      cout << "Distribution parameter of beta distribution " << distribution_params_[0] << " "  << distribution_params_[1] << " " << distribution_params_[2] << " " << distribution_params_[3] << endl;
      break;
    case INPAR::MLMC::pdf_lognormal:
      marginal_pdf_=lognormal;
      // Calculate mean of lognormal distribution based mu_N and sigma_N
      distribution_params_.push_back(exp(distribution_params_[0]+0.5*pow(distribution_params_[1],2)));
      // also calc variance and sigma
      //(exp(s^2 )-1) * exp(2m + s^2)
      sigma_0_ =sqrt((exp(pow(distribution_params_[1],2))-1)*exp(2*distribution_params_[0]+pow(distribution_params_[1],2)));
      cout << "sigma_0 "<< sigma_0_ << endl;
      if (myrank_ == 0)
        cout << "Distribution parameter of lognormal distribution " << distribution_params_[0] << " "  << distribution_params_[1] << " " << distribution_params_[2] << endl;
      break;
    default:
      dserror("Unknown Marginal pdf");
  }

  // create Multidimesional array to store the values
  double mydim= dim_;
  // transform needed because pow does not like to get two ints
  int size_of_field= int(pow(M_,mydim));
  values_ = new double[size_of_field];
  // The StoPro will have a period of 2*pi / Deltakappa == 2*pi*N*d / 6.
  // We want this to be >= 200.
  // ceil:= return next largest integer
  //N_ = (int)ceil( 600.0 / ( pi_ * d_ ) );
// Heuristic: PSD is of `insignificant magnitude' for
   //   abs(kappa) <= 6/d



  //new formular

  //dkappa_=pi_/(2*N_);
  //dkappa_=pi_/(3*N_);

  //dkappa_ = 2*pi_/(mlmcp.get<double>("PERIODICITY"));
  if (myrank_ == 0)
  {
    cout << "Random Field Parameters "<< endl;
    cout << "Periodicity L: " << periodicity_ << endl;
    cout << "M: " << M_ << endl;
    cout << "N: " << N_ << endl;
    cout << "kappa_u: " << kappa_u_ << endl;
    cout << "dkappa: " << dkappa_ << endl;
    cout << "dx " << dx_ << endl;
  }
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
  switch(dim_)
  {
    case 3:
      CalcDiscretePSD3D();
      SimGaussRandomFieldFFT3D();
      for (int i= 0 ; i<10000; i++)
       {
        CreateNewPhaseAngles(seed_+i);
        //SimGaussRandomFieldFFT3D();
        //TranslateToNonGaussian();
        ofstream File;
        File.open("RFatPoint.txt",ios::app);
        // use at() to get an error massage just in case
        //File << setprecision (9) << values_[2323]<< endl;
        File << setprecision (9) << TESTSimGaussRandomField3D(5.0, 5.0, 5.0)<< endl;

        File.close();
        }
      dserror("You did it");
   //   EOF 3D testing
      break;
    case 2:
      CalcDiscretePSD();
      SimGaussRandomFieldFFT();
      // testing
            for (int i= 0 ; i<10000; i++)
             {
              CreateNewPhaseAngles(seed_+i);
              SimGaussRandomFieldFFT();
              TranslateToNonGaussian();
              ofstream File;
              File.open("RFatPoint.txt",ios::app);
              // use at() to get an error massage just in case
              File << setprecision (9) << values_[2323]<< endl;
              File.close();
              }
            dserror("You did it");
            //EOF 3D testing
      break;
    default:
      dserror("Dimension of random field must be 2 or 3, fix your input file");
      break;
  }
  WriteRandomFieldToFile();
  TranslateToNonGaussian();


}
void GenRandomField::CreateNewSample(unsigned int seed)
{
  CreateNewPhaseAngles(seed);
  switch(dim_)
    {
      case 3:
        SimGaussRandomFieldFFT3D();
        break;
      case 2:
        SimGaussRandomFieldFFT();;
        break;
      default:
        dserror("Dimension of random field must be 2 or 3, fix your input file");
        break;
    }
    TranslateToNonGaussian();
}


void GenRandomField::CreateNewPhaseAngles(unsigned int seed)
{

  // This defines a random number genrator
  boost::mt19937 mt;
  // Defines random number generator for numbers between 0 and 2pi
  boost::uniform_real<double> random( 0, 2*pi_ );

  // set seed of random number generator
  // same seed produces same string of random numbers
  mt.seed(seed);
  switch (dim_){

  case 3:
    Phi_0_.clear();
    Phi_1_.clear();
    Phi_2_.clear();
    Phi_3_.clear();

    for ( int k5 = 0; k5 < N_ * N_ * N_; ++k5 )
      {
        Phi_0_.push_back( random( mt ) );
        Phi_1_.push_back( random( mt ) );
        Phi_2_.push_back( random( mt ) );
        Phi_3_.push_back( random( mt ) );
      }
    break;
  case 2:
    Phi_0_.clear();
    Phi_1_.clear();
    for ( int k5 = 0; k5 < N_ * N_ ; ++k5 )
         {
           Phi_0_.push_back( random( mt ) );
           Phi_1_.push_back( random( mt ) );

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
  // just compute PSD
  cout << "sigma_0_" << sigma_0_ << endl;
  for (int j=0;j<N_;j++)
  {
    for (int k=0;k<N_;k++)
    {
      if(k==0||j==0)
      {
       discrete_PSD_.push_back(0.5*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2)))));
      }
      else
      {
      //discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/(2*sqrt(pi_)),2))-(pow(d_*k*dkappa_/(2*sqrt(pi_)),2)))));
      discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2)))));

      }
    }
  }
  // Write to file
  if (myrank_ == 0)
    {
      ofstream File;
      File.open("DiscretePSD.txt",ios::out);
      int size = int (pow(N_,2.0));
      for(int i=0;i<size;i++)
      {
        File << discrete_PSD_[i]<< endl;
      }
      File.close();
    }

  if(marginal_pdf_!=normal)
  {
    // compute underlying gaussian distribution based on shields2011
    //cout << "NO SPECTRAL MATHCING"<< endl;
   SpectralMatching();
  }
  else
  {
    if (myrank_ == 0)
      cout << " Nothing to do marginal pdf gaussian " << endl;
  }
}

void GenRandomField::CalcDiscretePSD3D()
{
  // just compute PSD
   cout << "sigma_0_" << sigma_0_ << endl;
  // just compute PSD
  for (int j=0;j<N_;j++)
  {
    for (int k=0;k<N_;k++)
    {
      for (int l=0;l<N_;l++)
      {
        if(k==0||j==0||l==0)
        {
          discrete_PSD_.push_back(0.25*(pow(sigma_0_,2)*pow(d_,3)/(pow((2*sqrt(pi_)),3))*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))-(pow(d_*l*dkappa_/2,2)) )));
        }
        else
        {
        //discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/(2*sqrt(pi_)),2))-(pow(d_*k*dkappa_/(2*sqrt(pi_)),2)))));
        discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,3)/(pow((2*sqrt(pi_)),3))*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))-(pow(d_*l*dkappa_/2,2)) )));
        }
      }
    }
  }


  if(marginal_pdf_!=normal)
  {
    // compute underlying gaussian distribution based on shields2011
   // cout << "NO SPECTRAL MATHCING"<< endl;

   SpectralMatching3D();
  }
  else
  {
    if (myrank_ == 0)
      cout << " Nothing to do marginal pdf gaussian " << endl;
  }
}


void GenRandomField::SimGaussRandomFieldFFT()
{
  double A; // store some stuff
  // store coefficients
  Teuchos::RCP<Teuchos::Array <complex<double> > > b1=  Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > b2=  Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));

  Teuchos::RCP<Teuchos::Array <complex<double> > > d1= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d2= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));

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
        (*b1)[k+M_*j]=0.0;
        (*b2)[k+M_*j]=0.0;
      }
      else
      {
        //A=sqrt(2*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
        A=sqrt(2*(discrete_PSD_[k+j*N_]*(pow(dkappa_,2))));
        real((*b1)[k+M_*j])=A*sqrt(2)*cos(Phi_0_[k+N_*j]);
        imag((*b1)[k+M_*j])= A*sqrt(2)*sin(Phi_0_[k+N_*j]);
        real((*b2)[k+M_*j])= A*sqrt(2)*cos(Phi_1_[k+N_*j]);
        imag((*b2)[k+M_*j])= A*sqrt(2)*sin(Phi_1_[k+N_*j]);
        // last 4 lines was M_ before
      }
    }
  }


  int rank = 1; /* not 2: we are computing 1d transforms 1d transforms of length M_ */
  int N_fftw = M_;
  int howmany = M_; // same here
  int idist = M_; //  the distance in memory between the first element of the first array and the  first element of the second array */
  int odist = M_;
  int istride =1;
  int ostride = 1; /* distance between two elements in the same row/column/rank */
  //int *inembed = n;
  //int *onembed = n;


  fftw_plan ifft_of_rows;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_collums;
  // Why we have to do this &((*b1)[0]))
  // the fftw requires the adress of the first element of a fftw_complex array. We can  get to a normal array by casting.
  // Since we do not want to use standard arrays with new double[size] but RCP::Teuchos::Arrays instead we have to do (*b1)[0]
  // to get the first element of the array and than &(..) to get the adress

  ifft_of_rows = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
                  NULL,
                  ostride,
                  odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);
  ifft_of_rows2 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b2)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d2)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);

  istride =M_;
  ostride=M_;
  idist=1;
  odist=1;

  ifft_of_collums = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d2)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);

  fftw_execute(ifft_of_rows);
  fftw_execute(ifft_of_rows2);
  //complex<double> scaling (M_,M_);
  // transpose d1
  for (int k=0;k<M_*M_;k++)
  {
    (*d2)[k]=conj((*d2)[k]);
    (*d1)[k]=(*d1)[k]+(*d2)[k];
  }
  fftw_execute(ifft_of_collums);

  for(int i=0;i<M_*M_;i++)
  {
    values_[i]=real((*d2)[i]);
  }
  // free memory
  fftw_destroy_plan(ifft_of_rows);
  fftw_destroy_plan(ifft_of_rows2);
  fftw_destroy_plan(ifft_of_collums);
}
void GenRandomField::SimGaussRandomFieldFFT3D()
{
  double A; // store some stuff
  // store coefficients
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > b1= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > b2= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > b3= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > b4= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > d1= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > d2= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d3= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d4= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //cout << "classic line " << __LINE__ << endl;
  Teuchos::RCP<Teuchos::Array <complex<double> > > d5= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d6= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d7= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  // define complex i
  complex<double> i_comp (0,1);
  for (int j=0;j<M_;j++)
  {
    for (int k=0;k<M_;k++)
    {
      for (int l=0;l<M_;l++)
      {

        //A=sqrt(2*(pow(sigma_0_,2)*d_/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
        // sort entries row major style
        // set first elements to zero
        if(k==0||j==0||l==0||j>(N_-2)||k>(N_-2)||l>(N_-2))
        {
          (*b1)[l+M_*(k+M_*j)]=0.0;
          (*b2)[l+M_*(k+M_*j)]=0.0;
          (*b3)[l+M_*(k+M_*j)]=0.0;
          (*b4)[l+M_*(k+M_*j)]=0.0;
        }
        else
        {
          //A=sqrt(2*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
          A=sqrt(2*(discrete_PSD_[l+N_*(k+N_*j)]*(pow(dkappa_,3))));
          real((*b1)[l+M_*(k+M_*j)])= A*sqrt(2)*cos(Phi_0_[l+N_*(k+N_*j)]);
          imag((*b1)[l+M_*(k+M_*j)])= A*sqrt(2)*sin(Phi_0_[l+N_*(k+N_*j)]);
          real((*b2)[l+M_*(k+M_*j)])= A*sqrt(2)*cos(Phi_1_[l+N_*(k+N_*j)]);
          imag((*b2)[l+M_*(k+M_*j)])= A*sqrt(2)*sin(Phi_1_[l+N_*(k+N_*j)]);
          real((*b3)[l+M_*(k+M_*j)])= A*sqrt(2)*cos(Phi_2_[l+N_*(k+N_*j)]);
          imag((*b3)[l+M_*(k+M_*j)])= A*sqrt(2)*sin(Phi_2_[l+N_*(k+N_*j)]);
          real((*b4)[l+M_*(k+M_*j)])= A*sqrt(2)*cos(Phi_3_[l+N_*(k+N_*j)]);
          imag((*b4)[l+M_*(k+M_*j)])= A*sqrt(2)*sin(Phi_3_[l+N_*(k+N_*j)]);
      }
      }
    }
  }
  // set up FFTW Plans
  int rank = 1; /* not 2: we are computing 1d transforms  1d transforms of length M_*/
  int N_fftw = M_;
  int howmany = M_*M_; // same here
  int idist = M_; //  the distance in memory between the first element of the first array and the  first element of the second array */
  int odist = M_;
  int istride =1;
  int ostride = 1; /* distance between two consecutive elements */
  // four plans for the rows
  fftw_plan ifft_of_rows1;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_rows3;
  fftw_plan ifft_of_rows4;
  // two for the collumns
  fftw_plan ifft_of_collums1;
  fftw_plan ifft_of_collums2;
  // and one for the rank ( third dim)
  fftw_plan ifft_of_rank;

  ifft_of_rows1 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
                  NULL,
                  ostride,
                  odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);
  ifft_of_rows2 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b2)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d2)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
  ifft_of_rows3 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b3)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d3)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
   ifft_of_rows4 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b4)[0]))),
        NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d4)[0]))),
                    NULL,
                    ostride,
                    odist,
                    FFTW_BACKWARD,
                    FFTW_ESTIMATE);
   // start here we need to transpose the arrays first
   // stride and dist are the same for collumns because arrays are transposed before FFTs
   istride =1;
   ostride=1;
   idist=M_;
   odist=M_;
   ifft_of_collums1= fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
         NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d5)[0]))),
                     NULL,
                     ostride,
                     odist,
                     FFTW_BACKWARD,
                     FFTW_ESTIMATE);
   ifft_of_collums2= fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*d3)[0]))),
            NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d6)[0]))),
                        NULL,
                        ostride,
                        odist,
                        FFTW_BACKWARD,
                        FFTW_ESTIMATE);
   // and now the rank
   istride =M_*M_;
   ostride=M_*M_;
   idist=1;
   odist=1;
   ifft_of_rank= fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*d5)[0]))),
               NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d7)[0]))),
                           NULL,
                           ostride,
                           odist,
                           FFTW_BACKWARD,
                           FFTW_ESTIMATE);
   fftw_execute(ifft_of_rows1);
   fftw_execute(ifft_of_rows2);
   fftw_execute(ifft_of_rows3);
   fftw_execute(ifft_of_rows4);
   for (int k=0;k<M_*M_*M_;k++)
    {
      (*d2)[k]=conj((*d2)[k]);
      (*d4)[k]=conj((*d4)[k]);
      (*d2)[k]=(*d1)[k]+(*d2)[k];
      (*d4)[k]=(*d3)[k]+(*d4)[k];
      //old
//      (*d2)[k]=conj((*d2)[k]);
//      (*d4)[k]=conj((*d4)[k]);
//      (*d1)[k]=(*d1)[k]+(*d2)[k];
//      (*d3)[k]=(*d3)[k]+(*d4)[k];

    }
   // start here we need to transpose the arrays first
   // We need cannot do a pure 1D decomposition of 3D FFT
   // transpose
    for (int j=0;j<M_;j++)
       {
         for (int k=0;k<M_;k++)
         {
           for (int l=0;l<M_;l++)
           {
             real((*d1)[l+(M_)*(k+(M_)*j)])=real((*d2)[k+(M_)*(l+(M_)*j)]);
             imag((*d1)[l+(M_)*(k+(M_)*j)])=imag((*d2)[k+(M_)*(l+(M_)*j)]);
             real((*d3)[l+(M_)*(k+(M_)*j)])=real((*d4)[k+(M_)*(l+(M_)*j)]);
             imag((*d3)[l+(M_)*(k+(M_)*j)])=imag((*d4)[k+(M_)*(l+(M_)*j)]);
           }
         }
       }
   fftw_execute(ifft_of_collums1);
   fftw_execute(ifft_of_collums2);
   for (int j=0;j<M_;j++)
   {
     for (int k=0;k<M_;k++)
     {
       for (int l=0;l<M_;l++)
       {
         real((*d1)[l+(M_)*(k+(M_)*j)])=real((*d5)[k+(M_)*(l+(M_)*j)]);
         imag((*d1)[l+(M_)*(k+(M_)*j)])=imag((*d5)[k+(M_)*(l+(M_)*j)]);
         real((*d3)[l+(M_)*(k+(M_)*j)])=real((*d6)[k+(M_)*(l+(M_)*j)]);
         imag((*d3)[l+(M_)*(k+(M_)*j)])=imag((*d6)[k+(M_)*(l+(M_)*j)]);
       }
     }
   }
   // and back here

   for (int k=0;k<M_*M_*M_;k++)
       {
         (*d6)[k]=conj((*d3)[k]);
         (*d5)[k]=(*d1)[k]+(*d6)[k];
         // old
//         (*d6)[k]=conj((*d6)[k]);
//         (*d5)[k]=(*d5)[k]+(*d6)[k];
       }
   fftw_execute(ifft_of_rank);
   for(int i=0;i<M_*M_*M_;i++)
     {
       values_[i]=real((*d7)[i]);
       //cout << "d7 " << real((*d7)[i])<< endl;
     }
   fftw_destroy_plan(ifft_of_rows1);
   fftw_destroy_plan(ifft_of_rows2);
   fftw_destroy_plan(ifft_of_rows3);
   fftw_destroy_plan(ifft_of_rows4);
   fftw_destroy_plan(ifft_of_collums1);
   fftw_destroy_plan(ifft_of_collums2);
   fftw_destroy_plan(ifft_of_rank);
}
double GenRandomField::TESTSimGaussRandomField3D(double x, double y, double z)
{
  //double A; // store some stuff
  // store grid
  Teuchos::RCP<Teuchos::Array <complex<double> > > grid= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  //for( int g=0;g<M_*M_*M_;g++)
  {
    //grid
  }
  double result = 0;
  for (int j=0;j<N_;j++)
    {
      for (int k=0;k<N_;k++)
      {
        for (int l=0;l<N_;l++)
        {
          // set first elements to zero
         // if(k==0||j==0||l==0||j>(N_-2)||k>(N_-2)||l>(N_-2))
          {

          }
          //else
          {
            //A=sqrt(2*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
            //A=sqrt(2*(discrete_PSD_[l+N_*(k+N_*j)]*(pow(dkappa_,3))));

            result += sqrt( 2*discrete_PSD_[l+N_*(k+N_*j)]* pow( dkappa_, 3 ) ) *
                                           (cos( (l) * dkappa_ * x  + (k) * dkappa_ * y + (j) * dkappa_ * z + Phi_0_[l+N_*(k+N_*j)])+
                                               cos( (l) * dkappa_ * x  + (k) * dkappa_ * y - (j) * dkappa_ * z + Phi_1_[l+N_*(k+N_*j)])+
                                               cos( (l) * dkappa_ * x  - (k) * dkappa_ * y - (j) * dkappa_ * z + Phi_2_[l+N_*(k+N_*j)])+
                                               cos( (l) * dkappa_ * x  - (k) * dkappa_ * y + (j) * dkappa_ * z + Phi_3_ [l+N_*(k+N_*j)] )
                                               );

        }
        }
      }
    }
  return sqrt( 2 ) * result;


}

  //
void GenRandomField::ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret)
{
  // root bounding Box
  vector<double> maxrbb;
//  maxrbb.push_back(10.0e-19);
//  maxrbb.push_back(10.0e-19);
//  maxrbb.push_back(10.0e-19);
  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);
 // maxrbb.push_back(10.0e-19);
  vector<double> minrbb;
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);


//  bb_max_.push_back(10.0e-19);
//  bb_max_.push_back(10.0e-19);
//  bb_max_.push_back(10.0e-19);
  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);

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

  if (myrank_ == 0)
  {
    cout << "min " << bb_min_[0] << " "<< bb_min_[1]  << " "<< bb_min_[2] << endl;
    cout << "max " << bb_max_[0] << " "<< bb_max_[1]  << " "<< bb_max_[2] << endl;
  }

}
double GenRandomField::EvalFieldAtLocation(vector<double> location, bool writetofile, bool output)
{
  int index_x;
  int index_y;
  int index_z;

  // Compute indices
  index_x=int(floor((location[0]-bb_min_[0])/dx_));
  index_y=int(floor((location[1]-bb_min_[1])/dx_));
  // HACK for 2D art_aorta_case SET z to y
   if (myrank_ == 0&& output && dim_==2)
   {
     cout << "hack in use" << endl;
   }
   if (dim_==2)
     index_y=int(floor((location[2]-bb_min_[2])/dx_));
  index_z=int(floor((location[2]-bb_min_[2])/dx_));
  // check index
  if (index_x>M_||index_y>M_||index_z>M_)
    dserror("Index out of bounds");
  if (writetofile && myrank_==0 )
  {
    ofstream File;
    File.open("RFatPoint.txt",ios::app);
    // use at() to get an error massage just in case
    File << setprecision (9) << values_[index_x+M_*index_y]<< endl;
    File.close();
  }
  if (dim_==2)
  {
    return values_[index_x+M_*index_y];
  }
  else
    return values_[index_x+M_*(index_y+M_*index_z)];

  cout << "end of writing field" << endl;


}
// Translate Gaussian to nonGaussian process based on Mircea Grigoriu's translation process
// theory
void GenRandomField::TranslateToNonGaussian()
{
  double dim = dim_;
  // check wether pdf is gaussian
  switch(marginal_pdf_)
  {
    case normal:
      //cout << RED_LIGHT << "WARNING: Target marginal PDF is gaussian so nothing to do here"<< END_COLOR << endl;
    break;

    case beta:
    {
      dserror("fix this function");
      normal_distribution<>  my_norm(0,sigma_0_);
      cout << "sigma_ul_g_cur_it_" << sigma_ul_g_cur_it_ << endl;
      // init a prefixed beta distribution
      //double a_beta = 4.0;
      //double b_beta = 2.0;
      beta_distribution<> my_beta(distribution_params_[0],distribution_params_[1]);
      //translate

      for(int i=0;i<(pow(M_,dim));i++)
      {
        //cout << "value[i] " << values_[i] << endl;
        values_[i]=quantile(my_beta,cdf(my_norm, values_[i]))*distribution_params_[3]+distribution_params_[2];;
      }
    }
    break;

    case lognormal:
    {
      //normal_distribution<>  my_norm2(0,sigma_ul_g_cur_it_);
      boost::accumulators::accumulator_set<double, stats<boost::accumulators::tag::variance> > acc;
      // loop over i < pow(M_,dim_) so we dont need seperate functions for each dim
      for(int i=0;i<(pow(M_,dim));i++)
      {
        acc(values_[i]);
      }
      cout << "mean(acc) "<< mean(acc) << endl;
      cout << "var "<< sqrt(variance(acc)) << endl;
      normal_distribution<>  my_norm2(0,sqrt(variance(acc)));
      //cout << "sigma_ul_g_cur_it_" << sigma_ul_g_cur_it_ << endl;
      // init params for logn distribution based on shields paper
      //double mu_bar=1.8;
      //double sigma_N = sqrt(log(1+1/(pow(1.8,2))));
      //double mu_N = log(mu_bar)-pow(sigma_N,2)/2;
      //lognormal_distribution<>  my_lognorm(mu_N,sigma_N);
      lognormal_distribution<>  my_lognorm(distribution_params_[0],distribution_params_[1]);
      cout << "distribution_params_[0]" << distribution_params_[0] << endl;
      cout << "distribution_params_[1]" << distribution_params_[1] << endl;
      // Try to calc sigma here
      //double sigma_logn=sqrt(my_lognorm.variance());
      //double sigma_logn= sqrt(variance(my_lognorm));
      //cout << "sigma_logn " << sigma_logn << endl;
      //dserror("Stop rigth here");
      //accumulator_set<int, stats<tag::variance> > acc2;
      for(int i=0;i<(pow(M_,dim));i++)
      {
        values_[i]=quantile(my_lognorm,cdf(my_norm2, values_[i]));
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
  //error to target psd init with 100 %
  double psd_error=100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  vector<double> PSD_ng_target;
  vector<double> PSD_ng(N_*N_,0.0);
  vector<double> PSD_ul_g(N_*N_,0.0);
  vector<double> rho(2*N_*2*N_,0.0);

  PSD_ng_target=discrete_PSD_;
  // do i need to set this zero??
  PSD_ng_target[0]=0.0;
  // calc sigma form discrete PSD

  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_autocorr= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr_ng= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2));

  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ul_g_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_PSD_ng_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ng_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >( N_*2*N_*2));

  for (int j=0;j<N_*2;j++)
  {
    for (int k=0;k<N_*2;k++)
    {
      // sort entries ro w major style
      if(j>(N_-1)||k>(N_-1))
      {
        real((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
        imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
      }
      else
      {
        real((*PSD_ul_g_complex)[k+N_*2*j])=discrete_PSD_[k+j*N_];
        imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
        if(k==0||j==0)
        {
          // comment later on
          //real((*PSD_ul_g_complex)[k+N_*2*j])=real((*PSD_ul_g_complex)[k+N_*2*j])*0.5;
          //imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
        }
        else if (k==0&&j==0)
        {
          //real((*PSD_ul_g_complex)[k+N_*2*j])=0;
          //real((*PSD_ul_g_complex)[k+N_*2*j])=real((*PSD_ul_g_complex)[k+N_*2*j])*0.25;
          //imag((*PSD_ul_g_complex)[k+N_*2*j])=0;
        }
      }
    }
  }
//  fftw_plan test1;
//  fftw_plan test2;
//  test1= fftw_plan_dft_2d(2*N_, 2*N_, (reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))),
//       (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))),
//                                   1, FFTW_ESTIMATE);
//
//   test2= fftw_plan_dft_2d(2*N_, 2*N_,(reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
//         (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))),
//                                     1, FFTW_ESTIMATE);
  // TWO DIM FFTS for spectral matching
  int rank = 1; /* not 2: we are computing 1d transforms int n[] = {1024}; 1d transforms of length 10 */
  int N_fftw = 2*N_;
  int howmany = 2*N_; // same here
  int idist = 2*N_;
  int odist = 2*N_;
  int istride =1;
  int ostride = 1; /* distance between two elements in the same row/column/rank */


  fftw_plan ifft_of_rows_of_psd;
  fftw_plan ifft_of_columns_of_psd;
  // ifft for autocorr
  fftw_plan ifft_of_rows_of_autocorr_ng;
  fftw_plan ifft_of_columns_of_autocorr_ng;

  ifft_of_rows_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))),
       NULL,
                   ostride,
                   odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);

  ifft_of_rows_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))),
            NULL,
                      ostride,
                      odist,
                     FFTW_FORWARD,
                     FFTW_ESTIMATE);


  istride =2*N_;
  ostride=2*N_;
  idist=1;
  odist=1;
  ifft_of_columns_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
  ifft_of_columns_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))),
         NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))),
                      NULL,
                      ostride,
                      odist,
                      FFTW_FORWARD,
                      FFTW_ESTIMATE);

  sigma_ul_g_cur_it_=sigma_0_;


  //while error > 0.5%
  int i =0;
  do
  {
    if(i!=0)// set new psd_ul_g
    {
      for (int j=0;j<N_*2;j++)
      {
        for (int k=0;k<N_*2;k++)
        {
          if(j>(N_-1)||k>(N_-1))
          {
            real((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
            imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
          }
          else
          {
            real((*PSD_ul_g_complex)[k+N_*2*j])=PSD_ul_g[k+j*N_];
            imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
            if(k==0||j==0)
            {
              //real((*PSD_ul_g_complex)[k+N_*2*j])=real((*PSD_ul_g_complex)[k+N_*2*j])*0.5;
            }
          }
        }
      }
      double sigma_ul_g_cur_it_helper_=0;
      for (int g=0; g<4*N_*N_;g++)
      {
        sigma_ul_g_cur_it_helper_+=real((*PSD_ul_g_complex)[g])*pow(dkappa_,2);
      }
      sigma_ul_g_cur_it_helper_=sigma_ul_g_cur_it_helper_-real((*PSD_ul_g_complex)[0])*pow(dkappa_,2)*0.0;
      sigma_ul_g_cur_it_=sqrt(4*sigma_ul_g_cur_it_helper_);
      cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << endl;
    }
    cout << "PSD( 132) "<< (*PSD_ul_g_complex)[280] << endl;
    //fftw_execute(test1);
    fftw_execute(ifft_of_rows_of_psd);
    fftw_execute(ifft_of_columns_of_psd);
    double scaling_fac= 2*pi_;
    // Check is this supposed to be kappa_cutof (probably yes)
     scaling_fac=dkappa_*N_;

    // loop over vectorlength
    for(int k=0;k<2*N_*2*N_;k++)
    {
      rho[k] = real((*autocorr)[k])*2*pow((scaling_fac),2)/(2*N_*2*N_*(pow(sigma_0_,2)));
        //lets go for +_ 3 * sigma here
      real((*autocorr_ng)[k])=Integrate(-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,rho[k]);
      // maybe intergration is different for beta
    }
    fftw_execute(ifft_of_rows_of_autocorr_ng);
    fftw_execute(ifft_of_columns_of_autocorr_ng);
    cout << "PSD ng ( 132) "<< real((*PSD_ng_complex)[280])/(pow(scaling_fac,2))*2 << endl;

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
         PSD_ng[k+j*N_]=real((*PSD_ng_complex)[k+j*2*N_])/(pow(scaling_fac,2));
         PSD_ul_g[k+j*N_]=real((*PSD_ul_g_complex)[k+j*2*N_]);
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


    //ofstream File2;
    //File2.open("Psd_coplex.txt",ios::app);
    //for (int j=0;j<2*N_;j++)
     // {
      //for (int k=0;k<2*N_;k++)
       // {
        //File2 << k << " " << j << " " << real(PSD_ng_complex[k+j*2*N_])/(4*pow(pi_,2)) << endl;
        //}
     //}
     //File2.close();
     // compute error based on equation(19) from shield2011
    error_numerator=0.0;
    error_denominator=0.0;
    for(int g=0;g<N_*N_;g++)
    {
      error_numerator+=pow((PSD_ng[g]-PSD_ng_target[g]),2);
      error_denominator+=pow((PSD_ng_target[g]),2);
    }
    psd_error=100*sqrt(error_numerator/error_denominator);
    if (myrank_ == 0)
      cout<< "Error to target PSD: " << psd_error << endl;
    // increase counter
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while(psd_error >0.5);



  // free memory
  fftw_destroy_plan(ifft_of_columns_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_columns_of_autocorr_ng);


  // Write PSD_ul_g PSD_ng_target and PSD_ng to a file
  // Dimension is 128* 128

  for(int h=0;h<N_*N_;h++)
  {
    if(PSD_ng[h]>10e-10)// change that
      discrete_PSD_[h]=PSD_ul_g[h];
    else
      //remove all the very small entries to get rid of the wiggles
      discrete_PSD_[h]=0.0;
  }
  // Write to file
   if (myrank_ == 0)
     {
       ofstream File;
       File.open("DiscretePSDTranslated.txt",ios::out);
       int size = int (pow(N_,2.0));
       for(int i=0;i<size;i++)
       {
         File << discrete_PSD_[i]<< endl;
       }
       File.close();
     }
  //ofstream File;
  //File.open("Psd.txt",ios::app);
 // for (int j=0;j<N_;j++)
 // {
 //   for (int k=0;k<N_;k++)
  //  {
  //    File << k << " " << j << " " << PSD_ng_target[k+j*N_] << " " << PSD_ng[k+j*N_]<< " "  << PSD_ul_g[k+j*N_] << endl;
  //  }
  //}
  //File.close();
  cout<< "Spectral Matching done "<< endl;
}
// Transform PSD of underlying gauusian process
void GenRandomField::SpectralMatching3D()
{
  //error to target psd init with 100 %
  double psd_error=100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  vector<double> PSD_ng_target;
  vector<double> PSD_ng(N_*N_*N_,0.0);
  vector<double> PSD_ul_g(N_*N_*N_,0.0);
  vector<double> rho(2*N_*2*N_*2*N_,0.0);

  PSD_ng_target=discrete_PSD_;
  // do i need to set this zero??
  PSD_ng_target[0]=0.0;
  // calc sigma form discrete PSD

  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr= Teuchos::rcp( new Teuchos::Array<complex<double> >               ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_autocorr= Teuchos::rcp( new Teuchos::Array<complex<double> >        ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_autocorr2= Teuchos::rcp( new Teuchos::Array<complex<double> >       ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr_ng= Teuchos::rcp( new Teuchos::Array<complex<double> >            ( N_*2*N_*2*N_*2,0.0));

  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ul_g_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >       ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_PSD_ng_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >  ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > almost_PSD_ng_complex2= Teuchos::rcp( new Teuchos::Array<complex<double> > ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ng_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >         ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > temp= Teuchos::rcp( new Teuchos::Array<complex<double> >         ( N_*2*N_*2*N_*2,0.0));

  for (int j=0;j<N_*2;j++)
  {
    for (int k=0;k<N_*2;k++)
    {
      for (int l=0;l<N_*2;l++)
      {
        // sort entries ro w major style
        if(j>(N_-1)||k>(N_-1)||l>(N_-1))
        {
          //l+M_*(k+M_*j)
          real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
          imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
        }
        else
        {
          real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=discrete_PSD_[l+(N_)*(k+(N_)*j)];
          imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
        }
      }
    }
  }

   int rank = 1; /* not 3: we are computing 1d transforms int n[] = {1024} 1d transforms of length 10 */
   int N_fftw = 2*N_;
   int howmany = 2*N_*2*N_; // same here
   int idist = 2*N_;// the distance in memory between the first element  of the first array and the first element of the second array */
   int odist = 2*N_;
   int istride =1;
   int ostride = 1; /* distance between two elements in the same row/column/rank */

  fftw_plan ifft_of_rows_of_psd;
  fftw_plan ifft_of_columns_of_psd;
  fftw_plan ifft_of_rank_of_psd;
  // ifft for autocorr
  fftw_plan ifft_of_rows_of_autocorr_ng;
  fftw_plan ifft_of_columns_of_autocorr_ng;
  fftw_plan ifft_of_rank_of_autocorr_ng;


  ifft_of_rows_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))),
         NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))),
         NULL,
                     ostride,
                     odist,
                    FFTW_BACKWARD,
                    FFTW_ESTIMATE);


    ifft_of_rows_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
        NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))),
              NULL,
                        ostride,
                        odist,
                        FFTW_BACKWARD,
                       FFTW_ESTIMATE);

    // now the same for collumns
    idist=2*N_;
    odist=2*N_;// the distance in memory between the first element  of the first array and the first element of the second array */
    istride = 1;/* distance between two elements in the same row/column/rank */
    ostride = 1; /* distance between two elements in the same row/column/rank */
    howmany = 2*N_*2*N_;// just 2*N_ we call this multiple times
    rank = 1; /* not 3: we are computing 1d transforms int n[] = {1024} 1d transforms of length 10 */


    ifft_of_columns_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*temp)[0]))),
         NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*almost_autocorr2)[0]))),
                     NULL,
                     ostride,
                     odist,
                     FFTW_BACKWARD,
                     FFTW_ESTIMATE);

    ifft_of_columns_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
           NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex2)[0]))),
                        NULL,
                        ostride,
                        odist,
                        FFTW_BACKWARD,
                        FFTW_ESTIMATE);

      // and now the rank (aka third dim of the array)
      howmany = 2*N_*2*N_;
      idist=1;// the distance in memory between the first element  of the first array and the first element of the second array */
      odist=1;
      istride =2*N_*2*N_;
      ostride=2*N_*2*N_; /* distance between two elements in the same row/column/rank */
      ifft_of_rank_of_psd = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))),
           NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))),
                       NULL,
                       ostride,
                       odist,
                       FFTW_BACKWARD,
                       FFTW_ESTIMATE);

      ifft_of_rank_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))),
             NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))),
                          NULL,
                          ostride,
                          odist,
                          FFTW_BACKWARD,
                          FFTW_ESTIMATE);

  sigma_ul_g_cur_it_=sigma_0_;
  //while error > 0.5%
  int i =0;
  do
  {
    if(i!=0)// set new psd_ul_g
    {
      for (int j=0;j<N_*2;j++)
      {
        for (int k=0;k<N_*2;k++)
        {
          for (int l=0;l<N_*2;l++)
          {
            if(j>(N_-1)||k>(N_-1)||l>(N_-1))
            {
              real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
              imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
            }
            else
            {
              real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=PSD_ul_g[l+(N_)*(k+(N_)*j)];
              imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
              if(k==0||j==0||l==0)
              {
                //real((*PSD_ul_g_complex)[k+N_*2*j])=real((*PSD_ul_g_complex)[k+N_*2*j])*0.25;
              }
            }
          }
        }
      }
    }
    cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << endl;


    cout << "PSD( 16644) "<< (*PSD_ul_g_complex)[10+(2*N_)*(10+(2*N_)*2)] << endl;
    cout << "PSD( 16644) "<< (*PSD_ul_g_complex)[16644] << endl;
    fftw_execute(ifft_of_rows_of_psd);

    for (int j=0;j<N_*2;j++)
       {
         for (int k=0;k<N_*2;k++)
         {
           for (int l=0;l<N_*2;l++)
           {
             real((*temp)[l+(2*N_)*(k+(2*N_)*j)])=real((*autocorr)[k+(2*N_)*(l+(2*N_)*j)]);
             imag((*temp)[l+(2*N_)*(k+(2*N_)*j)])=imag((*autocorr)[k+(2*N_)*(l+(2*N_)*j)]);
           }
         }
       }
  fftw_execute(ifft_of_columns_of_psd);
  // and transpose back
  // transpose
      for (int j=0;j<N_*2;j++)
         {
           for (int k=0;k<N_*2;k++)
           {
             for (int l=0;l<N_*2;l++)
             {
               real((*almost_autocorr)[l+(2*N_)*(k+(2*N_)*j)])=real((*almost_autocorr2)[k+(2*N_)*(l+(2*N_)*j)]);
               // if you forget the imaginary part it will mess with your scaling
               imag((*almost_autocorr)[l+(2*N_)*(k+(2*N_)*j)])=imag((*almost_autocorr2)[k+(2*N_)*(l+(2*N_)*j)]);
             }
           }
         }


    fftw_execute(ifft_of_rank_of_psd);
    double scaling_fac= dkappa_*N_;

    for(int k=0;k<2*N_*2*N_*2*N_;k++)
    {
      //scaling with sigma to get autocorrelation function
      // Factor 2 to here at the back is essential (altough not quite sure were it comes from)
      rho[k] = real((*autocorr)[k])*2*pow((scaling_fac),3)/(2*N_*2*N_*2*N_*(pow(sigma_0_,2)));
      real((*autocorr_ng)[k])=Integrate(-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,rho[k]);
      imag((*autocorr_ng)[k])=0.0;
      // The followign lines are good for testing if the FFT works correctly
      //rho[k] = real((*autocorr)[k])*pow((scaling_fac),3)/(2*N_*2*N_*2*N_);
      //if (k==16644)
      //  cout << "rho[]" << rho[16644]  << endl;
     // real((*autocorr_ng)[k])=rho[k] ;
    }
   fftw_execute(ifft_of_rows_of_autocorr_ng);
   // transpose
  for (int j=0;j<N_*2;j++)
     {
       for (int k=0;k<N_*2;k++)
       {
         for (int l=0;l<N_*2;l++)
         {
           real((*autocorr_ng)[l+(2*N_)*(k+(2*N_)*j)])=real((*PSD_ng_complex)[k+(2*N_)*(l+(2*N_)*j)]);
           imag((*autocorr_ng)[l+(2*N_)*(k+(2*N_)*j)])=imag((*PSD_ng_complex)[k+(2*N_)*(l+(2*N_)*j)]);
         }
       }
     }
   fftw_execute(ifft_of_columns_of_autocorr_ng);
   // transpose back
    for (int j=0;j<N_*2;j++)
       {
         for (int k=0;k<N_*2;k++)
         {
           for (int l=0;l<N_*2;l++)
           {
             real((*almost_PSD_ng_complex)[l+(2*N_)*(k+(2*N_)*j)])=real((*almost_PSD_ng_complex2)[k+(2*N_)*(l+(2*N_)*j)]);
             imag((*almost_PSD_ng_complex)[l+(2*N_)*(k+(2*N_)*j)])=imag((*almost_PSD_ng_complex2)[k+(2*N_)*(l+(2*N_)*j)]);
           }
         }
       }


     fftw_execute(ifft_of_rank_of_autocorr_ng);
     // cout << "PSD after fft ( 10 10 10) "<< real((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*2<< "imag "<< imag((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8.  << endl;
     // cout << "PSD after fft ( 16644) "<< real((*PSD_ng_complex)[16644])/pow((scaling_fac),3)*2 << endl;
    for (int j=0;j<N_*2;j++)
    {
      for (int k=0;k<N_*2;k++)
      {
        for (int l=0;l<N_*2;l++)
        {
          if(j>(N_-1)||k>(N_-1)||l>(N_-1))
          {}
          else
          {
            PSD_ng[l+(N_)*(k+(N_)*j)]=real((*PSD_ng_complex)[l+(2*N_)*(k+(2*N_)*j)])/(pow(scaling_fac,3));
            PSD_ul_g[l+(N_)*(k+(N_)*j)]=real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)]);
          }
        }
      }
    }

    PSD_ng[0]=0.0;
    for(int k=0;k<N_*N_*N_;k++)
    {
      if(PSD_ng[k]>10e-10)
      {
        PSD_ul_g[k]=pow(PSD_ng_target[k]/PSD_ng[k],1.4)*PSD_ul_g[k];
      cout << "psd_ul_g" << PSD_ul_g[k] <<  endl;
      }
      // do not set to zero because if once zero you'll never get it non-zero again
      else
      {
        //cout<<  " PSD ul g " << PSD_ul_g[k] << endl;
        PSD_ul_g[k]=10e-10;
      }
    }

     // compute error based on equation(19) from shield2011
    error_numerator=0.0;
    error_denominator=0.0;
    for(int g=0;g<N_*N_*N_;g++)
    {
      error_numerator+=pow((PSD_ng[g]-PSD_ng_target[g]),2);
      error_denominator+=pow((PSD_ng_target[g]),2);
    }
    psd_error=100*sqrt(error_numerator/error_denominator);
    if (myrank_ == 0)
      cout<< "Error to target PSD: " << psd_error << endl;
    // increase counter
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while(psd_error >0.5);
  // free memory
  fftw_destroy_plan(ifft_of_rows_of_psd);
  fftw_destroy_plan(ifft_of_columns_of_psd);
  fftw_destroy_plan(ifft_of_rank_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_columns_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_rank_of_autocorr_ng);
}

// Routine to calculate
// Transform PSD of underlying gauusian process using direct 3DFFT routine provided by FFTW
void GenRandomField::SpectralMatching3D3D()
{
  //error to target psd init with 100 %
  double psd_error=100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  vector<double> PSD_ng_target;
  vector<double> PSD_ng(N_*N_*N_,0.0);
  vector<double> PSD_ul_g(N_*N_*N_,0.0);
  vector<double> rho(2*N_*2*N_*2*N_,0.0);

  PSD_ng_target=discrete_PSD_;
  // do i need to set this zero??
  PSD_ng_target[0]=0.0;
  // calc sigma form discrete PSD

  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr= Teuchos::rcp( new Teuchos::Array<complex<double> >               ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > autocorr_ng= Teuchos::rcp( new Teuchos::Array<complex<double> >            ( N_*2*N_*2*N_*2,0.0));

  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ul_g_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >       ( N_*2*N_*2*N_*2,0.0));
  Teuchos::RCP<Teuchos::Array <complex<double> > > PSD_ng_complex= Teuchos::rcp( new Teuchos::Array<complex<double> >         ( N_*2*N_*2*N_*2,0.0));

  for (int j=0;j<N_*2;j++)
  {
    for (int k=0;k<N_*2;k++)
    {
      for (int l=0;l<N_*2;l++)
      {
        // sort entries ro w major style
        if(j>(N_-1)||k>(N_-1)||l>(N_-1))
        {
          //l+M_*(k+M_*j)
          real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
          imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
        }
        else
        {
          real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=discrete_PSD_[l+(N_)*(k+(N_)*j)];
          imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
        }
      }
    }
  }
  fftw_plan fft_of_psd;
  fftw_plan fft_of_autocorr_ng;
  //sign, can be either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1),
  fft_of_psd= fftw_plan_dft_3d(2*N_, 2*N_,2*N_, (reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))),
      (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))),
      1, FFTW_ESTIMATE);

  fft_of_autocorr_ng= fftw_plan_dft_3d(2*N_, 2*N_, 2*N_, (reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
           (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))),
                                       1, FFTW_ESTIMATE);
  sigma_ul_g_cur_it_=sigma_0_;
  //while error > 0.5%
  int i =0;
  do
  {
    if(i!=0)// set new psd_ul_g
    {
      for (int j=0;j<N_*2;j++)
      {
        for (int k=0;k<N_*2;k++)
        {
          for (int l=0;l<N_*2;l++)
          {
            if(j>(N_-1)||k>(N_-1)||l>(N_-1))
            {
              real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
              imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
            }
            else
            {
              real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=PSD_ul_g[l+(N_)*(k+(N_)*j)];
              imag((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)])=0.0;
              if(k==0||j==0||l==0)
              {
                real((*PSD_ul_g_complex)[k+N_*2*j])=real((*PSD_ul_g_complex)[k+N_*2*j])*0.25;
              }
            }
          }
        }
      }
    }

    fftw_execute(fft_of_psd);
    double scaling_fac= dkappa_*N_;
    // Check is this supposed to be kappa_cutof (probably yes)
    //scaling_fac=;

    for(int k=0;k<2*N_*2*N_*2*N_;k++)
    {
      // Factor 2 to here at the back is essential (altough not quite sure were it comes from)
      rho[k] = real((*autocorr)[k])*2*pow((scaling_fac),3)/(2*N_*2*N_*2*N_*(pow(sigma_0_,2)));
      real((*autocorr_ng)[k])=Integrate(-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,-3*sigma_ul_g_cur_it_,3*sigma_ul_g_cur_it_,rho[k]);
      // The followign lines are good for testing if the FFT works correctly
      //rho[k] = real((*autocorr)[k])*pow((scaling_fac),3)/(2*N_*2*N_*2*N_);
      //rho[k] = real((*autocorr)[k])*2*pow((scaling_fac),3)/(2*N_*2*N_*2*N_*(pow(sigma_0_,2)));
      // if (k==16644)
      //  cout << "rho[]" << rho[16644]  << endl;
      ///real((*autocorr_ng)[k])=rho[k] ;
    }
    fftw_execute(fft_of_autocorr_ng);
    // cout << "PSD after fft ( 10 10 10) "<< real((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8. << "imag "<< imag((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8.  << endl;
    //cout << "PSD after fft ( 16644) "<< real((*PSD_ng_complex)[16644])/pow((scaling_fac),3)*8 << endl;

    for (int j=0;j<N_*2;j++)
    {
      for (int k=0;k<N_*2;k++)
      {
        for (int l=0;l<N_*2;l++)
        {
          if(j>(N_-1)||k>(N_-1)||l>(N_-1))
          {}
          else
          {
            PSD_ng[l+(N_)*(k+(N_)*j)]=real((*PSD_ng_complex)[l+(2*N_)*(k+(2*N_)*j)])/(pow(scaling_fac,3));
            PSD_ul_g[l+(N_)*(k+(N_)*j)]=real((*PSD_ul_g_complex)[l+(2*N_)*(k+(2*N_)*j)]);
          }
        }
      }
    }

    PSD_ng[0]=0.0;

    for(int k=0;k<N_*N_*N_;k++)
    {
      if(PSD_ng[k]>10e-10)
        PSD_ul_g[k]=pow(PSD_ng_target[k]/PSD_ng[k],1.4)*PSD_ul_g[k];
      // do not set to zero because if once zero you'll never get it non-zero again
      else
      {
        PSD_ul_g[k]=10e-10;
      }
    }
    // compute error based on equation(19) from shield2011
    error_numerator=0.0;
    error_denominator=0.0;
    for(int g=0;g<N_*N_*N_;g++)
    {
      error_numerator+=pow((PSD_ng[g]-PSD_ng_target[g]),2);
      error_denominator+=pow((PSD_ng_target[g]),2);
    }
    psd_error=100*sqrt(error_numerator/error_denominator);
    if (myrank_ == 0)
      cout<< "Error to target PSD: " << psd_error << endl;
    // increase counter
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while(psd_error >0.5);
  // free memory
  fftw_destroy_plan(fft_of_psd);
  fftw_destroy_plan(fft_of_autocorr_ng);
}
double GenRandomField::Integrate(double xmin, double xmax, double ymin, double ymax, double rho)
{
  // get trillios gausspoints with high order
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
  double result = 0.0;
  //double mu_bar = 1.8;
  //double sigma_N = sqrt(log(1+1/(pow(1.8,2))));
  //double mu_N = log(mu_bar)-pow(sigma_N,2)/2;
  //beta_distribution<>
  //lognormal_distribution<>  my_lognorm(mu_N,sigma_N);
  normal_distribution<> my_normal(0,sigma_ul_g_cur_it_);
  //get parametres
  switch(marginal_pdf_)
  {
    case lognormal:
    {
      //lognormal_distribution<>  my_lognorm(mu_N,sigma_N);
      lognormal_distribution<>  my_lognorm(distribution_params_[0],distribution_params_[1]);

      result=quantile(my_lognorm,(cdf(my_normal,(argument_x))))*quantile(my_lognorm,(cdf(my_normal,argument_y)))*
      (1/(2*pi_*pow(sigma_ul_g_cur_it_,2)*sqrt(1-pow(rho,2))))*exp(-(pow((argument_x),2)+pow((argument_y),2)-2*rho*(argument_x)*(argument_y))
      /(2*pow(sigma_ul_g_cur_it_,2)*(1-pow(rho,2))));

    }
    break;
    case beta:
    {
      dserror("fix this function");
      beta_distribution<>  my_beta(distribution_params_[0],distribution_params_[1]);
      result=(quantile(my_beta,(cdf(my_normal,(argument_x))))*distribution_params_[3]+distribution_params_[2])*(quantile(my_beta,(cdf(my_normal,argument_y)))*distribution_params_[3]+distribution_params_[2])*
       (1/(2*pi_*pow(sigma_0_,2)*sqrt(1-pow(rho,2))))*exp(-(pow((argument_x),2)+pow((argument_y),2)-2*rho*(argument_x)*(argument_y))
       /(2*pow(sigma_0_,2)*(1-pow(rho,2))));
       // cout << " quantile " << quantile(my_beta,0.9) << endl;
        //cout << "distribution_params_[0],distribution_params_[1] " << distribution_params_[0] << " " << distribution_params_[1]<< endl;
    }
    break;
    default:
    dserror("Only Beta and Lognorm distribution supported so far fix your input file");
    break;
  }
  // calc function value

   //Testing do with to gauss distributions
  //result=quantile(my_normal,(cdf(my_normal,(argument_x))))*quantile(my_normal,(cdf(my_normal,argument_y)))*
    //    (1/(2*pi_*pow(sigma_0_,2)*sqrt(1-pow(rho,2))))*exp(-(pow((argument_x),2)+pow((argument_y),2)-2*rho*(argument_x)*(argument_y))
      //      /(2*pow(sigma_0_,2)*(1-pow(rho,2))));

  return result;
}

// Write Random Field to file
void GenRandomField::WriteRandomFieldToFile()
{
  if (myrank_ == 0)
  {
    ofstream File;
    File.open("RandomField.txt",ios::out);
   int size = int (pow(M_,double(dim_)));
    for(int i=0;i<size;i++)
    {
      File << values_[i]<< endl;
    }
    File.close();
  }
}
// Compute PSD from current sample
void GenRandomField::GetPSDFromSample(Teuchos::RCP<Teuchos::Array <double> > sample_psd)
{
  //check wether sample_psd has correct size
  if(sample_psd->length()!=M_*M_)
    dserror("Sizemismatch");

  Teuchos::RCP<Teuchos::Array <complex<double> > > b1=
        Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));

  // define complex i
  complex<double> i_comp (0,1);

  for (int j=0;j<M_*M_;j++)
  {
    real((*b1)[j])= values_[j];
    imag((*b1)[j])= 0.0;
  }
  // allocate output arrays
  Teuchos::RCP<Teuchos::Array <complex<double> > > d1=
        Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));
  Teuchos::RCP<Teuchos::Array <complex<double> > > d2=
         Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_));

  fftw_plan fft_of_rows;
  fftw_plan fft_of_collums;


  int rank = 1; /* not 2: we are computing 1d transforms */
  /* 1d transforms of length M_ */
  int N_fftw = M_;
  int howmany = M_; // same here
  int idist = M_;
  int odist = M_;
  int istride =1;
  int ostride = 1; /* distance between two elements in the same column */

  fft_of_rows = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
                  NULL,
                  ostride,
                  odist,
                  FFTW_FORWARD,
                  FFTW_ESTIMATE);

  istride =M_;
  ostride=M_;
  idist=1;
  odist=1;
  fft_of_collums = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(&((*d2)[0]))),
                   NULL,
                   ostride,
                   odist,
                   FFTW_FORWARD,
                   FFTW_ESTIMATE);

  fftw_execute(fft_of_rows);
  fftw_execute(fft_of_collums);

  // move values into class variable
  for(int i=0;i<M_*M_;i++)
  {
    (*sample_psd)[i]=1/(pow(dkappa_,2))*1/(M_*M_)*1/(M_*M_)*pow(abs((*d2)[i]),2);
  }

  fftw_destroy_plan(fft_of_rows);
  fftw_destroy_plan(fft_of_collums);
}
// Compute PSD from current sample
void GenRandomField::GetPSDFromSample3D(Teuchos::RCP<Teuchos::Array <double> > sample_psd)
{
  //check wether sample_psd has correct size
  if(sample_psd->length()!=M_*M_*M_)
    dserror("Sizemismatch");

  Teuchos::RCP<Teuchos::Array <complex<double> > > b1=Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  // allocate output arrays
  Teuchos::RCP<Teuchos::Array <complex<double> > > d1= Teuchos::rcp( new Teuchos::Array<complex<double> >( M_*M_*M_));
  // define complex i
  complex<double> i_comp (0,1);

  for (int j=0;j<M_*M_*M_;j++)
  {
    real((*b1)[j])= values_[j];
    imag((*b1)[j])= 0.0;
  }
  fftw_plan fft_of_rf;
  //sign, can be either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1),
    fft_of_rf= fftw_plan_dft_3d(M_, M_,M_, (reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
        (reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
        -1, FFTW_ESTIMATE);

  fftw_execute(fft_of_rf);
  // move values into class variable
  for(int i=0;i<M_*M_*M_;i++)
  {
    (*sample_psd)[i]=1/(pow(dkappa_,3))*1/(M_*M_*M_)*1/(M_*M_*M_)*pow(abs((*d1)[i]),2);
  }

  fftw_destroy_plan(fft_of_rf);
}
// Write Random Field to file
void GenRandomField::WriteSamplePSDToFile(Teuchos::RCP<Teuchos::Array <double> > sample_psd)
{
  if (myrank_ == 0)
  {
    ofstream File;
    File.open("SamplePSD.txt",ios::out);
   int size = int (pow(M_,double(dim_)));
    for(int i=0;i<size;i++)
    {
      File << (*sample_psd)[i]<< endl;
    }
    File.close();
  }
}

#endif
