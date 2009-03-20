/*!----------------------------------------------------------------------
\file turbulence_statistics_bfs.cpp

\brief calculate pressures, mean velocity values and fluctuations for
turbulent flow over a backward-facing step

<pre>
o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - 21 lines in x2-direction
  - lines along upper and lower wall

o loop nodes closest to / on lines

o values on lines are averaged in time over all steps between two
  outputs

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "turbulence_statistics_bfs.H"

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines

  o Allocate distributed vector for squares
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfs::TurbulenceStatisticsBfs(
  RefCountPtr<DRT::Discretization> actdis,
  ParameterList&                   params)
  :
  discret_(actdis),
  params_ (params)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap,true);

  toggleu_      = LINALG::CreateVector(*dofrowmap,true);
  togglev_      = LINALG::CreateVector(*dofrowmap,true);
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);

  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates for in x1- and x2-direction
  //----------------------------------------------------------------------
  x1coordinates_ = rcp(new vector<double> );
  x2coordinates_ = rcp(new vector<double> );

  // the criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1avcoords;
  set<double,LineSortCriterion> x2avcoords;

  // loop nodes and build sets of lines in x1- and x2-direction
  // accessible on this proc
  // For x1-direction: consider upper wall
  // and assume no change in discretization behind the step
  // For x2-direction: consider vertical line at x1=0
  // and assume no change in discretization behind the step
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

    if (node->X()[1]<0.0820002 && node->X()[1]>0.0819998)
      x1avcoords.insert(node->X()[0]);
    if (node->X()[0]<2e-9 && node->X()[0]>-2e-9)
      x2avcoords.insert(node->X()[1]);

    if (x3min_>node->X()[2]) x3min_=node->X()[2];
    if (x3max_<node->X()[2]) x3max_=node->X()[2];
  }

  // communicate x3mins and x3maxs
  double min;
  discret_->Comm().MinAll(&x3min_,&min,1);
  x3min_=min;

  double max;
  discret_->Comm().MaxAll(&x3max_,&max,1);
  x3max_=max;

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------
  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();

    vector<char> sblock;
    vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
           x1line!=x1avcoords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(sblock,*x1line);
      }
#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblock.size();

      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);

      rblock.clear();

      // receive from predecessor
      frompid=(myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblock,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#else
      // dummy communication
      rblock.clear();
      rblock=sblock;
#endif

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        vector<double> coordsvec;

        coordsvec.clear();

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x1avcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x2line=x2avcoords.begin();
           x2line!=x2avcoords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(sblock,*x2line);
      }
#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblock.size();

      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);

      rblock.clear();

      // receive from predecessor
      frompid=(myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblock,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#else
      // dummy communication
      rblock.clear();
      rblock=sblock;
#endif

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        vector<double> coordsvec;

        coordsvec.clear();

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2avcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    x1coordinates_ = rcp(new vector<double> );
    x2coordinates_ = rcp(new vector<double> );

    for(set<double,LineSortCriterion>::iterator coord1=x1avcoords.begin();
        coord1!=x1avcoords.end();
        ++coord1)
    {
      x1coordinates_->push_back(*coord1);
    }

    for(set<double,LineSortCriterion>::iterator coord2=x2avcoords.begin();
        coord2!=x2avcoords.end();
        ++coord2)
    {
      x2coordinates_->push_back(*coord2);
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  numx1coor_ = x1coordinates_->size();
  numx2coor_ = x2coordinates_->size();

  // x1-direction
  x1sump_ =  rcp(new Epetra_SerialDenseMatrix);
  x1sump_->Reshape(2,numx1coor_);

  x1sumtauw_ =  rcp(new Epetra_SerialDenseMatrix);
  x1sumtauw_->Reshape(2,numx1coor_);

  // the following vectors are only necessary for low-Mach-number flow
  x1sumrho_ =  rcp(new Epetra_SerialDenseMatrix);
  x1sumrho_->Reshape(2,numx1coor_);

  x1sumT_ =  rcp(new Epetra_SerialDenseMatrix);
  x1sumT_->Reshape(2,numx1coor_);

  // x2-direction
  // first-order moments
  x2sumu_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumu_->Reshape(numx1coor_,numx2coor_);

  x2sumv_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumv_->Reshape(numx1coor_,numx2coor_);

  x2sumw_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumw_->Reshape(numx1coor_,numx2coor_);

  x2sump_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sump_->Reshape(numx1coor_,numx2coor_);

  // second-order moments
  x2sumsqu_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqu_->Reshape(numx1coor_,numx2coor_);

  x2sumsqv_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqv_->Reshape(numx1coor_,numx2coor_);

  x2sumsqw_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqw_->Reshape(numx1coor_,numx2coor_);

  x2sumsqp_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqp_->Reshape(numx1coor_,numx2coor_);

  x2sumuv_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumuv_->Reshape(numx1coor_,numx2coor_);

  x2sumuw_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumuw_->Reshape(numx1coor_,numx2coor_);

  x2sumvw_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumvw_->Reshape(numx1coor_,numx2coor_);

  // mean and rms of subgrid viscosity
  x2sumsv_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsv_->Reshape(numx1coor_,numx2coor_);

  x2sumsqsv_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqsv_->Reshape(numx1coor_,numx2coor_);

  // the following vectors are only necessary for low-Mach-number flow
  // first-order moments
  x2sumrho_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumrho_->Reshape(numx1coor_,numx2coor_);

  x2sumT_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumT_->Reshape(numx1coor_,numx2coor_);

  // second-order moments
  x2sumsqrho_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqrho_->Reshape(numx1coor_,numx2coor_);

  x2sumsqT_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumsqT_->Reshape(numx1coor_,numx2coor_);

  x2sumrhou_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumrhou_->Reshape(numx1coor_,numx2coor_);

  x2sumrhouT_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumrhouT_->Reshape(numx1coor_,numx2coor_);

  x2sumrhov_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumrhov_->Reshape(numx1coor_,numx2coor_);

  x2sumrhovT_ =  rcp(new Epetra_SerialDenseMatrix);
  x2sumrhovT_->Reshape(numx1coor_,numx2coor_);

  // set number of samples to zero
  numsamp_ = 0;

  return;
}// TurbulenceStatisticsBfs::TurbulenceStatisticsBfs

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfs::~TurbulenceStatisticsBfs()
{
  return;
}// TurbulenceStatisticsBfs::~TurbulenceStatisticsBfs()


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::DoTimeSample(
Teuchos::RefCountPtr<Epetra_Vector> velnp,
Teuchos::RefCountPtr<Epetra_Vector> subgrvisc,
Epetra_Vector &                     force
)
{
  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // loop nodes in x1-direction
  //----------------------------------------------------------------------
  for (vector<double>::iterator x1line=x1coordinates_->begin();
       x1line!=x1coordinates_->end();
       ++x1line)
  {
    x1nodnum++;
    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (vector<double>::iterator x2line=x2coordinates_->begin();
         x2line!=x2coordinates_->end();
         ++x2line)
    {
      x2nodnum++;

      // skip non-existing area before step
      if (*x1line<-2e-9 and *x2line<-2e-9) continue;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes=0;

      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if ((node->X()[0]<(*x1line+2e-9) and node->X()[0]>(*x1line-2e-9)) and
            (node->X()[1]<(*x2line+2e-9) and node->X()[1]>(*x2line-2e-9)))
        {
          vector<int> dof = discret_->Dof(node);
          double      one = 1.0;

          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
          togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
          togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs=0;

      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure and subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_,&u);
        velnp->Dot(*togglev_,&v);
        velnp->Dot(*togglew_,&w);
        velnp->Dot(*togglep_,&p);

        double sv;
        subgrvisc->Dot(*toggleu_,&sv);

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        //----------------------------------------------------------------------
        double usm=u/countnodesonallprocs;
        double vsm=v/countnodesonallprocs;
        double wsm=w/countnodesonallprocs;
        double psm=p/countnodesonallprocs;
        double svsm=sv/countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum,x2nodnum)+=usm;
        (*x2sumv_)(x1nodnum,x2nodnum)+=vsm;
        (*x2sumw_)(x1nodnum,x2nodnum)+=wsm;
        (*x2sump_)(x1nodnum,x2nodnum)+=psm;
        (*x2sumsv_)(x1nodnum,x2nodnum)+=svsm;

        (*x2sumsqu_)(x1nodnum,x2nodnum)+=usm*usm;
        (*x2sumsqv_)(x1nodnum,x2nodnum)+=vsm*vsm;
        (*x2sumsqw_)(x1nodnum,x2nodnum)+=wsm*wsm;
        (*x2sumsqp_)(x1nodnum,x2nodnum)+=psm*psm;
        (*x2sumsqsv_)(x1nodnum,x2nodnum)+=svsm*svsm;

        (*x2sumuv_)(x1nodnum,x2nodnum)+=usm*vsm;
        (*x2sumuw_)(x1nodnum,x2nodnum)+=usm*wsm;
        (*x2sumvw_)(x1nodnum,x2nodnum)+=vsm*wsm;

        // values at lower and upper wall
        if (x2nodnum == 0)
        {
          double tauw=0.0;
          for(int rr=0;rr<(*toggleu_).MyLength();++rr)
          {
            tauw += force[rr]*(*toggleu_)[rr];
          }
          double tauwsm=tauw/countnodesonallprocs;

          (*x1sump_)(0,x1nodnum)+=psm;
          (*x1sumtauw_)(0,x1nodnum)+=tauwsm;
        }
        else if (x2nodnum == (numx2coor_-1))
        {
          double tauw=0.0;
          for(int rr=0;rr<(*toggleu_).MyLength();++rr)
          {
            tauw += force[rr]*(*toggleu_)[rr];
          }
          double tauwsm=tauw/countnodesonallprocs;

          (*x1sump_)(1,x1nodnum)+=psm;
          (*x1sumtauw_)(1,x1nodnum)+=tauwsm;
        }
      }
    }
  }

  return;
}// TurbulenceStatisticsBfs::DoTimeSample

//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::DoLomaTimeSample(
Teuchos::RefCountPtr<Epetra_Vector> velnp,
Teuchos::RefCountPtr<Epetra_Vector> vedenp,
Teuchos::RefCountPtr<Epetra_Vector> subgrvisc,
Epetra_Vector &                     force,
const double                        eosfac)
{
  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // loop nodes in x1-direction
  //----------------------------------------------------------------------
  for (vector<double>::iterator x1line=x1coordinates_->begin();
       x1line!=x1coordinates_->end();
       ++x1line)
  {
    x1nodnum++;
    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (vector<double>::iterator x2line=x2coordinates_->begin();
         x2line!=x2coordinates_->end();
         ++x2line)
    {
      x2nodnum++;

      // skip non-existing area before step
      if (*x1line<-2e-9 and *x2line<-2e-9) continue;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes=0;

      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if ((node->X()[0]<(*x1line+2e-9) and node->X()[0]>(*x1line-2e-9)) and
            (node->X()[1]<(*x2line+2e-9) and node->X()[1]>(*x2line-2e-9)))
        {
          vector<int> dof = discret_->Dof(node);
          double      one = 1.0;

          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
          togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
          togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs=0;

      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs-=1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure, density, temperature, and
        // subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_,&u);
        velnp->Dot(*togglev_,&v);
        velnp->Dot(*togglew_,&w);
        velnp->Dot(*togglep_,&p);

        double sv;
        subgrvisc->Dot(*toggleu_,&sv);

        double rho;
        vedenp->Dot(*togglep_,&rho);

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        //----------------------------------------------------------------------
        double usm=u/countnodesonallprocs;
        double vsm=v/countnodesonallprocs;
        double wsm=w/countnodesonallprocs;
        double psm=p/countnodesonallprocs;
        double svsm=sv/countnodesonallprocs;
        double rhosm=rho/countnodesonallprocs;
        // compute temperature: T = eosfac/rho
        double Tsm=eosfac/rhosm;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum,x2nodnum)+=usm;
        (*x2sumv_)(x1nodnum,x2nodnum)+=vsm;
        (*x2sumw_)(x1nodnum,x2nodnum)+=wsm;
        (*x2sump_)(x1nodnum,x2nodnum)+=psm;
        (*x2sumsv_)(x1nodnum,x2nodnum)+=svsm;

        (*x2sumT_)(x1nodnum,x2nodnum)+=Tsm;
        (*x2sumrho_)(x1nodnum,x2nodnum)+=rhosm;

        (*x2sumsqu_)(x1nodnum,x2nodnum)+=usm*usm;
        (*x2sumsqv_)(x1nodnum,x2nodnum)+=vsm*vsm;
        (*x2sumsqw_)(x1nodnum,x2nodnum)+=wsm*wsm;
        (*x2sumsqp_)(x1nodnum,x2nodnum)+=psm*psm;
        (*x2sumsqsv_)(x1nodnum,x2nodnum)+=svsm*svsm;

        (*x2sumsqT_)(x1nodnum,x2nodnum)+=Tsm*Tsm;
        (*x2sumsqrho_)(x1nodnum,x2nodnum)+=rhosm*rhosm;

        (*x2sumuv_)(x1nodnum,x2nodnum)+=usm*vsm;
        (*x2sumuw_)(x1nodnum,x2nodnum)+=usm*wsm;
        (*x2sumvw_)(x1nodnum,x2nodnum)+=vsm*wsm;

        (*x2sumrhou_)(x1nodnum,x2nodnum)+=rhosm*usm;
        (*x2sumrhouT_)(x1nodnum,x2nodnum)+=rhosm*usm*Tsm;
        (*x2sumrhov_)(x1nodnum,x2nodnum)+=rhosm*vsm;
        (*x2sumrhovT_)(x1nodnum,x2nodnum)+=rhosm*vsm*Tsm;

        // values at lower and upper wall
        if (x2nodnum == 0)
        {
          double tauw = 0.0;
          for(int rr=0;rr<(*toggleu_).MyLength();++rr)
          {
            tauw += force[rr]*(*toggleu_)[rr];
          }
          double tauwsm = tauw/countnodesonallprocs;

          (*x1sump_)(0,x1nodnum)+=psm;
          (*x1sumtauw_)(0,x1nodnum)+=tauwsm;

          (*x1sumrho_)(0,x1nodnum)+=rhosm;
          (*x1sumT_)(0,x1nodnum)+=Tsm;
        }
        else if (x2nodnum == (numx2coor_-1))
        {
          double tauw = 0.0;
          for(int rr=0;rr<(*toggleu_).MyLength();++rr)
          {
            tauw += force[rr]*(*toggleu_)[rr];
          }
          double tauwsm = tauw/countnodesonallprocs;

          (*x1sump_)(1,x1nodnum)+=psm;
          (*x1sumtauw_)(1,x1nodnum)+=tauwsm;

          (*x1sumrho_)(1,x1nodnum)+=rhosm;
          (*x1sumT_)(1,x1nodnum)+=Tsm;
        }
      }
    }
  }

  return;
}// TurbulenceStatisticsBfc::DoLomaTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent flow over a backward-facing step (first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";
    (*log) << scientific;

    (*log) << "\n";
    (*log) << "# lower wall behind step and complete upper wall\n";
    (*log) << "#     x1";
    (*log) << "        lw-pmean       lw-tauw       lw-utau      uw-pmean       uw-tauw       uw-utau\n";

    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      double lwx1p    = (*x1sump_)(0,i)/numsamp_;
      double lwx1tauw = (*x1sumtauw_)(0,i)/numsamp_;
      double lwx1utau = sqrt(lwx1tauw);

      double uwx1p    = (*x1sump_)(1,i)/numsamp_;
      double uwx1tauw = (*x1sumtauw_)(1,i)/numsamp_;
      double uwx1utau = sqrt(uwx1tauw);

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x1coordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << lwx1p;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1tauw;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1utau;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1p;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1tauw;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1utau;
      (*log) << "\n";
    }

    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      (*log) << "\n";
      (*log) << "# line in x2-direction at x1 = " << setw(11) << setprecision(4) << (*x1coordinates_)[i] << "\n";
      (*log) << "#     x2";
      (*log) << "           umean         vmean         wmean         pmean        svmean";
      (*log) << "         urms          vrms          wrms          prms         svrms";
      (*log) << "          u'v'          u'w'          v'w'\n";

      for(unsigned j=0; j<x2coordinates_->size(); ++j)
      {
        double x2u  = (*x2sumu_)(i,j)/numsamp_;
        double x2v  = (*x2sumv_)(i,j)/numsamp_;
        double x2w  = (*x2sumw_)(i,j)/numsamp_;
        double x2p  = (*x2sump_)(i,j)/numsamp_;
        double x2sv = (*x2sumsv_)(i,j)/numsamp_;

        double x2urms  = sqrt((*x2sumsqu_)(i,j)/numsamp_-x2u*x2u);
        double x2vrms  = sqrt((*x2sumsqv_)(i,j)/numsamp_-x2v*x2v);
        double x2wrms  = sqrt((*x2sumsqw_)(i,j)/numsamp_-x2w*x2w);
        double x2prms  = sqrt((*x2sumsqp_)(i,j)/numsamp_-x2p*x2p);
        double x2svrms = sqrt((*x2sumsqsv_)(i,j)/numsamp_-x2sv*x2sv);

        double x2uv   = (*x2sumuv_)(i,j)/numsamp_-x2u*x2v;
        double x2uw   = (*x2sumuw_)(i,j)/numsamp_-x2u*x2w;
        double x2vw   = (*x2sumvw_)(i,j)/numsamp_-x2v*x2w;

        (*log) <<  " "  << setw(11) << setprecision(4) << (*x2coordinates_)[j];
        (*log) << "   " << setw(11) << setprecision(4) << x2u;
        (*log) << "   " << setw(11) << setprecision(4) << x2v;
        (*log) << "   " << setw(11) << setprecision(4) << x2w;
        (*log) << "   " << setw(11) << setprecision(4) << x2p;
        (*log) << "   " << setw(11) << setprecision(4) << x2sv;
        (*log) << "   " << setw(11) << setprecision(4) << x2urms;
        (*log) << "   " << setw(11) << setprecision(4) << x2vrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2wrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2prms;
        (*log) << "   " << setw(11) << setprecision(4) << x2svrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2uv;
        (*log) << "   " << setw(11) << setprecision(4) << x2uw;
        (*log) << "   " << setw(11) << setprecision(4) << x2vw;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsBfs::DumpStatistics

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::DumpLomaStatistics(int          step,
                                                      const double eosfac)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".loma_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent flow over a backward-facing step at low Mach number (first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";
    (*log) << scientific;

    (*log) << "\n";
    (*log) << "# lower wall behind step and complete upper wall\n";
    (*log) << "#     x1";
    (*log) << "        lw-pmean       lw-tauw       lw-utau    lw-rhomean      lw-Tmean      uw-pmean       uw-tauw       uw-utau    uw-rhomean      uw-Tmean\n";

    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      double lwx1p    = (*x1sump_)(0,i)/numsamp_;
      double lwx1tauw = (*x1sumtauw_)(0,i)/numsamp_;
      double lwx1utau = sqrt(lwx1tauw);

      double lwx1rho  = (*x1sumrho_)(0,i)/numsamp_;
      double lwx1T    = (*x1sumT_)(0,i)/numsamp_;

      double uwx1p   = (*x1sump_)(1,i)/numsamp_;
      double uwx1tauw = (*x1sumtauw_)(1,i)/numsamp_;
      double uwx1utau = sqrt(uwx1tauw);

      double uwx1rho  = (*x1sumrho_)(1,i)/numsamp_;
      double uwx1T    = (*x1sumT_)(1,i)/numsamp_;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x1coordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << lwx1p;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1tauw;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1utau;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1rho;
      (*log) << "   " << setw(11) << setprecision(4) << lwx1T;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1p;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1tauw;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1utau;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1rho;
      (*log) << "   " << setw(11) << setprecision(4) << uwx1T;
      (*log) << "\n";
    }

    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      (*log) << "\n";
      (*log) << "# line in x2-direction at x1 = " << setw(11) << setprecision(4) << (*x1coordinates_)[i] << "\n";
      (*log) << "#     x2";
      (*log) << "           umean         vmean         wmean         pmean        svmean       rhomean         Tmean      rhoumean     rhouTmean      rhovmean     rhovTmean";
      (*log) << "         urms          vrms          wrms          prms         svrms        rhorms          Trms";
      (*log) << "          u'v'          u'w'          v'w'\n";

      for(unsigned j=0; j<x2coordinates_->size(); ++j)
      {
        double x2u  = (*x2sumu_)(i,j)/numsamp_;
        double x2v  = (*x2sumv_)(i,j)/numsamp_;
        double x2w  = (*x2sumw_)(i,j)/numsamp_;
        double x2p  = (*x2sump_)(i,j)/numsamp_;
        double x2sv = (*x2sumsv_)(i,j)/numsamp_;

        double x2rho   = (*x2sumrho_)(i,j)/numsamp_;
        double x2T     = (*x2sumT_)(i,j)/numsamp_;
        double x2rhou  = (*x2sumrhou_)(i,j)/numsamp_;
        double x2rhouT = (*x2sumrhouT_)(i,j)/numsamp_;
        double x2rhov  = (*x2sumrhov_)(i,j)/numsamp_;
        double x2rhovT = (*x2sumrhovT_)(i,j)/numsamp_;

        double x2urms  = sqrt((*x2sumsqu_)(i,j)/numsamp_-x2u*x2u);
        double x2vrms  = sqrt((*x2sumsqv_)(i,j)/numsamp_-x2v*x2v);
        double x2wrms  = sqrt((*x2sumsqw_)(i,j)/numsamp_-x2w*x2w);
        double x2prms  = sqrt((*x2sumsqp_)(i,j)/numsamp_-x2p*x2p);
        double x2svrms = sqrt((*x2sumsqsv_)(i,j)/numsamp_-x2sv*x2sv);

        double x2rhorms = sqrt((*x2sumsqrho_)(i,j)/numsamp_-x2rho*x2rho);
        double x2Trms   = sqrt((*x2sumsqT_)(i,j)/numsamp_-x2T*x2T);

        double x2uv   = (*x2sumuv_)(i,j)/numsamp_-x2u*x2v;
        double x2uw   = (*x2sumuw_)(i,j)/numsamp_-x2u*x2w;
        double x2vw   = (*x2sumvw_)(i,j)/numsamp_-x2v*x2w;

        double x2rhouppTpp = x2rhouT-eosfac*x2rhou/x2rho;
        double x2rhovppTpp = x2rhovT-eosfac*x2rhov/x2rho;

        (*log) <<  " "  << setw(11) << setprecision(4) << (*x2coordinates_)[j];
        (*log) << "   " << setw(11) << setprecision(4) << x2u;
        (*log) << "   " << setw(11) << setprecision(4) << x2v;
        (*log) << "   " << setw(11) << setprecision(4) << x2w;
        (*log) << "   " << setw(11) << setprecision(4) << x2p;
        (*log) << "   " << setw(11) << setprecision(4) << x2sv;
        (*log) << "   " << setw(11) << setprecision(4) << x2rho;
        (*log) << "   " << setw(11) << setprecision(4) << x2T;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhou;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhouT;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhov;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhovT;
        (*log) << "   " << setw(11) << setprecision(4) << x2urms;
        (*log) << "   " << setw(11) << setprecision(4) << x2vrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2wrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2prms;
        (*log) << "   " << setw(11) << setprecision(4) << x2svrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhorms;
        (*log) << "   " << setw(11) << setprecision(4) << x2Trms;
        (*log) << "   " << setw(11) << setprecision(4) << x2uv;
        (*log) << "   " << setw(11) << setprecision(4) << x2uw;
        (*log) << "   " << setw(11) << setprecision(4) << x2vw;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhouppTpp;
        (*log) << "   " << setw(11) << setprecision(4) << x2rhovppTpp;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsBfs::DumpLomaStatistics


#endif /* CCADISCRET       */
