
#ifdef STKADAPTIVE

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../stk_lib/stk_discret.H"

#include "fluid_implicit.H"
#include "fluid_problem.H"

#include "../linalg/linalg_solver.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../stk_refine/stk_mesh.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::FLD::Problem::Problem()
  : meta_( Teuchos::rcp( new MetaMesh ) )
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FLD::Problem::Setup()
{
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis( genprob.numff, 0 );
  if ( not actdis->HaveDofs() )
  {
    actdis->FillComplete();
  }

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->FluidSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  dis_ = Teuchos::rcp( new STK::Discretization( actdis->Comm() ) );

  fluid_ = Teuchos::rcp( new STK::FLD::Fluid( *dis_, solver ) );

  // setup mesh part definitions

  dis_->MetaSetup( meta_, *actdis );

  // declare fields

  fluid_->declare_fields( meta_->MetaData() );

  // done with meta data

  meta_->Commit();

  // create uniform mesh

  mesh_ = Teuchos::rcp( new Mesh( *meta_, MPI_COMM_WORLD ) );

  // setup mesh

  dis_->MeshSetup( mesh_, *actdis );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FLD::Problem::Execute()
{
  Setup();

#if 0
  Teuchos::RCP<IO::DiscretizationWriter> output = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);
  output->NewStep(0,0.0);
  output->WriteElementData();
#endif

  dis_->AdaptMesh( std::vector<stk::mesh::EntityKey>(),
                   std::vector<stk::mesh::EntityKey>() );

  fluid_->Integrate();

  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis( genprob.numff, 0 );

  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(actdis->Comm());
}

#endif
