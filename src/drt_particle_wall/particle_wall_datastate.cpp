/*---------------------------------------------------------------------------*/
/*! \file
\brief wall data state container for particle wall handler
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_wall_datastate.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDataState::WallDataState(const Teuchos::ParameterList& params) : params_(params)
{
  // empty constructor
}

void PARTICLEWALL::WallDataState::Init(const Teuchos::RCP<DRT::Discretization> walldiscretization)
{
  // set wall discretization
  walldiscretization_ = walldiscretization;

  // get flags defining considered states of particle wall
  bool ismoving = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_MOVING");
  bool isloaded = DRT::INPUT::IntegralValue<int>(params_, "PARTICLE_WALL_LOADED");

  // set current dof row and column map
  curr_dof_row_map_ = Teuchos::rcp(new Epetra_Map(*walldiscretization_->DofRowMap()));

  // create states needed for moving walls
  if (ismoving)
  {
    disp_row_ = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_), true);
    disp_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap()), true);
    disp_row_last_transfer_ = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_), true);
    vel_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap()), true);
    acc_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap()), true);
  }

  // create states needed for loaded walls
  if (isloaded)
  {
    force_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap()), true);
  }
}

void PARTICLEWALL::WallDataState::Setup()
{
  // nothing to do
}

void PARTICLEWALL::WallDataState::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEWALL::WallDataState::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEWALL::WallDataState::CheckForCorrectMaps()
{
  if (disp_row_ != Teuchos::null)
    if (not disp_row_->Map().SameAs(*walldiscretization_->DofRowMap()))
      dserror("map of state 'disp_row_' corrupt!");

  if (disp_col_ != Teuchos::null)
    if (not disp_col_->Map().SameAs(*walldiscretization_->DofColMap()))
      dserror("map of state 'disp_col_' corrupt!");

  if (disp_row_last_transfer_ != Teuchos::null)
    if (not disp_row_last_transfer_->Map().SameAs(*walldiscretization_->DofRowMap()))
      dserror("map of state 'disp_row_last_transfer_' corrupt!");

  if (vel_col_ != Teuchos::null)
    if (not vel_col_->Map().SameAs(*walldiscretization_->DofColMap()))
      dserror("map of state 'vel_col_' corrupt!");

  if (acc_col_ != Teuchos::null)
    if (not acc_col_->Map().SameAs(*walldiscretization_->DofColMap()))
      dserror("map of state 'acc_col_' corrupt!");

  if (force_col_ != Teuchos::null)
    if (not force_col_->Map().SameAs(*walldiscretization_->DofColMap()))
      dserror("map of state 'force_col_' corrupt!");
}

void PARTICLEWALL::WallDataState::UpdateMapsOfStateVectors()
{
  if (disp_row_ != Teuchos::null and disp_col_ != Teuchos::null)
  {
    // export row map based displacement vector
    Teuchos::RCP<Epetra_Vector> temp = disp_row_;
    disp_row_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofRowMap(), true));
    LINALG::Export(*temp, *disp_row_);

    // update column map based displacement vector
    disp_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap(), true));
    LINALG::Export(*disp_row_, *disp_col_);

    // store displacements after last transfer
    disp_row_last_transfer_ = Teuchos::rcp(new Epetra_Vector(*disp_row_));
  }

  if (vel_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*vel_col_, *temp);
    // export old row map based vector to new column map based vector
    vel_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap(), true));
    LINALG::Export(*temp, *vel_col_);
  }

  if (acc_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*acc_col_, *temp);
    // export old row map based vector to new column map based vector
    acc_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap(), true));
    LINALG::Export(*temp, *acc_col_);
  }

  if (force_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*force_col_, *temp);
    // export old row map based vector to new column map based vector
    force_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization_->DofColMap(), true));
    LINALG::Export(*temp, *force_col_);
  }

  // set new dof row map
  curr_dof_row_map_ = Teuchos::rcp(new Epetra_Map(*walldiscretization_->DofRowMap()));
}
