/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for parallel grid generator functionality

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <Epetra_MpiComm.h>
#include "globalproblem.H"
#include "discret.H"
#include "matpar_material.H"
#include "matpar_bundle.H"
#include "gridgenerator.H"
#include "io_pstream.H"

namespace
{
  void CreateMaterialInGlobalProblem()
  {
    const auto mat_stvenant = Teuchos::rcp(new MAT::PAR::Material(
        1, INPAR::MAT::MaterialType::m_stvenant, "MAT_Struct_StVenantKirchhoff"));

    mat_stvenant->Add("YOUNG", 1.0);
    mat_stvenant->Add("NUE", 0.1);
    mat_stvenant->Add("DENS", 2.0);
    mat_stvenant->Add("THEXPANS", 1.0);

    DRT::Problem::Instance()->Materials()->Insert(1, mat_stvenant);
  }

  class GridGeneratorTest : public ::testing::Test
  {
   public:
    GridGeneratorTest()
    {
      inputData_.bottom_corner_point_ = std::array<double, 3>{-1.0, -2.0, -3.0};
      inputData_.top_corner_point_ = std::array<double, 3>{2.5, 3.5, 4.5};
      inputData_.interval_ = std::array<int, 3>{5, 10, 15};
      inputData_.node_gid_of_first_new_node_ = 17;
    };

   protected:
    void SetUp() override
    {
      CreateMaterialInGlobalProblem();
      comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      IO::cout.setup(false, false, false, IO::standard, comm_, 0, 0, "dummyFilePrefix");
      testdis_ = Teuchos::rcp(new DRT::Discretization("dummy", comm_));
    }

    void TearDown() override { IO::cout.close(); }

   public:
    DRT::GRIDGENERATOR::RectangularCuboidInputs inputData_{};
    Teuchos::RCP<DRT::Discretization> testdis_;
    Teuchos::RCP<Epetra_Comm> comm_;
  };

  TEST_F(GridGeneratorTest, TestGridGeneratorWithHex27Elements)
  {
    inputData_.elementtype_ = "SOLIDH27";
    inputData_.distype_ = "HEX27";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear";

    DRT::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const double* nodePosition = lastNode->X();

    if (comm_->MyPID() == 0)
    {
      EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
      EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
      EXPECT_NEAR(nodePosition[2], -0.5, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 2541);
      EXPECT_EQ(testdis_->NumMyRowElements(), 250);
      EXPECT_EQ(testdis_->NumMyColNodes(), 3003);
      EXPECT_EQ(testdis_->NumMyColElements(), 300);
      EXPECT_EQ(lastNode->Id(), 2557);
    }
    else if (comm_->MyPID() == 1)
    {
      EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
      EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
      EXPECT_NEAR(nodePosition[2], 2.0, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 2310);
      EXPECT_EQ(testdis_->NumMyRowElements(), 250);
      EXPECT_EQ(testdis_->NumMyColNodes(), 3003);
      EXPECT_EQ(testdis_->NumMyColElements(), 300);
      EXPECT_EQ(lastNode->Id(), 4867);
    }
    else if (comm_->MyPID() == 2)
    {
      EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
      EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
      EXPECT_NEAR(nodePosition[2], 4.5, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 2310);
      EXPECT_EQ(testdis_->NumMyRowElements(), 250);
      EXPECT_EQ(testdis_->NumMyColNodes(), 2541);
      EXPECT_EQ(testdis_->NumMyColElements(), 250);
      EXPECT_EQ(lastNode->Id(), 7177);
    }
  }

  TEST_F(GridGeneratorTest, TestGridGeneratorWithWedge6Elements)
  {
    inputData_.elementtype_ = "SOLIDW6";
    inputData_.distype_ = "WEDGE6";
    inputData_.elearguments_ = "MAT 1 KINEM nonlinear";
    inputData_.autopartition_ = true;

    DRT::GRIDGENERATOR::CreateRectangularCuboidDiscretization(*testdis_, inputData_, true);

    testdis_->FillComplete(false, false, false);

    DRT::Node* lastNode = testdis_->lRowNode(testdis_->NumMyRowNodes() - 1);
    const double* nodePosition = lastNode->X();

    if (comm_->MyPID() == 0)
    {
      EXPECT_NEAR(nodePosition[0], -0.3, 1e-14);
      EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
      EXPECT_NEAR(nodePosition[2], 2.0, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 352);
      EXPECT_EQ(testdis_->NumMyRowElements(), 511);
      EXPECT_EQ(testdis_->NumMyColNodes(), 467);
      EXPECT_EQ(testdis_->NumMyColElements(), 596);
      EXPECT_EQ(lastNode->Id(), 4859);
    }
    else if (comm_->MyPID() == 1)
    {
      EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
      EXPECT_NEAR(nodePosition[1], 0.75, 1e-14);
      EXPECT_NEAR(nodePosition[2], 2.0, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 335);
      EXPECT_EQ(testdis_->NumMyRowElements(), 519);
      EXPECT_EQ(testdis_->NumMyColNodes(), 465);
      EXPECT_EQ(testdis_->NumMyColElements(), 590);
      EXPECT_EQ(lastNode->Id(), 4757);
    }
    else if (comm_->MyPID() == 2)
    {
      EXPECT_NEAR(nodePosition[0], 2.5, 1e-14);
      EXPECT_NEAR(nodePosition[1], 3.5, 1e-14);
      EXPECT_NEAR(nodePosition[2], 4.5, 1e-14);
      EXPECT_EQ(testdis_->NumMyRowNodes(), 369);
      EXPECT_EQ(testdis_->NumMyRowElements(), 470);
      EXPECT_EQ(testdis_->NumMyColNodes(), 456);
      EXPECT_EQ(testdis_->NumMyColElements(), 570);
      EXPECT_EQ(lastNode->Id(), 7177);
    }
  }

}  // namespace
