/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/
#include "baci_cut_meshintersection.hpp"
#include "baci_cut_options.hpp"
#include "baci_cut_side.hpp"
#include "baci_cut_tetmeshintersection.hpp"
#include "baci_cut_volumecell.hpp"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_hex8_twintri()
{
  CORE::GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    CORE::LINALG::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.5;
    tri3_xyze(1, 0) = 0.0;
    tri3_xyze(2, 0) = 1.0;
    tri3_xyze(0, 1) = 0.5;
    tri3_xyze(1, 1) = 1.0;
    tri3_xyze(2, 1) = 0.0;
    tri3_xyze(0, 2) = 0.25;
    tri3_xyze(1, 2) = 1.0;
    tri3_xyze(2, 2) = 1.0;
    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, CORE::FE::CellType::tri3);
  }
  {
    CORE::LINALG::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.5;
    tri3_xyze(1, 0) = 0.0;
    tri3_xyze(2, 0) = 1.0;
    tri3_xyze(0, 1) = 0.4;
    tri3_xyze(1, 1) = 0.0;
    tri3_xyze(2, 1) = 0.0;
    tri3_xyze(0, 2) = 0.5;
    tri3_xyze(1, 2) = 1.0;
    tri3_xyze(2, 2) = 0.0;
    nids.clear();
    nids.push_back(11);
    nids.push_back(14);
    nids.push_back(12);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, CORE::FE::CellType::tri3);
  }

  CORE::LINALG::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, CORE::FE::CellType::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8_twinQuad()
{
  CORE::GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.1;
    quad4_xyze(1, 0) = 0.02;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 1.0;
    quad4_xyze(1, 1) = 0.02;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 1.0;
    quad4_xyze(1, 2) = 0.02;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.1;
    quad4_xyze(1, 3) = 0.02;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }
  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.1;
    quad4_xyze(1, 0) = 0.02;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 0.1;
    quad4_xyze(1, 1) = 0.02;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 0.1;
    quad4_xyze(1, 2) = 1.0;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.1;
    quad4_xyze(1, 3) = 1.0;
    quad4_xyze(2, 3) = 0.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(14);
    nids.push_back(15);
    nids.push_back(16);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }

  CORE::LINALG::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, CORE::FE::CellType::hex8);

  intersection.Status();

  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8_chairCut()
{
  CORE::GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.01;
    quad4_xyze(1, 0) = 0.0;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 0.02;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 0.02;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.01;
    quad4_xyze(1, 3) = 0.0;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }
  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.02;
    quad4_xyze(1, 0) = 0.45;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 1.0;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 1.0;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.02;
    quad4_xyze(1, 3) = 0.45;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(12);
    nids.push_back(15);
    nids.push_back(16);
    nids.push_back(13);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }

  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.0;
    quad4_xyze(1, 0) = 0.55;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 0.0;
    quad4_xyze(1, 1) = 0.55;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 0.8;
    quad4_xyze(1, 2) = 0.55;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.8;
    quad4_xyze(1, 3) = 0.55;
    quad4_xyze(2, 3) = 0.0;

    nids.clear();
    nids.push_back(17);
    nids.push_back(18);
    nids.push_back(19);
    nids.push_back(20);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }

  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.95;
    quad4_xyze(1, 0) = 1.0;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = 0.8;
    quad4_xyze(1, 1) = 0.55;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 0.8;
    quad4_xyze(1, 2) = 0.55;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 0.95;
    quad4_xyze(1, 3) = 1.0;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(21);
    nids.push_back(20);
    nids.push_back(19);
    nids.push_back(22);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }

  CORE::LINALG::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, CORE::FE::CellType::hex8);

  intersection.Status();

  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8_VCut()
{
  CORE::GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.5;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.5;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = -0.5;
    quad4_xyze(1, 2) = 1.5;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = -0.5;
    quad4_xyze(1, 3) = 1.5;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }
  {
    CORE::LINALG::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.9;
    quad4_xyze(1, 0) = 1.5;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.9;
    quad4_xyze(1, 1) = 1.5;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.5;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.5;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(16);
    nids.push_back(15);
    nids.push_back(12);
    nids.push_back(11);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, CORE::FE::CellType::quad4);
  }

  CORE::LINALG::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, CORE::FE::CellType::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
