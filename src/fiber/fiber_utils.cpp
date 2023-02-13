/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of utility functions for fiber interpolation

\level 3
*/
/*----------------------------------------------------------------------*/


#include "fiber_utils.H"
#include "fiber_node.H"
#include "fiber_nodal_fiber_holder.H"
#include <cmath>

template <DRT::Element::DiscretizationType distype>
void DRT::FIBER::UTILS::ProjectFibersToGaussPoints(DRT::Node** nodes,
    const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        1>>& shapefcts,
    FIBER::NodalFiberHolder& gpFiberHolder)
{
  // number of nodes per element
  constexpr int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  std::array<DRT::FIBER::FiberNode*, nen> fiberNodes;
  std::vector<std::array<std::array<double, 3>, nen>> fibers;
  std::map<DRT::FIBER::CoordinateSystemDirection, std::array<std::array<double, 3>, nen>>
      coordinateSystemDirections;
  std::map<DRT::FIBER::AngleType, std::array<double, nen>> angles;

  for (int inode = 0; inode < nen; ++inode)
  {
    fiberNodes[inode] = dynamic_cast<DRT::FIBER::FiberNode*>(nodes[inode]);

    if (fiberNodes[inode] == nullptr)
    {
      dserror("At least one node of the element does not provide fibers.");
    }

    for (const auto& pair : fiberNodes[inode]->CoordinateSystemDirections())
    {
      coordinateSystemDirections[pair.first][inode] = pair.second;
    }

    for (std::size_t fiberid = 0; fiberid < fiberNodes[inode]->Fibers().size(); ++fiberid)
    {
      if (inode == 0)
      {
        // The first node has to create the array of fibers for each node
        fibers.emplace_back(std::array<std::array<double, 3>, nen>());
      }
      fibers[fiberid][inode] = fiberNodes[inode]->Fibers()[fiberid];
    }

    for (const auto& pair : fiberNodes[inode]->Angles())
    {
      angles[pair.first][inode] = pair.second;
    }
  }


  // project fibers
  for (const auto& fiber : fibers)
  {
    std::vector<LINALG::Matrix<3, 1>> gpQuantity;
    ProjectQuantityWithShapeFunctions<distype, 3>(fiber, shapefcts, gpQuantity);

    // normalize fiber vectors
    for (LINALG::Matrix<3, 1>& quantity : gpQuantity)
    {
      quantity.Scale(1.0 / quantity.Norm2());
    }

    // add quantity to the fiber holder
    gpFiberHolder.AddFiber(gpQuantity);
  }


  // project coordinate system directions
  for (const auto& pair : coordinateSystemDirections)
  {
    const DRT::FIBER::CoordinateSystemDirection type = pair.first;
    std::vector<LINALG::Matrix<3, 1>> gpQuantity;
    ProjectQuantityWithShapeFunctions<distype, 3>(pair.second, shapefcts, gpQuantity);

    // normalize fiber vectors
    for (LINALG::Matrix<3, 1>& quantity : gpQuantity)
    {
      quantity.Scale(1.0 / quantity.Norm2());
    }

    // add quantity to the fiber holder
    gpFiberHolder.SetCoordinateSystemDirection(type, gpQuantity);
  }

  // project angles
  for (const auto& pair : angles)
  {
    const DRT::FIBER::AngleType type = pair.first;
    std::vector<double> gpAngle;
    ProjectQuantityWithShapeFunctions<distype>(pair.second, shapefcts, gpAngle);

    // add quantity to the fiber holder
    gpFiberHolder.SetAngle(type, gpAngle);
  }


  // orthogonalize coordinate system
  if (gpFiberHolder.ContainsCoordinateSystemDirection(CoordinateSystemDirection::Circular) &&
      gpFiberHolder.ContainsCoordinateSystemDirection(CoordinateSystemDirection::Tangential))
  {
    const std::vector<LINALG::Matrix<3, 1>>& cir =
        gpFiberHolder.GetCoordinateSystemDirection(CoordinateSystemDirection::Circular);
    std::vector<LINALG::Matrix<3, 1>>& tan =
        gpFiberHolder.GetCoordinateSystemDirectionMutual(CoordinateSystemDirection::Tangential);

    // orthogonalize tangential vector, preserve circular direction
    for (std::size_t gp = 0; gp < tan.size(); ++gp)
    {
      double tancir = tan[gp].Dot(cir[gp]);
      tan[gp].Update(-tancir, cir[gp], 1.0);
      tan[gp].Scale(1.0 / tan[gp].Norm2());
    }

    // orthogonalize radial vector, preserve circular and tangential direction
    if (gpFiberHolder.ContainsCoordinateSystemDirection(CoordinateSystemDirection::Radial))
    {
      std::vector<LINALG::Matrix<3, 1>>& rad =
          gpFiberHolder.GetCoordinateSystemDirectionMutual(CoordinateSystemDirection::Radial);
      for (std::size_t gp = 0; gp < tan.size(); ++gp)
      {
        double radcir = rad[gp].Dot(cir[gp]);
        double radtan = rad[gp].Dot(tan[gp]);
        // double
        rad[gp].Update(-radcir, cir[gp], -radtan, tan[gp], 1.0);
        rad[gp].Scale(1.0 / rad[gp].Norm2());
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype, std::size_t dim>
void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions(
    const std::array<std::array<double, dim>,
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        1>>& shapefcts,
    std::vector<LINALG::Matrix<dim, 1>>& quantityProjected)
{
  quantityProjected.resize(shapefcts.size());

  // number of nodes per element
  constexpr std::size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  constexpr std::size_t nsd = DRT::UTILS::DisTypeToDim<distype>::dim;

  // number of gauss points
  const std::size_t ngp = shapefcts.size();

  for (LINALG::Matrix<dim, 1>& quantity : quantityProjected)
  {
    quantity.Clear();
  }

  for (std::size_t i = 0; i < nsd; ++i)
  {
    for (std::size_t q = 0; q < ngp; ++q)
    {
      for (std::size_t j = 0; j < nen; ++j)
      {
        quantityProjected[q](i) += shapefcts[q](j) * quantity[j][i];
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions(
    const std::array<double, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        1>>& shapefcts,
    std::vector<double>& quantityProjected)
{
  quantityProjected.resize(shapefcts.size());

  // number of nodes per element
  constexpr int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // number of gauss points
  const int ngp = shapefcts.size();

  for (double& quantity : quantityProjected)
  {
    quantity = 0.0;
  }

  for (int q = 0; q < ngp; ++q)
  {
    for (int j = 0; j < nen; ++j)
    {
      quantityProjected[q] += shapefcts[q](j) * quantity[j];
    }
  }
}

template <std::size_t dim>
void DRT::FIBER::UTILS::SetupCardiacFibers(
    const NodalFiberHolder& fibers, std::vector<LINALG::Matrix<dim, 1>>& f)
{
  if (fibers.FibersSize() > 0)
  {
    const std::vector<LINALG::Matrix<3, 1>>& fib = fibers.GetFiber(0);
    f.resize(fib.size());
    for (std::size_t gp = 0; gp < fib.size(); ++gp)
    {
      for (std::size_t i = 0; i < dim; ++i)
      {
        f[gp](i) = fib[gp](i);
      }
    }
  }
  else if (fibers.ContainsCoordinateSystemDirection(CoordinateSystemDirection::Circular) &&
           fibers.ContainsCoordinateSystemDirection(CoordinateSystemDirection::Tangential))
  {
    const std::vector<LINALG::Matrix<3, 1>>& cir =
        fibers.GetCoordinateSystemDirection(CoordinateSystemDirection::Circular);
    const std::vector<LINALG::Matrix<3, 1>>& tan =
        fibers.GetCoordinateSystemDirection(CoordinateSystemDirection::Tangential);
    const std::vector<double>& helix = fibers.GetAngle(AngleType::Helix);
    const std::vector<double>& transverse = fibers.GetAngle(AngleType::Transverse);
    f.resize(cir.size());

    double deg2rad = M_PI / 180.;
    for (unsigned int gp = 0; gp < cir.size(); ++gp)
    {
      LINALG::Matrix<3, 1> rad(false);
      rad.CrossProduct(cir[gp], tan[gp]);

      double tmp1 = cos(helix[gp] * deg2rad) * cos(transverse[gp] * deg2rad);
      double tmp2 = sin(helix[gp] * deg2rad) * cos(transverse[gp] * deg2rad);
      double tmp3 = sin(transverse[gp] * deg2rad);

      for (unsigned int i = 0; i < 3; ++i)
      {
        f[gp](i) = tmp1 * cir[gp](i, 0) + tmp2 * tan[gp](i, 0) + tmp3 * rad(i, 0);
      }
      f[gp].Scale(1.0 / f[gp].Norm2());
    }
  }
  else
  {
    dserror("You have to specify either FIBER1 or CIR, TAN, HELIX and TRANS");
  }
}

template <DRT::Element::DiscretizationType distype>
bool DRT::FIBER::UTILS::HaveNodalFibers(DRT::Node** nodes)
{
  constexpr int numberOfNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // check whether node can be casted to FiberNode
  auto nodeHasFibers = [](const Node* n)
  { return dynamic_cast<const DRT::FIBER::FiberNode*>(n) != nullptr; };
  return std::all_of(nodes, &nodes[numberOfNodes], nodeHasFibers);
}

template void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions<DRT::Element::tet4>(
    const std::array<double,
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement, 1>>& shapefcts,
    std::vector<double>& quantityProjected);
template void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions<DRT::Element::tet10>(
    const std::array<double,
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement, 1>>& shapefcts,
    std::vector<double>& quantityProjected);
template void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions<DRT::Element::hex8>(
    const std::array<double,
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement, 1>>& shapefcts,
    std::vector<double>& quantityProjected);
template void DRT::FIBER::UTILS::ProjectQuantityWithShapeFunctions<DRT::Element::tri3>(
    const std::array<double,
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement>
        quantity,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>& shapefcts,
    std::vector<double>& quantityProjected);

template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::tet4>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::tet10>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::hex8>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::tri3>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::quad4>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::quad9>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad9>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::hex27>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::pyramid5>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::pyramid5>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);
template void DRT::FIBER::UTILS::ProjectFibersToGaussPoints<DRT::Element::nurbs9>(DRT::Node**,
    const std::vector<LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::nurbs9>::numNodePerElement, 1>>&,
    FIBER::NodalFiberHolder&);

template bool DRT::FIBER::UTILS::HaveNodalFibers<DRT::Element::tet4>(DRT::Node** nodes);
template bool DRT::FIBER::UTILS::HaveNodalFibers<DRT::Element::tet10>(DRT::Node** nodes);
template bool DRT::FIBER::UTILS::HaveNodalFibers<DRT::Element::hex8>(DRT::Node** nodes);

template void DRT::FIBER::UTILS::SetupCardiacFibers<3>(
    const NodalFiberHolder& fibers, std::vector<LINALG::Matrix<3, 1>>& f);
template void DRT::FIBER::UTILS::SetupCardiacFibers<2>(
    const NodalFiberHolder& fibers, std::vector<LINALG::Matrix<2, 1>>& f);
