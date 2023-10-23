#ifndef BACI_SCATRA_ELE_CALC_NO_PHYSICS_FWD_HPP
#define BACI_SCATRA_ELE_CALC_NO_PHYSICS_FWD_HPP

/*----------------------------------------------------------------------*/
/*! \file

\brief forward declarations for scatra_ele_calc_no_physics classes

\level 2


*/

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::DiscretizationType::nurbs27>;
#endif
