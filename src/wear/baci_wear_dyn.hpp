/*---------------------------------------------------------------------*/
/*! \file
\brief Control routine for structure with ale problems (wear).

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_WEAR_DYN_HPP
#define BACI_WEAR_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

/* entry point for the solution of structural ale wear problems */
void wear_dyn_drt(int restart);

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // WEAR_DYN_H
