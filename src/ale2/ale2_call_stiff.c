
/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_keku' which calculates the stiffness
matrix at one integration point for a 2d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates usual stiffness matrix in total lagrangian formulation 

<pre>                                                              mn 06/02
This routine calculates usual stiffness matrix in total lagrangian
formulation.

</pre>
\param **s   DOUBLE    (o)  element stiffness matrix 
\param bs  **DOUBLE    (i)  derivative operator
\param d   **DOUBLE    (i)  constitutive matrix
\param fac   DOUBLE    (i)  integration factor
\param nd    INT       (i)  total number degrees of freedom of element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; caled by: ale2_static_ke(), ale3_static_ke()

*----------------------------------------------------------------------*/
void ale2_keku(
    DOUBLE  **s, 
    DOUBLE  **bs, 
    DOUBLE  **d, 
    DOUBLE    fac, 
    INT       nd,
    INT       neps
    )
{
  INT            i, j, k, l, m;
  DOUBLE         dum;
  DOUBLE         db[24];


#ifdef DEBUG 
  dstrc_enter("ale2_keku");
#endif

  for (j=0; j<nd; j++)
  {
    for (k=0; k<neps; k++)
    {
      db[k] = 0.0;
      for (l=0; l<neps; l++)
      {
        db[k] = db[k] + d[k][l]*bs[l][j]*fac ;
      }
    }
    for (i=0; i<nd; i++)
    {
      dum = 0.0;
      for (m=0; m<neps; m++)
      {
        dum = dum + bs[m][i]*db[m] ;
      }
      s[i][j] = s[i][j] + dum ;
    }
  }

#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of ale2_keku */

#endif
/*! @} (documentation module close)*/
