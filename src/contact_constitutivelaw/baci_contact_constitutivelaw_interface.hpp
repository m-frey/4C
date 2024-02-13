/*---------------------------------------------------------------------*/
/*! \file

\brief Contact interface for contact using micro information like the roughness profile

\level 1


*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_CONSTITUTIVELAW_INTERFACE_HPP
#define BACI_CONTACT_CONSTITUTIVELAW_INTERFACE_HPP

#include "baci_config.hpp"

#include "baci_contact_interface.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace MORTAR
{
  class InterfacedataContainer;
}

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class ConstitutiveLaw;
  }

  class ConstitutivelawInterface : public Interface
  {
   public:
    /*!
    \brief Standard constructor creating empty contact interface

    This initializes the employed shape function set for lagrangian mutlipliers
    to a specific setting. Throughout the evaluation process, this set will be employed
    for the field of lagrangian multipliers.

    \param idata_ptr (in): data container
    \param id (in): Unique interface id
    \param comm (in): A communicator object
    \param dim (in): Global problem dimension
    \param icontact (in): Global contact parameter list
    \param selfcontact (in): Flag for self contact status
    */
    ConstitutivelawInterface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData,
        const int id, const Epetra_Comm& comm, const int dim,
        const Teuchos::ParameterList& icontact, bool selfcontact,
        const int contactconstitutivelawid);

    /**
     * \brief Assemble gap-computed lagrange multipliers and nodal linlambda derivatives into nodal
     * quantities using the Macauley bracket
     *
     * When dealing with penalty methods, the lagrange multipliers are not independent variables
     * anymore. Instead, they can be computed in terms of the weighted gap and the penalty
     * parameter. This is done here so every node stores the correct lm and thus we integrate
     * smoothly into the overlaying algorithm.
     *
     * Additionally, we use the performed loop over all nodes to store the nodal derivlambda_j
     * matrix right there.
     *
     * As a result, the function notifies the calling routine if any negative gap was detected
     * and thus wether the interface is in contact or not. In consequence, after calling this
     * routine from within the penalty strategy object, the contact status is known at a global
     * level.
     *
     * Note: To be able to perform this computation, weighted gaps and normals have to be available
     * within every node! Since this computation is done via Interface::Evaluate() in the Integrator
     * class, these corresponding methods have to be called before AssembleMacauley()!
     *
     * \params[in/out] localisincontact true if at least one node is in contact
     * \params[in/out] localactivesetchange true if the active set changed
     *
     */
    void AssembleRegNormalForces(bool& localisincontact, bool& localactivesetchange) override;

    /**
     * \brief Throws an error since frictional contact is not yet tested
     */
    void AssembleRegTangentForcesPenalty() override;

    /** \brief return the multi-scale constitutive law used for the contact containing information
     i.e. on the micro roughness
     */
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> GetConstitutiveContactLaw()
    {
      return coconstlaw_;
    }

   private:
    /** \brief multi-scale constitutive law used for the contact containing information
     i.e. on the micro roughness
     */
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> coconstlaw_;
  };  // class Interface
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_CONSTITUTIVELAW_INTERFACE_H
