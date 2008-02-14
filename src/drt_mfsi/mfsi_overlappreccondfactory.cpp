#ifdef CCADISCRET

#include "mfsi_overlappreccondfactory.H"
#include "mfsi_overlappreccondoperator.H"

#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_InverseLinearOperator.hpp>

#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::OverlappingPCFactory::OverlappingPCFactory(
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  structure,
  //Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  interface,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  fluid,
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  ale
  )
  : structure_(structure),
    //interface_(interface),
    fluid_(fluid),
    ale_(ale)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MFSI::OverlappingPCFactory::isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc ) const
{
  dserror("MFSI::PreconditionerFactory::isCompatible() not implemented");
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Thyra::PreconditionerBase<double> > MFSI::OverlappingPCFactory::createPrec() const
{
  // construct the default preconditioner that will be initialized with an
  // appropiate operator
  return rcp(new Thyra::DefaultPreconditioner<double>());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlappingPCFactory::initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
                                                Thyra::PreconditionerBase<double> *precOp,
                                                const Thyra::ESupportSolveUse supportSolveUse) const
{
  Teuchos::RCP< const Thyra::LinearOpBase<double,double> > fsiOp = fwdOpSrc->getOp();
  Teuchos::RCP< const Thyra::DefaultBlockedLinearOp<double> > blockFsiOp =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultBlockedLinearOp<double> >(fsiOp);

#if 0

  // very simple block Jacobian preconditioner

  Thyra::ConstLinearOperator<double> structInnerOp = blockFsiOp->getBlock(0,0);
  Thyra::ConstLinearOperator<double> fluidInnerOp  = blockFsiOp->getBlock(1,1);
  Thyra::ConstLinearOperator<double> aleInnerOp    = blockFsiOp->getBlock(2,2);

  // Build inverse operators

  Thyra::ConstLinearOperator<double> invstruct = inverse(*structure_,structInnerOp,Thyra::IGNORE_SOLVE_FAILURE);
  Thyra::ConstLinearOperator<double> invfluid  = inverse(*fluid_,    fluidInnerOp, Thyra::IGNORE_SOLVE_FAILURE);
  Thyra::ConstLinearOperator<double> invale    = inverse(*ale_,      aleInnerOp,   Thyra::IGNORE_SOLVE_FAILURE);

  // set special MFSI block preconditioner object to Thyra Wrapper

  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > M =
    Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());
  M->beginBlockFill(3,3);

  M->setBlock(0,0,invstruct.constPtr());
  M->setBlock(1,1,invfluid.constPtr());
  M->setBlock(2,2,invale.constPtr());

  M->endBlockFill();

#else

  // constuct the preconditioner operator

  Teuchos::RCP<OverlappingPCOperator> M =
    Teuchos::rcp(new OverlappingPCOperator(structure_,fluid_,ale_,blockFsiOp));

#endif

  // stupid type conversations
  Teuchos::RCP<Thyra::LinearOpBase<double> > myM = M;
  Thyra::LinearOperator<double> constM = myM;

  Thyra::DefaultPreconditioner<double>* defaultPrec =
    &Teuchos::dyn_cast<Thyra::DefaultPreconditioner<double> >(*precOp);

  (*defaultPrec).initializeRight(constM.constPtr());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlappingPCFactory::uninitializePrec(Thyra::PreconditionerBase<double> *prec,
                                                  Teuchos::RCP<const Thyra::LinearOpSourceBase<double> >  *fwdOpSrc,
                                                  Thyra::ESupportSolveUse *supportSolveUse) const
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::OverlappingPCFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>&)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> MFSI::OverlappingPCFactory::getParameterList()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> MFSI::OverlappingPCFactory::unsetParameterList()
{
  return Teuchos::null;
}


#endif
