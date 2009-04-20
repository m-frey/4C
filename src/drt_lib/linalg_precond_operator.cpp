/*!----------------------------------------------------------------------
\file  linalg_precond_operator.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "linalg_precond_operator.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::LinalgPrecondOperator(
  Teuchos::RCP<Epetra_Operator> precond,
  bool                          project) :
  project_(project),
  precond_(precond)
{
  return;
} // LINALG::LinalgPrecondOperator::LinalgPrecondOperator

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::LinalgPrecondOperator::~LinalgPrecondOperator()
{
  return;
} // LINALG::LinalgPrecondOperator::~LinalgPrecondOperator

/* --------------------------------------------------------------------
                    (Modified) ApplyInverse call
   -------------------------------------------------------------------- */
int LINALG::LinalgPrecondOperator::ApplyInverse(
  const Epetra_MultiVector &X, 
  Epetra_MultiVector       &Y
  ) const
{ 

  int ierr=0;

  // if necessary, project out matrix kernel
  if(project_)
  {
    // Apply the inverse preconditioner to get new basis vector for the
    // Krylov space
    ierr=precond_->ApplyInverse(X,Y);

    // check for vectors for matrix kernel and weighted basis mean vector
    if(c_ == Teuchos::null || w_ == Teuchos::null)
    {
      dserror("no c_ and w_ supplied");
    }

    // loop all solution vectors
    for(int sv=0;sv<Y.NumVectors();++sv)
    {
      // loop all basis vectors of kernel and orthogonalize against them
      for(int mm=0;mm<c_->NumVectors();++mm)
      {
        // loop all weight vectors 
        for(int rr=0;rr<w_->NumVectors();++rr)
        {
          /*
                   T
                  w * c
          */
          double wTc=0.0;

          ((*c_)(mm))->Dot(*((*w_)(rr)),&wTc);
          
          if(fabs(wTc)<1e-14)
          {
            dserror("weight vector must not be orthogonal to c");
          }

          /*
                   T
                  c * Y
          */
          double cTY=0.0;

          ((*c_)(mm))->Dot(*(Y(sv)),&cTY);

          /*
                                  T
                       T         c * Y
                      P Y = Y - ------- * w
                                  T
                                 w * c
          */
          (Y(sv))->Update(-cTY/wTc,*((*w_)(rr)),1.0);
          
        } // loop all weight vectors
      } // loop kernel basis vectors
    } // loop all solution vectors
  }
  else
  {
    ierr=precond_->ApplyInverse(X,Y);
  }

  return(ierr);
} // LINALG::LinalgPrecondOperator::ApplyInverse


#endif // CCADISCRET
