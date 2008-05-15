/*!
\file dof_distribution_switcher.cpp

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

//#include <blitz/array.h>
#include "xfem.H"
#include "dof_management.H"
#include "dof_distribution_switcher.H"
//#include "xdofmapcreation.H"
//#include "enrichment_utils.H"
#include "dofkey.H"
//#include "../io/gmsh.H"
#include "../drt_lib/drt_discret.H"
//#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
//#include "../drt_lib/linalg_systemmatrix.H"


using namespace std;

static XFEM::Enrichment genAlternativeEnrichment(
        const int                    gnodeid,
        const XFEM::PHYSICS::Field   oldphysvar,
        const RCP<XFEM::DofManager>  dofman
        )
{
    std::set<XFEM::FieldEnr> fieldset = dofman->getNodeDofSet(gnodeid);
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenriter = fieldset.begin(); fieldenriter != fieldset.end(); ++fieldenriter)
    {
        const XFEM::PHYSICS::Field physvar = fieldenriter->getField();
        if (oldphysvar == physvar)
        {
            return fieldenriter->getEnrichment();
            break;
        }
        else
        {
            dserror("bug!");
            exit(1);
        }
    }
    dserror("bug!");
    exit(1);
    //return XFEM::Enrichment();
}

void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistribution(
        RCP<Epetra_Vector>&             vector
        ) const
{
    // create new vector with new number of dofs 
    RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);
    
    if (vector == null)
    {
        cout << "created new vector with all zeros" << endl;
    }
    else
    {
        const RCP<Epetra_Vector> oldVector = vector;
        const Epetra_BlockMap& oldmap = oldVector->Map();
//        cout << "olddofrowmap_" << endl;
//        cout << (olddofrowmap_) << endl;
//        cout << "newdofrowmap_" << endl;
//        cout << (newdofrowmap_) << endl;
        
        if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");
        
        
        for (DofPosMap::const_iterator newdof = newNodalDofDistrib_.begin();
                                       newdof != newNodalDofDistrib_.end();
                                       ++newdof)
        {
            const XFEM::DofKey<XFEM::onNode> newdofkey = newdof->first;
            const int newdofpos = newdof->second;
            
            DofPosMap::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
            if (olddof != oldNodalDofDistrib_.end())  // if dofkey has existed before, use old value
            {
                //const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
                const int olddofpos = olddof->second;
                //cout << "init to old value" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
            }
            else // if dofkey has not been existed before, initialize to zero
            {
                //cout << "init to zero" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
            }
        }

        for (DofPosMap::const_iterator olddof = oldNodalDofDistrib_.begin();
                                       olddof != oldNodalDofDistrib_.end();
                                       ++olddof)
        {
            const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
            const int olddofpos = olddof->second;
            const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();
            
            // try to find successor
            DofPosMap::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
            if (newdof == newNodalDofDistrib_.end())  // no successor
            {
                dserror("bug: the interface is not moving at the moment");
                const int gnodeid = olddofkey.getGid();
                const BlitzVec3 actpos(toBlitzArray(ih_->xfemdis()->gNode(gnodeid)->X()));
                const XFEM::FieldEnr oldfieldenr = olddofkey.getFieldEnr();
                const XFEM::Enrichment oldenr = oldfieldenr.getEnrichment();
                const double enrval = oldenr.EnrValue(actpos, ih_->cutterdis(), Enrichment::approachUnknown);
                
                // create alternative dofkey
                const XFEM::Enrichment altenr(genAlternativeEnrichment(gnodeid, oldphysvar, dofman_));
                
                // find dof position of alternative key
                const XFEM::FieldEnr altfieldenr(oldfieldenr.getField(), altenr);
                const XFEM::DofKey<XFEM::onNode> altdofkey(gnodeid, altfieldenr);
                const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;
                if (newdofpos < 0)
                    dserror("bug!");
                
                // add old value to already existing values
                (*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofpos];
            }
        }
    }
    // set vector to zero or initialized vector
    vector = newVector;
}






#endif  // #ifdef CCADISCRET
