/*!
\file dof_management_element.cpp

\brief provides the element dofmanager class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "dof_management_element.H"

#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager() :
  numElemDof_(0),
  nodalDofSet_()
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
    const DRT::Element& ele,
    const map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
) :
  nodalDofSet_(nodalDofSet),
  DisTypePerElementField_(element_ansatz)
{
  ComputeDependendInfo(ele, nodalDofSet, enrfieldset, element_ansatz);
 
  return;
}
    
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::ElementDofManager::ComputeDependendInfo(
    const DRT::Element& ele,
    const map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
)
{
  // count number of dofs for each node
  map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp;
  for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp)
  {
    const int gid = tmp->first;
    nodalNumDof_[gid] = tmp->second.size();
  }
  
  // set number of parameters per field to zero
  for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp)
  {
    const std::set<XFEM::FieldEnr> enrfieldset = tmp->second;
    for (set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      numParamsPerField_[field] = 0;
      paramsLocalEntries_[field] = vector<int>();
    }
  }
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    numParamsPerField_[field] = 0;
    paramsLocalEntries_[field] = vector<int>();
  }
      
      
  unique_enrichments_.clear();
  // count number of parameters per field
  // define local position of unknown by looping first over nodes and then over its unknowns!
  int dofcounter = 0;
  const int* nodeids = ele.NodeIds();
  for (int inode=0; inode<ele.NumNode(); ++inode)
  {
    const int gid = nodeids[inode];
    map<int, const set <XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
    if (entry == nodalDofSet_.end())
      dserror("impossible ;-)");
    const std::set<XFEM::FieldEnr> enrfieldset = entry->second;
    
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      numParamsPerField_[field] += 1;
      paramsLocalEntries_[field].push_back(dofcounter);
      unique_enrichments_.insert(enrfield->getEnrichment());
      dofcounter++;
    }
  }
  // loop now over element dofs
  // we first loop over the fields and then over the params
  numElemDof_ = 0;
  enrichedFieldperPhysField_.clear();
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    //cout << physVarToString(field) << endl;
    std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator schnack = element_ansatz.find(field);
    if (schnack == element_ansatz.end())
    {
      cout << XFEM::PHYSICS::physVarToString(field) << endl;
      dserror("field not found -> bug");
    }
    
    enrichedFieldperPhysField_[field].insert(*enrfield);
    
    const DRT::Element::DiscretizationType eledofdistype = schnack->second;
    const int numparam = DRT::UTILS::getNumberOfElementNodes(eledofdistype);
    for (int inode=0; inode<numparam; ++inode)
    {
      numElemDof_ +=1;
      numParamsPerField_[field] += 1;
      paramsLocalEntries_[field].push_back(dofcounter);
      unique_enrichments_.insert(enrfield->getEnrichment());
      dofcounter++;
    }
  }
  
  return;
}
    
/*----------------------------------------------------------------------*
 |  construct element dof manager                               ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
    const DRT::Element&  ele,
    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    const XFEM::DofManager& dofman
) :
  DisTypePerElementField_(element_ansatz)
{
  // nodal dofs for ele
  nodalDofSet_.clear();
  for (int inode = 0; inode < ele.NumNode(); ++inode)
  {
    const int gid = ele.NodeIds()[inode];
    nodalDofSet_.insert(make_pair(gid,dofman.getNodeDofSet(gid)));
  }
  
  // element dofs for ele
  std::set<XFEM::FieldEnr> enrfieldset(dofman.getElementDofSet(ele.Id()));
  
  ComputeDependendInfo(ele, nodalDofSet_, enrfieldset, element_ansatz);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(const ElementDofManager&)
{
  dserror("no copying");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::~ElementDofManager()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string XFEM::ElementDofManager::toString() const
{
  std::stringstream s;
  map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp;
  for (tmp = nodalDofSet_.begin(); tmp != nodalDofSet_.end(); ++tmp)
  {
    const int gid = tmp->first;
    const set <XFEM::FieldEnr> actset = tmp->second;
    for ( std::set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
    {
      s << "Node: " << gid << ", " << var->toString() << endl;
    };
  };
  return s.str();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::AssemblyType XFEM::CheckForStandardEnrichmentsOnly(
    const ElementDofManager&   eleDofManager,
    const int                  numnode,
    const int*                 nodeids
)
{
  // find out whether we can use standard assembly or need xfem assembly
  XFEM::AssemblyType assembly_type = XFEM::standard_assembly;
  for (int inode = 0; inode < numnode; ++inode)
  {
    if (assembly_type == XFEM::xfem_assembly)
    {
      break;
    }
    const int gid = nodeids[inode];
    const std::set<XFEM::FieldEnr>& fields = eleDofManager.FieldEnrSetPerNode(gid);
    if (fields.size() != 4)
    {
      assembly_type = XFEM::xfem_assembly;
      break;
    };
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fields.begin(); fieldenr != fields.end(); ++fieldenr)
    {
      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        assembly_type = XFEM::xfem_assembly;
        break;
      };
    };
  };
  const int eledof = eleDofManager.NumDofPerElement();
  if (eledof != 0)
  {
    assembly_type = XFEM::xfem_assembly;
  }
  
  return assembly_type;
}

    
#endif  // #ifdef CCADISCRET
