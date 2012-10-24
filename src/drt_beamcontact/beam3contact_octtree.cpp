/*!----------------------------------------------------------------------
\file beam3contact.cpp
\brief Octtree for beam contact search

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*----------------------------------------------------------------------*/

#include "beam3contact_octtree.H"
#include "beam3contact.H"
#include "beam3contact_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_Time.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <map>
#include <math.h>

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/11|
 *----------------------------------------------------------------------*/
Beam3ContactOctTree::Beam3ContactOctTree(ParameterList& params, DRT::Discretization& discret,DRT::Discretization& searchdis, const int& dofoffset):
discret_(discret),
searchdis_(searchdis),
basisnodes_(discret.NumGlobalNodes()),
dofoffset_(dofoffset)
{
  // octree specs
  // extrusion factor
  extrusionfactor_ = params.get<double>("BEAMS_EXTFAC", 1.05);
  // extrusion factor
  radialextrusion_ = params.get<double>("BEAMS_RADFAC", 1.05);
  // max tree depth
  maxtreedepth_ = params.get<int>("BEAMS_TREEDEPTH", 6);
  // max number of bounding boxes per leaf octant
  minbboxesinoctant_ = params.get<int>("BEAMS_BOXESINOCT",8);

  // set flag signaling the existence of periodic boundary conditions
  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // retrieve the dimensions of the periodic boundary box
  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(statmechparams,"PERIODLENGTH"));
    std::string word;
    char* input;
    while (PL >> word)
      periodlength_->push_back(std::strtod(word.c_str(), &input));
  }
  if((int)periodlength_->size()<3)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());

  if(periodlength_->at(0)>1e-12)
    periodicBC_ = true;
  else
    periodicBC_ = false;

  // determine bounding box type
  string boundingbox = params.get<string>("BEAMS_OCTREE","None");
  if(boundingbox == "octree_axisaligned")
  {
    if(!discret_.Comm().MyPID())
      cout<<"Search routine:\nOctree + Axis Aligned BBs"<<endl;
    boundingbox_ = Beam3ContactOctTree::axisaligned;
  }
  else if(boundingbox == "octree_cylorient")
  {
    if(!discret_.Comm().MyPID())
      cout<<"Search routine:\nOctree + Cylindrical Oriented BBs"<<endl;
    boundingbox_ = Beam3ContactOctTree::cyloriented;
  }
  else if(boundingbox == "octree_spherical")
  {
    if(!discret_.Comm().MyPID())
      cout<<"Search routine:\nOctree + Spherical BBs"<<endl;
    boundingbox_ = Beam3ContactOctTree::spherical;
  }
  else
    dserror("No Octree declared in your Input file!");

  if(!discret_.Comm().MyPID())
    cout<<"max. tree depth        = "<<maxtreedepth_<<"\nmax. BB per octant     = "<<minbboxesinoctant_<<"\nextrusion factor       = "<<extrusionfactor_<<endl;

  // get line conditions
  bbox2line_ = rcp(new Epetra_Vector(*(searchdis_.NodeColMap())));
  bbox2line_->PutScalar(-1.0);
  std::vector<DRT::Condition*> lines;
  discret_.GetCondition("FilamentNumber", lines);

  if((int)lines.size()==0)
    dserror("For octree-based search,define line conditions in input file section FILAMENT NUMBERS.");

  for(int i=0; i<(int)lines.size(); i++)
    for(int j=0; j<(int)lines[i]->Nodes()->size(); j++)
      (*bbox2line_)[searchdis_.NodeColMap()->LID( lines[i]->Nodes()->at(j))] = lines[i]->GetInt("Filament Number");

  return;
}


/*----------------------------------------------------------------------*
 |  calls the almighty Octtree (public)                      meier 01/11|
 *----------------------------------------------------------------------*/
vector<RCP<Beam3contact> > Beam3ContactOctTree::OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, int step)
{
#ifdef OCTREEDEBUG
  double t_start = Teuchos::Time::wallTime();
#endif
  // initialize class vectors
  InitializeOctreeSearch();
  // build axis aligned bounding boxes
  CreateBoundingBoxes(currentpositions);
  // call recursive octtree build
  // clear vector for assigning bounding boxes to octants to be on the safe side before (re)assigning bounding boxes
  bool bboxesfound = locateAll();
  // intersection checks
  vector<RCP<Beam3contact> > contactpairs;
  if(bboxesfound)
  {
    BoundingBoxIntersection(currentpositions, &contactpairs);
    // output
    OctreeOutput(contactpairs , step);
  }
  else
    contactpairs.clear();
#ifdef OCTREEDEBUG
  if(!discret_.Comm().MyPID())
    cout<<"Octree Search time:\t\t"<<Teuchos::Time::wallTime()-t_start<<" seconds"<<endl;
#endif
  return contactpairs;
}// OctTreeSearch()
/*----------------------------------------------------------------------*
 |  Return the octants to which this bounding box belongs               |
 |  (public)                                               mueller 01/11|
 *----------------------------------------------------------------------*/
std::vector<int> Beam3ContactOctTree::InWhichOctantLies(const int& thisBBoxID)
{
  std::vector<int> octants(bbox2octant_->NumVectors(),-1);
  int bboxcolid = searchdis_.ElementColMap()->LID(thisBBoxID);
  for(int i=0; i<bbox2octant_->NumVectors(); i++)
    octants[i] = (int)(*bbox2octant_)[i][bboxcolid];
  return octants;
}
/*----------------------------------------------------------------------*
 |  Intersect the bounding boxes of a certain octant with a given       |
 |  bounding box (public)                                  mueller 01/11|
 *----------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectBBoxesWith(Epetra_SerialDenseMatrix& nodecoords, Epetra_SerialDenseMatrix& nodeLID)
{
  /* notes:
   * 1) do not apply this before having constructed the octree. This is merely a query tool
   * 2)"boxid" does not necessarily coincide with the bounding box we are going to intersect with the other boxes
   * in the octant. The reason: The bounding box may actually not exist.
   * Of course, if it does exist, "boxid" will be the id of the bounding box we actually want to check the
   * other boxes against. However, if the bounding box is merely a hypothetical construct (i.e. there is no actual beam element), then
   * we have to give a box id that does exist in order to find the correct octant. Ideally, that means that "boxid" should be the ID of
   * a bounding box which is a direct neighbor of our (hypothetical) bounding box.
   * 3) nodecoords are the coordinates of the nodes of the (non-)existing element.
   */
  bool intersection = false;

  // determine bounding box limits
  RCP<Epetra_SerialDenseMatrix> bboxlimits = rcp(new Epetra_SerialDenseMatrix(1,1));

  // build bounding box according to given type
  switch(boundingbox_)
  {
    case Beam3ContactOctTree::axisaligned:
      CreateAABB(nodecoords,0,bboxlimits);
    break;
    case Beam3ContactOctTree::cyloriented:
      CreateCOBB(nodecoords,0,bboxlimits);
    break;
    default: dserror("No or an invalid Octree type was chosen. Check your input file!");
  }

  // retrieve octants in which the bounding box with ID thisBBoxID is located
  std::vector<std::vector<int> > octants;
  octants.clear();
  // get the octants for two bounding boxe (element) GIDs adjacent to each given node LID
  for(int i=0; i<nodeLID.M(); i++)
    octants.push_back(InWhichOctantLies(searchdis_.lColNode((int)nodeLID(i,0))->Elements()[0]->Id()));

  // intersection of given bounding box with all other bounding boxes in given octant
  for(int ibox=0; ibox<(int)octants.size(); ibox++)
  {
    for(int oct=0; oct<(int)octants[ibox].size(); oct++)
    {
      if(octants[ibox][oct]!=-1)
      {
        for(int i=0; i<bboxesinoctants_->NumVectors(); i++)
        {
          // take only values of existing bounding boxes and not the filler values (-9)
          if((int)(*bboxesinoctants_)[i][octants[ibox][oct]]>-0.9)
          {
            // get the second bounding box ID
            std::vector<int> bboxinoct(2,-1);
            bboxinoct[0] = (int)(*bboxesinoctants_)[i][octants[ibox][oct]];
            /*check for adjacent nodes: if there are adjacent nodes, then, of course, there
             * the intersection test will turn out positive. We skip those cases.*/
            // note: bounding box IDs are equal to element GIDs
            bool sharednode = false;
            for(int j=0; j<searchdis_.gElement(bboxinoct[0])->NumNode(); j++)
            {
              for(int k=0; k<(int)nodeLID.M(); k++)
              {
                if(searchdis_.NodeColMap()->LID(searchdis_.gElement(bboxinoct[0])->NodeIds()[j])==(int)nodeLID(k,0))
                {
                  sharednode = true;
                  break;
                }
              }
              if(sharednode)
                break;
            }
            // apply different bounding box intersection schemes
            if(!sharednode)
            {
              switch(boundingbox_)
              {
                case Beam3ContactOctTree::axisaligned:
                  intersection = IntersectionAABB(bboxinoct, bboxlimits);
                break;
                case Beam3ContactOctTree::cyloriented:
                  intersection = IntersectionCOBB(bboxinoct, bboxlimits);
                break;
                default: dserror("No or an invalid Octree type was chosen. Check your input file!");
              }
            }

            if(intersection)
              break;
          }
          else // loop reached the first bogus value (-9)
            break;
        }
      }
      else
        break;
      if(intersection)
        break;
    }
    if(intersection)
      break;
  }

  return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Output of octants, bounding boxes and contact pairs (public)       mueller 01/12 |
 *----------------------------------------------------------------------------------.*/
void Beam3ContactOctTree::OctreeOutput(std::vector<RCP<Beam3contact> >& cpairs, int step)
{
  if(!discret_.Comm().MyPID() && step!=-1)
  {
    // active contact pairs
    if((int)cpairs.size()>0)
    {
      //Print ContactPairs to .dat-file and plot with Matlab....................
      std::ostringstream filename;
      if(step!=-2)
        filename << "ContactPairs"<<std::setw(6) << setfill('0') << step <<".dat";
      else
        filename << "ContactPairsInit.dat" <<endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int i=0;i<(int)cpairs.size();i++)
        myfile << (cpairs[i]->Element1())->Id() <<"  "<< (cpairs[i]->Element2())->Id() <<endl;
      fprintf(fp, myfile.str().c_str());
      fclose(fp);
    }
    // octant limits output
    if((int)octreelimits_.size()>0)
    {
      std::ostringstream filename;
      if(step!=-2)
        filename << "OctreeLimits"<<std::setw(6) << setfill('0') << step <<".dat";
      else
        filename << "OctreeLimitsInit.dat"<<endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<(int)octreelimits_.size(); u++)
      {
        for (int v=0; v<(int)octreelimits_[u].M(); v++)
          myfile <<scientific<<octreelimits_[u](v)<<" ";
        myfile <<endl;
      }
      // root box
      for(int u=0; u<(int)rootbox_.M(); u++)
        myfile<<scientific<<rootbox_(u)<<" ";
      myfile<<endl;
      fprintf(fp, myfile.str().c_str());
      fclose(fp);

#ifdef OCTREEDEBUG
      for (int u=0; u<(int)octreelimits_.size(); u++)
        for (int v=0; v<(int)octreelimits_[u].M(); v++)
          if(v%2==0 && octreelimits_[u](v)<rootbox_(v) && fabs(octreelimits_[u](v)-rootbox_(v))>1e-8)
            dserror("Octant minimum %4.10f below root box minimum %4.10f",octreelimits_[u](v),rootbox_(v));
          else if(v%2==1 && octreelimits_[u](v)>rootbox_(v) && fabs(octreelimits_[u](v)-rootbox_(v))>1e-8)
            dserror("Octant maximum %4.10f above root box minimum %4.10f",octreelimits_[u](v),rootbox_(v));
#endif
    }
    // bounding box coords output
    if(allbboxes_!=Teuchos::null)
    {
      std::ostringstream filename;
      if(step!=-2)
        filename << "BoundingBoxCoords"<<std::setw(6) << setfill('0') << step <<".dat";
      else
        filename << "BoundingBoxCoordsInit.dat"<<endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<allbboxes_->MyLength(); u++)
      {
        for (int v=0; v<allbboxes_->NumVectors(); v++)
          myfile <<scientific<<setprecision(10)<<(*allbboxes_)[v][u] <<" ";
        myfile <<endl;
      }
      fprintf(fp, myfile.str().c_str());
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)               mueller 01/12|
 |  Initialize octree class vectors and specifications                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::InitializeOctreeSearch()
{
#ifdef OCTREEDEBUG
  if(!discret_.Comm().MyPID())
    cout<<"Searchdis: "<<searchdis_.ElementColMap()->NumMyElements()<<", Probdis: "<<discret_.NumGlobalElements()<<endl;
#endif
  // mapping bounding boxes to octants with -1.0 for empty with 4 columns (max number of octants a single BB can belong to)
  bbox2octant_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),4));
  bbox2octant_->PutScalar(-1.0);
  // number of shifts across volume boundaries in case of periodic boundary conditions (for intersection optimization)
  if(periodicBC_)
    numshifts_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap()),true));

  // determine radius factor by looking at the absolute mean variance of a bounding box (not quite sure...)
  //beam diameter
  diameter_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap())));
  for(int i=0; i<searchdis_.ElementColMap()->NumMyElements(); i++)
  {
    DRT::Element* beamelement = searchdis_.lColElement(i);
    const DRT::ElementType & eot = beamelement->ElementType();

    if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      (*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(beamelement))->Izz()) / M_PI));
    if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
      (*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3ii*>(beamelement))->Izz()) / M_PI));
    // feasibility check
    if ((*diameter_)[i] <= 0.0) dserror("ERROR: Did not receive feasible element radius.");
  }
  // storage of Bounding Boxes
  // (components 0,...,5 contain bounding box limits)
  // (components 6,...,23 contain bounding box limits in case of periodic boundary conditions):
  // a box may b subject to a boundary shift up to 3 times -> 4 segments -> 24 values + 1 bounding box ID
  allbboxes_ = Teuchos::null;
  if(periodicBC_)
    allbboxes_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),4*6+1, true));
  else
    allbboxes_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),7, true));

  return;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)                 meier 01/11|
 |  generates bounding boxes extended with factor 1.05                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions)
{
#ifdef MEASURETIME
  double t_AABB = Teuchos::Time::wallTime();
#endif
  // Get Nodes from discretization....................
  // build bounding boxes according to input parameter
  for (int elecolid=0; elecolid<searchdis_.ElementColMap()->NumMyElements(); elecolid++)
  {
    int elegid = searchdis_.ElementColMap()->GID(elecolid);
    // only do stuff for Row Elements
    if(searchdis_.ElementRowMap()->LID(elegid)>-1)
    {
      // get the element with local ID (LID) elecolid
      DRT::Element* element = searchdis_.lColElement(elecolid);

      // vector for the global IDs (GID) of element
      std::vector<int> nodelids(2,0);
      for(int i=0; i<(int)nodelids.size(); i++)
      {
        int gid = element->Nodes()[i]->Id();
        nodelids[i] = searchdis_.NodeColMap()->LID(gid);
      }

      //store nodal positions into matrix coords
      LINALG::SerialDenseMatrix coord(3,2,true);
      for(int i=0; i<(int)nodelids.size(); i++)
      {
        const map<int, LINALG::Matrix<3,1> >::const_iterator nodepos = currentpositions.find(nodelids.at(i));
        for(int j=0; j<coord.M(); j++)
          coord(j,i) = (nodepos->second)(j);
      }

      // build bounding box according to given type
      switch(boundingbox_)
      {
        case Beam3ContactOctTree::axisaligned:
          CreateAABB(coord, elecolid);
        break;
        case Beam3ContactOctTree::cyloriented:
          CreateCOBB(coord, elecolid);
        break;
        case Beam3ContactOctTree::spherical:
          CreateSPBB(coord, elecolid);
        break;
        default: dserror("No or an invalid Octree type was chosen. Check your input file!");
      }
    }
  } //end for-loop which goes through all elements

  // communication of findings
  Epetra_MultiVector allbboxesrow(*(searchdis_.ElementRowMap()),allbboxes_->NumVectors(),true);
  CommunicateMultiVector(allbboxesrow, *allbboxes_);

  if(periodicBC_)
  {
    Epetra_Vector numshiftsrow(*(searchdis_.ElementRowMap()),true);
    CommunicateVector(numshiftsrow, *numshifts_);
  }

#ifdef MEASURETIME
  double bbgentimelocal = Teuchos::Time::wallTime() - t_AABB;
  double bbgentimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&bbgentimelocal, &bbgentimeglobal, 1);

  if(!searchdis_.Comm().MyPID())
    cout << "\n\nBBox creation time:\t\t" << bbgentimeglobal<< " seconds"<<endl;
#endif

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create an Axis Aligned Bounding Box   (private)                           mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateAABB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits seperately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.

  // factor by which the box is extruded in each dimension (->input file parameter??)
  double extrusionfactor = extrusionfactor_;
  if(bboxlimits!=Teuchos::null)
    extrusionfactor = 1.0;

  //Decision if periodic boundary condition is applied....................
  //number of spatial dimensions
  const int ndim = 3;

  if(elecolid<0 || elecolid>=searchdis_.ElementColMap()->NumMyElements())
    dserror("Given Element Column Map ID is %d !", elecolid);

  int elegid = searchdis_.ElementColMap()->GID(elecolid);

  /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
   * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
   * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
   * at PeriodLength);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
   * 0 close to boundary at PeriodLength);*/
  LINALG::Matrix<3,1> cut;
  cut.Clear();

  /* In order to determine the correct vector "dir" of the visualization at the boundaries,
   * a copy of "coord" with adjustments in the proper places is introduced
   * in unshift, always the second node lies outside of the volume*/
  LINALG::SerialDenseMatrix unshift(coord.M(), coord.N());

  // Compute "cut"-matrix (only in case of periodic BCs)
  // number of overall shifts
  int numshifts = 0;
  // dof at which the bounding box segment is shifted (used in case of a single shift)
  int shiftdof = -1;
  if(periodicBC_)
  {
    // We have to first make sure that the given coordinates lie within the boundary volume. Otherwise,
    // the bounding box creation does not work properly. Why do we have to do this? During a Nweton step,
    // the element standing behind this bounding box might be displaced out of the simulated volume.
    // In order for this method to work, we need the first node to be within the volume, the second one
    // potentially outside.
    // shift into volume
    for(int dof=0; dof<coord.M(); dof++)
    {
      for(int node=0; node<coord.N(); node++)
      {
        if(coord(dof,node)>periodlength_->at(dof))
          coord(dof,node) -= periodlength_->at(dof);
        else if(coord(dof,node)<0.0)
          coord(dof,node) += periodlength_->at(dof);
      }
    }
    // shift second node outside of volume if bounding box was cut before
    for (int dof=0; dof<ndim; dof++)
    {
      // initialize unshift with coord values
      unshift(dof,0) = coord(dof,0);
      unshift(dof,1) = coord(dof,1);
      if (fabs(coord(dof,1)-periodlength_->at(dof)-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
      {
        cut(dof) = 1.0;
        shiftdof = dof;
        unshift(dof,1) -= periodlength_->at(dof);
        numshifts++;
      }
      if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
      {
        cut(dof) = 2.0;
        shiftdof = dof;
        unshift(dof,1) += periodlength_->at(dof);
        numshifts++;
      }
    }
    if(bboxlimits!=Teuchos::null)
      bboxlimits = rcp(new Epetra_SerialDenseMatrix((numshifts+1)*6,1));
    else
      (*numshifts_)[elecolid] = numshifts; // store number of shifts
  }
  else
  {
    if(bboxlimits!=Teuchos::null)
      bboxlimits = rcp(new Epetra_SerialDenseMatrix(6,1));
  }
  /* take action according to number of shifts
   * this may seem not too elegant, but consider that among the cut elements the majority is only shifted once.
   * A single shift can be performed with much less computational effort than multiple shifts since in the case of
   * multiple shifts, one has to reiterate the dof-wise shifts, determine the position of the nodes after each shift
   * and calculate the coordinates of the found segments.
   */
  double bboxdiameter = (*diameter_)[elecolid];;
  if(bboxlimits!=Teuchos::null)
    bboxdiameter = (*diameter_)[searchdis_.ElementColMap()->NumMyElements()-1];

  // standard unshifted bounding box
  switch(numshifts)
  {
    case 0:
    {
      //do normal process with nodecoords0 and nodecoords1.
      //Calculate Center Point of AABB
      LINALG::Matrix<3,1> midpoint;
      for(int i=0; i<(int)midpoint.M(); i++)
        midpoint(i) = 0.5*(coord(i,0) + coord(i,1));

      //Calculate edgelength of AABB
      LINALG::Matrix<3,1> edgelength;
      for(int i=0; i<(int)edgelength.M(); i++)
        edgelength(i) = fabs(coord(i,1) - coord(i,0));

      //Check for edgelength of AABB
      for(int i=0; i<(int)edgelength.M(); i++)
        if (edgelength(i)<bboxdiameter)
          edgelength(i) = bboxdiameter;

      // Calculate limits of AABB with extrusion around midpoint
      if(bboxlimits!=Teuchos::null)
      {
        for(int i=0; i<6; i++)
          if(i%2==0)
            (*bboxlimits)(i,0) = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
          else if(i%2==1)
            (*bboxlimits)(i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
      }
      else
      {
        for(int i=0; i<6; i++)
        {
          if(i%2==0)
            (*allbboxes_)[i][elecolid] = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
          else if(i%2==1)
            (*allbboxes_)[i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
        }
      }
    }
    break;
    default: // broken bounding boxes due to periodic BCs
    {
      // directional vector
      LINALG::Matrix<3, 1> dir;
      for (int dof = 0; dof < ndim; dof++)
        dir(dof) = unshift(dof,1) - unshift(dof,0);
      dir.Scale(1.0/dir.Norm2());

      /* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
       * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
       * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
      LINALG::Matrix<3,2> lambdaorder;
      lambdaorder.PutScalar(1e6);
      // collect lambdas
      for(int dof=0; dof<(int)lambdaorder.M(); dof++)
      {
        switch((int)cut(dof))
        {
          case 1:
          {
            lambdaorder(dof,0) = -coord(dof, 0) / dir(dof);
            lambdaorder(dof,1) = dof;
          }
          break;
          case 2:
          {
            lambdaorder(dof,0) = (periodlength_->at(dof) - coord(dof,0)) / dir(dof);
            lambdaorder(dof,1) = dof;
          }
          break;
          default:
          {
            lambdaorder(dof,1) = dof;
          }
        }
      }

      // sort the lambdas (ascending values) and indices accordingly
      // in case of multiple shifts
      if(numshifts>1)
      {
        for(int i=0; i<(int)lambdaorder.M()-1; i++)
          for(int j=i+1; j<(int)lambdaorder.M(); j++)
            if(lambdaorder(j,0)<lambdaorder(i,0))
            {
              double temp = lambdaorder(i,0);
              int tempindex = (int)lambdaorder(i,1);
              lambdaorder(i,0) = lambdaorder(j,0);
              lambdaorder(i,1) = lambdaorder(j,1);
              lambdaorder(j,0) = temp;
              lambdaorder(j,1) = tempindex;
            }
      }
      else  // for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
        for(int i=0; i<(int)lambdaorder.N(); i++)
          lambdaorder(0,i) = lambdaorder(shiftdof,i);

      // calculate segment lambdas
      for(int dof=numshifts-1; dof>0; dof--)
        lambdaorder(dof,0) -= lambdaorder(dof-1,0);

      // the idea is to gradually shift the matrix "unshift" back into the volume and, while doing so,
      // calculate the segments except for the last one
      // determine closest boundary component-wise
      for(int shift=0; shift<numshifts; shift++)
      {
        //second point
        for(int i=0 ;i<unshift.M(); i++)
          unshift(i,1) = unshift(i,0) + lambdaorder(shift,0)*dir(i);

        //Calculate Center Point and edgelength of AABB
        LINALG::Matrix<3,1> midpoint;
        LINALG::Matrix<3,1> edgelength;
        for(int i=0; i<(int)midpoint.M(); i++)
        {
          midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
          edgelength(i) = unshift(i,1) - unshift(i,0);
          //Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
          for(int i=0; i<(int)edgelength.M(); i++)
            if (edgelength(i)<bboxdiameter)
              edgelength(i) = bboxdiameter;
        }
        // Calculate limits of AABB of the current segment (which definitely lies in the volume) with extrusion around midpoint
        if(bboxlimits!=Teuchos::null)
        {
          for(int i=0; i<6; i++)
          {
            if(i%2==0)
              (*bboxlimits)(shift*6+i,0) = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
            else if(i%2==1)
              (*bboxlimits)(shift*6+i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
          }
        }
        else
        {
          for(int i=0; i<6; i++)
          {
            if(i%2==0)
              (*allbboxes_)[shift*6+i][elecolid] = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
            else if(i%2==1)
              (*allbboxes_)[shift*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
          }
        }

        int currshift = (int)lambdaorder(shift,1);
        // shift the coordinates of the second point
        if(cut(currshift)==1.0)
          unshift(currshift,1) += periodlength_->at(currshift);
        else if(cut(currshift)==2.0)
          unshift(currshift,1) -= periodlength_->at(currshift);
        // make second point the first and calculate new second point in the next iteration!
        for(int i=0; i<unshift.M(); i++)
          unshift(i,0) = unshift(i,1);
      }

      // the last segment
      for(int dof=0; dof<unshift.M(); dof++)
        unshift(dof,1) = coord(dof,1);

      //Calculate Center Point and edgelength of AABB
      LINALG::Matrix<3,1> midpoint;
      LINALG::Matrix<3,1> edgelength;
      for(int i=0; i<(int)midpoint.M(); i++)
      {
        midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
        edgelength(i) = unshift(i,1) - unshift(i,0);
        //Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
        if (edgelength(i)<bboxdiameter)
          edgelength(i) = bboxdiameter;
      }

      // limits of the last bounding box
      if(bboxlimits!=Teuchos::null)
      {
        for(int i=0; i<6; i++)
        {
          if(i%2==0)
            (*bboxlimits)(numshifts*6+i,0) = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
          else if(i%2==1)
            (*bboxlimits)(numshifts*6+i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
        }
      }
      else
      {
        for(int i=0; i<6; i++)
        {
          if(i%2==0)
            (*allbboxes_)[numshifts*6+i][elecolid] = midpoint(i/2) - 0.5*edgelength(i/2)*extrusionfactor;
          else if(i%2==1)
            (*allbboxes_)[numshifts*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
        }
      }
    }
  }

  if(bboxlimits==Teuchos::null)
  {
    //fill up the rest of the 24 values with bogus values
    if(periodicBC_ && numshifts<3)
      for(int i=allbboxes_->NumVectors()-1; i>(numshifts+1)*6-1; i--)
        (*allbboxes_)[i][elecolid] = -1e9;

     // store GID (=box number)
     (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;
  }

  //Bring coordinates in case of periodic boundary condition in right order ( "-1e9 signals the bogus value from above)
  //[xmin xmax ymin ymax zmin zmax ...]
  if (periodicBC_ && numshifts>0)
  {
    double minimum =0.0;    double maximum =0.0;

    if(bboxlimits!=Teuchos::null)
    {
      for(int i=6; i<(bboxlimits->M())/2;i++)
      {
        minimum = min((*bboxlimits)(2*i,0),(*bboxlimits)(2*i+1,0));
        maximum = max((*bboxlimits)(2*i,0),(*bboxlimits)(2*i+1,0));
        (*bboxlimits)(2*i,0) = minimum;    (*bboxlimits)(2*i+1,0) = maximum;
      }
    }
    else
    {
      for(int i=6; i<(allbboxes_->NumVectors()-1)/2;i++)
      {
        // leave loop at first bogus entry
        if((2*i)%6==0 && (*allbboxes_)[2*i][elecolid]==-1e9)
          break;
        minimum = min((*allbboxes_)[2*i][elecolid],(*allbboxes_)[2*i+1][elecolid]);
        maximum = max((*allbboxes_)[2*i][elecolid],(*allbboxes_)[2*i+1][elecolid]);
        //cout << minimum << endl;
        (*allbboxes_)[2*i][elecolid] = minimum;    (*allbboxes_)[2*i+1][elecolid] = maximum;
      }// end of correct
    }
  } // end of normal AAABB

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box   (private)                   mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateCOBB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits seperately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.

  double extrusionfactor = extrusionfactor_;
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    extrusionfactor = 1.0;
  const int ndim = 3;
  int elegid = searchdis_.ElementColMap()->GID(elecolid);
  int numshifts = 0;
  int shiftdof = -1;
  LINALG::Matrix<3,1> cut;
  cut.Clear();
  LINALG::SerialDenseMatrix unshift(coord.M(),coord.N());

  if(periodicBC_)
  {
    // shift into volume
    for(int dof=0; dof<coord.M(); dof++)
    {
      for(int node=0; node<coord.N(); node++)
      {
        if(coord(dof,node)>periodlength_->at(dof))
          coord(dof,node) -= periodlength_->at(dof);
        else if(coord(dof,node)<0.0)
          coord(dof,node) += periodlength_->at(dof);
      }
    }

    // shift out second node (if possible)
    for (int dof=0; dof<ndim; dof++)
    {
      // initialize unshift with coord values
      unshift(dof,0) = coord(dof,0);
      unshift(dof,1) = coord(dof,1);
      if (fabs(coord(dof,1)-periodlength_->at(dof)-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
      {
        cut(dof) = 1.0;
        shiftdof = dof;
        unshift(dof,1) -= periodlength_->at(dof);
        numshifts++;
      }
      if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
      {
        cut(dof) = 2.0;
        shiftdof = dof;
        unshift(dof,1) += periodlength_->at(dof);
        numshifts++;
      }
    }
    if(bboxlimits!=Teuchos::null)
      bboxlimits = rcp(new Epetra_SerialDenseMatrix((numshifts+1)*6,1));
    else
      (*numshifts_)[elecolid] = numshifts;
  }
  else
  {
    unshift = coord;
    if(bboxlimits!=Teuchos::null)
      bboxlimits = rcp(new Epetra_SerialDenseMatrix(6,1));
  }

  //directional vector
  LINALG::Matrix<3,1> dir;
  for(int dof=0; dof<(int)dir.M(); dof++)
    dir(dof) = unshift(dof,1) - unshift(dof,0);

  switch(numshifts)
  {
    case 0:
    {
      dir.Scale(extrusionfactor);

      if(bboxlimits!=Teuchos::null)
      {
        for(int dof=0; dof<unshift.M(); dof++)
        {
          (*bboxlimits)(dof,0) = unshift(dof,1)-dir(dof);
          (*bboxlimits)(dof+3,0) = unshift(dof,0)+dir(dof);
        }
      }
      else
      {
        for(int dof=0; dof<unshift.M(); dof++)
        {
          (*allbboxes_)[dof][elecolid] = unshift(dof,1)-dir(dof);
          (*allbboxes_)[dof+3][elecolid] = unshift(dof,0)+dir(dof);
        }
      }
    }
    break;
    default: // broken bounding boxes due to periodic BCs
    {
      /* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
       * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
       * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
      dir.Scale(1.0/dir.Norm2());
      LINALG::Matrix<3,2> lambdaorder;
      lambdaorder.PutScalar(1e6);
      // collect lambdas
      for(int dof=0; dof<(int)lambdaorder.M(); dof++)
      {
        switch((int)cut(dof))
        {
          case 1:
          {
            lambdaorder(dof,0) = -unshift(dof, 0) / dir(dof);
            lambdaorder(dof,1) = dof;
          }
          break;
          case 2:
          {
            lambdaorder(dof,0) = (periodlength_->at(dof) - unshift(dof,0)) / dir(dof);
            lambdaorder(dof,1) = dof;
          }
          break;
          default:
          {
            lambdaorder(dof,1) = dof;
          }
        }
      }
      // sort the lambdas (ascending values) and indices accordingly
      // in case of multiple shifts
      if(numshifts>1)
      {
        for(int j=0; j<(int)lambdaorder.M()-1; j++)
          for(int k=j+1; k<(int)lambdaorder.M(); k++)
            if(lambdaorder(k,0)<lambdaorder(j,0))
            {
              double temp = lambdaorder(j,0);
              int tempindex = (int)lambdaorder(j,1);
              lambdaorder(j,0) = lambdaorder(k,0);
              lambdaorder(j,1) = lambdaorder(k,1);
              lambdaorder(k,0) = temp;
              lambdaorder(k,1) = tempindex;
            }
        // calculate segment lambdas
        for(int dof=numshifts-1; dof>0; dof--)
          lambdaorder(dof,0) -= lambdaorder(dof-1,0);
      }
      else // for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
      {
        for(int n=0; n<(int)lambdaorder.N(); n++)
        {
          double tmp = lambdaorder(shiftdof,n);
          lambdaorder(0,n) = tmp;
        }
      }

      for(int shift=0; shift<numshifts; shift++)
      {
        //second point
        for(int dof=0 ;dof<unshift.M(); dof++)
          unshift(dof,1) = unshift(dof,0) + lambdaorder(shift,0)*dir(dof);
        // Calculate limits of the bounding box segment (convenient because lambdas are segment lengths)
        if(bboxlimits!=Teuchos::null)
        {
          for(int dof=0; dof<unshift.M(); dof++)
          {
             (*bboxlimits)(shift*6+dof,0) = unshift(dof,1)-lambdaorder(shift,0)*extrusionfactor*dir(dof);
             (*bboxlimits)(shift*6+dof+3,0) = unshift(dof,0)+lambdaorder(shift,0)*extrusionfactor*dir(dof);
          }
        }
        else
        {
          for(int dof=0; dof<unshift.M(); dof++)
          {
             (*allbboxes_)[shift*6+dof][elecolid] = unshift(dof,1)-lambdaorder(shift,0)*extrusionfactor*dir(dof);
             (*allbboxes_)[shift*6+dof+3][elecolid] = unshift(dof,0)+lambdaorder(shift,0)*extrusionfactor*dir(dof);
          }
        }
        int currshift = (int)lambdaorder(shift,1);
        if(cut(currshift)==1.0)
          unshift(currshift,1) += periodlength_->at(currshift);
        else if(cut(currshift)==2.0)
          unshift(currshift,1) -= periodlength_->at(currshift);
        for(int dof=0; dof<unshift.M(); dof++)
          unshift(dof,0) = unshift(dof,1);
      }

      // the last segment
      double llastseg = 0.0;
      for(int dof=0; dof<unshift.M(); dof++)
      {
        unshift(dof,1) = coord(dof,1);
        llastseg += (unshift(dof,1)-unshift(dof,0))*(unshift(dof,1)-unshift(dof,0));
      }
      llastseg  = sqrt(llastseg);

      // limits of the last bounding box
      if(bboxlimits!=Teuchos::null)
      {
        for(int dof=0; dof<unshift.M(); dof++)
        {
          (*bboxlimits)(numshifts*6+dof,0) = unshift(dof,1)-llastseg*extrusionfactor*dir(dof);
          (*bboxlimits)(numshifts*6+dof+3,0) = unshift(dof,0)+llastseg*extrusionfactor*dir(dof);
        }
      }
      else
      {
        for(int dof=0; dof<unshift.M(); dof++)
        {
          (*allbboxes_)[numshifts*6+dof][elecolid] = unshift(dof,1)-llastseg*extrusionfactor*dir(dof);
          (*allbboxes_)[numshifts*6+dof+3][elecolid] = unshift(dof,0)+llastseg*extrusionfactor*dir(dof);
        }
      }
    }
  }

  // fill all latter entries  except for the last one (->ID) with bogus values (in case of periodic BCs)
  if(bboxlimits==Teuchos::null)
  {
    if(periodicBC_ && numshifts<3)
      for(int i=allbboxes_->NumVectors()-1; i>(numshifts+1)*6-1; i--)
        (*allbboxes_)[i][elecolid] = -1e9;
    // last entry: element GID
    (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;
  };
  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box   (private)                   mueller 1/12|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateSPBB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  if(bboxlimits!=Teuchos::null)
  {
    bboxlimits = rcp(new Epetra_SerialDenseMatrix(3,1));
    for(int dof=0; dof<coord.M(); dof++)
      (*bboxlimits)(dof,0) = coord(dof,0);
  }
  else
  {
    for(int dof=0; dof<coord.M(); dof++)
      (*allbboxes_)[dof][elecolid] = coord(dof,0);
  }
  return;
}
/*-----------------------------------------------------------------------------------------*
 |  locateAll function (private); Recursive division of a 3-dimensional set.    meier 02/11|
 |  [Oct_ID, element-ID] = locateAll( &allbboxes_)                                          |
 |  Performs recursive tree-like division of a set of Bounding Boxes.                      |
 |  N0 is maximum permissible number of "counted" boxes in the leaf octant.                |
 |  Returns vector IND of the same size as rows of allbboxes_ showing which region each     |
 |  box of a set belongs to; binary matrices BX, BY, BZ where each row shows               |
 |  "binary address" of each region are written to .dat- files each.                       |
 *----------------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::locateAll()
{
#ifdef MEASURETIME
  double t_octree = Teuchos::Time::wallTime();
#endif
  // get the root box
  rootbox_ = GetRootBox();

  // Convert Epetra_MultiVector allbboxes_ to vector(vector<double>)
  std::vector<std::vector<double> > allbboxesstdvec(allbboxes_->MyLength(), std::vector<double>(allbboxes_->NumVectors(),0.0));
  for(int i=0; i < allbboxes_->MyLength(); i++)
    for(int j=0; j<allbboxes_->NumVectors(); j++)
      allbboxesstdvec[i][j] = (*allbboxes_)[j][i];

  //initial tree depth value (will be incremented with each recursive call of locateBox()
  int treedepth = 0;

  // Parameters and Initialization
  std::vector<std::vector<int> > bboxesinoctants;
  bboxesinoctants.clear();
  octreelimits_.clear();
  // Recursively construct octree; Proc 0 only (parallel computing impossible)
  if(searchdis_.Comm().MyPID()==0)
    locateBox(allbboxesstdvec, rootbox_, octreelimits_, bboxesinoctants, treedepth);

  Epetra_MultiVector bbox2octantrow(*(searchdis_.ElementRowMap()),4,true);
  CommunicateMultiVector(bbox2octantrow, *bbox2octant_, true);

  //determine maximum depth of OctreeMap
  int maxdepthlocal = 0;
  int bboxlengthlocal = 0;
  if(discret_.Comm().MyPID()==0)
  {
    bboxlengthlocal = (int)bboxesinoctants.size();
    for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
      if((int)bboxesinoctants[i].size()>maxdepthlocal)
        maxdepthlocal = (int)bboxesinoctants[i].size();
  }

  int maxdepthglobal = 0;
  int bboxlengthglobal = 0;
  discret_.Comm().MaxAll(&maxdepthlocal, &maxdepthglobal, 1);
  discret_.Comm().MaxAll(&bboxlengthlocal, &bboxlengthglobal, 1);

  /* build temporary, fully overlapping map and row map for octree
   * Note: maxdepthglobal does not occur for a converging Newton iteration. Yet, in some cases, when
   * encountering divergence for the Newton scheme, this might happen.
   * In biopolymer network simulations, this setting is not unlikely and unavoidable. a maximum depth of
   * 0 means, there are no bounding boxes/elements in any octants. Hence, we will not detect any contact and
   * therefore skip the rest of the octree algorithm.*/
  if(maxdepthglobal>0)
  {
    // create octree maps
    std::vector<int> gids;
    for (int i=0 ; i<bboxlengthglobal; i++ )
      gids.push_back(i);
    // crosslinker column and row map
    Epetra_Map octtreerowmap((int)gids.size(), 0, discret_.Comm());
    Epetra_Map octtreemap(-1, (int)gids.size(), &gids[0], 0, discret_.Comm());

    // build Epetra_MultiVectors which hold the BBs of the OctreeMap; for communication
    bboxesinoctants_ = rcp(new Epetra_MultiVector(octtreemap,maxdepthglobal));
    Epetra_MultiVector bboxinoctrow(octtreerowmap,maxdepthglobal, true);

    // fill bboxinoct for Proc 0
    if(searchdis_.Comm().MyPID()==0)
    {
      bboxesinoctants_->PutScalar(-9.0);
      for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
        for(int j=0; j<(int)bboxesinoctants[i].size(); j++)
          (*bboxesinoctants_)[j][i] = bboxesinoctants[i][j];
    }

    // Communication
    CommunicateMultiVector(bboxinoctrow, *bboxesinoctants_,true);

#ifdef OCTREEDEBUG
    std::ostringstream filename;
    if(!discret_.Comm().MyPID())
    {
      filename << "BBinOct.dat";
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<bboxesinoctants_->MyLength(); u++)
      {
        for (int v=0; v<bboxesinoctants_->NumVectors(); v++)
          myfile <<scientific<<(*bboxesinoctants_)[v][u] <<" ";
        myfile <<endl;
      }
      fprintf(fp, myfile.str().c_str());
      fclose(fp);

      cout<<"bboxesinoctants_ : "<<bboxesinoctants_->MyLength()<<"x"<<bboxesinoctants_->NumVectors()<<endl;
    }
#endif
#ifdef MEASURETIME
     if(!searchdis_.Comm().MyPID())
       cout << "\nOctree building time:\t\t" << Teuchos::Time::wallTime() - t_octree<< " seconds" << endl;
#endif
     return true;
  }
  else
  {
    return false;
  }
}// end of method locateAll


/*-----------------------------------------------------------------------------------------*
 |  locateBox function (private);                                               meier 02/11|
 |  [ind, bx, by, bz] = locateBox(&allbboxes, lim, n0)                                     |
 |  Primitive for locateAll                                                                |
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateBox(std::vector<std::vector<double> >& allbboxesstdvec,
                                    LINALG::Matrix<6,1>& lim,
                                    std::vector<LINALG::Matrix<6,1> >& OctreeLimits,
                                    std::vector<std::vector<int> >& bboxesinoctants,
                                    int& treedepth)
{
  // Divide further
  double extrusionfactor = max(extrusionfactor_,radialextrusion_);
  // Center of octant
  LINALG::Matrix<3,1> center;
  // edge length vector of the suboctants
  LINALG::Matrix<3,1> newedgelength;
  for(int i=0; i<(int)center.M(); i++)
  {
    center(i) = (lim(2*i)+lim(2*i+1))/2.0;
    newedgelength(i) = fabs(lim(2*i+1)-(lim)(2*i))/2.0;
  }
  std::vector<LINALG::Matrix<6,1> > limits;
  limits.clear();
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        LINALG::Matrix<6,1> sublim;
        sublim(0) = center(0) + (i-1)*newedgelength(0);
        sublim(1) = center(0) +  i   *newedgelength(0);
        sublim(2) = center(1) + (j-1)*newedgelength(1);
        sublim(3) = center(1) +  j   *newedgelength(1);
        sublim(4) = center(2) + (k-1)*newedgelength(2);
        sublim(5) = center(2) +  k   *newedgelength(2);

        limits.push_back(sublim);
      }

  /* Decision to which child box belongs....................
  *
  *           5 ======================== 7
  *           //|                       /||
  *          // |                      //||
  *         //  |                     // ||
  *        //   |                    //  ||
  *       //    |                   //   ||
  *      //     |                  //    ||
  *     //      |                 //     ||
  *    1 ========================= 3     ||
  *    ||       |                ||      ||
  *    ||       |                ||      ||
  *    ||       |      o (center)||      ||
  *    ||      4 ----------------||------ 6
  *    ||      /                 ||     //
  *    ||     /                  ||    //
  *    ||    /                   ||   //
  *    ||   /                    ||  //
  *    ||  /                     || //      y  z
  *    || /                      ||//       | /
  *    ||/                       ||/        |/
  *    0 ========================= 2        ---> x
  *
  */
  Teuchos::RCP<LINALG::Matrix<3,1> > octcenter = Teuchos::null;

  //Goes through all suboctants
  for( int oct=0; oct<8; oct++)
  {
    // Define temporary vector of same size as current allbboxesstdvec
    std::vector<std::vector<double> > bboxsubset;
    bboxsubset.clear();

    // we need the octant centers when applying cylindrical bounding boxes.
    if(boundingbox_==Beam3ContactOctTree::cyloriented)
    {
      octcenter = Teuchos::rcp(new LINALG::Matrix<3,1>);
      for(int i=0; i<(int)octcenter->M(); i++)
        (*octcenter)(i) = (limits[oct](2*i)+limits[oct](2*i+1))/2.0;
    }

    if(periodicBC_)
    {
      for( int i=0; i<(int)allbboxesstdvec.size(); i++)
      {
        // flag for a bounding box located in the octant or so close to it that its cylindrical hull intersects with the octant
        bool inoctant = false;
        // a bounding box is at maximum divided into 4 subsegments due to periodic boundary conditions
        for(int isub=0; isub<4; isub++)
        {
          // 1) remember: the gid of the bounding box (=element gid) is at the last position
          // 2) loop over the limits of the current octant and check if the current bounding box lies within this octant.
          // 3) Then, check componentwise and leave after first "hit"
          if(allbboxesstdvec[i][6*isub]!=-1e9)
          {
            switch(boundingbox_)
            {
              case Beam3ContactOctTree::axisaligned:
              {
                if(!((limits[oct](0) >= allbboxesstdvec[i][6*isub+1]) || (allbboxesstdvec[i][6*isub] >= limits[oct](1)) || (limits[oct](2) >= allbboxesstdvec[i][6*isub+3]) || (allbboxesstdvec[i][6*isub+2] >= limits[oct](3)) || (limits[oct](4) >= allbboxesstdvec[i][6*isub+5]) || (allbboxesstdvec[i][6*isub+4] >= limits[oct](5))))
                {
                  bboxsubset.push_back(allbboxesstdvec[i]);
                  inoctant = true;
                }
              }
              break;
              case Beam3ContactOctTree::cyloriented:
              {
                // loop over end points of the bounding box
                for(int j=0; j<2; j++)
                {
                  // Idea: The largest absolute component value of the directional vector v from octant center to
                  // bounding box end point position indicates the octant face which is intersected first by the line with
                  // direction v. Octant faces are each parallel to one of the global spatial directions.

                  // component value of the directional vector from octant center to bounding box end point position
                  double vmax = allbboxesstdvec[i][6*isub+3*j]-(*octcenter)(0);
                  // distance between bounding j-th bounding box end point and oct-th octant center
                  double d = vmax * vmax;
                  // index for the maximum absolute value of the directional vector
                  int kmax = 0;
                  for(int k=1; k<(int)octcenter->M(); k++)
                  {
                    d += (allbboxesstdvec[i][6*isub+3*j+k]-(*octcenter)(k))*(allbboxesstdvec[i][6*isub+3*j+k]-(*octcenter)(k));
                    if(fabs(allbboxesstdvec[i][6*isub+3*j+k]-(*octcenter)(k))>fabs(vmax))
                    {
                      vmax = allbboxesstdvec[i][6*isub+3*j+k]-(*octcenter)(k);
                      kmax = k;
                    }
                  }
                  d = sqrt(d);

                  double boxradius = 0.5*extrusionfactor*(*diameter_)[searchdis_.ElementColMap()->LID((int)allbboxesstdvec[i][(int)allbboxesstdvec[i].size()-1])];

                  if(d<=(0.5*newedgelength(kmax)*sqrt(3.0))+boxradius)
                  {
                    // unit vector component
                    vmax /= d;
                    // normal component
                    double normal = -1.0;
                    // note: it's always n!=0
                    if(vmax<0)
                      normal = 1.0;
                    // segment length from octant center to intersection the line with directional vector v and the closest octant face
                    double lambda = -0.5*newedgelength(kmax)/(vmax*normal);

                    // 2 cases: end point is in the octant or it is outside but its cylindrical hull intersects with the octant
                    if(lambda>=d || (d>lambda && d-lambda<=boxradius))
                    {
                      inoctant = true;
                      bboxsubset.push_back(allbboxesstdvec[i]);
                      // Since we found a bounding box end point to lie in the octant, we do not need to investigate further
                      break; //j-loop
                    }
                  }
                }
              }
              break;
              default: dserror("No or an invalid Octree type was chosen. Check your input file!");
            }

            if(inoctant) // isub-loop
              break;
          }
          else
            break; // isub-loop
        } // end of isub-loop
      } // end of for-loop which goes through all elements of input
    }
    else // standard procedure without periodic boundary conditions
    {
      for( int i=0; i<(int)allbboxesstdvec.size(); i++)
      {
        // Processes colums indices 1 to 6
        // 2)loop over the limits of the current octant and check if the current bounding box lies within this octant.
        // 3)Then, check componentwise and leave after first "hit"
        switch(boundingbox_)
        {
          case Beam3ContactOctTree::axisaligned:
          {
            if(!((limits[oct](0) >= allbboxesstdvec[i][1]) || (allbboxesstdvec[i][0] >= limits[oct](1)) || (limits[oct](2) >= allbboxesstdvec[i][3]) || (allbboxesstdvec[i][2] >= limits[oct](3)) || (limits[oct](4) >= allbboxesstdvec[i][5]) || (allbboxesstdvec[i][4] >= limits[oct](5))))
              bboxsubset.push_back(allbboxesstdvec[i]);
          }
          break;
          case Beam3ContactOctTree::cyloriented:
          {
            for(int j=0; j<2; j++)
            {
              int kmax = 0;
              double vmax = allbboxesstdvec[i][3*j] - (*octcenter)(0);
              double d = vmax*vmax;
              for(int k=1; k<3; k++)
              {
                d += (allbboxesstdvec[i][3*j+k] - (*octcenter)(k))*(allbboxesstdvec[i][3*j+k] - (*octcenter)(k));
                if(fabs(allbboxesstdvec[i][3*j+k] - (*octcenter)(k))>fabs(vmax))
                {
                  vmax = allbboxesstdvec[i][3*j+k] - (*octcenter)(k);
                  kmax = k;
                }
              }
              d = sqrt(d);

              double boxradius = 0.5*extrusionfactor*(*diameter_)[searchdis_.ElementColMap()->LID((int)allbboxesstdvec[i][(int)allbboxesstdvec[i].size()-1])];
              if(d<=newedgelength(kmax)/2.0*sqrt(3.0)+boxradius)
              {
                double normal = -1.0;
                if(vmax<0)
                  normal = 1.0;

                vmax /= d;
                double lambda = - newedgelength(kmax)/ (2.0*vmax*normal);

                if(lambda>=d || (d>lambda && d-lambda<=boxradius))
                {
                  bboxsubset.push_back(allbboxesstdvec[i]);
                  break;
                }
              }
            }
          }
          break;
          default: dserror("No or an invalid Octree type was chosen. Check your input file!");
        }
      }// end of for-loop which goes through all elements of input
    }

    // current tree depth
    int currtreedepth = treedepth+1;
    // Check for further recursion by checking number of boxes in octant (first criterion)....................
    int N = (int)bboxsubset.size();

    //If to divide further, let LocateBox call itself with updated inputs
    if (N > minbboxesinoctant_ && currtreedepth < maxtreedepth_-1)
      locateBox(bboxsubset, limits[oct], OctreeLimits, bboxesinoctants, currtreedepth);
    else
    {
      // no further discretization of the volume because either the maximal tree depth or the minimal number of bounding
      // boxes per octant has been reached
      // this vector holds the IDs of the bounding boxes in this octant
      if(N>0)
      {
        std::vector<int> boxids;
        boxids.clear();
        //Push back Limits of suboctants to OctreeLimits
        OctreeLimits.push_back(limits[oct]);

        for (int m = 0; m < (int)bboxsubset.size(); m++)
        {
          // note: the Bounding Box ID is the last column entry of the m-th entry vector bboxsubset
          boxids.push_back((int)bboxsubset[m][(int)bboxsubset[m].size()-1]);
          // assign current octant number to the bounding box
          for(int n=0; n<bbox2octant_->NumVectors(); n++)
            if((*bbox2octant_)[n][searchdis_.ElementColMap()->LID((int)boxids.size()-1)]<-0.9)
            {
              (*bbox2octant_)[n][searchdis_.ElementColMap()->LID((int)boxids.size()-1)] = (double)bboxesinoctants.size();
              break; //leave after finding first empty slot
            }
        }
        // add bounding box IDs of this octant to the global vector
        bboxesinoctants.push_back(boxids);
      }
    }
  }// end of loop which goes through all suboctants
  return;
} // end of method locateBox

/*----------------------------------------------------------------------------------*
 |  Calculate limits of the root box   (public)                        mueller 11/11|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<6,1> Beam3ContactOctTree::GetRootBox()
{
  LINALG::Matrix<6,1> lim;
  // if peLriodic BCs are applied
  if(periodicBC_)
  {
    Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
    for(int i=0; i<(int)lim.M(); i++)
    {
      if(i%2==0)
        lim(i)= 0.0;
      else
        lim(i) =  periodlength_->at((i-1)/2);
    }
  }
  else // standard procedure to find root box limits
  {
    // initialize
    for(int i=0; i<(int)lim.M(); i++)
         if(i%2==0)
           lim(i) = 1e13;
         else
           lim(i) = -1e13;

    switch(boundingbox_)
    {
      case Beam3ContactOctTree::axisaligned:
      {
        // loop over allbboxes_ and determine the extremes
        for(int i=0; i<allbboxes_->MyLength(); i++)
          for(int j=0; j<allbboxes_->NumVectors()-1; j++)
            if(j%2==0 && (*allbboxes_)[j][i]<lim(j)) // mins
              lim(j) = (*allbboxes_)[j][i];
            else if(j%2!=0 && (*allbboxes_)[j][i]>lim(j)) // maxes
              lim(j) = (*allbboxes_)[j][i];
        // determine bounds for cubic root box
        LINALG::Matrix<3,1> maxdist;
        for(int i=0; i<(int)maxdist.M(); i++)
          maxdist(i) = fabs(lim(2*i)+lim(2*i+1) / 2.0 - lim(2*i));
        for(int i=0; i<(int)maxdist.M(); i++)
        {
          lim(2*i) = lim(2*i)+lim(2*i+1) / 2.0 - maxdist.MaxValue()*extrusionfactor_;
          lim(2*i+1) = lim(2*i) + 2.0 * maxdist.MaxValue()*extrusionfactor_;
        }
      }
      break;
      case Beam3ContactOctTree::cyloriented:
      {
        for(int i=0; i<allbboxes_->MyLength(); i++) // loop boxes
        {
          for(int j=0; j<allbboxes_->NumVectors()-1; j++)
          {
            // x
            if(j%3==0)
            {
              // min
              if((*allbboxes_)[j][i]<lim(0))
                lim(0) = (*allbboxes_)[j][i];
              // max
              else if((*allbboxes_)[j][i]>lim(1))
                lim(1) = (*allbboxes_)[j][i];
            }
            // y
            else if(j%3-1==0)
            {
              // min
              if((*allbboxes_)[j][i]<lim(2))
                lim(2) = (*allbboxes_)[j][i];
              // max
              else if((*allbboxes_)[j][i]>lim(3))
                lim(3) = (*allbboxes_)[j][i];
            }
            // z
            else if(j%3-2==0)
            {
              // min
              if((*allbboxes_)[j][i]<lim(4))
                lim(4) = (*allbboxes_)[j][i];
              // max
              else if((*allbboxes_)[j][i]>lim(5))
                lim(5) = (*allbboxes_)[j][i];
            }
          }
        }
        // determine bounds for cubic root box
        LINALG::Matrix<3,1> maxdist;
        for(int i=0; i<(int)maxdist.M(); i++)
          maxdist(i) = fabs(lim(2*i)+lim(2*i+1) / 2.0 - lim(2*i));
        for(int i=0; i<(int)maxdist.M(); i++)
        {
          lim(2*i) = lim(2*i)+lim(2*i+1) / 2.0 - maxdist.MaxValue()*extrusionfactor_;
          lim(2*i+1) = lim(2*i) + 2.0 * maxdist.MaxValue()*extrusionfactor_;
        }
      }
      break;
      case Beam3ContactOctTree::spherical:
      {

      }
      break;
      default: dserror("selected bounding box typ is not implemented!");
    }
  }
  return lim;
}

/*-----------------------------------------------------------------------------------*
 |  Bounding Box Intersection function (private)                          meier 01/11|
 |  Intersects Bounding Boxes in same line of map OctreeMap                          |
 |  Gives back vector of intersection pairs                                          |
 *----------------------------------------------------------------------------------*/
void Beam3ContactOctTree::BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
                                                  vector<RCP<Beam3contact> >* contactpairs)
{
#ifdef MEASURETIME
  double t_search = Teuchos::Time::wallTime();
#endif
  // Build contact pair Map
  std::map<int, std::vector<int> > contactpairmap;
  // create contact pair vector, redundant on all Procs; including redundant pairs
  //for-loop lines of map
  for (int i=0 ; i<bboxesinoctants_->MyLength(); i++ )
  {
    //for-loop index first box
    for(int j=0; j<bboxesinoctants_->NumVectors(); j++)
    {
      // first box ID
      std::vector<int> bboxIDs(2,0);
      bboxIDs[0] = (int)(*bboxesinoctants_)[j][i];

      //for-loop second box
      for(int k=j+1; k<bboxesinoctants_->NumVectors(); k++)
      {
        bboxIDs[1] = (int)(*bboxesinoctants_)[k][i];

        // exclude element pairs sharing one node
        // contact flag
        bool considerpair = false;
        // only consider existing bounding boxes, i.e. no dummy entries "-9.0"
        if(bboxIDs[0]>-1 && bboxIDs[1]>-1)
        {
          considerpair = true;
          DRT::Element* element1 = searchdis_.gElement(bboxIDs[0]);
          DRT::Element* element2 = searchdis_.gElement(bboxIDs[1]);

          for(int k=0; k<element1->NumNode(); k++)
          {
            for(int l=0; l<element2->NumNode(); l++)
            {
              if((*bbox2line_)[searchdis_.NodeColMap()->LID(element1->NodeIds()[k])]==(*bbox2line_)[searchdis_.NodeColMap()->LID(element2->NodeIds()[l])])
              {
                considerpair = false;
                break;
              }
            }
            if(!considerpair)
              break;
          }
        }

        if (considerpair)
        {
          //cout<<"IDs: "<<bboxIDs[0]<<", "<< bboxIDs[1]<<endl;
          // apply different bounding box intersection schemes
          bool intersection = false;
          switch(boundingbox_)
          {
            case Beam3ContactOctTree::axisaligned:
              intersection = IntersectionAABB(bboxIDs);
            break;
            case Beam3ContactOctTree::cyloriented:
              intersection = IntersectionCOBB(bboxIDs);
            break;
            case Beam3ContactOctTree::spherical:
              intersection = IntersectionSPBB(bboxIDs);
            break;
            default: dserror("No or an invalid Octree type was chosen. Check your input file!");
          }

          if (intersection)
          {
            // note: creation of unique "first" entries in map, attention: IDs identical to crosslinker GIDs!!
            int mapfirst = (bboxIDs[0] + 1)*basisnodes_ + bboxIDs[1];
            contactpairmap.insert ( pair<int, vector<int> > (mapfirst, bboxIDs));
          }
        }
      }
    }
  }

  // build Pair Vector from contactpairmap
  std::map<int, std::vector<int> >::iterator it;
  int counter = 0;
  for ( it=contactpairmap.begin() ; it !=contactpairmap.end(); it++ )
  {
    counter++;
    //if(!discret_.Comm().MyPID())
      //cout << scientific << (*it).first <<"  "<< ((*it).second)[0]<<" "<< ((*it).second)[1]<<endl;
    int collid1 = searchdis_.ElementColMap()->LID(((*it).second)[0]);
    int collid2 = searchdis_.ElementColMap()->LID(((*it).second)[1]);

    DRT::Element* tempele1 = searchdis_.lColElement(collid1);
    DRT::Element* tempele2 = searchdis_.lColElement(collid2);

    // matrices to store nodal coordinates
    Epetra_SerialDenseMatrix ele1pos(3,tempele1->NumNode());
    Epetra_SerialDenseMatrix ele2pos(3,tempele2->NumNode());

    // store nodal coordinates of element 1
    for (int m=0;m<tempele1->NumNode();++m)
    {
     int tempGID = (tempele1->NodeIds())[m];
     LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
     for(int n=0;n<3;n++) ele1pos(n,m) = temppos(n);
    }

    // store nodal coordinates of element 2
    for (int m=0;m<tempele2->NumNode();++m)
    {
     int tempGID = (tempele2->NodeIds())[m];
     LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
     for(int n=0;n<3;n++) ele2pos(n,m) = temppos(n);
    }

    // add to pair vector
    contactpairs->push_back(rcp (new Beam3contact(discret_,searchdis_,dofoffset_,tempele1,tempele2,ele1pos,ele2pos)));
  }
  //if(!discret_.Comm().MyPID())
    //cout<<"number of boxes: "<<counter<<endl;

#ifdef MEASURETIME
  double isectimelocal = Teuchos::Time::wallTime() - t_search;
  double isectimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&isectimelocal, &isectimeglobal, 1);
  discret_.Comm().Barrier();
  if(!searchdis_.Comm().MyPID())
    cout << "Intersection time:\t\t" << isectimeglobal << " seconds"<<endl;
#endif

  return;
}//end of method BoundingBoxIntersection()

/*-----------------------------------------------------------------------------------*
 |  Axis Aligned Bounding Box Intersection function when both bounding boxes         |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionAABB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  /* Why have bboxlimits seperately? In certain cases, it is required to intersect hypothetical bounding boxes
   * (i.e. without an existing element) with bounding boxes of existing elements. Then, the second bounding box ID
   * is not relevant anymore since it does not exist in the octree. bboxlimits takes over the part of defining the
   * limits of the (hypothetical) bounding box.*/

  bool intersection = false;
  // translate box / element GIDs to ElementColMap()-LIDs
  // note: GID and ColumnMap LID are usually the same except for crosslinker elements from statistical mechanics
  int entry1 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int entry2 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  //Initialization....................
  double a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax;
  double b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax;

  if(periodicBC_)
  {
    int numshifts = 0;
    if(bboxlimits!=Teuchos::null)
      numshifts = (bboxlimits->M()%6)-1;
    else
      numshifts = (int)(*numshifts_)[entry2];
    // note: n shifts means n+1 segments
    for(int i=0; i<(int)(*numshifts_)[entry1]+1; i++)
    {
      for(int j=0; j<numshifts+1; j++)
      {
        //Intersection Test
        a_xmin=(*allbboxes_)[i*6][entry1];     a_xmax=(*allbboxes_)[i*6+1][entry1];
        a_ymin=(*allbboxes_)[i*6+2][entry1];   a_ymax=(*allbboxes_)[i*6+3][entry1];
        a_zmin=(*allbboxes_)[i*6+4][entry1];   a_zmax=(*allbboxes_)[i*6+5][entry1];

        if(bboxlimits!=Teuchos::null)
        {
          b_xmin=(*bboxlimits)(j*6,0);     b_xmax=(*bboxlimits)(j*6+1,0);
          b_ymin=(*bboxlimits)(j*6+2,0);   b_ymax=(*bboxlimits)(j*6+3,0);
          b_zmin=(*bboxlimits)(j*6+4,0);   b_zmax=(*bboxlimits)(j*6+5,0);
        }
        else
        {
          b_xmin=(*allbboxes_)[j*6][entry2];     b_xmax=(*allbboxes_)[j*6+1][entry2];
          b_ymin=(*allbboxes_)[j*6+2][entry2];   b_ymax=(*allbboxes_)[j*6+3][entry2];
          b_zmin=(*allbboxes_)[j*6+4][entry2];   b_zmax=(*allbboxes_)[j*6+5][entry2];
        }

        // if intersection exists, return true
        if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
        {
          intersection = true;
          break;
        }
      }
      if(intersection)
        break;
    }
  }
  else  // standard procedure without periodic boundary conditions
  {
    //Intersection Test
    a_xmin=(*allbboxes_)[0][entry1];   a_xmax=(*allbboxes_)[1][entry1];
    a_ymin=(*allbboxes_)[2][entry1];   a_ymax=(*allbboxes_)[3][entry1];
    a_zmin=(*allbboxes_)[4][entry1];   a_zmax=(*allbboxes_)[5][entry1];

    if(bboxlimits!=Teuchos::null)
    {
      b_xmin=(*bboxlimits)(0,0);   b_xmax=(*bboxlimits)(1,0);
      b_ymin=(*bboxlimits)(2,0);   b_ymax=(*bboxlimits)(3,0);
      b_zmin=(*bboxlimits)(4,0);   b_zmax=(*bboxlimits)(5,0);
    }
    else
    {
      b_xmin=(*allbboxes_)[0][entry2];   b_xmax=(*allbboxes_)[1][entry2];
      b_ymin=(*allbboxes_)[2][entry2];   b_ymax=(*allbboxes_)[3][entry2];
      b_zmin=(*allbboxes_)[4][entry2];   b_zmax=(*allbboxes_)[5][entry2];
    }
    // if intersection exists, return true
    if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
      intersection = true;
  }

  return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Cylindrical Oriented Bounding Box Intersection function when both bounding boxes |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionCOBB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  /* intersection test by calculating the distance between the two bounding box center lines
   * and comparing it to the respective diameters of the beams*/
  bool intersection = false;
  int bboxid0 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int bboxid1 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  // In case of a hypothetical BB, simply take the last beam element's diameter (does the job for now)
  double bbox1diameter = 0.0;
  if(bboxlimits!=Teuchos::null)
    bbox1diameter = (*diameter_)[diameter_->MyLength()-1];
  else
    bbox1diameter = (*diameter_)[bboxid1];

  // A heuristic value (for now). It allows us to detect contact in advance by enlarging the beam radius.
  double radiusextrusion = radialextrusion_;
  if(bboxlimits!=Teuchos::null)
    radiusextrusion = 1.0;

  if(periodicBC_)
  {
    int numshifts = 0;
    if(bboxlimits!=Teuchos::null)
      numshifts = (bboxlimits->M()%6)-1;
    else
      numshifts = (int)(*numshifts_)[bboxid1];

    for(int i=0; i<(int)(*numshifts_)[bboxid0]+1; i++)
    {
      for(int j=0; j<numshifts+1; j++)
      {
        // first points and directional vectors of the bounding boxes
        LINALG::Matrix<3,1> v;
        LINALG::Matrix<3,1> w;

        // angle between the bounding boxes
        double alpha = 0.0;
        for(int k=0; k<(int)v.M(); k++) // first BB point and direction
          v(k) = (*allbboxes_)[i*6+k+3][bboxid0]-(*allbboxes_)[i*6+k][bboxid0];
        if(bboxlimits!=Teuchos::null) // second BB point and direction
        {
          for(int k=0; k<(int)v.M(); k++)
            w(k) = (*bboxlimits)(j*6+k+3,0)-(*bboxlimits)(j*6+k,0);
        }
        else
        {
          for(int k=0; k<(int)v.M(); k++)
            w(k) = (*allbboxes_)[j*6+k+3][bboxid1]-(*allbboxes_)[j*6+k][bboxid1];
        }
        alpha = acos(v.Dot(w)/(v.Norm2()*w.Norm2()));

        // non-parallel case
        /* note: We distinguish between a parallel and a non-parallel case because of the singularity
         * in the calculation of the binormal due to the cross product in the denominator*/
        if(alpha>1e-10)
        {
          // first points of both BBs
          LINALG::Matrix<3,1> x;
          LINALG::Matrix<3,1> y;
          for(int k=0; k<(int)v.M(); k++)
            x(k) = (*allbboxes_)[i*6+k][bboxid0];
          if(bboxlimits!=Teuchos::null)
          {
            for(int k=0; k<(int)v.M(); k++)
              y(k) = (*bboxlimits)(j*6+k,0);
          }
          else
          {
            for(int k=0; k<(int)v.M(); k++)
              y(k) = (*allbboxes_)[j*6+k][bboxid1];
          }

          //note: d = abs(dot(y-x,n))
          LINALG::Matrix<3,1> yminusx = y;
          yminusx -= x;

          // binormal vector
          LINALG::Matrix<3,1> n;
          n(0) = v(1)*w(2)-v(2)*w(1);
          n(1) = v(2)*w(0)-v(0)*w(2);
          n(2) = v(0)*w(1)-v(1)*w(0);
          n.Scale(1.0/n.Norm2());

          int index0 = 1;
          int index1 = 2;
          for(int k=0; k<(int)n.M(); k++)
          {
            if(n(k)>1e-12)
              break;
            else
              index0 = (k+1)%3;
          }
          index1 = (index0+1)%3;

          // 1. distance criterium
          double d = yminusx.Dot(n);

          if (fabs(d)<=radiusextrusion*((*diameter_)[bboxid0]+bbox1diameter)/2.0)
          {
            // 2. Do the two bounding boxes actually intersect?
            double lbb0 = v.Norm2();
            double lbb1 = w.Norm2();
            v.Scale(1.0/lbb0);
            w.Scale(1.0/lbb1);
            // shifting the point on the second line by d*n facilitates the calculation of the mu and lambda (segment lengths)
            for(int k=0; k<(int)y.M(); k++)
              y(k) = y(k) - d * n(k);
            // line-wise check of the segment lengths
            double mu = (v(index1)*(y(index0)-x(index0))-v(index0)*(y(index1)-x(index1)))/(v(index0)*w(index1)-v(index1)*w(index0));
            if(mu>=0 && mu<=lbb0)
            {
              double lambda = (y(index0)-x(index0)+w(index0)*mu)/v(index0);
              if(lambda>=0 && lambda<=lbb1)
              {
                intersection = true;
                break; // leave j-loop
              }
            }
          }
        }
        else // parallel case -> d = abs(cross(v0,(x1-x0))/abs(v0)
        {
          LINALG::Matrix<3,1> x;
          LINALG::Matrix<3,1> v;
          LINALG::Matrix<3,1> yminusx;
          for(int k=0; k<(int)x.M();k++)
          {
            x(k) = (*allbboxes_)[i*6+k][bboxid0];
            v(k) = (*allbboxes_)[i*6+k+3][bboxid0]-x(k);
          }
          if(bboxlimits!=Teuchos::null)
          {
            for(int k=0; k<(int)x.M();k++)
              yminusx(k) = (*bboxlimits)(j*6+k,0)-x(k);
          }
          else
          {
            for(int k=0; k<(int)x.M();k++)
              yminusx(k) = (*allbboxes_)[j*6+k][bboxid1]-x(k);
          }

          double phi = acos(fabs(v.Dot(yminusx)/(v.Norm2()*yminusx.Norm2())));
          double d = yminusx.Norm2()*sin(phi);

          if(d<radiusextrusion*((*diameter_)[bboxid0]+bbox1diameter)/2.0)
          {
            // distance between first point of first BB and second point of second BB
            double d2 = 0.0;
            // length of first and second BB
            double l0 = v.Norm2();
            double l1 = 0.0;
            if(bboxlimits!=Teuchos::null)
            {
              for(int k=0; k<(int)x.M(); k++)
              {
                d2 += ((*bboxlimits)(j*6+k,0) - x(k))*((*bboxlimits)(j*6+k,0) - x(k));
                l1 += ((*bboxlimits)(j*6+k+3,0) - (*bboxlimits)(j*6+k,0))*((*bboxlimits)(j*6+k+3,0) - (*bboxlimits)(j*6+k,0));
              }
            }
            else
            {
              for(int k=0; k<(int)x.M(); k++)
              {
                d2 += ((*allbboxes_)[j*6+k+3][bboxid1] - x(k))*((*allbboxes_)[j*6+k+3][bboxid1] - x(k));
                l1 += ((*allbboxes_)[j*6+k+3][bboxid1] - (*allbboxes_)[j*6+k][bboxid1])*((*allbboxes_)[j*6+k+3][bboxid1] - (*allbboxes_)[j*6+k][bboxid1]);
              }
            }
            d2 = sqrt(d2);
            l1 = sqrt(l1);
            if(d2<=l0+l1)
            {
              intersection = true;
              break; // (leave j-loop)
            }
          }
        }
      }
      if(intersection)
        break; // (leave i-loop)
    }
  }
  else  // standard procedure without periodic boundary conditions
  {
    // first points and directional vectors of the bounding boxes
    LINALG::Matrix<3,1> v;
    LINALG::Matrix<3,1> w;

    // angle between the bounding boxes
    double alpha = 0.0;
    for(int k=0; k<(int)v.M(); k++) // first BB
      v(k) = (*allbboxes_)[k+3][bboxid0]-(*allbboxes_)[k][bboxid0];
    if(bboxlimits!=Teuchos::null) //second BB
    {
      for(int k=0; k<(int)v.M(); k++)
        w(k) = (*bboxlimits)(k+3,0)-(*bboxlimits)(k,0);
    }
    else
    {
      for(int k=0; k<(int)v.M(); k++)
        w(k) = (*allbboxes_)[k+3][bboxid1]-(*allbboxes_)[k][bboxid1];
    }

    alpha = acos(v.Dot(w)/(v.Norm2()*w.Norm2()));
    // non-parallel case
    if(alpha>1e-10)
    {
      LINALG::Matrix<3,1> x;
      LINALG::Matrix<3,1> y;

      for(int k=0; k<(int)v.M(); k++)
        x(k) = (*allbboxes_)[k][bboxid0];
      if(bboxlimits!=Teuchos::null)
      {
        for(int k=0; k<(int)v.M(); k++)
          y(k) = (*bboxlimits)(k,0);
      }
      else
      {
        for(int k=0; k<(int)v.M(); k++)
          y(k) = (*allbboxes_)[k][bboxid1];
      }

      //note: d = abs(dot(y-x,n))
      LINALG::Matrix<3,1> yminusx = y;
      yminusx -= x;

      // binormal vector
      LINALG::Matrix<3,1> n;
      n(0) = v(1)*w(2)-v(2)*w(1);
      n(1) = v(2)*w(0)-v(0)*w(2);
      n(2) = v(0)*w(1)-v(1)*w(0);
      n.Scale(1.0/n.Norm2());

      int index0 = 1;
      int index1 = 2;
      for(int k=0; k<(int)n.M(); k++)
      {
        if(n(k)>1e-12)
          break;
        else
          index0 = (k+1)%3;
      }
      index1 = (index0+1)%3;

      // 1. distance criterium
      double d = yminusx.Dot(n);

      if (fabs(d)<=radiusextrusion*((*diameter_)[bboxid0]+bbox1diameter)/2.0)
      {
        // 2. Do the two bounding boxes actually intersect?
        double lbb0 = v.Norm2();
        double lbb1 = w.Norm2();
        v.Scale(1.0/lbb0);
        w.Scale(1.0/lbb1);
        // shifting the point on the second line by d*n facilitates the calculation of the mu and lambda (segment lengths)
        for(int k=0; k<(int)y.M(); k++)
          y(k) = y(k) - d * n(k);
        // line-wise check of the segment lengths
        double mu = (v(index1)*(y(index0)-x(index0))-v(index0)*(y(index1)-x(index1)))/(v(index0)*w(index1)-v(index1)*w(index0));
        if(mu>=0 && mu<=lbb0)
        {
          double lambda = (y(index0)-x(index0)+w(index0)*mu)/v(index0);
          if(lambda>=0 && lambda<=lbb1)
            intersection = true;
        }
      }
    }
    else
    {
      LINALG::Matrix<3,1> x;
      LINALG::Matrix<3,1> v;
      LINALG::Matrix<3,1> yminusx;
      for(int k=0; k<(int)x.M();k++)
      {
        x(k) = (*allbboxes_)[k][bboxid0];
        v(k) = (*allbboxes_)[k+3][bboxid0]-x(k);
      }
      if(bboxlimits!=Teuchos::null)
      {
        for(int k=0; k<(int)x.M();k++)
          yminusx(k) = (*bboxlimits)(k,0)-x(k);
      }
      else
      {
        for(int k=0; k<(int)x.M();k++)
          yminusx(k) = (*allbboxes_)[k][bboxid1]-x(k);
      }

      double phi = acos(fabs(v.Dot(yminusx)/(v.Norm2()*yminusx.Norm2())));
      double d = yminusx.Norm2()*sin(phi);

      if(d<radiusextrusion*((*diameter_)[bboxid0]+bbox1diameter)/2.0)
      {
        // distance between first point of first BB and second point of second BB
        double d2 = 0.0;
        // length of first and second BB
        double l0 = v.Norm2();
        double l1 = 0.0;
        if(bboxlimits!=Teuchos::null)
        {
          for(int k=0; k<(int)x.M(); k++)
          {
            d2 += ((*bboxlimits)(k,0) - x(k))*((*bboxlimits)(k,0) - x(k));
            l1 += ((*bboxlimits)(k+3,0) - (*bboxlimits)(k,0))*((*bboxlimits)(k+3,0) - (*bboxlimits)(k,0));
          }
        }
        else
        {
          for(int k=0; k<(int)x.M(); k++)
          {
            d2 += ((*allbboxes_)[k+3][bboxid1] - x(k))*((*allbboxes_)[k+3][bboxid1] - x(k));
            l1 += ((*allbboxes_)[k+3][bboxid1] - (*allbboxes_)[k][bboxid1])*((*allbboxes_)[k+3][bboxid1] - (*allbboxes_)[k][bboxid1]);
          }
        }
        d2 = sqrt(d2);
        l1 = sqrt(l1);
        if(d2<=l0+l1)
          intersection = true;
      }
    }
  }

  return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Spherical Bounding Box Intersection function when both bounding boxes           |
 |  for linkers                                       (private)        mueller 01/12|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionSPBB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  bool intersection = false;
  int bboxid0 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int bboxid1 = searchdis_.ElementColMap()->LID(bboxIDs[1]);
  LINALG::Matrix<3,1> v;
  double radiusextrusion = 1.1;

  double bbox1diameter = 0.0;
  if(bboxlimits!=Teuchos::null)
    bbox1diameter = (*diameter_)[diameter_->MyLength()-1];
  else
    bbox1diameter = (*diameter_)[bboxid1];

  if(bboxlimits!=Teuchos::null)
  {
    for (int i=0; i<(int)v.M(); i++)
      v(i)=(*bboxlimits)(i,0)-(*allbboxes_)[i][bboxid0];
  }
  else
  {
    for (int i=0; i<(int)v.M(); i++)
      v(i)=(*allbboxes_)[i][bboxid1]-(*allbboxes_)[i][bboxid0];
  }

  double d=1e9;
  d=v.Norm2();

  if (d<radiusextrusion*((*diameter_)[bboxid0]+bbox1diameter)/2.0)
    intersection=true;

  return intersection;
}

/*-----------------------------------------------------------------------*
 | communicate Vector to all Processors                    mueller 11/11 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::CommunicateVector(Epetra_Vector& InVec, Epetra_Vector& OutVec, bool zerofy, bool doexport, bool doimport)
{
  /* zerofy InVec at the beginning of each search except for Proc 0
   * for subsequent export and reimport. This way, we guarantee redundant information
   * on all processors. */

  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_.Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::CommunicateMultiVector(Epetra_MultiVector& InVec, Epetra_MultiVector& OutVec, bool zerofy, bool doexport, bool doimport)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_.Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
  return;
}

