/*----------------------------------------------------------------------*/
/*!
\file post_drt_monitor.cpp

\brief monitoring filter for one data

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/foerster
            089 - 289-15262
</pre>

*/
/*----------------------------------------------------------------------*/

/*!
\addtogroup Monitoring
*//*! @{ (documentation module open)*/
#ifdef CCADISCRET

#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_discret.H"

#include "post_drt_monitor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MonWriter::MonWriter(
  PostProblem& problem,
  string& infieldtype,
  int node
  )
  : myrank_(problem.comm()->MyPID()) // get my processor id
{
  // determine the owner of the node
  nodeowner_ = false;
  
  int numdis = problem.num_discr();
  //std::string fieldtype = "";
  // loop over all available discretizations
  for (int i=0; i<numdis; ++i)
  {
    PostField* field = problem.get_discretization(i);
    if (field->name()==infieldtype)
    {
      // pointer (rcp) to actual discretisation
      Teuchos::RCP< DRT::Discretization > mydiscrete = field->discretization();
      // store, if this node belongs to me
      if (mydiscrete->HaveGlobalNode(node))
      {
        nodeowner_ = mydiscrete->HaveGlobalNode(node);
      }
    }
  }// end loop over dis
    
  //ensure that we really found exactly one node owner
  {
    int localnodeowner = (int) nodeowner_;
    int numnodeowner = 0;
    (problem.comm())->SumAll(&localnodeowner,&numnodeowner,1);
    if ((myrank_==0) and (numnodeowner==0))
      dserror("Could not find node %d",node);
    if ((myrank_==0) and (numnodeowner>1))
      dserror("Found more than one owner of node %d: %d",node,numnodeowner);
  }

  return;
}


/*----------------------------------------------------------------------*/
void MonWriter::WriteMonFile(
  PostProblem& problem,
  string& infieldtype,
  int node
  )
{
  // create my output file
  std::string filename = problem.outname() + ".mon";
  std::ofstream outfile;
  if (nodeowner_)
  {
    outfile.open(filename.c_str());
  }
  //int numdis = problem.num_discr();

  // get pointer to discretisation of actual field
  PostField* field = GetFieldPtr(problem);
  if (field == NULL) dserror("Could not obtain field");

  CheckInfieldType(infieldtype);

  // pointer (rcp) to actual discretisation
  Teuchos::RCP< DRT::Discretization > mydiscrete = field->discretization();
  // space dimension of the problem
  int dim = problem.num_dim();

  // get actual results of total problem
  PostResult result = PostResult(field);

  // compute offset = datamap.MinAllGID() - field->discretization()->DofRowMap()->MinAllGID().
  // Note that datamap can only be computed in WriteResult(...), which is pure virtual on
  // this level. Hence offset is split up into two parts!
  // First part:
  const int offset1 = - field->discretization()->DofRowMap()->MinAllGID();

  // global nodal dof numbers
  std::vector<int> gdof;

  if (nodeowner_)
  {
    // test, if this node belongs to me
    bool ismynode = mydiscrete->HaveGlobalNode(node);
    if (!ismynode) // if this node does not belong to this field ( or proc, but we should be serial)
      FieldError(node);

    // pointer to my actual node
    const DRT::Node* mynode = mydiscrete->gNode(node);

    // global nodal dof numbers
    gdof = mydiscrete->Dof(mynode);
    // set some dummy values
    for(unsigned i=0;i < gdof.size();i++)
    {
      gdof[i]+=offset1;
    }

    // write header
    WriteHeader(outfile);
    outfile << node << "\n";
    outfile << "# control information: nodal coordinates   ";
    outfile << "x = " << mynode->X()[0] << "    ";
    outfile << "y = " << mynode->X()[1] << "    ";
    if (dim > 2) outfile << "z = " << mynode->X()[2];
    outfile << "\n";
    outfile << "#\n";

    WriteTableHead(outfile,dim);
  }
  else // this proc is not the node owner
  {
    // set some dummy values
    for(int i=0; i < dim+1; ++i)
    {
      gdof.push_back(-1);
    }
  }

  // this is a loop over all time steps that should be written
  // writing step size is considered
  if (nodeowner_)
  {
    while(result.next_result())
      WriteResult(outfile,result,gdof,dim);
  }

  // close file
  if (outfile.is_open())
    outfile.close();
}

/*----------------------------------------------------------------------*/
void MonWriter::WriteMonStressFile(
  PostProblem& problem,
  string& infieldtype,
  string stresstype,
  int node
  )
{
  // stop it now
  if ( (stresstype != "none") and (stresstype != "ndxyz") )
    dserror("Cannot deal with requested stress output type: %s", stresstype.c_str());
  
  // write stress
  if (stresstype != "none")
  { 
    // file name
    const std::string filename = problem.outname() + ".stress.mon";

    // define kind of stresses
    std::vector<std::string> groupnames;
    groupnames.push_back("gauss_cauchy_stresses_xyz");
    groupnames.push_back("gauss_2PK_stresses_xyz");

    // write it, now
    WriteMonStrFile(filename, problem, infieldtype, "stress", stresstype, groupnames, node);  
  }

  return;
}

/*----------------------------------------------------------------------*/
void MonWriter::WriteMonStrainFile(
  PostProblem& problem,
  string& infieldtype,
  string straintype,
  int node
  )
{
  // stop it now
  if ( (straintype != "none") and (straintype != "ndxyz") )
    dserror("Cannot deal with requested strain output type: %s", straintype.c_str());

  if (straintype != "none")
  {
    // output file name
    const std::string filename = problem.outname() + ".strain.mon";

    // define kind of strains
    std::vector<std::string> groupnames;
    groupnames.push_back("gauss_GL_strains_xyz");
    groupnames.push_back("gauss_EA_strains_xyz");

    // write, now
    WriteMonStrFile(filename, problem, infieldtype, "strain", straintype, groupnames, node);
  }

  return;
}

/*----------------------------------------------------------------------*/
void MonWriter::WriteMonStrFile(
  const string& filename,
  PostProblem& problem,
  string& infieldtype,
  const string strname,
  const string strtype,
  std::vector<std::string> groupnames,
  int node
  )
{
  // create my output file
  std::ofstream outfile;
  if (nodeowner_)
  {
    outfile.open(filename.c_str());
  }
  //int numdis = problem.num_discr();

  // get pointer to discretisation of actual field
  PostField* field = GetFieldPtr(problem);
  if (field == NULL) dserror("Could not obtain field");

  CheckInfieldType(infieldtype);

  // pointer (rcp) to actual discretisation
  Teuchos::RCP< DRT::Discretization > mydiscrete = field->discretization();
  // space dimension of the problem
  const int dim = problem.num_dim();

  // get actual results of total problem
  PostResult result = PostResult(field);

  // global nodal dof numbers
  std::vector<int> gdof;

  // compute offset = datamap.MinAllGID() - field->discretization()->DofRowMap()->MinAllGID().
  // Note that datamap can only be compute in WriteResult(...), which is pure virtual on
  // this level. Hence offset is split up into two parts!
  // First part:
  const int offset1 = - field->discretization()->DofRowMap()->MinAllGID();

  if (nodeowner_)
  {
    // test, if this node belongs to me
    bool ismynode = mydiscrete->HaveGlobalNode(node);
    if (!ismynode) // if this node does not belong to this field ( or proc, but we should be seriell)
      FieldError(node);

    // pointer to my actual node
    const DRT::Node* mynode = mydiscrete->gNode(node);

    // global nodal dof numbers
    gdof = mydiscrete->Dof(mynode);
    // set some dummy values
    for(unsigned i=0;i < gdof.size();i++)
    {
      gdof[i] += offset1;
    }

    // write header
    WriteHeader(outfile);
    outfile << node << "\n";
    outfile << "# control information: nodal coordinates   ";
    outfile << "x = " << mynode->X()[0] << "    ";
    outfile << "y = " << mynode->X()[1] << "    ";
    if (dim > 2) outfile << "z = " << mynode->X()[2];
    outfile << "\n";
    outfile << "#\n";

    WriteStrTableHead(outfile,strname,strtype,dim);
  }
  else // this proc is not the node owner
  {
    // set some dummy values
    for(int i=0; i < dim+1; ++i)
    {
      gdof.push_back(-1);
    }
  }

  // This is a loop over all possible stress or strain modes (called groupnames).
  // The call is handed to _all_ processors, because the extrapolation of the
  // stresses/strains from Gauss points to nodes is done by DRT::Discretization
  // utilising an assembly call. The assembly is parallel and thus all processors
  // have to be incoporated --- at least I think so.
  // (culpit: bborn, 07/09)
  for (std::vector<std::string>::iterator gn=groupnames.begin(); gn!=groupnames.end(); ++gn)
    WriteStrResults(outfile,problem,result,gdof,dim,strtype,*gn,node);

  if (outfile.is_open())
    outfile.close();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FieldMonWriter::GetFieldPtr(PostProblem& problem)
{
  return problem.get_discretization(0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "fluid")
    cout << "\nPure fluid problem, field option other than fluid has been ignored!\n\n";
}

/*----------------------------------------------------------------------*/
void FluidMonWriter::FieldError(int node)
{
  dserror("Node %i does not belong to fluid field!",node);
}

/*----------------------------------------------------------------------*/
void FluidMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# fluid problem, writing nodal data of node ";
}

/*----------------------------------------------------------------------*/
void FluidMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     u_x      u_y      p\n";
    break;
  case 3:
   outfile << "# step   time     u_x      u_y      u_z      p\n";
   break;
  default:
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

/*----------------------------------------------------------------------*/
void FluidMonWriter::WriteResult(
  ofstream& outfile,
  PostResult& result,
  std::vector<int>& gdof,
  int dim
  )
{
  // get actual result vector
  Teuchos::RCP< Epetra_Vector > resvec = result.read_result("velnp");
  const Epetra_BlockMap& velmap = resvec->Map();
  // do output of general time step data
  outfile << right << std::setw(10) << result.step();
  outfile << right << std::setw(16) << scientific << result.time();

  //compute second part of offset
  int offset2 = velmap.MinAllGID();

  // do output for velocity and pressure
  for(unsigned i=0; i < gdof.size(); ++i)
  {
    const int lid = velmap.LID(gdof[i]+offset2);
    outfile << right << std::setw(16) << scientific << (*resvec)[lid];
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "structure")
    cout << "\nPure structural problem, field option other than structure has been ignored!\n\n";
}

/*----------------------------------------------------------------------*/
void StructMonWriter::FieldError(int node)
{;
  dserror("Node %i does not belong to structure field!", node);
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# structure problem, writing nodal data of node ";
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "#"
            << std::right << std::setw(9) << "step"
            << std::right << std::setw(16) << "time"
            << std::right << std::setw(16) << "d_x"
            << std::right << std::setw(16) << "d_y"
            << std::right << std::setw(16) << "v_x"
            << std::right << std::setw(16) << "v_y"
            << std::right << std::setw(16) << "a_x"
            << std::right << std::setw(16) << "a_y"
            << std::endl;
    break;
  case 3:
    outfile << "#"
            << std::right << std::setw(9) << "step"
            << std::right << std::setw(16) << "time"
            << std::right << std::setw(16) << "d_x"
            << std::right << std::setw(16) << "d_y"
            << std::right << std::setw(16) << "d_z"
            << std::right << std::setw(16) << "v_x"
            << std::right << std::setw(16) << "v_y"
            << std::right << std::setw(16) << "v_z"
            << std::right << std::setw(16) << "a_x"
            << std::right << std::setw(16) << "a_y"
            << std::right << std::setw(16) << "a_z"
            << std::right << std::setw(16) << "p"
            << std::endl;
    break;
  default:
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteResult(
  ofstream& outfile,
  PostResult& result,
  std::vector<int>& gdof,
  int dim
  )
{
  // write front

  // do output of general time step data
  outfile << std::right << std::setw(10) << result.step();
  outfile << std::right << std::setw(16) << std::scientific << result.time();

  // check dimensions
  unsigned noddof = 0;
  if (dim == 2)
  {
    noddof = gdof.size();
  }
  else if (dim == 3)
  {
    if (gdof.size() == (unsigned)dim)  // ordinary case: 3 displ DOFs
      noddof = (unsigned)dim;
    else if (gdof.size() == (unsigned)dim+1)  // displacement+pressure: 3+1 DOFs
      noddof = (unsigned)dim;
    else  // eg. shell with displacement+rotation: 3+3 DOFs
      noddof = gdof.size();
  }

  // displacement

  // get actual result vector displacement
  Teuchos::RCP<Epetra_Vector> resvec = result.read_result("displacement");
  const Epetra_BlockMap& dispmap = resvec->Map();

  // compute second part of offset
  int offset2 = dispmap.MinAllGID();

  // do output of displacement
  for(unsigned i=0; i < noddof; ++i)
  {
    const int lid = dispmap.LID(gdof[i]+offset2);
    if (lid == -1) dserror("illegal gid %d at %d!",gdof[i],i);
    outfile << std::right << std::setw(16) << std::scientific << (*resvec)[lid];
  }

  // velocity

  // get actual result vector velocity
  resvec = result.read_result("velocity");
  const Epetra_BlockMap& velmap = resvec->Map();

  // compute second part of offset
  offset2 = velmap.MinAllGID();

  // do output of velocity
  for(unsigned i=0; i <  noddof; ++i)
  {
    const int lid = velmap.LID(gdof[i]+offset2);
    if (lid == -1) dserror("illegal gid %d at %d!",gdof[i],i);
    outfile << std::right << std::setw(16) << std::scientific << (*resvec)[lid];
  }

  // acceleration

  // get actual result vector acceleration
  resvec = result.read_result("acceleration");
  const Epetra_BlockMap& accmap = resvec->Map();

  //compute second part of offset
  offset2 = accmap.MinAllGID();

  // do output for acceleration
  for(unsigned i=0; i <  noddof; ++i)
  {
    const int lid = accmap.LID(gdof[i]+offset2);
    if (lid == -1) dserror("illegal gid %d at %d!",gdof[i],i);
    outfile << std::right << std::setw(16) << std::scientific << (*resvec)[lid];
  }

  // pressure
  if (gdof.size() == (unsigned)dim+1)
  {
    // get actual result vector displacement/pressure
    resvec = result.read_result("displacement");
    const Epetra_BlockMap& pressmap = resvec->Map();

    // compute second part of offset
    offset2 = pressmap.MinAllGID();

    // do output of pressure
    {
      const unsigned i = (unsigned)dim;
      const int lid = pressmap.LID(gdof[i]+offset2);
      if (lid == -1) dserror("illegal gid %d at %d!",gdof[i],i);
      outfile << std::right << std::setw(16) << std::scientific << (*resvec)[lid];
    }
  }

  outfile << "\n";
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteStrTableHead(
  ofstream& outfile,
  const string strname,
  const string strtype,
  const int dim
  )
{
  switch (dim)
  {
  case 2:
    outfile << "#"
            << std::right << std::setw(9) << "step"
            << std::right << std::setw(16) << "time"
            << std::right << std::setw(16) << strname+"_xx"
            << std::right << std::setw(16) << strname+"_yy"
            << std::right << std::setw(16) << strname+"_xy"
            << std::endl;
    break;
  case 3:
    outfile << "#"
            << std::right << std::setw(9) << "step"
            << std::right << std::setw(16) << "time"
            << std::right << std::setw(16) << strname+"_xx"
            << std::right << std::setw(16) << strname+"_yy"
            << std::right << std::setw(16) << strname+"_zz"
            << std::right << std::setw(16) << strname+"_xy"
            << std::right << std::setw(16) << strname+"_yz"
            << std::right << std::setw(16) << strname+"_zx"
            << std::endl;
    break;
  default:
    dserror("Number of dimensions in space differs from 2 and 3!");
  }

  return;
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteStrResults(
  ofstream& outfile,
  PostProblem& problem,
  PostResult& result,
  std::vector<int>& gdof,
  int dim,
  string strtype,
  string groupname,
  const int node
  )
{
  result.next_result();  // needed
  if (map_has_map(result.group(), groupname.c_str()))
  {
    // strings
    std::string name;
    std::string out;
    if (groupname == "gauss_2PK_stresses_xyz")
    {
      name = "nodal_2PK_stresses_xyz";
      out = "2nd Piola-Kirchhoff stresses";
    }
    else if (groupname == "gauss_cauchy_stresses_xyz")
    {
      name = "nodal_cauchy_stresses_xyz";
      out = "Cauchy stresses";
    }
    else if (groupname == "gauss_GL_strains_xyz")
    {
      name = "nodal_GL_strains_xyz";
      out = "Green-Lagrange strains";
    }
    else if (groupname == "gauss_EA_strains_xyz")
    {
      name = "nodal_EA_strains_xyz";
      out = "Euler-Almansi strains";
    }
    else
    {
      dserror("trying to write something that is not a stress or a strain");
      exit(1);
    }

    // get pointer to discretisation of actual field
    PostField* field = GetFieldPtr(problem);
    // base file name
    const std::string basename = problem.outname();

    // open file
    const std::string filename = basename + "_"+ field->name() + "."+ name;
    cout << "reading from " << filename << endl;
    std::ofstream file;
    int startfilepos = 0;
    if (myrank_ == 0)
    {
      file.open(filename.c_str());
      startfilepos = file.tellp(); // file position should be zero, but we stay flexible
    }

    if (myrank_ == 0)
      cout << "writing node-based " << out << endl;

    // DOFs at node
    int numdf = 0;
    if (dim == 3)
      numdf = 6;
    else if (dim == 2)
      numdf = 3;
    else
      dserror("Cannot handle dimension %d", dim);

    // this is a loop over all time steps that should be written
    // bottom control here, because first set has been read already
    do {
      WriteStrResult(outfile,field,result,groupname,name,numdf,node);
    } while (result.next_result());

    // close result file
    if (file.is_open())
      file.close();
  }

  return;
}

/*----------------------------------------------------------------------*/
void StructMonWriter::WriteStrResult(
  ofstream& outfile,
  PostField*& field,
  PostResult& result,
  const string groupname,
  const string name,
  const int numdf,
  const int node
  ) const
{
  // get stresses/strains at Gauss points
  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > data
    = result.read_result_serialdensematrix(groupname);
  // discretisation (once more)
  const Teuchos::RCP<DRT::Discretization> dis = field->discretization();

  // extrapolate stresses/strains to nodes
  // and assemble them in two global vectors
  Teuchos::ParameterList p;
  p.set("action", "postprocess_stress");
  p.set("stresstype", "ndxyz");
  p.set("gpstressmap", data);
  Teuchos::RCP<Epetra_Vector> normal_stresses = Teuchos::rcp(new Epetra_Vector(*(dis->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> shear_stresses = Teuchos::rcp(new Epetra_Vector(*(dis->DofRowMap())));
  dis->Evaluate(p,Teuchos::null,Teuchos::null,normal_stresses,shear_stresses,Teuchos::null);

  // average stresses/strains and print to file
  if (nodeowner_)
  {
    const DRT::Node* lnode = dis->gNode(node);
    const std::vector<int> lnodedofs = dis->Dof(lnode);
    const int adjele = lnode->NumElement();
    std::vector<double> nodal_stresses;
    if (numdf == 6)
    {
      if (lnodedofs.size() < 3)
        dserror("Too few DOFs at node of interest");
      nodal_stresses.push_back((*normal_stresses)[lnodedofs[0]]/adjele);
      nodal_stresses.push_back((*normal_stresses)[lnodedofs[1]]/adjele);
      nodal_stresses.push_back((*normal_stresses)[lnodedofs[2]]/adjele);
      nodal_stresses.push_back((*shear_stresses)[lnodedofs[0]]/adjele);
      nodal_stresses.push_back((*shear_stresses)[lnodedofs[1]]/adjele);
      nodal_stresses.push_back((*shear_stresses)[lnodedofs[2]]/adjele);      
    }
    else if (numdf == 3)
    {
      if (lnodedofs.size() < 2)
        dserror("Too few DOFs at node of interest");
      nodal_stresses.push_back((*normal_stresses)[lnodedofs[0]]/adjele);
      nodal_stresses.push_back((*normal_stresses)[lnodedofs[1]]/adjele);
      nodal_stresses.push_back((*shear_stresses)[lnodedofs[0]]/adjele);
    }
    else
    {
      dserror("Don't know what to do with %d DOFs per node", numdf);
    }
    
    // print to file
    outfile << std::right << std::setw(10) << result.step();
    outfile << std::right << std::setw(16) << std::scientific << result.time();
    for (std::vector<double>::iterator ns=nodal_stresses.begin(); ns!=nodal_stresses.end(); ++ns)
      outfile << std::right << std::setw(16) << std::scientific << *ns;
    outfile << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "ale")
    cout << "\nPure ALE problem, field option other than ale has been ignored!\n\n";
}

/*----------------------------------------------------------------------*/
void AleMonWriter::FieldError(int node)
{
  dserror("Node %i does not belong to ALE field!", node);
}

/*----------------------------------------------------------------------*/
void AleMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# ALE problem, writing nodal data of node ";
}

/*----------------------------------------------------------------------*/
void AleMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     d_x      d_y\n";
    break;
  case 3:
   outfile << "# step   time     d_x      d_y      d_z\n";
   break;
  default:
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

/*----------------------------------------------------------------------*/
void AleMonWriter::WriteResult(
  ofstream& outfile,
  PostResult& result,
  std::vector<int>& gdof,
  int dim
  )
{
  // get actual result vector for displacement
  Teuchos::RCP< Epetra_Vector > resvec = result.read_result("displacement");
  const Epetra_BlockMap& dispmap = resvec->Map();
  // do output of general time step data
  outfile << right << std::setw(10) << result.step();
  outfile << right << std::setw(16) << scientific << result.time();

  //compute second part of offset
  int offset2 = dispmap.MinAllGID();

  // do output for velocity and pressure
  for(unsigned i=0; i < gdof.size()-1; ++i)
  {
    const int lid = dispmap.LID(gdof[i]+offset2);
    outfile << right << std::setw(16) << scientific << (*resvec)[lid];
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiFluidMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(1);
  if (myfield->name() != "fluid")
    dserror("Fieldtype of field 1 is not fluid.");
  return myfield;
}

/*----------------------------------------------------------------------*/
void FsiFluidMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of fluid node ";
}

/*----------------------------------------------------------------------*/
void FsiFluidMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     d_x      d_y      u_x      u_y      p\n";
    break;
  case 3:
   outfile << "# step   time     d_x      d_y      d_z     u_x      u_y      u_z      p\n";
   break;
  default:
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

/*----------------------------------------------------------------------*/
void FsiFluidMonWriter::WriteResult(
  ofstream& outfile,
  PostResult& result,
  std::vector<int>& gdof,
  int dim
  )
{
  // get actual result vector for displacement
  Teuchos::RCP< Epetra_Vector > resvec = result.read_result("dispnp");
  const Epetra_BlockMap& dispmap = resvec->Map();
  // do output of general time step data
  outfile << right << std::setw(10) << result.step();
  outfile << right << std::setw(16) << scientific << result.time();

  //compute second part of offset
  int offset2 = dispmap.MinAllGID();

  for(unsigned i=0; i < gdof.size()-1; ++i)
  {
    const int lid = dispmap.LID(gdof[i]+offset2);
    outfile << right << std::setw(16) << scientific << (*resvec)[lid];
  }


  // get actual result vector for velocity
  resvec = result.read_result("velnp");
  const Epetra_BlockMap& velmap = resvec->Map();

  //compute second part of offset
  offset2 = velmap.MinAllGID();

  // do output for velocity and pressure
  for(unsigned i=0; i < gdof.size(); ++i)
  {
    const int lid = velmap.LID(gdof[i]+offset2);
    outfile << right << std::setw(16) << scientific << (*resvec)[lid];
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiStructMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(0);
  if (myfield->name() != "structure")
    dserror("Fieldtype of field 1 is not structure.");
  return myfield;
}

/*----------------------------------------------------------------------*/
void FsiStructMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of structure node ";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiAleMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(1);
  if (myfield->name() != "fluid")
    dserror("Fieldtype of field 1 is not fluid.");
  return myfield;
}

/*----------------------------------------------------------------------*/
void FsiAleMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of ALE node ";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*!
 * \brief filter main routine for monitoring filter
 *
 * Write ASCII file of one nodes history
 *
 * Note: Works in seriell version only! Requires to read one instance of the discretisation!!
 *
 * \author chfoe
 * \date 11/07
 */
int main(int argc, char** argv)
{
  // command line processor to deal with arguments
  Teuchos::CommandLineProcessor my_comlinproc;
  my_comlinproc.setDocString("Post DRT monitoring filter\n\nwrite nodal result data of specified node into outfile.mon");

  bool required = true;
  /* Set an additional integer command line option which is the global node Id
     of the node you're interested in */
  int node = 0;
  my_comlinproc.setOption("node", &node, "Global node number",required);
  /* Set a std::string command line option */
  std::string infieldtype = "fluid";
  my_comlinproc.setOption("field", &infieldtype, "Field to which output node belongs (fluid, structure, ale)");

  // my post processing problem itself
  PostProblem problem(my_comlinproc,argc,argv);


  switch (problem.Problemtype())
  {
    case prb_fsi:
    {
      if(infieldtype == "fluid")
      {
        FsiFluidMonWriter mymonwriter(problem,infieldtype,node);
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else if(infieldtype == "structure")
      {
        FsiStructMonWriter mymonwriter(problem,infieldtype,node);
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else if(infieldtype == "ale")
      {
        dserror("There is no ALE output. Displacements of fluid nodes can be printed.");
        FsiAleMonWriter mymonwriter(problem,infieldtype,node);
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else
      {
        dserror("handling for monitoring of this fieldtype not yet implemented");
      }
      break;
    }
    case prb_structure:
    {
      StructMonWriter mymonwriter(problem,infieldtype,node);
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      mymonwriter.WriteMonStressFile(problem,infieldtype,problem.stresstype(),node);
      mymonwriter.WriteMonStrainFile(problem,infieldtype,problem.straintype(),node);
      break;
    }
    case prb_fluid:
    {
      FluidMonWriter mymonwriter(problem,infieldtype,node);
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      break;
    }
    case prb_ale:
    {
      AleMonWriter mymonwriter(problem,infieldtype,node);
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      break;
    }
    default:
    {
      dserror("problem type %d not yet supported", problem.Problemtype());
    }
  }

  return 0;
}
#endif
/*! @} (documentation module close)*/
