/*!----------------------------------------------------------------------
\file designnode.cpp
\brief A node that is part of a CAD design description

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "designnode.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignNode::DesignNode(int id, const double* coords) :
Node(id,coords)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignNode::DesignNode(const CCADISCRETIZATION::DesignNode& old) :
Node(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignNode::~DesignNode()
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignNode* CCADISCRETIZATION::DesignNode::Clone() const
{
  CCADISCRETIZATION::DesignNode* newnode = new CCADISCRETIZATION::DesignNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::DesignNode::Print(ostream& os) const
{
  Node::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::DesignNode::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  int basesize=0;
  const char* basedata = Node::Pack(basesize);
  
  size = 
  sizeint  + // holds size 
  basesize + // holds basedata
  0;         // continue to add stuff here
  
  char* data = new char[size];
  int position=0;
  
  // add size
  AddtoPack(position,data,size);    
  // add basedata
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;  

  if (position != size)
  {
    cout << "CCADISCRETIZATION::DesignNode::Pack:\n"
         << "Mismatch in size of data " << size << " <-> " << position << endl
         << __FILE__ << ":" << __LINE__ << endl;
    exit(EXIT_FAILURE);
  }
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::DesignNode::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size=0;
  ExtractfromPack(position,data,size);
  
  // extract base class
  int basesize = Size(&data[position]);
  Node::Unpack(&data[position]);
  position += basesize;
  
  if (position != size)
  {
    cout << "CCADISCRETIZATION::DesignNode::Unpack:\n"
         << "Mismatch in size of data " << size << " <-> " << position << endl
         << __FILE__ << ":" << __LINE__ << endl;
    exit(EXIT_FAILURE);
  }
  return true;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
