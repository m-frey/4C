
#include "cut_element.H"
#include "cut_volumecell.H"


void GEO::CUT::Node::RegisterCuts()
{
  if ( Position()==Point::oncutsurface )
  {
    for ( plain_edge_set::iterator i=edges_.begin(); i!=edges_.end(); ++i )
    {
      Edge * e = *i;
      point_->AddEdge( e );
    }
  }
}


void GEO::CUT::Node::AssignNodalCellSet( std::vector<plain_volumecell_set> & ele_vc_sets,
                                         std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets)
{
    for( std::vector<plain_volumecell_set>::iterator s=ele_vc_sets.begin();
         s!=ele_vc_sets.end();
         s++)
    {
        plain_volumecell_set  cell_set = *s;

        for ( plain_volumecell_set::const_iterator i=cell_set.begin(); i!=cell_set.end(); ++i )
        {
            VolumeCell * cell = *i;

            // if at least one cell of this cell_set contains the point, then the whole cell_set contains the point
            if ( cell->Contains( point() ) )
            {
                nodal_cell_sets[this].push_back( cell_set );

                // the rest of cells in this set has not to be checked for this node
                break; // finish the cell_set loop, breaks the inner for loop!
            }
        }
    }
}


void GEO::CUT::Node::FindDOFSets( bool include_inner )
{
  const plain_element_set & elements = Elements();

  std::map<Node*, plain_volumecell_set > nodal_cells;
  plain_volumecell_set cells;

  for ( plain_element_set::const_iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    {
      const plain_volumecell_set & element_cells = e->VolumeCells();
      if ( include_inner )
      {
        std::copy( element_cells.begin(), element_cells.end(), std::inserter( cells, cells.begin() ) );
      }
      else
      {
        for ( plain_volumecell_set::const_iterator i=element_cells.begin(); i!=element_cells.end(); ++i )
        {
          VolumeCell * vc = *i;
          if ( vc->Position()==Point::outside )
          {
            cells.insert( vc );
          }
        }
      }

      const std::vector<Node*> & nodes = e->Nodes();
      for ( std::vector<Node*>::const_iterator i=nodes.begin();
            i!=nodes.end();
            ++i )
      {
        Node * n = *i;
        for ( plain_volumecell_set::const_iterator i=element_cells.begin(); i!=element_cells.end(); ++i )
        {
          VolumeCell * cell = *i;

          if ( not include_inner and cell->Position()!=Point::outside )
            continue;

          if ( cell->Contains( n->point() ) )
          {
            nodal_cells[n].insert( cell );
          }
        }
      }
    }
  }

//   std::cout << "id=" << Id()
//             << " #ele=" << elements.size()
//             << " #cells=" << cells.size()
//             << " #nodal=" << nodal_cells.size()
//     ;

  // First, get the nodal cells that make up the first dofset. In most cases
  // this loop has one pass only. But if the node is cut, there will be more
  // than one set of cells that are attached to this node.

  plain_volumecell_set done;

  BuildDOFCellSets( point(), cells, nodal_cells[this], done );

  nodal_cells.erase( this );

  for ( std::map<Node*, plain_volumecell_set >::iterator i=nodal_cells.begin();
        i!=nodal_cells.end();
        ++i )
  {
    Node * n = i->first;
    plain_volumecell_set & cellset = i->second;
    BuildDOFCellSets( n->point(), cells, cellset, done );
  }

  // do any remaining internal volumes that are not connected to any node
  BuildDOFCellSets( NULL, cells, cells, done );

//   std::cout << " #dofsets: " << dofsets_.size() << " [";
//   for ( unsigned i=0; i<dofsets_.size(); ++i )
//   {
//     std::cout << " " << dofsets_[i].size();
//   }
//   std::cout << " ]\n";
}


void GEO::CUT::Node::FindDOFSetsNEW( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets,
                                     std::vector<plain_volumecell_set> & cell_sets)
{

    // finally: fill dof_cellsets_

    // do the connection between elements
    plain_volumecell_set done;
    plain_volumecell_set cells;

    // get cell_sets as a plain_volume_set
    for(std::vector<plain_volumecell_set>::iterator i=cell_sets.begin(); i!=cell_sets.end(); i++)
    {
       std::copy((*i).begin(), (*i).end(), std::inserter( cells, cells.begin() ) );
    }



    // First, get the nodal cells that make up the first dofset. In most cases
    // this loop has one pass only. But if the node is cut, there will be more
    // than one set of cells that are attached to this node.

    BuildDOFCellSets( point(), cell_sets, cells, nodal_cell_sets[this], done );

    nodal_cell_sets.erase( this );


    for ( std::map<Node*, std::vector<plain_volumecell_set> >::iterator i=nodal_cell_sets.begin();
          i!=nodal_cell_sets.end();
          ++i )
    {
        Node * n = i->first;

        std::vector<plain_volumecell_set> & cellset = i->second;
        BuildDOFCellSets( n->point(), cell_sets, cells, cellset, done );
    }

    // do any remaining internal volumes that are not connected to any node
    BuildDOFCellSets( NULL, cell_sets, cells, cell_sets, done );



}

void GEO::CUT::Node::BuildDOFCellSets( Point * p,
                                       const std::vector<plain_volumecell_set> & cell_sets,
                                       const plain_volumecell_set & cells,
                                       const std::vector<plain_volumecell_set> & nodal_cell_sets,
                                       plain_volumecell_set & done )
{

    for( std::vector<plain_volumecell_set>::const_iterator s=nodal_cell_sets.begin(); s!=nodal_cell_sets.end(); s++)
    {
        plain_volumecell_set nodal_cells = *s;

        for ( plain_volumecell_set::const_iterator i=nodal_cells.begin();
              i!=nodal_cells.end();
              ++i )
        {

            VolumeCell * cell = *i;
            if ( done.count( cell )==0 )
            {
                plain_volumecell_set connected;
                // REMARK: here use the version without! elements check:
                // here we build cell sets within one global element with vcs of subelements
                // maybe the vcs of one subelement are not connected within one subelement, but within one global element,
                // therefore more than one vc of one subelements may be connected.
                cell->Neighbors( p, cells, done, connected);

                if ( connected.size()>0 )
                {
                    std::set<plain_volumecell_set> connected_sets;

                    int count=0;
                    // find all cells of connected in cell_sets and add the corresponding cell_sets
                    for(plain_volumecell_set::iterator c=connected.begin(); c!= connected.end(); c++)
                    {
                        VolumeCell* connected_cell=*c;
                        count++;

                        int cell_it = 0;
                        for(std::vector<plain_volumecell_set>::const_iterator i=cell_sets.begin();
                            i!=cell_sets.end();
                            i++ )
                        {
                            cell_it++;

                            // contains the current cell
                            if((*i).count( connected_cell ) > 0)
                            {
                                connected_sets.insert(*i);
                            }
                        }
                    }

                    dof_cellsets_.push_back(connected_sets);
                    std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );

                }
            }
        }
    }


//	plain_volumecell_set collected_nodal_cell_sets;
//
//	for( std::vector<plain_volumecell_set>::const_iterator s=nodal_cell_sets.begin(); s!=nodal_cell_sets.end(); s++)
//	{
//		plain_volumecell_set nodal_cells = *s;
//
//		  for ( plain_volumecell_set::const_iterator i=nodal_cells.begin();
//		        i!=nodal_cells.end();
//		        ++i )
//		  {
//			  collected_nodal_cell_sets.insert(*i);
//		  }
//	}
//
//	for(plain_volumecell_set::iterator i=collected_nodal_cell_sets.begin(); i!= collected_nodal_cell_sets.end(); i++)
//	{
//	    VolumeCell * cell = *i;
//	    if ( done.count( cell )==0 )
//	    {
//	      plain_volumecell_set connected;
//	      plain_element_set elements;
//	      cell->Neighbors( p, cells, done, connected, elements );
//
//	      if ( connected.size()>0 )
//	      {
//if( Id() == 6) cout << "connected.size() " << connected.size() << endl;
//	    	  std::set<plain_volumecell_set> connected_sets;
//
//	    	  int count=0;
//	    	  // find all cells of connected in cell_sets and add the corresponding cell_sets
//	    	  for(plain_volumecell_set::iterator c=connected.begin(); c!= connected.end(); c++)
//	    	  {
//	    		  VolumeCell* connected_cell=*c;
//	    		  count++;
//
//	    		  int cell_it = 0;
//	    		  for(std::vector<plain_volumecell_set>::const_iterator i=cell_sets.begin();
//		    			i!=cell_sets.end();
//		    			i++ )
//		    	  {
//	    			  cell_it++;
//
//		    		  // contains the current cell
//		    		  if((*i).count( connected_cell ) > 0)
//		    		  {
//
//		    			  connected_sets.insert(*i);
//		    		  }
//		    	  }
//	    	  }
//
//	    	  if(Id()==6) cout << "connected_sets.size() " << connected_sets.size() << endl;
//
//	    	  if(Id()==6)  	cout<< "new dof_cellset created" << endl;
//	    	    dof_cellsets_.push_back(connected_sets);
//
//
//	        std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );
//	      }
//	    }
//	}



}

void GEO::CUT::Node::BuildDOFCellSets( Point * p,
                                       const plain_volumecell_set & cells,
                                       const plain_volumecell_set & nodal_cells,
                                       plain_volumecell_set & done )
{
  for ( plain_volumecell_set::const_iterator i=nodal_cells.begin();
        i!=nodal_cells.end();
        ++i )
  {
    VolumeCell * cell = *i;
    if ( done.count( cell )==0 )
    {
      plain_volumecell_set connected;
      plain_element_set elements;
      cell->Neighbors( p, cells, done, connected, elements );

      if ( connected.size()>0 )
      {
        dofsets_.push_back( connected );

        std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );
      }
    }
  }
}

int GEO::CUT::Node::DofSetNumber( VolumeCell * cell )
{
  int dofset = -1;
  for ( unsigned i=0; i<dofsets_.size(); ++i )
  {
    plain_volumecell_set & cells = dofsets_[i];
    if ( cells.count( cell ) > 0 )
    {
      if ( dofset==-1 )
      {
        dofset = i;
      }
      else
      {
        throw std::runtime_error( "volume dofset not unique" );
      }
    }
  }
  if ( dofset==-1 )
  {
    cout << "dofset not found for node " << this->Id() << endl;
    throw std::runtime_error( "volume dofset not found" );
  }
  return dofset;
}

int GEO::CUT::Node::DofSetNumberNEW( plain_volumecell_set & cells )
{
  int dofset = -1;

  // find the first cell of cells, this is only a volume cell of a subelement
  if(cells.size() == 0) dserror( "cells is empty");

//  VolumeCell* cell = cells[0];
  VolumeCell* cell = *(cells.begin());

  for ( unsigned i=0; i<dof_cellsets_.size(); ++i ) // loop over sets
  {
    std::set< plain_volumecell_set >& cellsets = dof_cellsets_[i];

    for(std::set<plain_volumecell_set>::iterator j = cellsets.begin();
    		j!=cellsets.end();
    		j++)
    {
    	plain_volumecell_set set = *j;


        if ( set.count( cell ) > 0 )
        {
          if ( dofset==-1 )
          {
            dofset = i;
          }
          else
          {
            cout << "node: " << Id() << endl;
            cout << "first dofset id: " << dofset << endl;
            cout << "new dofset id: " << i << endl;
            cell->Print(cout);
            throw std::runtime_error( "volume dofset not unique" );
          }
        }
    }
  }
  if ( dofset==-1 )
  {
    cout << "dofset not found for node " << this->Id() << endl;
//    throw std::runtime_error( "volume dofset not found" );
  }
  return dofset;
}


#if 0
int GEO::CUT::Node::NumDofSets( bool include_inner )
{
  if ( include_inner )
  {
    return DofSets().size();
  }
  else
  {
    int numdofsets = 0;
    for ( std::vector<plain_volumecell_set >::iterator i=dofsets_.begin();
          i!=dofsets_.end();
          ++i )
    {
      plain_volumecell_set & cells = *i;
      GEO::CUT::Point::PointPosition position = GEO::CUT::Point::undecided;
      for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        VolumeCell * c = *i;
        GEO::CUT::Point::PointPosition cp = c->Position();
        if ( cp == GEO::CUT::Point::undecided )
        {
          throw std::runtime_error( "undecided volume cell position" );
        }
        if ( position!=GEO::CUT::Point::undecided and position!=cp )
        {
          throw std::runtime_error( "mixed volume cell set" );
        }
        position = cp;
      }
      if ( position==GEO::CUT::Point::outside )
      {
        numdofsets += 1;
      }
    }
    return numdofsets;
  }
}
#endif
