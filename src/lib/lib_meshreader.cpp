/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality for reading nodes

\level 0

*/
/*---------------------------------------------------------------------*/

#include "lib_discret.H"
#include "lib_meshreader.H"
#include "lib_domainreader.H"
#include "lib_elementreader.H"
#include "lib_globalproblem.H"
#include "lib_inputreader.H"
#include "nurbs_discret_control_point.H"
#include "immersed_problem_immersed_node.H"
#include "fiber_node.H"
#include "io_pstream.H"

#include <string>

namespace DRT::INPUT
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  MeshReader::MeshReader(const DRT::INPUT::DatFileReader& reader, std::string sectionname)
      : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
  {
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
      const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
      const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
      const std::string* geofilepath)
  {
    switch (geometrysource)
    {
      case INPAR::geometry_full:
      {
        std::string fullsectionname("--" + sectionname + " ELEMENTS");
        Teuchos::RCP<ElementReader> er =
            Teuchos::rcp(new DRT::INPUT::ElementReader(dis, reader, fullsectionname, elementtypes));
        element_readers_.push_back(er);
        break;
      }
      case INPAR::geometry_box:
      {
        std::string fullsectionname("--" + sectionname + " DOMAIN");
        Teuchos::RCP<DomainReader> dr =
            Teuchos::rcp(new DRT::INPUT::DomainReader(dis, reader, fullsectionname));
        domain_readers_.push_back(dr);
        break;
      }
      case INPAR::geometry_file:
      {
        dserror("Unfortunately not yet implemented, but feel free ...");
        break;
      }
      default:
        dserror("Unknown geometry source");
        break;
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
      const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
      const INPAR::GeometryType geometrysource, const std::string* geofilepath)
  {
    std::set<std::string> dummy;
    AddAdvancedReader(dis, reader, sectionname, dummy, geometrysource, geofilepath);
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<DRT::Discretization>> MeshReader::FindDisNode(int global_node_id)
  {
    std::vector<Teuchos::RCP<DRT::Discretization>> list_of_discretizations;
    for (const auto& element_reader : element_readers_)
      if (element_reader->HasNode(global_node_id))
        list_of_discretizations.emplace_back(element_reader->MyDis());

    return list_of_discretizations;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::ReadAndPartition()
  {
    // We need to track the max global node ID to offset node numbering and for sanity checks
    int max_node_id = 0;

    ReadMeshFromDatFile(max_node_id);
    CreateInlineMesh(max_node_id);

    // last check if there are enough nodes
    node_reader->ThrowIfNotEnoughNodes(max_node_id);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::ReadMeshFromDatFile(int& max_node_id)
  {
    TEUCHOS_FUNC_TIME_MONITOR("MeshReader::ReadMeshFromDatFile");

    // read and partition element information
    for (const auto& element_reader : element_readers_) element_reader->ReadAndPartition();

    // TODO: This should go into the node reader
    {
      // Check if there are any nodes to be read. If not, leave right away.
      const int numnodes = reader_.ExcludedSectionLength(sectionname_);

      if (numnodes == 0) return;
      const int myrank = comm_->MyPID();

      // We will read the nodes block wise. We will use one block per processor
      // so the number of blocks is numproc
      // OR number of blocks is numnodes if less nodes than procs are read in
      // determine a rough blocksize
      int number_of_blocks = std::min(comm_->NumProc(), numnodes);
      int blocksize = std::max(numnodes / number_of_blocks, 1);

      // an upper limit for blocksize
      const int maxblocksize = 200000;

      if (blocksize > maxblocksize)
      {
        // without an additional increase of number_of_blocks by 1 the last block size
        // could reach a maximum value of (2*maxblocksize)-1, potentially
        // violating the intended upper limit!
        number_of_blocks = 1 + static_cast<int>(numnodes / maxblocksize);
        blocksize = maxblocksize;
      }

      // open input file at the right position
      // note that stream is valid on proc 0 only!
      const std::string inputfile_name = reader_.MyInputfileName();
      std::ifstream file;
      if (myrank == 0)
      {
        file.open(inputfile_name.c_str());
        file.seekg(reader_.ExcludedSectionPosition(sectionname_));
      }
      std::string tmp;
      std::string tmp2;

      // note that the last block is special....
      int filecount = 0;
      for (int block = 0; block < number_of_blocks; ++block)
      {
        if (myrank == 0)
        {
          int block_counter = 0;
          for (; file; ++filecount)
          {
            file >> tmp;

            if (tmp == "NODE")
            {
              double coords[6];
              int nodeid;
              // read in the node coordinates
              file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
              // store current position of file reader
              int length = file.tellg();
              file >> tmp2;
              nodeid--;
              max_node_id = std::max(max_node_id, nodeid) + 1;
              std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(nodeid);
              if (tmp2 != "ROTANGLE")  // Common (Boltzmann) Nodes with 3 DoFs
              {
                // go back with file reader in order to make the expression in tmp2 available to
                // tmp in the next iteration step
                file.seekg(length);
                for (unsigned i = 0; i < diss.size(); ++i)
                {
                  // create node and add to discretization
                  Teuchos::RCP<DRT::Node> node =
                      Teuchos::rcp(new DRT::Node(nodeid, coords, myrank));
                  diss[i]->AddNode(node);
                }
              }
              else  // Cosserat Nodes with 6 DoFs
              {
                // read in the node ritations in case of cosserat nodes
                file >> coords[3] >> coords[4] >> coords[5];
                for (unsigned i = 0; i < diss.size(); ++i)
                {
                  // create node and add to discretization
                  Teuchos::RCP<DRT::Node> node =
                      Teuchos::rcp(new DRT::Node(nodeid, coords, myrank, true));
                  diss[i]->AddNode(node);
                }
              }
              ++block_counter;
              if (block != number_of_blocks - 1)  // last block takes all the rest
                if (block_counter == blocksize)   // block is full
                {
                  ++filecount;
                  break;
                }
            }
            // this is a specialized node for immersed problems
            else if (tmp == "INODE")
            {
              double coords[6];
              int nodeid;
              // read in the node coordinates
              file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
              // store current position of file reader
              int length = file.tellg();
              file >> tmp2;
              nodeid--;
              max_node_id = std::max(max_node_id, nodeid) + 1;
              std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(nodeid);

              if (tmp2 != "ROTANGLE")  // Common (Boltzmann) Nodes with 3 DoFs
              {
                // go back with file reader in order to make the expression in tmp2 available to
                // tmp in the next iteration step
                file.seekg(length);
                for (unsigned i = 0; i < diss.size(); ++i)
                {
                  // create node and add to discretization
                  Teuchos::RCP<DRT::Node> node =
                      Teuchos::rcp(new IMMERSED::ImmersedNode(nodeid, coords, myrank));
                  diss[i]->AddNode(node);
                }
              }
              else  // Cosserat Nodes with 6 DoFs
              {
                dserror("no valid immersed node definition");
              }
              ++block_counter;
              if (block != number_of_blocks - 1)  // last block takes all the rest
                if (block_counter == blocksize)   // block is full
                {
                  ++filecount;
                  break;
                }
            }
            // this node is a Nurbs control point
            else if (tmp == "CP")
            {
              // read control points for isogeometric analysis (Nurbs)
              double coords[3];
              double weight;

              int cpid;
              file >> cpid >> tmp >> coords[0] >> coords[1] >> coords[2] >> weight;
              cpid--;
              max_node_id = std::max(max_node_id, cpid) + 1;
              if (cpid != filecount)
                dserror("Reading of control points failed: They must be numbered consecutive!!");
              if (tmp != "COORD") dserror("failed to read control point %d", cpid);
              std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(cpid);
              for (unsigned i = 0; i < diss.size(); ++i)
              {
                Teuchos::RCP<DRT::Discretization> dis = diss[i];
                // create node/control point and add to discretization
                Teuchos::RCP<DRT::NURBS::ControlPoint> node =
                    Teuchos::rcp(new DRT::NURBS::ControlPoint(cpid, coords, weight, myrank));
                dis->AddNode(node);
              }
              ++block_counter;
              if (block != number_of_blocks - 1)  // last block takes all the rest
                if (block_counter == blocksize)   // block is full
                {
                  ++filecount;
                  break;
                }
            }
            // this is a special node with additional fiber information
            else if (tmp == "FNODE")
            {
              enum class FiberType
              {
                Unknown,
                Angle,
                Fiber,
                CosyDirection
              };

              // read fiber node
              std::array<double, 3> coords = {0.0, 0.0, 0.0};
              std::map<FIBER::CoordinateSystemDirection, std::array<double, 3>> cosyDirections;
              std::vector<std::array<double, 3>> fibers;
              std::map<FIBER::AngleType, double> angles;

              int nodeid;
              // read in the node coordinates and fiber direction
              file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
              nodeid--;
              max_node_id = std::max(max_node_id, nodeid) + 1;

              while (true)
              {
                // store current position of file reader
                std::ifstream::pos_type length = file.tellg();
                // try to read new fiber direction or coordinate system
                file >> tmp2;

                DRT::FIBER::CoordinateSystemDirection coordinateSystemDirection;
                DRT::FIBER::AngleType angleType;
                FiberType type = FiberType::Unknown;

                if (tmp2 == "FIBER" + std::to_string(1 + fibers.size()))
                {
                  type = FiberType::Fiber;
                }
                else if (tmp2 == "CIR")
                {
                  coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Circular;
                  type = FiberType::CosyDirection;
                }
                else if (tmp2 == "TAN")
                {
                  coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Tangential;
                  type = FiberType::CosyDirection;
                }
                else if (tmp2 == "RAD")
                {
                  coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Radial;
                  type = FiberType::CosyDirection;
                }
                else if (tmp2 == "HELIX")
                {
                  angleType = DRT::FIBER::AngleType::Helix;
                  type = FiberType::Angle;
                }
                else if (tmp2 == "TRANS")
                {
                  angleType = DRT::FIBER::AngleType::Transverse;
                  type = FiberType::Angle;
                }
                else
                {
                  // No more fiber information. Jump to last position.
                  file.seekg(length);
                  break;
                }

                // add fiber / angle to the map
                switch (type)
                {
                  case FiberType::Unknown:
                  {
                    dserror(
                        "Unknown fiber node attribute. Numbered fibers must be in order, i.e. "
                        "FIBER1, FIBER2, ...");
                  }
                  case FiberType::Angle:
                  {
                    file >> angles[angleType];
                    break;
                  }
                  case FiberType::Fiber:
                  {
                    std::array<double, 3> fiber_components;
                    file >> fiber_components[0] >> fiber_components[1] >> fiber_components[2];
                    fibers.emplace_back(fiber_components);
                    break;
                  }
                  case FiberType::CosyDirection:
                  {
                    file >> cosyDirections[coordinateSystemDirection][0] >>
                        cosyDirections[coordinateSystemDirection][1] >>
                        cosyDirections[coordinateSystemDirection][2];
                    break;
                  }
                  default:
                    dserror("Unknown number of components");
                }
              }

              // add fiber information to node
              std::vector<Teuchos::RCP<DRT::Discretization>> discretizations = FindDisNode(nodeid);
              for (auto& dis : discretizations)
              {
                auto node = Teuchos::rcp(new DRT::FIBER::FiberNode(
                    nodeid, coords, cosyDirections, fibers, angles, myrank));
                dis->AddNode(node);
              }

              ++block_counter;
              if (block != number_of_blocks - 1)  // last block takes all the rest
              {
                if (block_counter == blocksize)  // block is full
                {
                  ++filecount;
                  break;
                }
              }
            }
            else if (tmp.find("--") == 0)
              break;
            else
              dserror("unexpected word '%s'", tmp.c_str());
          }
        }

        // export block of nodes to other processors as reflected in rownodes, changes ownership of
        // nodes
        for (const auto& element_reader : element_readers_)
          element_reader->MyDis()->ProcZeroDistributeNodesToAll(*element_reader->MyRowNodes());
      }
    }

    // last thing to do here is to produce nodal ghosting/overlap
    for (const auto& element_reader : element_readers_)
    {
      element_reader->MyDis()->ExportColumnNodes(*element_reader->MyColNodes());
      element_reader->Complete();
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::CreateInlineMesh(int& max_node_id)
  {
    for (const auto& domain_reader : domain_readers_)
    {
      // communicate node offset to all procs
      int local_max_node_id = max_node_id;
      comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

      domain_reader->CreatePartitionedMesh(max_node_id);
      domain_reader->Complete();
      max_node_id = domain_reader->MyDis()->NodeRowMap()->MaxAllGID() + 1;
    }
  }
}  // namespace DRT::INPUT
