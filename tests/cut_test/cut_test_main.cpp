
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_element.H"
#include "cut_test_utils.H"

//#include <boost/program_options.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <map>
#include <string>
#include <sstream>

void test_hex8_simple();
void test_tet4_simple();
void test_pyramid5_simple();
void test_wedge6_simple();
void test_hex8_fullside();
void test_hex8_diagonal();
void test_hex8_tet4();
void test_hex8_hex8();
void test_hex8_touch();
void test_hex8_touch2();
void test_hex8_schraeg();
void test_hex8_tet4_touch();
void test_hex8_tet4_touch2();
void test_hex8_mesh();
void test_hex8_double();
void test_hex8_multiple();
void test_hex8_bad1();
void test_hex8_bad2();
void test_hex8_bad3();
void test_hex8_bad4();
void test_hex8_wedge6();
void test_hex8_quad4_touch();
void test_hex8_quad4_touch2();
void test_hex8_quad4_touch3();
void test_hex8_quad4_cut();
void test_hex8_quad4_gedreht();
void test_hex8_hex8_durchstoss();
void test_hex8_hex8_onside();
void test_hex8_quad4_schnitt();
void test_hex8_quad4_touch4();
void test_hex8_quad4_touch5();
void test_hex8_quad4_touch6();
void test_hex8_quad4_touch7();
void test_quad4_quad4_simple();
void test_hex8_quad4_mesh();

void test_hex27_quad9_simple();
void test_hex20_quad9_simple();
void test_hex20_quad9_moved();
void test_tet10_quad9_simple();
void test_tet10_quad9_moved();

int main( int argc, char ** argv )
{
  typedef void ( *testfunct )();

  std::map<std::string, testfunct> functable;

  functable["hex8_simple"] = test_hex8_simple;
  functable["tet4_simple"] = test_tet4_simple;
  functable["pyramid5_simple"] = test_pyramid5_simple;
  functable["wedge6_simple"] = test_wedge6_simple;
  functable["hex8_diagonal"] = test_hex8_diagonal;
  functable["hex8_fullside"] = test_hex8_fullside;
  functable["hex8_hex8"] = test_hex8_hex8;
  functable["hex8_tet4"] = test_hex8_tet4;
  functable["hex8_touch"] = test_hex8_touch;
  functable["hex8_touch2"] = test_hex8_touch2;
  functable["hex8_schraeg"] = test_hex8_schraeg;
  functable["hex8_tet4_touch"] = test_hex8_tet4_touch;
  functable["hex8_tet4_touch2"] = test_hex8_tet4_touch2;
  functable["hex8_mesh"] = test_hex8_mesh;
  functable["hex8_double"] = test_hex8_double;
  functable["hex8_multiple"] = test_hex8_multiple;
  functable["hex8_bad1"] = test_hex8_bad1;
  functable["hex8_bad2"] = test_hex8_bad2;
  functable["hex8_bad3"] = test_hex8_bad3;
  functable["hex8_bad4"] = test_hex8_bad4;
  functable["hex8_wedge6"] = test_hex8_wedge6;
  functable["hex8_quad4_touch"] = test_hex8_quad4_touch;
  functable["hex8_quad4_touch2"] = test_hex8_quad4_touch2;
  functable["hex8_quad4_touch3"] = test_hex8_quad4_touch3;
  functable["hex8_quad4_cut"] = test_hex8_quad4_cut;
  functable["hex8_quad4_gedreht"] = test_hex8_quad4_gedreht;
  functable["hex8_hex8_durchstoss"] = test_hex8_hex8_durchstoss;
  functable["hex8_hex8_onside"] = test_hex8_hex8_onside;
  //functable["hex8_quad4_schnitt"] = test_hex8_quad4_schnitt;
  functable["hex8_quad4_touch4"] = test_hex8_quad4_touch4;
  functable["hex8_quad4_touch5"] = test_hex8_quad4_touch5;
  functable["hex8_quad4_touch6"] = test_hex8_quad4_touch6;
  //functable["hex8_quad4_touch7"] = test_hex8_quad4_touch7;
  functable["hex8_quad4_mesh"] = test_hex8_quad4_mesh;

  functable["hex27_quad9_simple"] = test_hex27_quad9_simple;
  functable["hex20_quad9_simple"] = test_hex20_quad9_simple;
  functable["hex20_quad9_moved"] = test_hex20_quad9_moved;
  functable["tet10_quad9_simple"] = test_tet10_quad9_simple;
  functable["tet10_quad9_moved"] = test_tet10_quad9_moved;

  Teuchos::CommandLineProcessor clp( false );

  std::string indent = "\t\t\t\t\t";
  std::stringstream doc;
  doc << "Available tests:\n"
      << indent << "(all)\n";
  for ( std::map<std::string, testfunct>::iterator i=functable.begin(); i!=functable.end(); ++i )
  {
    const std::string & name = i->first;
    doc << indent << name << "\n";
  }

  std::string testname = "(all)";
  clp.setOption( "test", &testname, doc.str().c_str() );

  switch ( clp.parse( argc, argv ) )
  {
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    break;
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
    return 0;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    std::cerr << argv[0] << ": unrecognized option\n";
    return 1;
  }

  if ( testname == "(all)" )
  {
    for ( std::map<std::string, testfunct>::iterator i=functable.begin(); i!=functable.end(); ++i )
    {
      ( *i->second )();
    }
  }
  else
  {
    std::map<std::string, testfunct>::iterator i = functable.find( testname );
    if ( i==functable.end() )
    {
      std::cerr << argv[0] << ": test '" << testname << "' not found\n";
      return 1;
    }
    else
    {
      ( *i->second )();
    }
  }

  return 0;
}
