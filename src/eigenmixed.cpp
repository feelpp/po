#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/dh.hpp>

using namespace Feel;
using namespace Feel::vf;

class EigenProblem
:
public Simget
{
  typedef Simget super;
public:
  typedef Mesh<Simplex<3 > > mesh_type;

  typedef bases<RaviartThomas<0>, Nedelec<0,NedelecKind::NED1> > basis_type;
  typedef FunctionSpace<mesh_type, basis_type > space_type;
  typedef boost::shared_ptr<space_type> space_ptrtype;

  void run();
};


void
EigenProblem::run()
{
  if ( Environment::worldComm().isMasterRank() )
  {
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute EigenProblem\n";
  }

  auto mesh = loadMesh( new mesh_type );

  auto Xh = space_type::New( mesh );
  auto U = Xh->element();
  auto w = U.element<0>();
  auto u = U.element<1>();

  auto a = form2( _test=Xh, _trial=Xh );
  a = integrate( _range=elements(mesh),
  		 _expr=inner(idt(w),id(w))
  		 - inner(curlt(u),id(w))
  		 + inner(idt(w),id(u)) );

  auto b = form2( _test=Xh, _trial=Xh );
  b = integrate( _range=elements(mesh), _expr=inner(idt(u),id(u)) );

  auto l = form1( _test=Xh );
  a += on(_range=boundaryfaces(mesh), _rhs=l, _element=w, _expr=zero<3,1>() );

  auto modes = veigs( _formA=a, _formB=b );

  auto e = exporter(mesh);
  auto i = 0;

  for( auto const& mode : modes ) {
    auto g = mode.second.element<1>();
    auto cg = mode.second.element<0>();
    auto lambda = mode.first;

    auto n = normL2(_range=boundaryfaces(mesh),_expr=inner(idv(g),N()) );
    auto cn = normL2(_range=boundaryfaces(mesh),_expr=inner(idv(cg),N()) );

    auto d = normL2(_range=elements(mesh),_expr=divv(g) );

    auto err = normL2(_range=elements(mesh), _expr=idv(cg)-curlv(g) );
    auto errEigen = normL2(_range=elements(mesh),_expr=curlv(g)-lambda*idv(g) );
    auto errEigen2 = normL2(_range=elements(mesh), _expr=idv(g)-lambda*idv(g) );

    if ( Environment::worldComm().isMasterRank() ) {
      std::cout << "normale(u) : " << n << std::endl;
      std::cout << "normale(w) : " << cn << std::endl;
      std::cout << "divergence(u) : " << d << std::endl;
      std::cout << "e(w - curl(u)) : " << err << std::endl;
      std::cout << "e(curl(u) - l u) : " << errEigen << std::endl;
      std::cout << "e(w - l u) : " << errEigen2 << std::endl;
    }

    e->add( (boost::format("mode-%1%")%i).str(), g);
  }

  e->save();
}


int main(int argc, char **argv) 
{
  Environment env( _argc=argc, _argv=argv,
		   _desc=feel_options(),
		   _about=about(_name="po_eigenmixt",
				_author="Romain Hild",
				_email="romain.hild@plasticomnium.com") );

  Application app;
  
  app.add( new EigenProblem() );
  app.run();
}
