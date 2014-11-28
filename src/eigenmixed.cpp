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

  typedef bases<RaviartThomas<0>, Nedelec<0,NedelecKind::NED1>, Lagrange<2,Scalar > > basis_type;
  typedef FunctionSpace<mesh_type, basis_type > space_type;

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

  auto epsilon = doption("parameters.epsilon");
  auto gamma = doption("parameters.gamma");

  auto mesh = loadMesh( new mesh_type );

  auto Xh = space_type::New( mesh );

  auto U = Xh->element();

#ifdef PO_UHCURL
  // if u is in h(curl) and w in h(div)
  auto u = U.element<1>();
  auto w = U.element<0>();
  auto p = U.element<2>();

  auto a = form2( _test=Xh, _trial=Xh );
  a = integrate( elements(mesh),
		 - trans(curlt(u))*id(w) - trans(idt(u))*trans(grad(p))
  		 + inner(idt(w),id(w)) + trans(idt(w))*curl(u)
		 + gradt(p)*id(u)
		 + epsilon*inner(idt(u),id(u))
		 + gamma*inner(idt(p),id(p))
  		 );

  a += integrate( boundaryfaces( mesh ),
       		        trans(idt(u))*cross(id(w),N()) );

#else
  // if u is in h(div) and w in h(curl)
  auto u = U.element<0>();
  auto w = U.element<1>();
  auto p = U.element<2>();

  auto a = form2( _test=Xh, _trial=Xh );
  a = integrate( elements(mesh),
		 - trans(idt(u))*curl(w) + divt(u)*id(p)
  		 + inner(idt(w),id(w)) + trans(curlt(w))*id(u)
		 - idt(p)*div(u)
		 + epsilon*inner(idt(u),id(u))
		 + gamma*inner(idt(p),id(p))
  		 );

  a += integrate( boundaryfaces( mesh ),
		  -trans(idt(w))*cross(id(u),N()) );

#endif

  auto b = form2( _test=Xh, _trial=Xh );
  b = integrate( elements(mesh),
  		           inner(idt(u),id(u))
		              //  + gamma*inner(idt(w),id(w))
		              //  + gamma*inner(idt(p),id(p))
  		               );

  auto l = form1( _test=Xh );
  if( boption("bcVec") )
      a += on(_range=boundaryfaces(mesh), _rhs=l, _element=w, _expr=zero<3,1>() );
  else
      a += on(_range=boundaryfaces(mesh), _rhs=l, _element=w, _expr=cst(0.) );


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
    i++;
  }

  e->save();
}


int main(int argc, char **argv)
{
  po::options_description appoptions( "" );
	appoptions.add_options()
    ( "bcVec", po::value<bool>()->default_value( true ), "Dirichlet is vectorial" );

  Environment env( _argc=argc, _argv=argv,
		   _desc=appoptions,
		   _desc_lib=feel_options(),
		   _about=about(_name=
#ifdef PO_UHCURL
				"po_eigenmixed_curl",
#else
				"po_eigenmixed_div",
#endif
				_author="Romain Hild",
				_email="romain.hild@plasticomnium.com") );

  Application app;

  app.add( new EigenProblem() );
  app.run();
}
