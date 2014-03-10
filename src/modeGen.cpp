/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <boost/mpi/timer.hpp>
//#include <feel/feelalg/backend.hpp>
//#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feeldiscr/region.hpp>
//#include <feel/feelpoly/im.hpp>
//#include <feel/feelpoly/polynomialset.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


//template<int Dim>
class EigenProblem
:
public Simget
{
    typedef Simget super;
public:
    static const uint16_type Order = 2;
    static const uint16_type Dim = 3;
    /// [typedef]
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef bases<Lagrange<Order, Scalar>, Lagrange<Order, Scalar>, Lagrange<Order, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    /// [typedef]
    typedef bases<Lagrange<Order,Vectorial>> vbasis_type;
    typedef FunctionSpace<mesh_type, vbasis_type> vspace_type;
    typedef Exporter<mesh_type> export_type;
    //typedef double value_type;
    //typedef Backend<value_type> backend_type;
    //typedef boost::shared_ptr<backend_type> backend_ptrtype;
    //typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //typedef typename backend_type::vector_type vector_type;
    //typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //typedef boost::shared_ptr<space_type> space_ptrtype;
    //typedef typename space_type::element_type element_type;
    //typedef typename vspace_type::element_type element_type;

    void run();
private:

}; // EigenProblem

//template<int Dim> const uint16_type EigenProblem<Dim>::Order;
const uint16_type EigenProblem::Order;
const uint16_type EigenProblem::Dim;

//template<int Dim>
void
//EigenProblem<Dim>::run()
EigenProblem::run()
{
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "Execute EigenProblem<" << Dim << ">\n";
    }
    Environment::changeRepository( boost::format( "po/%1%/%2%D-P%3%/" )
                                  % this->about().appName()
                                  % Dim
                                  % Order );

    boost::mpi::timer t;

    /// [mesh]
    int nev = option(_name="solvereigen.nev").template as<int>();
    int ncv = option(_name="solvereigen.ncv").template as<int>();

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "nev= " << nev <<std::endl;
        std::cout << "ncv= " << ncv <<std::endl;
    }

    auto mesh = loadMesh(_mesh = new mesh_type );
    /// [mesh]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "mesh= " << t.elapsed() <<std::endl;
    }

    auto Vh = vspace_type::New( mesh );
    auto W = Vh->element();
    /// [space]
    auto Xh = space_type::New( mesh );
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto u3 = U.template element<2>();
    auto V = Xh->element();
    auto v1 = V.template element<0>();
    auto v2 = V.template element<1>();
    auto v3 = V.template element<2>();
    /// [space]

    auto l = form1( _test=Xh );
    /// [bilinear]
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ), (dyt(u3)-dzt(u2)) * (dy(v3)-dz(v2))
                   + (dzt(u1)-dxt(u3)) * (dz(v1)-dx(v3))
                   + (dxt(u2)-dyt(u1)) * (dx(v2)-dy(v1))
                   + (dxt(u1)+dyt(u2)+dzt(u3)) * (dx(u1)+dy(u2)+dz(u3)));
    /// [bilinear]

    /// [boundary]
    a += on(_range=markedfaces(mesh, 1), _rhs=l, _element=u3, _expr=cst(0.));
    a += on(_range=markedfaces(mesh, 2), _rhs=l, _element=u3, _expr=cst(0.));
    a += on(_range=markedfaces(mesh, 3), _rhs=l, _element=u1, _expr=cst(0.));
    a += on(_range=markedfaces(mesh, 3), _rhs=l, _element=u2, _expr=cst(0.));
    /// [boundary]

    /// [rhs]
    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), idt( u1 )*id( v1 )
                   + idt( u2 )*id( v2 )
                   + idt( u3 )*id( v3 ) );
    /// [rhs]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "matrices= " << t.elapsed() <<std::endl;
    }

    SolverEigen<double>::eigenmodes_type modes;
    /// [modes]
    modes= eigs( _matrixA=a.matrixPtr(),
         _matrixB=b.matrixPtr(),
         _nev=nev,
         _ncv=ncv,
         _transform=SINVERT,
         _spectrum=SMALLEST_MAGNITUDE,
         _verbose = true );
    /// [modes]

    auto e =  exporter( _mesh=mesh );
    auto femodes = std::vector<decltype( Xh->element() )>( modes.size(), Xh->element() );
    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";
        LOG(INFO) << "nev " << nev << " ncv " << ncv << " proc " << Environment::numberOfProcessors() << " timer " << t.elapsed() << std::endl;

        int i = 0;
        for( auto const& mode : modes )
        {
            std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second << ") ";
            femodes[i] = *mode.second.get<2>();
            /// [project]
            W = vf::project(_space=Vh, _range=elements(mesh),
                            _expr=vec(idv(femodes[i].template element<0>()),
                                      idv(femodes[i].template element<1>()),
                                      idv(femodes[i].template element<2>()) ) );
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), W );
            /// [project]
            std::cout << " div u = " << integrate(elements(mesh), divv(W)).evaluate()(0,0) << std::endl;
            ++i;
        }
        e->save();
    }
}

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                    _desc=feel_options(),
                    _about=about(_name="po_modeGen",
                                 _author="Romain Hild",
                                 _email="romain.hild@plasticomnium.com") );

    Application app;

    //    app.add( new EigenProblem<3>() );
    app.add( new EigenProblem() );
    app.run();
}





