/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <boost/mpi/timer.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

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
    typedef bases<Lagrange<Order, Vectorial>, Lagrange<Order, Scalar> > mlbasis_type;
    typedef FunctionSpace<mesh_type, mlbasis_type> mlspace_type;

    typedef Exporter<mesh_type> export_type;

    void run();
private:

}; // EigenProblem

const uint16_type EigenProblem::Order;
const uint16_type EigenProblem::Dim;

void
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

    /// [option]
    int nev = option(_name="solvereigen.nev").template as<int>();
    int ncv = option(_name="solvereigen.ncv").template as<int>();
    /// [option]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "nev= " << nev <<std::endl;
        std::cout << "ncv= " << ncv <<std::endl;
    }
    /// [mesh]
    auto mesh = loadMesh(_mesh = new mesh_type );
    /// [mesh]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "maillage = " << t.elapsed() << std::endl;
    }

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
                   + (dxt(u1)+dyt(u2)+dzt(u3)) * (dx(u1)+dy(u2)+dz(u3)) );
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
        std::cout << "matrices = " << t.elapsed() << std::endl;
    }

    /// [modes]
    SolverEigen<double>::eigenmodes_type modes;
    modes=
    eigs( _matrixA=a.matrixPtr(),
         _matrixB=b.matrixPtr(),
         _nev=nev,
         _ncv=ncv,
         _transform=SINVERT,
         _spectrum=SMALLEST_MAGNITUDE,
         _verbose = true );
    /// [modes]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "modes = " << t.elapsed() << std::endl;
    }

    auto e =  exporter( _mesh=mesh );
    auto femodes = std::vector<decltype( Xh->element() )>( modes.size(), Xh->element() );

    auto Vh = Pchv<Order>( mesh );
    auto g = std::vector<decltype( Vh->element() )>( modes.size(), Vh->element() );
    /*
    auto g0 = std::vector<decltype( Vh->element() )>( modes.size(), Vh->element() );
    auto vG0 = Vh->element();
    auto a2 = form2( _test=Vh, _trial=Vh );
    auto l2 = form1( _test=Vh );

    auto MLh = mlspace_type::New( mesh );
    auto psi = std::vector<decltype( MLh->element() )>( modes.size(), MLh->element() );
    auto VPsi = MLh->element();
    auot a3 = form2( _test=MLh, _trial=MLh );
    auto l3 = form1( _trial=MLh );
    */
    if ( !modes.empty() )
    {
        std::cout << "nev " << nev << ", ncv " << ncv << ", proc " << Environment::numberOfProcessors() << ", time " << t.elapsed() << std::endl;

        int i = 0;
        for( auto const& mode : modes )
        {
            femodes[i] = *mode.second.get<2>();
            /// [project]
            g[i] = vf::project(_space=Vh, _range=elements(mesh),
                            _expr=vec(idv(femodes[i].template element<0>()),
                                      idv(femodes[i].template element<1>()),
                                      idv(femodes[i].template element<2>()) ) );
            /// [project]
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), g[i] );

            std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second.get<1>() << ") ";
            std::cout << " div u = " << integrate(elements(mesh), divv(g[i])).evaluate()(0,0) << std::endl;
            //std::cout << "(g" << i << ",g" << i-1 << ") = " << integrate(elements(mesh), trans(idv(g[i]))*idv(g[i-1]) ).evaluate()(0,0) << std::endl;

            //a2 = integrate( _range=elements(mesh),
            //                _expr=divt(g0[i])*div(vG0) );
            //a2 +=integrate( _range=elements(mesh),
            //                _expr=gradt(g0[i])*grad(vG0) );
            //l2 = integrate( _range=elements(mesh),
            //                _expr=trans(gradv(g[i]))*grad(vG0) );
            //a2 += on( _range=boundaryfaces(mesh),
            //          _element=g0[i], _rhs=l2, _expr=cst(0.) );
            //a2.solve( _rhs=l2, _solution=g0[i] );
            //e->add( (boost::format("g0-%1%" ) % i ).str(), g0[i] );

            //auto uPsi = psi[i].template element<0>();
            //auto lambdaPsi = psi[i].template element<1>();
            //auto vPsi = VPsi.template element<0>();
            //auto nuPsi = VPsi.template element<1>();
            //a3 = integrate( elements(mesh), gradt(uPsi)*trans(grad(vPsi)) + id(vPsi)*idt(lambdaPsi) + idt(uPsi)*id(nuPsi) );
            //l3 = integrate( elements(mesh), divv(g0[i])*id(vPsi) );
            //a3.solve( -rhs=l3, _solution=psi[i] );

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
                    _about=about(_name="po_mode_gen_curl",
                                 _author="Romain Hild",
                                 _email="romain.hild@plasticomnium.com") );

    Application app;

    app.add( new EigenProblem() );
    app.run();
}





