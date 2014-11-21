/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 18 Nov 2014
 
 Copyright (C) 2014 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <feel/feel.hpp>
#include <feel/feeldiscr/ned1h.hpp>

/** use Feel namespace */
using namespace Feel;


template<int Dim, int Order>
class EigenProblem
:
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;
    // [Nh]
    typedef Mesh<Simplex<Dim>> mesh_type;
    typedef meta::Ned1h<mesh_type,1>::type space_type;
    typedef meta::Ned1h<mesh_type,1>::ptrtype space_ptrtype;
    // [Nh]
    // [P1ch]
    typedef Mesh<Simplex<Dim-1,1,Dim>> boundary_mesh_type;
    typedef meta::P1ch<boundary_mesh_type,1>::type scalar_space_type;
    typedef meta::P1ch<boundary_mesh_type,1>::ptrtype scalar_space_ptrtype;
    // [P1ch]

    void run();
private:

}; // EigenProblem

template<int Dim, int Order>
void
EigenProblem<Dim, Order>::run()
{
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "Execute EigenProblem<" << Dim << ">\n";
    }

    Environment::changeRepository( boost::format( "eigen/%1%/%2%D-P%3%/h_%4%/" )
                                   % this->about().appName()
                                   % Dim
                                   % Order
                                   % option(_name="gmsh.hsize").template as<double>() );

    auto mesh = loadMesh(_mesh = new mesh_type );

    auto Xh = space_type::New( mesh );
    // [boundarymesh]
    auto boundarymesh = createSubmesh(mesh, boundaryfaces(mesh));
    auto Vh = scalar_space_type::New( boundarymesh );
    // [boundarymesh]
    
    auto u = Xh->element();
    auto v = Xh->element();
    auto l = form1( _test=Xh );
    // [forms]
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();
    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    // [forms]

    //  now we need to compute C, the change of basis to handle the
    //  decomposition of the field u as a hcurl function + gradient of a H1
    //  scalar field
    auto ctilde  = form2(_test=Xh,_trial=XhxVh );
    // we have computed the pattern (non-zero entries) of the matrix
    auto Ctilde = ctilde.matrixPtr();

    // now we need to compute the set of dof indices of Xh and Xh that are
    // involved in the decomposition: we need to remove the set of dof which are
    // on the boundary (edges on the boundary)
    //for( auto dof: )
    
    // then extract the submatrix of the matrix Ctilde
    //auto C = Ctilde.extract( rows, cols );
    
    // now we fill C

    // next step is to build C^T A C and C^T B C: in feature/operator we have a
    // operator framework that allows to define such objects
#if 0
    auto opCtrans = op(C,"C",true);
    auto opC = op(C,"C");
    auto Atilde = compose(opCtrans,compose(op(A,"A"),opC) );
    auto Btilde = compose(opCtrans,compose(op(B,"B"),opC) );
#endif
    
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "number of eigenvalues computed= " << option(_name="solvereigen.nev").template as<int>() <<std::endl;
        std::cout << "number of eigenvalues for convergence= " << option(_name="solvereigen.ncv").template as<int>() <<std::endl;
    }

    auto modes= veigs( _formA=a, _formB=b );

    auto e =  exporter( _mesh=mesh );

    if ( e->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";
        int i = 0;
        for( auto const& mode: modes )
        {
            auto norml2_div = normL2(_range=elements(mesh), _expr=divv(mode.second));
            if ( Environment::isMasterRank() )
            {
                std::cout << "||div(u_" << i << ")||_0 = " << norml2_div << "\n";
            }
            e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), mode.second );
            i++;
        }

        e->save();
        LOG(INFO) << "exportResults done\n";
    }

}

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="ge_curl",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    //app.add( new EigenProblem<2,2>() );
    app.add( new EigenProblem<3,1>() );
    app.run();
}

