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

enum EdgeType {
    EDGE_INTERIOR = 0, // edge in the interior
    EDGE_BOUNDARY, // edge on boundary
    EDGE_BOUNDARY_VERTEX_1, // edge touches boundary with one vertex
    EDGE_BOUNDARY_VERTEX_2 // edge touches boundary with two vertices (but still is in the interior)
};

struct DofEdgeInfo
{
    size_type index;
    size_type sign;
    EdgeType type;
    size_type dof_vertex_id1;
    size_type dof_vertex_id2;

};

inline
std::ostream&
operator<<( std::ostream& os, DofEdgeInfo const& dei  )
{
    os << "-----------Dof-Edge-Info------------\n"
         << "index       : " << dei.index << "\n"
         << "sign        : " << dei.sign << "\n"
         << "type        : " << dei.type << "\n"
         << "id vertex 1 : " << dei.dof_vertex_id1 << "\n"
         << "id vertex 2 : " << dei.dof_vertex_id2 << "\n";
    return os;
}

inline
std::ostream&
operator<<( std::ostream& os, std::vector<size_type> const& vec  )
{
    for( auto v : vec )
        os << v << " ";
    os << std::endl;
    return os;
}

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
    typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    // typedef typename meta::Ned1h<mesh_type,1>::type space_type;
    // typedef typename meta::Ned1h<mesh_type,1>::ptrtype space_ptrtype;
    // [Nh]
    // [P1ch]
    typedef Mesh<Simplex<Dim-1,1,Dim>> boundary_mesh_type;
    typedef typename meta::Pch<boundary_mesh_type,1>::type scalar_space_type;
    typedef typename meta::Pch<boundary_mesh_type,1>::ptrtype scalar_space_ptrtype;
    // [P1ch]

    //typedef bases<Nedelec<0,NedelecKind::NED1>, Lagrange<1, Scalar> > mixed_basis_type;
    // // mesh_type || boundary_mesh_type
    //typedef FunctionSpace<mesh_type, basis_type > space_type;
    //typedef boost::shared_ptr<space_type> space_ptrtype;

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

    Environment::changeRepository( boost::format( "%1%/h_%2%/" )
                                   % this->about().appName()
                                   % option(_name="gmsh.hsize").template as<double>() );

    auto mesh = loadMesh(_mesh = new mesh_type );
    if ( Environment::isMasterRank() ){
        std::cout << " number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << " number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << " number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << " number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << " number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " number of local vertices : " << mesh->numLocalVertices() << std::endl;
    }

    auto Vh = space_type::New( mesh );
    auto Sh = Pch<1>( mesh );
    if ( Environment::isMasterRank() ){
        std::cout << "[Vh] number of dof             : "
                  << Vh->nDof() << std::endl;
        std::cout << "[Vh] number of dof per proc    : "
                  << Vh->nLocalDof() << std::endl;
        std::cout << "[Vh] number of local dof       : "
                  << Vh->dof()->nLocalDof() << std::endl;
        std::cout << "[Vh] total number of local dof : "
                  << Vh->dof()->nLocalDof()*mesh->numGlobalElements() << std::endl;
        std::cout << "[Sh] number of dof             : "
                  << Sh->nDof() << std::endl;
        std::cout << "[Vh] number of dof per proc    : "
                  << Sh->nLocalDof() << std::endl;
        std::cout << "[Sh] number of local dof       : "
                  << Sh->dof()->nLocalDof() << std::endl;
        std::cout << "[Sh] total number of local dof : "
                  << Sh->dof()->nLocalDof()*mesh->numGlobalElements() << std::endl;
    }

    auto u = Vh->element();
    auto v = Vh->element();

    // [forms]
    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();

    auto b = form2( _test=Vh, _trial=Vh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    // [forms]


    //  now we need to compute C, the change of basis to handle the
    //  decomposition of the field u as a hcurl function + gradient of a H1
    //  scalar field

    // now we need to compute the set of dof indices of Xh that are
    // involved in the decomposition: we need to remove the set of dof which are
    // on the boundary (edges on the boundary)

    std::vector<bool> doneVh( Vh->nLocalDof(), false );
    std::vector<DofEdgeInfo> dof_edge_info( Vh->nLocalDof(), {invalid_uint16_type_value,1,EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value} );
    std::vector<size_type> internalIndexes;

    for( int i = 0; i < Vh->nDof(); ++i )
        for( auto const& dof : Vh->dof()->globalDof( i ) ) {
    // auto dofbegin = Vh->dof()->globalDof().first;
    // auto dofend = Vh->dof()->globalDof().second;

    // for( auto dofit = dofbegin; dofit != dofend; ++dofit ) {
    //     auto const& dof = *dofit;
        if ( !doneVh[ dof.first.index() ] ) {
            auto const& edge = mesh->element(dof.second.elementId()).edge(dof.second.localDof());
            dof_edge_info[dof.first.index()].index = dof.first.index();
            dof_edge_info[dof.first.index()].sign = dof.first.sign();

            auto const& pt1 = edge.point( 0 );
            auto const& pt2 = edge.point( 1 );
            size_type dofid1 = invalid_uint16_type_value;
            size_type dofid2 = invalid_uint16_type_value;

            for ( uint16_type i = 0; i < mesh->numLocalVertices(); ++i ) {
                if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt1.id() ) {
                    dofid1 = Sh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                }
                if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt2.id() ) {
                    dofid2 = Sh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                }
            }

            if ( edge.isOnBoundary() ) {
                dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY;
                dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
            }
            if ( !edge.isOnBoundary() ) {
                // if the edge is not on the boundary, we keep its index to extract the colomn
                internalIndexes.push_back(dof.first.index());
                //both points touch the boundary
                if ( pt1.isOnBoundary() && pt2.isOnBoundary() ) {
                    dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_2;
                    dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                    dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                }
                // one of the end points touch the boundary
                else if ( pt1.isOnBoundary()  ) {
                    dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_1;
                    dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                    dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                }
                else if ( pt2.isOnBoundary()  ) {
                    dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_1;
                    dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid2;
                    dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id1";
                }
                else {
                    dof_edge_info[dof.first.index()].type = EDGE_INTERIOR;
                    dof_edge_info[dof.first.index()].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                }
            }
        }
        doneVh[dof.first.index()] = true;
    }

    // now we fill C
    auto cInternal = backend()->newMatrix(_test=Vh,_trial=Vh);
    auto cBoundary = backend()->newMatrix(_test=Vh,_trial=Sh );

    // works only with 1 proc
    for( int i = 0; i < Vh->nDof(); ++i ) {
        cInternal->set(i,i,1);
        auto dei = dof_edge_info[i];
        if( dei.type == EDGE_BOUNDARY || dei.type == EDGE_BOUNDARY_VERTEX_2 ) {
            cBoundary->set(i,dei.dof_vertex_id1, 1*dei.sign);
            cBoundary->set(i,dei.dof_vertex_id2, -1*dei.sign);
        }
        else if( dei.type == EDGE_BOUNDARY_VERTEX_1 ) {
            cBoundary->set(i,dei.dof_vertex_id1, 1*dei.sign);
        }
    }

    cInternal->close();
    cBoundary->close();


    // then extract the submatrix of the matrix Ctilde

    // We keep the index of the vertex on the boundary
    std::vector<size_type> boundaryIndexes;
    std::vector<bool> doneSh( Sh->nLocalDof(), false );

    // for( int i = 0; i < Sh->nDof(); ++i ) {
    //     for( auto const& dof : Sh->dof()->globalDof( i ) ) {
    //         if ( !doneSh[ dof.first.index() ] ) {
    //             auto const& vertex = mesh->element(dof.second.elementId()).vertices(dof.second.localDof();
    //             if( vertex.isOnBoundary() ) {
    //                 boundaryIndexes.push_back( Vh->nDof()+dof.first.index() );
    //             }
    //         }
    //         doneSh[dof.first.index()] = true;
    //     }
    // }

    internalIndexes.insert( internalIndexes.end(), boundaryIndexes.begin(), boundaryIndexes.end() );
    std::cout << internalIndexes;


    std::vector<size_type> rows( Vh->nDof() );
    for( int i = 0; i < Vh->nDof(); ++i ) {
        rows[i] = i;
    }

    BlocksBaseSparseMatrix<double> cBlock(1,2);
    cBlock(0,0) = cInternal;
    cBlock(0,1) = cBoundary;
    auto cTilde = backend()->newBlockMatrix(_block=cBlock, _copy_values=true);
    cTilde->close();

    auto C = backend()->newMatrix(Vh->nDof(), internalIndexes.size(), Vh->nDof(), internalIndexes.size() );
    cTilde->createSubmatrix( *C, rows, internalIndexes);

    matA->printMatlab("a.m");
    matB->printMatlab("b.m");
    C->printMatlab("c.m");

    std::cout << "c : " << C->size1() << "x" << C->size2() << std::endl;


    // next step is to build C^T A C and C^T B C: in feature/operator we have a
    // operator framework that allows to define such objects


#if 0
    auto opCtrans = op(C,"C",true);
    auto opC = op(C,"C");
    auto Atilde = compose(opCtrans,compose(op(A,"A"),opC) );
    auto Btilde = compose(opCtrans,compose(op(B,"B"),opC) );
#endif

    // if ( Environment::worldComm().isMasterRank() )
    // {
    //     std::cout << "number of eigenvalues computed= " << option(_name="solvereigen.nev").template as<int>() <<std::endl;
    //     std::cout << "number of eigenvalues for convergence= " << option(_name="solvereigen.ncv").template as<int>() <<std::endl;
    // }

    // auto modes= veigs( _formA=a, _formB=b );

    // auto e =  exporter( _mesh=mesh );

    // if ( e->doExport() )
    // {
    //     LOG(INFO) << "exportResults starts\n";
    //     int i = 0;
    //     for( auto const& mode: modes )
    //     {
    //         auto norml2_div = normL2(_range=elements(mesh), _expr=divv(mode.second));
    //         if ( Environment::isMasterRank() )
    //         {
    //             std::cout << "||div(u_" << i << ")||_0 = " << norml2_div << "\n";
    //         }
    //         e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), mode.second );
    //         i++;
    //     }

    //     e->save();
    //     LOG(INFO) << "exportResults done\n";
    // }

}

int
main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description basischangeoptions( "basischange options" );
    basischangeoptions.add_options()
        ("dof", po::value<int>()->default_value( 0 ), "global dof id" );

    Environment env( _argc=argc, _argv=argv,
                     _desc=basischangeoptions,
                     _desc_lib=feel_options(),
                     _about=about(_name="po_basischange",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    //app.add( new EigenProblem<2,2>() );
    app.add( new EigenProblem<3,1>() );
    app.run();
}

