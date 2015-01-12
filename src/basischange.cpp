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
//#include <feel/feel.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <boost/mpi/timer.hpp>

/** use Feel namespace */
using namespace Feel;

enum EdgeType {
    EDGE_INTERIOR = 0, // edge in the interior
    EDGE_BOUNDARY, // edge on boundary
    EDGE_BOUNDARY_VERTEX_1, // edge touches boundary with starting vertex
    EDGE_BOUNDARY_VERTEX_2 // edge touches boundary with ending vertex
    EDGE_BOUNDARY_VERTEX_3 // edge touches boundary with two vertices (but still is in the interior)
};

struct DofEdgeInfo
{
    size_type index;
    int8_type sign1;
    int8_type sign2;

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
    typedef typename mesh_type::element_type::edge_permutation_type edge_permutation_type;
    
    typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    // [Nh]
    // typedef typename meta::Ned1h<mesh_type,1>::type space_type;
    // typedef typename meta::Ned1h<mesh_type,1>::ptrtype space_ptrtype;

    typedef Mesh<Simplex<Dim-1,1,Dim>> boundary_mesh_type;
    typedef typename meta::Pch<boundary_mesh_type,1>::type scalar_space_type;
    typedef typename meta::Pch<boundary_mesh_type,1>::ptrtype scalar_space_ptrtype;

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

    LOG(INFO) << "[timer] h = " << doption("gmsh.hsize") << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] h = " << doption("gmsh.hsize") << std::endl;


    boost::mpi::timer t;
    boost::mpi::timer total;

    // auto mesh = loadMesh(_mesh = new mesh_type );
    auto mesh = unitSphere();
    if ( Environment::isMasterRank() ){
        std::cout << " number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << " number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << " number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << " number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << " number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " number of local vertices : " << mesh->numLocalVertices() << std::endl;
    }

    LOG(INFO) << "[timer] mesh = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] mesh = " << t.elapsed() << " sec" << std::endl;
    t.restart();

    // [Nh2]
    auto Vh = space_type::New( mesh );
    // [Nh2]
    // [P1ch]
    auto Sh = Pch<1>( mesh );
    // [P1ch]

    LOG(INFO) << "[info] Vh dof = " << Vh->nDof() << std::endl;
    LOG(INFO) << "[info] Sh dof = " << Sh->nDof() << std::endl;
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

    LOG(INFO) << "[timer] spaces = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] spaces = " << t.elapsed() << " sec" << std::endl;
    t.restart();


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

    LOG(INFO) << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    t.restart();


    //  now we need to compute C, the change of basis to handle the
    //  decomposition of the field u as a hcurl function + gradient of a H1
    //  scalar field

    // now we need to compute the set of dof indices of Xh that are
    // involved in the decomposition: we need to remove the set of dof which are
    // on the boundary (edges on the boundary)

    std::vector<bool> doneVh( Vh->nLocalDof(), false );
    std::vector<bool> doneSh( Sh->nLocalDof(), false );
    DofEdgeInfo einfo_default {invalid_uint16_type_value,1,EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value};
    std::vector<DofEdgeInfo> dof_edge_info( Vh->nLocalDof(), einfo_default );
    std::vector<size_type> indexesToKeep;
    std::vector<size_type> boundaryIndexesToKeep;

    for( int i = 0; i < Vh->nDof(); ++i )
    {
        for( auto const& dof : Vh->dof()->globalDof( i ) )
        {
            // auto dofbegin = Vh->dof()->globalDof().first;
            // auto dofend = Vh->dof()->globalDof().second;

            // for( auto dofit = dofbegin; dofit != dofend; ++dofit ) {
            //     auto const& dof = *dofit;

            // std::cout << "first : " << dof.first.sign() << std::endl;

            if ( !doneVh[ dof.first.index() ] )
            {
                // we retrieve the edge corresponding to the dof
                auto const& edge = mesh->element(dof.second.elementId()).edge(dof.second.localDof());
                auto perm = mesh->element(dof.second.elementId()).edgePermutation(dof.second.localDof());
                dof_edge_info[dof.first.index()].index = dof.first.index();

                auto const& pt1 = edge.point( 0 );
                auto const& pt2 = edge.point( 1 );
                size_type dofid1 = invalid_uint16_type_value;
                size_type dofid2 = invalid_uint16_type_value;
                // if edge is oriented from 1 to 2 (globally) then the sign
                // is -1 (the gradient associated with dof1 is oriented from
                // 2 to 1)
                dof_edge_info[dof.first.index()].sign1 = (perm==edge_permutation_type::IDENTITY)?-1:1;
                dof_edge_info[dof.first.index()].sign2 = -dof_edge_info[dof.first.index()].sign1;
 
                // we retrieve the index of the dof of Sh
                // corresponding to the vertex of the edge
                for ( uint16_type i = 0; i < mesh->numLocalVertices(); ++i )
                {
                    if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt1.id() )
                    {
                        dofid1 = Sh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                    }
                    if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt2.id() )
                    {
                        dofid2 = Sh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                    }
                }

                if ( edge.isOnBoundary() )
                {
                    // if the points are on the boundary,
                    // we keep their indexes to extract the columns
                    if( !doneSh[ dofid1 ] )
                    {
                        boundaryIndexesToKeep.push_back( Vh->nDof()+dofid1 );
                        doneSh[ dofid1 ] = true;
                    }
                    if( !doneSh[ dofid2 ] )
                    {
                        boundaryIndexesToKeep.push_back( Vh->nDof()+dofid2 );
                        doneSh[ dofid2 ] = true;
                    }

                    dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY;
                    dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                    dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
                }
                if ( !edge.isOnBoundary() )
                {
                    // if the edge is not on the boundary, we keep its index
                    // to extract the column
                    indexesToKeep.push_back(dof.first.index());

                    //both points touch the boundary
                    if ( pt1.isOnBoundary() && pt2.isOnBoundary() )
                    {
                        dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_3;
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
                        CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                        CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                    }
                    // one of the end points touch the boundary
                    else if ( pt1.isOnBoundary()  )
                    {
                        dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_1;
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                        CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    }
                    else if ( pt2.isOnBoundary()  )
                    {
                        dof_edge_info[dof.first.index()].type = EDGE_BOUNDARY_VERTEX_2;
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = invalid_size_type_value;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
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

        LOG(INFO) << "[timer] indexes = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] indexes = " << t.elapsed() << " sec" << std::endl;
        t.restart();

        // std::cout << indexesToKeep;


        // now we fill C

        // [cIntBound]
        // cInternal correspond to the matrix of all the dofs of Vh
        auto cInternal = backend()->newMatrix(_test=Vh,_trial=Vh);
        // cBoundary correspond to the matrix of all the dofs of Sh
        auto cBoundary = backend()->newMatrix(_test=Vh,_trial=Sh );
        // [cIntBound]

        auto pm = ioption("pm");

        // [fill]
        // works only with 1 proc !!
        for( int i = 0; i < Vh->nDof(); ++i )
        {
            // the matrix cInternal is the identity
            // we'll remove the colons later
            cInternal->set(i,i,1);

            // the matrix cBoundary has +1 or -1 on the edges' vertex
            // which are on the boundary
            auto const& dei = dof_edge_info[i];
            switch( dei.type )
            {
            case EDGE_BOUNDARY:
            case EDGE_BOUNDARY_VERTEX_3:
            {
                cBoundary->set(i,dei.dof_vertex_id1, dei.sign1);
                cBoundary->set(i,dei.dof_vertex_id2, dei.sign2);
            }
            break;
            case EDGE_BOUNDARY_VERTEX_1:
            {
                cBoundary->set(i,dei.dof_vertex_id1, dei.sign1);
            };
            break;
            case EDGE_BOUNDARY_VERTEX_2:
            {
                cBoundary->set(i,dei.dof_vertex_id2, dei.sign2);
            };
            break;
        }
        // [fill]

        // we close the matrices so we can work on them
        cInternal->close();
        cBoundary->close();


        // we first create a matrix with all the dof of Vh and Sh
        // [ctilde]
        BlocksBaseSparseMatrix<double> cBlock(1,2);
        cBlock(0,0) = cInternal;
        cBlock(0,1) = cBoundary;
        auto cTilde = backend()->newBlockMatrix(_block=cBlock, _copy_values=true);
        cTilde->close();
        // [ctilde]

        LOG(INFO) << "[timer] fill = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] fill = " << t.elapsed() << " sec" << std::endl;
        t.restart();


        // then extract the submatrix of the matrix Ctilde

        // we keep all the rows of this matrix
        std::vector<size_type> rows( Vh->nDof() );
        std::iota( rows.begin(), rows.end(), 0 );

        // we need to sort the indexes of the colons for Sh
        std::sort(boundaryIndexesToKeep.begin(), boundaryIndexesToKeep.end() );
        indexesToKeep.insert(indexesToKeep.end(),
                             boundaryIndexesToKeep.begin(),
                             boundaryIndexesToKeep.end() );

        // we keep only the dofs corresponding to the edges not on the boundary
        // and the points on the boundary
        // [submatrix]
        auto C = backend()->newMatrix(Vh->nDof(), indexesToKeep.size(),
                                      Vh->nDof(), indexesToKeep.size() );
        cTilde->createSubmatrix( *C, rows, indexesToKeep);
        // [submatrix]

        LOG(INFO) << "[timer] submatrix = " << t.elapsed() << " sec" << std::endl;
        LOG(INFO) << "[timer] total = " << total.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() ) {
            std::cout << "[timer] submatrix = " << t.elapsed() << " sec"
                      << std::endl;
            std::cout << "[timer] total = " << total.elapsed() << " sec" << std::endl;
        }

        LOG(INFO) << "[info] Csize = " << C->size1() << " x " << C->size2()
                  << std::endl;
        if ( Environment::worldComm().isMasterRank() ) {
            std::cout << "c : " << C->size1() << " x " << C->size2() << std::endl;
        }

        matA->printMatlab("a.m");
        matB->printMatlab("b.m");
        C->printMatlab("c.m");

#if 0
        // some test on C

        cInternal->printMatlab("cInt.m");
        cBoundary->printMatlab("cBound.m");
        cTilde->printMatlab("cTilde.m");

        std::vector<size_type> indexesComp( Sh->nDof() );
        for( int i = 0; i < Sh->nDof(); ++i ) {
            indexesComp[i] = i;
        }
        for( int i = boundaryIndexesToKeep.size() - 1; i >= 0; --i ) {
            indexesComp.erase( indexesComp.begin() + boundaryIndexesToKeep[i] - Vh->nDof() );
        }

        auto cComp = backend()->newMatrix( Vh->nDof(), indexesComp.size(),
                                           Vh->nDof(), indexesComp.size() );
        cBoundary->createSubmatrix( *cComp, rows, indexesComp );

        auto cIntNorm = cInternal->l1Norm(); // each column should be 1
        auto cCompNorm = cComp->l1Norm(); // each column should be 0
        auto cTildeNorm = cTilde->linftyNorm(); // each row should be less than 3
        auto cNorm = C->linftyNorm(); // each row should be less than 3

        if ( Environment::worldComm().isMasterRank() ) {
            std::cout << "#bouIndex : " << boundaryIndexesToKeep.size() << std::endl
                      << "normL1 of cInt : " << cIntNorm << std::endl
                      << "normInf of cTilde : " << cTildeNorm << std::endl
                      << "normL1 of c complementary : " << cCompNorm << std::endl
                      << "normInf of c : " << cNorm << std::endl;
        }
#endif

        auto Ahat = backend()->newMatrix(indexesToKeep.size(), indexesToKeep.size(),indexesToKeep.size(), indexesToKeep.size() );
        auto Bhat = backend()->newMatrix(indexesToKeep.size(), indexesToKeep.size(),indexesToKeep.size(), indexesToKeep.size() );
        backend()->PtAP( matA, C, Ahat );
        backend()->PtAP( matB, C, Bhat );

        Ahat->printMatlab("Ahat.m");
        Bhat->printMatlab("Bhat.m");

        if ( Environment::worldComm().isMasterRank() ) {
            std::cout << "nev = " << ioption(_name="solvereigen.nev")
                      << "ncv= " << ioption(_name="solvereigen.ncv") << std::endl;
        }

        auto modes = eigs( _matrixA=Ahat, _matrixB=Bhat );

        // auto e =  exporter( _mesh=mesh );

        int i = 0;
        for( auto const& mode: modes ) {
            if ( Environment::isMasterRank() ) {
                std::cout << "eigenvalue " << i << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")" << std::endl;
            }
            // e->add( ( boost::format( "mode-%1%" ) % i ).str(), *mode.second.get<2>() );
            i++;
        }

        // e->save();
    }


    int
        main( int argc, char** argv )
    {
        using namespace Feel;

        po::options_description basischangeoptions( "basischange options" );
        basischangeoptions.add_options()
            ("pm", po::value<int>()->default_value( 1 ), "plus or minus one" );

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

