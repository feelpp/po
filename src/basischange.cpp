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
// #include <feel/feelvf/vf.hpp>
// #include <feel/feelalg/solvereigen.hpp>
// #include <feel/feelfilters/exporter.hpp>
#include <boost/mpi/timer.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

enum EdgeType {
    EDGE_INTERIOR = 0, // edge in the interior
    EDGE_BOUNDARY, // edge on boundary
    EDGE_BOUNDARY_VERTEX_1, // edge touches boundary with starting vertex
    EDGE_BOUNDARY_VERTEX_2, // edge touches boundary with ending vertex
    EDGE_BOUNDARY_VERTEX_3 // edge touches boundary with two vertices (but still is in the interior)
};

struct DofEdgeInfo
{
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

inline
std::ostream&
operator<<( std::ostream& os, std::set<size_type> const& se  )
{
    for( auto s : se )
        os << s << " ";
    os << std::endl;
    return os;
}

inline
std::ostream&
operator<<( std::ostream& os, std::vector<std::vector<size_type> > const& vec  )
{
    for( auto v : vec )
        os << "[ " << v << "]\n";
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
    typedef Mesh<Simplex<Dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::element_type::edge_permutation_type edge_permutation_type;

    // [Nh]
    typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type > space_type;
    // [Nh]
    typedef typename space_type::element_type element_type;

    typedef VectorPetsc<double> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

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

    mesh_ptrtype mesh;
    if( boption("useSphere") )
        mesh = unitSphere();
    else
        mesh = loadMesh(_mesh = new mesh_type );

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
        std::cout << "[Sh] number of dof per proc    : "
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
#if 1
    auto u = Vh->element();
    auto v = Vh->element();

    // [forms]
    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Vh, _trial=Vh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    matB->close();
    // [forms]


    LOG(INFO) << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "mat A size = " << matA->size1() << "x" << matA->size2() << std::endl;
        //std::cout << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    t.restart();
#endif

    //  now we need to compute C, the change of basis to handle the
    //  decomposition of the field u as a hcurl function + gradient of a H1
    //  scalar field

    // now we need to compute the set of dof indices of Xh that are
    // involved in the decomposition: we need to remove the set of dof which are
    // on the boundary (edges on the boundary)

    std::vector<bool> doneVh( Vh->nLocalDof(), false );
    std::vector<bool> doneSh( Sh->nLocalDof(), false );
    DofEdgeInfo einfo_default {1,-1,EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value};
    std::vector<DofEdgeInfo> dof_edge_info( Vh->nLocalDof(), einfo_default );
    std::vector<size_type> interiorIndexesToKeep;
    std::vector<size_type> boundaryIndexesToKeep;

    auto dofbegin = Vh->dof()->localDof().first;
    auto dofend = Vh->dof()->localDof().second;

    for( auto dofit = dofbegin; dofit != dofend; ++dofit )
    {
        auto const& dof = *dofit;

        auto index = dof.second.index();
        if ( !doneVh[ index ] )
        {
            auto const& eltLocalId = dof.first.elementId();
            auto const& elt = mesh->element(eltLocalId);
            auto const& edgeLocalId = dof.first.localDof();
            auto const& edge = elt.edge(edgeLocalId);
            auto perm = elt.edgePermutation( edgeLocalId );

            dof_edge_info[index].sign1 = (perm==edge_permutation_type::IDENTITY) ? -1 : 1;
            dof_edge_info[index].sign2 = -dof_edge_info[index].sign1;

            auto i1 = elt.eToP(edgeLocalId, 0);
            auto i2 = elt.eToP(edgeLocalId, 1);
            auto const& pt1 = elt.point( i1 );
            auto const& pt2 = elt.point( i2 );

            size_type dofid1 = invalid_uint16_type_value;
            size_type dofid2 = invalid_uint16_type_value;
            for ( uint16_type i = 0; i < mesh->numLocalVertices(); ++i )
            {
                // [dofLh]
                if ( mesh->element( eltLocalId ).point( i ).id() == pt1.id() )
                {
                    dofid1 = Sh->dof()->localToGlobal( eltLocalId, i, 0 ).index();
                    // [dofLh]
                }
                if ( mesh->element( eltLocalId ).point( i ).id() == pt2.id() )
                {
                    dofid2 = Sh->dof()->localToGlobal( eltLocalId, i, 0 ).index();
                }
            }

            if ( edge.isOnBoundary() )
            {
                // if the points are on the boundary,
                // we keep their indexes to extract the columns
                if( !doneSh[ dofid1 ] )
                {
                    boundaryIndexesToKeep.push_back( Vh->nLocalDof()+dofid1 );
                    doneSh[ dofid1 ] = true;
                }
                if( !doneSh[ dofid2 ] )
                {
                    boundaryIndexesToKeep.push_back( Vh->nLocalDof()+dofid2 );
                    doneSh[ dofid2 ] = true;
                }

                dof_edge_info[index].type = EDGE_BOUNDARY;
                dof_edge_info[index].dof_vertex_id1 = dofid1;
                dof_edge_info[index].dof_vertex_id2 = dofid2;
            }
            if ( !edge.isOnBoundary() )
            {
                // if the edge is not on the boundary, we keep its index
                // to extract the column
                interiorIndexesToKeep.push_back(index);

                //both points touch the boundary
                if ( pt1.isOnBoundary() && pt2.isOnBoundary() )
                {
                    dof_edge_info[index].type = EDGE_BOUNDARY_VERTEX_3;
                    dof_edge_info[index].dof_vertex_id1 = dofid1;
                    dof_edge_info[index].dof_vertex_id2 = dofid2;
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                }
                // the starting point touch the boundary
                else if ( pt1.isOnBoundary()  )
                {
                    dof_edge_info[index].type = EDGE_BOUNDARY_VERTEX_1;
                    dof_edge_info[index].dof_vertex_id1 = dofid1;
                    dof_edge_info[index].dof_vertex_id2 = invalid_size_type_value;
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                }
                // the ending point touch the boundary
                else if ( pt2.isOnBoundary()  )
                {
                    dof_edge_info[index].type = EDGE_BOUNDARY_VERTEX_2;
                    dof_edge_info[index].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[index].dof_vertex_id2 = dofid2;
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id1";
                }
                // the edge doesn't touch the boundary
                else {
                    dof_edge_info[index].type = EDGE_INTERIOR;
                    dof_edge_info[index].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[index].dof_vertex_id2 = invalid_size_type_value;
                }
            }

            if(Environment::worldComm().globalRank() == 0 )
            {
                std::cout << "[" << Environment::worldComm().globalRank() << "]\t" << index << std::endl
                          << "\t elt  : " << eltLocalId << std::endl
                          << "\t edge : " << edgeLocalId << std::endl
                          << "\t perm : " << (int)dof_edge_info[index].sign1 << std::endl
                          << "\t type : " << dof_edge_info[index].type << std::endl
                          << "\t id1  : " << dof_edge_info[index].dof_vertex_id1 << std::endl
                          << "\t id2  : " << dof_edge_info[index].dof_vertex_id2 << std::endl;
            }

            Environment::worldComm().barrier();

            if(Environment::worldComm().globalRank() == 1 )
            {
                std::cout << "[" << Environment::worldComm().globalRank() << "]\t" << index << std::endl
                          << "\t elt  : " << eltLocalId << std::endl
                          << "\t edge : " << edgeLocalId << std::endl
                          << "\t perm : " << (int)dof_edge_info[index].sign1 << std::endl
                          << "\t type : " << dof_edge_info[index].type << std::endl
                          << "\t id1  : " << dof_edge_info[index].dof_vertex_id1 << std::endl
                          << "\t id2  : " << dof_edge_info[index].dof_vertex_id2 << std::endl;
            }

            doneVh[ index ] = true;
        }
    }

    auto cInternal = backend()->newMatrix(_test=Vh,_trial=Vh);
    // cBoundary correspond to the matrix of all the dofs of Sh
    auto cBoundary = backend()->newMatrix(_test=Vh,_trial=Sh );
    // [cIntBound]

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "cInternal size = " << cInternal->size1() << "x" << cInternal->size2() << std::endl;
        std::cout << "cBoundary size = " << cBoundary->size1() << "x" << cBoundary->size2() << std::endl;
    }

    // [fill]
    // works only with 1 proc !!
    for( int i = 0; i < Vh->nLocalDof(); ++i )
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
        }
        break;
        case EDGE_INTERIOR:
            break;
        }
    }
    // [fill]

    // we close the matrices so we can work on them
    cInternal->close();
    cBoundary->close();

    cInternal->printMatlab("cInt.m");
    cBoundary->printMatlab("cBou.m");

    BlocksBaseSparseMatrix<double> cBlock(1,2);
    cBlock(0,0) = cInternal;
    cBlock(0,1) = cBoundary;
    auto cTilde = backend()->newBlockMatrix(_block=cBlock, _copy_values=true);
    cTilde->close();

    std::vector<size_type> rows( Vh->nLocalDof() );
    std::iota( rows.begin(), rows.end(), 0 );

    std::sort(boundaryIndexesToKeep.begin(), boundaryIndexesToKeep.end() );
    std::vector<size_type> indexesToKeep = interiorIndexesToKeep;
    indexesToKeep.insert(indexesToKeep.end(),
                         boundaryIndexesToKeep.begin(),
                         boundaryIndexesToKeep.end() );

    std::vector<std::vector<size_type> > gathered2dIndexesToKeep;
    mpi::all_gather(Environment::worldComm().globalComm(),
                    indexesToKeep,
                    gathered2dIndexesToKeep);

    std::cout << "gathered : " << gathered2dIndexesToKeep << std::endl;

    std::vector<size_type> gatheredIndexesToKeep = gathered2dIndexesToKeep[0];
    for( int i = 1; i < Environment::numberOfProcessors(); i++)
        gatheredIndexesToKeep.insert(gatheredIndexesToKeep.end(),
                                     gathered2dIndexesToKeep[i].begin(),
                                     gathered2dIndexesToKeep[i].end() );

    std::set<size_type> globalIndexesToKeep( gatheredIndexesToKeep.begin(), gatheredIndexesToKeep.end() );

    std::cout << "global : " << globalIndexesToKeep << std::endl;

#if 0
    auto C = backend()->newMatrix(Vh->nDof(), indexesToKeep.size(),
                                  Vh->nLocalDof(), indexesToKeep.size() );
    cTilde->createSubmatrix( *C, rows, indexesToKeep);
    C->close();

    std::cout << "c : " << C->size1() << " x " << C->size2() << std::endl;
    C->printMatlab("c.m");


    auto e =  exporter( _mesh=mesh );

    if( boption("eigs"))
    {
        auto Ahat = backend()->newMatrix(indexesToKeep.size(), indexesToKeep.size(), indexesToKeep.size(), indexesToKeep.size() );
        auto Bhat = backend()->newMatrix(indexesToKeep.size(), indexesToKeep.size(), indexesToKeep.size(), indexesToKeep.size() );
        backend()->PtAP( matA, C, Ahat );
        backend()->PtAP( matB, C, Bhat );
        Ahat->close();
        Bhat->close();

        if ( Environment::worldComm().isMasterRank() )
        {
            std::cout << "nev = " << ioption(_name="solvereigen.nev")
                      << "\t ncv= " << ioption(_name="solvereigen.ncv") << std::endl;
        }

        auto modes = eigs( _matrixA=Ahat,
                           _matrixB=Bhat,
                           _solver=(EigenSolverType)EigenMap[soption("solvereigen.solver")],
                           _problem=(EigenProblemType)EigenMap[soption("solvereigen.problem")],
                           _transform=(SpectralTransformType)EigenMap[soption("solvereigen.transform")],
                           _spectrum=(PositionOfSpectrum)EigenMap[soption("solvereigen.spectrum")]
                           );


        int i = 0;
        for( auto const& mode: modes )
        {
            if ( Environment::isMasterRank() )
            {
                std::cout << "eigenvalue " << i << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")" << std::endl;
            }
            auto tmpVec = backend()->newVector( Vh );
            C->multVector( mode.second.get<2>(), tmpVec);
            element_type tmp = Vh->element();
            tmp.add(*tmpVec);
            // element_type tmp = *tmpVec;
            tmp.close();
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), tmp );
            i++;
        }

    }

    e->save();
#endif
}


int
main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description basischangeoptions( "basischange options" );
    basischangeoptions.add_options()
        ("isPrinting", po::value<bool>()->default_value( false ), "print matrices")
        ("useSphere", po::value<bool>()->default_value( true ), "use sphere or other geo")
        ("eigs", po::value<bool>()->default_value( true ), "compute the eigen problem")
        ("load", po::value<bool>()->default_value( true ), "load eigen vectors and values")
        ("pm", po::value<int>()->default_value( 1 ), "plus or minus one" );

    Environment env( _argc=argc, _argv=argv,
                     _desc=basischangeoptions,
                     _desc_lib=feel_options(),
                     _about=about(_name="po_basischange",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    app.add( new EigenProblem<3,1>() );
    app.run();
}
