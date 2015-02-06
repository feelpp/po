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
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelfilters/exporter.hpp>
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

    typedef Nedelec<0, NedelecKind::NED1> basis_edge_type;
    typedef Lagrange<1,Scalar> basis_vertex_type;
    typedef bases<basis_edge_type, basis_vertex_type> basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;

    typedef FunctionSpace<mesh_type, bases<basis_edge_type> > space_edge_type;
    typedef typename space_edge_type::element_type element_edge_type;

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

    auto Xh = space_type::New( mesh );
    auto Nh = Xh->template functionSpace<0>();
    auto Lh = Xh->template functionSpace<1>();

    LOG(INFO) << "[info] Nh dof = " << Nh->nDof() << std::endl;
    LOG(INFO) << "[info] Lh dof = " << Lh->nDof() << std::endl;

    std::ofstream file;
    file.open( (boost::format( "info_%1%" ) % Environment::worldComm().globalRank() ).str() );
    file << "Xh : " << *Xh << std::endl
         << "Nh : " << *Nh << std::endl
         << "Lh : " << *Lh << std::endl << std::endl
         << "Nh : " << Nh->dof()->mapGlobalProcessToGlobalCluster()
         << "Lh : " << Lh->dof()->mapGlobalProcessToGlobalCluster() << std::endl;
    file.close();

    auto u = Nh->element();
    auto v = Nh->element();

    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Nh, _trial=Nh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    auto matB = b.matrixPtr();
    matB->close();

    std::vector<bool> doneNh( Nh->nLocalDof(), false );
    std::vector<bool> doneLh( Lh->nLocalDof(), false );
    DofEdgeInfo einfo_default {1,-1,EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value};
    std::vector<DofEdgeInfo> dof_edge_info( Nh->nLocalDof(), einfo_default );
    std::vector<size_type> interiorIndexesToKeep;
    std::vector<size_type> boundaryIndexesToKeep;

    auto dofbegin = Nh->dof()->localDof().first;
    auto dofend = Nh->dof()->localDof().second;

    for( auto dofit = dofbegin; dofit != dofend; ++dofit )
    {
        auto const& dof = *dofit;
    // for( int i = 0; i < Nh->nLocalDof(); ++i )
    //     for( auto const& dof : Nh->dof()->localDof( i ) )
    //     {

        auto index = dof.second.index();
        if ( !doneNh[ index ] )
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
                if ( mesh->element( eltLocalId ).point( i ).id() == pt1.id() )
                {
                    dofid1 = Lh->dof()->localToGlobal( eltLocalId, i, 0 ).index() + Nh->nLocalDofWithGhost();
                }
                if ( mesh->element( eltLocalId ).point( i ).id() == pt2.id() )
                {
                    dofid2 = Lh->dof()->localToGlobal( eltLocalId, i, 0 ).index() + Nh->nLocalDofWithGhost();
                }
            }

            if ( edge.isOnBoundary() )
            {
                // if the points are on the boundary,
                // we keep their indexes to extract the columns
                if( !doneLh[ dofid1 ] )
                {
                    boundaryIndexesToKeep.push_back( dofid1 );
                    doneLh[ dofid1 ] = true;
                }
                if( !doneLh[ dofid2 ] )
                {
                    boundaryIndexesToKeep.push_back( dofid2 );
                    doneLh[ dofid2 ] = true;
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
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                }
                // the edge doesn't touch the boundary
                else {
                    dof_edge_info[index].type = EDGE_INTERIOR;
                    dof_edge_info[index].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[index].dof_vertex_id2 = invalid_size_type_value;
                }
            }

            file.open( (boost::format( "info_%1%" ) % Environment::worldComm().globalRank() ).str(), std::ios::out | std::ios::app);
            file << "index   : " << index << std::endl
                 << "\t perm : " << (int)dof_edge_info[index].sign1 << std::endl
                 << "\t type : " << dof_edge_info[index].type << std::endl
                 << "\t id1  : " << dof_edge_info[index].dof_vertex_id1 << std::endl
                 << "\t id2  : " << dof_edge_info[index].dof_vertex_id2 << std::endl
                 << "\t id1_[x,y,z] : [" << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << "]" << std::endl
                 << "\t id2_[x,y,z] : [" << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << "]" << std::endl;
            file.close();

            doneNh[ index ] = true;
        }
    }

    auto cTilde = backend()->newMatrix(_test=Nh, _trial=Xh);

    for( int i = 0; i < Nh->nLocalDof(); ++i )
    {
        cTilde->set(i,i,1);

        auto const& dei = dof_edge_info[i];
        switch( dei.type )
        {
        case EDGE_BOUNDARY:
        case EDGE_BOUNDARY_VERTEX_3:
        {
            cTilde->set(i,dei.dof_vertex_id1, dei.sign1);
            cTilde->set(i,dei.dof_vertex_id2, dei.sign2);
        }
        break;
        case EDGE_BOUNDARY_VERTEX_1:
        {
            cTilde->set(i,dei.dof_vertex_id1, dei.sign1);
        };
        break;
        case EDGE_BOUNDARY_VERTEX_2:
        {
            cTilde->set(i,dei.dof_vertex_id2, dei.sign2);
        }
        break;
        case EDGE_INTERIOR:
            break;
        }
    }
    cTilde->close();

    auto mapPToC = cTilde->mapCol().mapGlobalProcessToGlobalCluster();

    std::vector<size_type> rows( Nh->nLocalDof() );
    std::iota( rows.begin(), rows.end(), 0 );

    std::vector<size_type> indexesToKeep = interiorIndexesToKeep;
    indexesToKeep.insert(indexesToKeep.end(),
                         boundaryIndexesToKeep.begin(),
                         boundaryIndexesToKeep.end() );
    std::sort(indexesToKeep.begin(), indexesToKeep.end() );

    auto C = cTilde->createSubMatrix(rows, indexesToKeep, false, boption("checkAndFixRange"));
    C->close();

    file.open( (boost::format( "info_%1%" ) % Environment::worldComm().globalRank() ).str(), std::ios::out | std::ios::app);
    file << "Nh : " << Nh->dof()->mapGlobalProcessToGlobalCluster()
         << "Nh clusterToPro : " << Nh->dof()->mapGlobalClusterToGlobalProcess()
         << "Lh : " << Lh->dof()->mapGlobalProcessToGlobalCluster()
         << "indexes to keep : " << indexesToKeep
         << "cTilde procToCluster : " << mapPToC
         << "c procToCluster : " << C->mapColPtr()->mapGlobalProcessToGlobalCluster()
         << "c clusterToProc : " << C->mapColPtr()->mapGlobalClusterToGlobalProcess() << std::endl;

    if( boption("isPrinting") )
    {
        if( Environment::worldComm().isMasterRank() )
            std::cout << "Start printing" << std::endl;

        matA->printMatlab("a.m");
        matB->printMatlab("b.m");
        C->printMatlab("c.m");
        //C->printMatlab();
    }

    if( boption("testC") )
    {
        auto fLh  = Lh->element();
        auto fNh = Nh->element();
        for( int i = 0; i < 50; i++)
        {
            auto tmpNh = backend()->newVector( Nh );
            auto tmpLh = backend()->newVector( Lh );

            auto alphaPrime = backend()->newVector( C->mapColPtr() );
            auto map = alphaPrime->mapPtr()->mapGlobalProcessToGlobalCluster();

            std::srand(std::time(NULL));

            int n = C->mapColPtr()->nDof();
            int stop = 0;
            if(n < 100)
                stop = n/2;
            else
                stop = 50;

            file << "[" << i << "] : ";

            for( int j = 0; j < stop; j++)
            {
                int r = rand() % n;
                auto it = std::find(map.begin(), map.end(), r);
                auto index = it-map.begin();

                if( it != map.end() )
                {
                    alphaPrime->set(index, 1);

                    auto col = indexesToKeep[index];
                    file << col << " ";
                    if( col < Nh->nLocalDofWithGhost() )
                    {
                        tmpNh->set(col, 1);
                    } else
                    {
                        int dofLocalLh = col - Nh->nLocalDofWithGhost();
                        tmpLh->set(dofLocalLh, 1);
                    }
                }
            }

            file << std::endl;

            tmpNh->close();
            fNh = *tmpNh;
            tmpLh->close();
            fLh = *tmpLh;

            auto f = vf::project( _space=Nh, _range=elements(mesh),
                                  _expr= idv(fNh) + trans(gradv(fLh)) );

            alphaPrime->close();
            auto alphaVec = backend()->newVector( Nh );
            C->multVector( alphaPrime, alphaVec);
            auto alpha = Nh->element();
            alphaVec->close();
            alpha = *alphaVec;

            auto erreur = normL2( elements(mesh), idv(f) - idv(alpha) );
            auto curln = normL2( boundaryfaces(mesh), trans(curlv(f))*N() );

            if ( Environment::worldComm().isMasterRank() )
                std::cout << i  << "\t erreur : " << erreur << "\t curl*n : " << curln << std::endl;
        }
    }

    file.close();

    if( boption("eigs"))
    {
        auto Ahat = backend()->newMatrix(C->mapRowPtr(), C->mapColPtr() );
        auto Bhat = backend()->newMatrix(C->mapRowPtr(), C->mapColPtr() );
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

        auto e =  exporter( _mesh=mesh );

        int i = 0;
        for( auto const& mode: modes )
        {
            if ( Environment::isMasterRank() )
                std::cout << "eigenvalue " << i << " = (" << boost::get<0>(mode.second) << "," <<  boost::get<1>(mode.second) << ")" << std::endl;

            auto tmpVec = backend()->newVector( Nh );
            C->multVector( boost::get<2>(mode.second), tmpVec);
            tmpVec->close();
            auto tmp = Nh->element();
            tmp = *tmpVec;
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), tmp );
            i++;
        }

        e->save();
    }
}


int
main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description basischangeoptions( "basischange options" );
    basischangeoptions.add_options()
        ("isPrinting", po::value<bool>()->default_value( false ), "print matrices")
        ("useSphere", po::value<bool>()->default_value( true ), "use sphere or other geo")
        ("eigs", po::value<bool>()->default_value( false ), "compute the eigen problem")
        ("load", po::value<bool>()->default_value( false ), "load eigen vectors and values")
        ("testC", po::value<bool>()->default_value( false ), "test matrix C")
        ("checkAndFixRange", po::value<bool>()->default_value( true ), "check and fix range for createSubMatrix" );

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
