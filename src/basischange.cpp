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
#include <feel/feelinfo.h>
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

template<typename T>
std::ostream&
operator<<( std::ostream& out,  const std::vector<T>& vec  )
{
    out << "[ ";
    for( auto v : vec )
        out << v << " ";
    out << "]" << std::endl;
    return out;
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

template <class InputIterator>
double
vectorMeans(InputIterator first, InputIterator last)
{
    double res = std::accumulate(first, last, 0.0);
    return res/(last-first);
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
    typedef boost::shared_ptr<space_edge_type> space_edge_ptrtype;
    typedef typename space_edge_type::element_type element_type;

    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    void run();
private:
    space_edge_ptrtype Vh;
    sparse_matrix_ptrtype C;

    std::vector<double> lambda;
    std::vector<element_type> g;

    void load();
    void loadMatlab();
    void save();

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

    Feel::fs::path mesh_name(soption("gmsh.filename"));

    Environment::changeRepository( boost::format( "%1%/%2%/h_%3%/" )
                                   % this->about().appName()
                                   % mesh_name.stem().string()
                                   % doption("gmsh.hsize") );


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

    auto Xh = space_type::New( mesh );
    auto Nh = Xh->template functionSpace<0>();
    auto Lh = Xh->template functionSpace<1>();

    LOG(INFO) << "[info] Nh dof = " << Nh->nDof() << std::endl;
    LOG(INFO) << "[info] Lh dof = " << Lh->nDof() << std::endl;

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

    LOG(INFO) << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] form = " << t.elapsed() << " sec" << std::endl;
    t.restart();

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
                    dofid1 = Lh->dof()->localToGlobal( eltLocalId, i, 0 ).index();
                }
                if ( mesh->element( eltLocalId ).point( i ).id() == pt2.id() )
                {
                    dofid2 = Lh->dof()->localToGlobal( eltLocalId, i, 0 ).index();
                }
            }

            if ( edge.isOnBoundary() )
            {
                // if the points are on the boundary,
                // we keep their indexes to extract the columns
                if( !doneLh[ dofid1 ] )
                {
                    boundaryIndexesToKeep.push_back( dofid1 + Nh->nLocalDofWithGhost() );
                    doneLh[ dofid1 ] = true;
                }
                if( !doneLh[ dofid2 ] )
                {
                    boundaryIndexesToKeep.push_back( dofid2 + Nh->nLocalDofWithGhost() );
                    doneLh[ dofid2 ] = true;
                }

                dof_edge_info[index].type = EDGE_BOUNDARY;
                dof_edge_info[index].dof_vertex_id1 = dofid1 + Nh->nLocalDofWithGhost();
                dof_edge_info[index].dof_vertex_id2 = dofid2 + Nh->nLocalDofWithGhost();
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
                    dof_edge_info[index].dof_vertex_id1 = dofid1 + Nh->nLocalDofWithGhost();
                    dof_edge_info[index].dof_vertex_id2 = dofid2 + Nh->nLocalDofWithGhost();
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                }
                // the starting point touch the boundary
                else if ( pt1.isOnBoundary()  )
                {
                    dof_edge_info[index].type = EDGE_BOUNDARY_VERTEX_1;
                    dof_edge_info[index].dof_vertex_id1 = dofid1 + Nh->nLocalDofWithGhost();
                    dof_edge_info[index].dof_vertex_id2 = invalid_size_type_value;
                    CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                }
                // the ending point touch the boundary
                else if ( pt2.isOnBoundary()  )
                {
                    dof_edge_info[index].type = EDGE_BOUNDARY_VERTEX_2;
                    dof_edge_info[index].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[index].dof_vertex_id2 = dofid2 + Nh->nLocalDofWithGhost();
                    CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                }
                // the edge doesn't touch the boundary
                else {
                    dof_edge_info[index].type = EDGE_INTERIOR;
                    dof_edge_info[index].dof_vertex_id1 = invalid_size_type_value;
                    dof_edge_info[index].dof_vertex_id2 = invalid_size_type_value;
                }
            }

            doneNh[ index ] = true;
        }
    }

    LOG(INFO) << "[timer] info = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] info = " << t.elapsed() << " sec" << std::endl;
    t.restart();

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

    std::vector<size_type> rows( Nh->nLocalDof() );
    std::iota( rows.begin(), rows.end(), 0 );

    std::vector<size_type> indexesToKeep = interiorIndexesToKeep;
    indexesToKeep.insert(indexesToKeep.end(),
                         boundaryIndexesToKeep.begin(),
                         boundaryIndexesToKeep.end() );
    std::sort(indexesToKeep.begin(), indexesToKeep.end() );

    auto C = cTilde->createSubMatrix(rows, indexesToKeep, false, boption("checkAndFixRange"));
    C->close();

    LOG(INFO) << "[timer] matrixC = " << t.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] matrixC = " << t.elapsed() << " sec" << std::endl;
    t.restart();

    if( boption("isPrinting") )
    {
        if( Environment::worldComm().isMasterRank() )
            std::cout << "Start printing" << std::endl;

        matA->printMatlab("a.m");
        matB->printMatlab("b.m");
        C->printMatlab("c.m");
    }

    if( boption("testC") )
    {
        std::ofstream file;
        file.open( (boost::format( "info_%1%" ) % Environment::worldComm().globalRank() ).str() );
        file << "Xh : " << *Xh << std::endl
             << "Nh : " << *Nh << std::endl
             << "Lh : " << *Lh << std::endl;


        auto fLh  = Lh->element();
        auto fNh = Nh->element();
        for( int i = 0; i < 50; i++)
        {
            file << "[" << i << "] : ";

            auto tmpNh = backend()->newVector( Nh );
            auto tmpLh = backend()->newVector( Lh );

            auto alphaPrime = backend()->newVector( C->mapColPtr() );
            auto map = alphaPrime->mapPtr()->mapGlobalProcessToGlobalCluster();

            std::srand(std::time(NULL));

            int n = C->mapColPtr()->nDof();
            int stop = n < 100 ? n/2 : 50;

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

        file.close();

        LOG(INFO) << "[timer] testC = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] testC = " << t.elapsed() << " sec" << std::endl;
        t.restart();
    }

    if( boption("eigs"))
    {
        auto Ahat = backend()->newMatrix(C->mapRowPtr(), C->mapColPtr() );
        auto Bhat = backend()->newMatrix(C->mapRowPtr(), C->mapColPtr() );
#if FEELPP_VERSION_GREATER_THAN(0,99,2)
        backend()->PtAP( matA, C, Ahat );
        backend()->PtAP( matB, C, Bhat );
#endif
        Ahat->close();
        Bhat->close();

        LOG(INFO) << "[timer] hat = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] hat = " << t.elapsed() << " sec" << std::endl;
        t.restart();

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

        LOG(INFO) << "[timer] " << modes.size() << " eigs = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] " << modes.size() << " eigs = " << t.elapsed() << " sec" << std::endl;
        t.restart();

        auto e =  exporter( _mesh=mesh );

        lambda = std::vector<double>(modes.size(), 0);
        g = std::vector<element_type>(modes.size(), Vh->element());
        int i = 0;
        for( auto const& mode: modes )
        {
            if ( Environment::isMasterRank() )
                std::cout << "eigenvalue " << i << " = (" << boost::get<0>(mode.second) << "," <<  boost::get<1>(mode.second) << ")" << std::endl;
            lambda[i] = boost::get<0>(mode.second);

            auto tmpVec = backend()->newVector( Nh );
            C->multVector( boost::get<2>(mode.second), tmpVec);
            tmpVec->close();
            g[i] = *tmpVec;
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), g[i] );
            i++;
        }

        e->save();
        save();

        LOG(INFO) << "[timer] export = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] export = " << t.elapsed() << " sec" << std::endl;
    }

    if( boption("load") )
        load();

    if( boption("loadMatlab") )
        loadMatlab();

    if( boption("load") || boption("loadMatlab") )
    {
        int nbMode = ioption("nbMode");
        for( int i = 0; i < nbMode; i++)
        {
            int multiplicity = ioption("multiplicity");

            auto ev = std::vector<double>(multiplicity, 0);
            auto di = std::vector<double>(multiplicity, 0);
            auto nor = std::vector<double>(multiplicity, 0);
            auto curln = std::vector<double>(multiplicity, 0);
            auto errp = std::vector<double>(multiplicity, 0);
            auto errm = std::vector<double>(multiplicity, 0);
            auto curl2n = std::vector<double>(multiplicity, 0);
            auto errC2 = std::vector<double>(multiplicity, 0);

            for( int j = 0; j < multiplicity; j++)
            {
                int ind = i*multiplicity + j;

                ev[j] = std::sqrt(lambda[ind]);
                di[j] = normL2( elements(mesh), divv(g[ind]));
                nor[j] = normL2( boundaryfaces(mesh), trans(idv(g[ind]))*N() );
                curln[j] = normL2( boundaryfaces(mesh), trans(curlv(g[ind]))*N() );
                errp[j] = normL2( elements(mesh), curlv(g[ind])-ev[j]*idv(g[ind]) );
                errm[j] = normL2( elements(mesh), curlv(g[ind])+ev[j]*idv(g[ind]) );

                auto tmpCurl = vf::project( _space=Vh, _range=elements(mesh),
                                            _expr=curlv(g[ind]) );
                curl2n[j] = normL2( boundaryfaces(mesh), trans(curlv(tmpCurl))*N() );
                errC2[j] = normL2( elements(mesh), curlv(tmpCurl)-lambda[ind]*idv(g[ind]) );

                std::cout << "mode " << ind << " : " << ev[j] << std::endl
                          << "\t\t |div    | = " << di[j] << std::endl
                          << "\t\t |normale| = " << nor[j] << std::endl
                          << "\t\t |curl*n | = " << curln[j] << std::endl
                          << "\t\t |erreur+| = " << errp[j] << std::endl
                          << "\t\t |erreur-| = " << errm[j] << std::endl
                          << "\t\t |curl2*n| = " << curl2n[j] << std::endl
                          << "\t\t |err2   | = " << errC2[j] << std::endl;
            }

            double lambda_av = vectorMeans(ev.begin(), ev.end());
            double di_av = vectorMeans(di.begin(),di.end());
            double nor_av = vectorMeans(nor.begin(),nor.end());
            double curln_av = vectorMeans(curln.begin(),curln.end());
            double errp_av = vectorMeans(errp.begin(),errp.end());
            double errm_av = vectorMeans(errm.begin(),errm.end());
            double curl2n_av = vectorMeans(curl2n.begin(),curl2n.end());
            double errC2_av = vectorMeans(errC2.begin(),errC2.end());

            std::cout << "********** average **********" << std::endl
                      << "\t\t lambda    = " << lambda_av << std::endl
                      << "\t\t |div    | = " << di_av << std::endl
                      << "\t\t |normale| = " << nor_av << std::endl
                      << "\t\t |curl*n | = " << curln_av << std::endl
                      << "\t\t |erreur+| = " << errp_av << std::endl
                      << "\t\t |erreur-| = " << errm_av << std::endl
                      << "\t\t |curl2*n| = " << curl2n_av << std::endl
                      << "\t\t |err2   | = " << errC2_av << std::endl;
        }

        LOG(INFO) << "[timer] testMode = " << t.elapsed() << " sec" << std::endl;
        if ( Environment::worldComm().isMasterRank() )
            std::cout << "[timer] testMode = " << t.elapsed() << " sec" << std::endl;
        t.restart();


#if 0
        if( boption("space") )
        {
            auto fVh = Vh->element();
            auto fSh = Sh->element();

            srand(std::time(NULL));

            for( int i = 0; i < interiorIndexesToKeep.size(); i++)
                if(rand() % 10 > 8)
                    fVh.set(interiorIndexesToKeep[i], 1);
            for( int i = 0; i < boundaryIndexesToKeep.size(); i++)
                if(rand() % 10 > 8)
                    fSh.set(boundaryIndexesToKeep[i] - Vh->nDof(), 1);

            auto f = vf::project( _space=Vh, _range=elements(mesh),
                                  _expr= idv(fVh) + trans(gradv(fSh)) );

            auto dif = normL2( elements(mesh), divv(f));
            auto norf = normL2( boundaryfaces(mesh), trans(idv(f))*N() );
            auto curlnf = normL2( boundaryfaces(mesh), trans(curlv(f))*N() );

            std::cout << "f" << std::endl
                      << "\t\t |div    | = " << dif << std::endl
                      << "\t\t |normale| = " << norf << std::endl
                      << "\t\t |curl*n | = " << curlnf << std::endl;


            auto coef = std::vector<double>(nbMode, 0);
            auto tmpProj = Vh->element();
            double errProj = 0.0;

            for( int i = 0; i < nbMode; i++)
            {
                coef[i] = integrate(elements(mesh), idv(f)*trans(idv(g[i])) ).evaluate()(0,0);
                tmpProj.add(coef[i],g[i]);
                errProj = normL2( elements(mesh), idv(f)-idv(tmpProj) );
                std::cout << "i = " << i << "\t c = " << coef[i] << "\terreur = " << errProj << std::endl;
            }

            LOG(INFO) << "[timer] space = " << t.elapsed() << " sec" << std::endl;
            if ( Environment::worldComm().isMasterRank() )
                std::cout << "[timer] space = " << t.elapsed() << " sec" << std::endl;
            t.restart();
        }
#endif
    }

    LOG(INFO) << "[timer] total = " << total.elapsed() << " sec" << std::endl;
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "[timer] total = " << total.elapsed() << " sec" << std::endl;
}

template<int Dim, int Order>
void
EigenProblem<Dim, Order>::save()
{
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);
    for( int i = 0; i < lambda.size(); i++)
    {
        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].save(_path=path);
        if ( Environment::worldComm().isMasterRank() )
            s << lambda[i] << std::endl;
    }
    if ( Environment::worldComm().isMasterRank() )
        s.close();
}

template<int Dim, int Order>
void
EigenProblem<Dim, Order>::load()
{
    int nbMode = ioption("nbMode");
    int multiplicity = ioption("multiplicity");
    int total = nbMode * multiplicity;
    lambda = std::vector<double>(total, 0);
    g = std::vector<element_type>(total, Vh->element());

    std::fstream s;
    s.open ("lambda", std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found" << std::endl;
        exit(0);
    }

    for( int i = 0; i < total && s.good(); i++ ){
        std::string path = (boost::format("mode-%1%")%i).str();
        g[i].load(_path=path);

        s >> lambda[i];
    }

    s.close();
}

template<int Dim, int Order>
void
EigenProblem<Dim, Order>::loadMatlab()
{
    // ASSERT 1 PROC

    int nbMode = ioption("nbMode");
    int multiplicity = ioption("multiplicity");
    int total = nbMode * multiplicity;
    lambda = std::vector<double>(total, 0);
    g = std::vector<element_type>(total, Vh->element());

    std::fstream s;
    std::fstream sv;
    s.open("lambda", std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found" << std::endl;
        exit(0);
    }

    for( int i = 0; i < total; i++)
    {
        s >> lambda[i];
        int k = i + 1;
        sv.open ( (boost::format( "vec_%1%" ) % k ).str(), std::fstream::in);
        if( !sv.is_open() )
        {
            std::cout << "Eigen vectors not found" << std::endl;
            exit(0);
        }

        vector_ptrtype v = backend()->newVector( C->mapColPtr() );
        double val;
        for( int j = 0; j < C->mapColPtr()->nDof() && sv.good(); j++ )
        {
            sv >> val;
            v->set(j, val);
        }
        sv.close();

        auto tmp = backend()->newVector( Vh );
        C->multVector( v, tmp);
        tmp->close();
        g[i] = *tmp;
    }
    s.close();

}

int
main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description basischangeoptions( "basischange options" );
    basischangeoptions.add_options()
        ("isPrinting", po::value<bool>()->default_value( false ), "print matrices")
        ("useSphere", po::value<bool>()->default_value( true ), "use sphere or other geo")
        ("checkAndFixRange", po::value<bool>()->default_value( true ), "check and fix range for createSubMatrix" )
        ("eigs", po::value<bool>()->default_value( true ), "compute the eigen problem")
        ("load", po::value<bool>()->default_value( false ), "load eigen vectors and values")
        ("loadMatlab", po::value<bool>()->default_value( false ), "load eigen vectors and values")
        ("nbMode", po::value<int>()->default_value( 1 ), "number of modes to load" )
        ("multiplicity", po::value<int>()->default_value( 1 ), "multiplicity of the modes to load" )
        ("space", po::value<bool>()->default_value( false ), "verify the eigen space")
        ("testC", po::value<bool>()->default_value( false ), "test matrix C");

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
