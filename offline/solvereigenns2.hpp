#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

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

template<typename FunctionSpaceType1, typename FunctionSpaceType2>
class SolverEigenNS2
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::element_type::edge_permutation_type edge_permutation_type;

    typedef FunctionSpaceType1 space_ptrtype;
    typedef typename space_ptrtype::element_type space_type;

    typedef FunctionSpaceType2 scalar_space_ptrtype;
    typedef typename scalar_space_ptrtype::element_type scalar_space_type;
    typedef typename scalar_space_type::element_type scalar_element_type;

    typedef SolverEigenNS2<space_ptrtype, scalar_space_ptrtype> solvereigenns2_type;
    typedef typename boost::shared_ptr<solvereigenns2_type> solvereigenns2_ptrtype;

    typedef typename space_type::template sub_functionspace<0>::type space_edge_type;
    typedef boost::shared_ptr<space_edge_type> space_edge_ptrtype;
    typedef typename space_edge_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<1>::type space_vertex_type;
    typedef boost::shared_ptr<space_vertex_type> space_vertex_ptrtype;

    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef std::tuple<value_type, element_type, scalar_element_type> eigentuple_type;
    typedef std::vector<eigentuple_type> eigenmodes_type;

    typedef std::vector<element_type> eigen0_type;

    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = boost::shared_ptr<exporter_type>;


    mesh_ptrtype mesh;
    space_ptrtype Xh;
    space_edge_ptrtype Nh;
    space_vertex_ptrtype Lh;
    scalar_space_ptrtype Sh;

    std::vector<DofEdgeInfo> dof_edge_info;
    std::vector<size_type> interiorIndexesToKeep;
    std::vector<size_type> boundaryIndexesToKeep;
    std::vector<size_type> indexesToKeep;

    sparse_matrix_ptrtype matA;
    sparse_matrix_ptrtype matB;
    sparse_matrix_ptrtype C;
    sparse_matrix_ptrtype aHat;
    sparse_matrix_ptrtype bHat;

    eigenmodes_type modes;

    eigen0_type modes0;

    double tol = 1.e-3;

    exporter_ptrtype e;

    void setForms();
    void setInfo();
    void setDofsToRemove();
    void setMatrices();
    void solveEigen();
    void print();
    void save();
    void load();

    void testEigs();

public:
    static solvereigenns2_ptrtype build(const mesh_ptrtype& mesh, const FunctionSpaceType1& Xh, const FunctionSpaceType2& Sh);
    eigenmodes_type solve();
};

template<typename T1, typename T2>
typename SolverEigenNS2<T1,T2>::solvereigenns2_ptrtype
SolverEigenNS2<T1,T2>::build(const mesh_ptrtype& mesh, const T1& Xh, const T2& Sh)
{
    solvereigenns2_ptrtype ap( new SolverEigenNS2<T1,T2> );
    ap->mesh = mesh;
    ap->Xh = Xh;
    ap->Nh = Xh->template functionSpace<0>();
    ap->Lh = Xh->template functionSpace<1>();
    ap->Sh = Sh;
    return ap;
}

template<typename T1, typename T2>
typename SolverEigenNS2<T1,T2>::eigenmodes_type
SolverEigenNS2<T1,T2>::solve()
{
    tic();

    if( boption("eigen.compute"))
    {
        if ( Environment::isMasterRank() && ioption("offline.verbose") > 2)
            std::cout << " ---------- compute eigenmodes ----------\n";

        setForms();

        setInfo();
        setDofsToRemove();

        setMatrices();

        if( boption("eigen.print") )
            print();

        solveEigen();

        save();
    }
    else
    {
        if ( Environment::isMasterRank() && ioption("offline.verbose") > 2)
            std::cout << " ---------- load eigenmodes ----------\n";

        load();
    }

    if( boption( "eigen.export") )
    {
        tic();
        e = exporter( _mesh=mesh, _name="eigen" );
        for( int i = 0; i < modes.size(); i += 10)
        {
            e->add( ( boost::format( "mode-%1%" ) % i ).str(), std::get<1>(modes[i]) );
            e->add( ( boost::format( "psi-%1%" ) % i ).str(), std::get<2>(modes[i]) );
        }
        e->save();
        toc("export", ioption("offline.verbose") > 2);
    }
    toc("eigen", ioption("offline.verbose") > 1 );

    if( boption("eigen.test") )
        testEigs();

    return modes;
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setForms()
{
    tic();
    // [forms]
    auto u = Nh->element();
    auto v = Nh->element();

    // penalizaion for g.n, default value 0
    auto gamma = doption("parameters.gamma");

    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    a += integrate( boundaryfaces(mesh), gamma*(trans(idt(u))*N())*(trans(id(v))*N()) );
    matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Nh, _trial=Nh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    matB = b.matrixPtr();
    matB->close();
    // [forms]
    toc("forms", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setInfo()
{
    tic();
    std ::vector<bool> doneNh( Nh->nLocalDof(), false );
    std::vector<bool> doneLh( Lh->nLocalDof(), false );
    DofEdgeInfo einfo_default {1,-1,EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value};
    dof_edge_info = std::vector<DofEdgeInfo>( Nh->nLocalDof(), einfo_default );
    interiorIndexesToKeep = std::vector<size_type>();
    boundaryIndexesToKeep = std::vector<size_type>();

    // [loop]
    auto dofbegin = Nh->dof()->localDof().first;
    auto dofend = Nh->dof()->localDof().second;

    for( auto dofit = dofbegin; dofit != dofend; ++dofit )
        // [loop]
    {
        auto const& dof = *dofit;

        auto index = dof.second.index();
        if ( !doneNh[ index ] )
        {
            // [perm]
            auto const& eltLocalId = dof.first.elementId();
            auto const& elt = mesh->element(eltLocalId);
            auto const& edgeLocalId = dof.first.localDof();
            auto const& edge = elt.edge(edgeLocalId);
            auto perm = elt.edgePermutation( edgeLocalId );
            // [perm]

            dof_edge_info[index].sign1 = (perm==edge_permutation_type::IDENTITY) ? -1 : 1;
            dof_edge_info[index].sign2 = -dof_edge_info[index].sign1;

            // [points]
            auto i1 = elt.eToP(edgeLocalId, 0);
            auto i2 = elt.eToP(edgeLocalId, 1);
            auto const& pt1 = elt.point( i1 );
            auto const& pt2 = elt.point( i2 );
            // [points]

            size_type dofid1 = invalid_uint16_type_value;
            size_type dofid2 = invalid_uint16_type_value;
            // [dofLh]
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
            // [dofLh]

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
    toc( "infos", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setDofsToRemove()
{
    tic();
    // [remove]
    auto markers = vsoption("eigen.marker-list");
    auto s = markers.size();
    auto dofsToRemove = std::vector<int>();

    for( auto it = markers.begin(); it != markers.end(); ++it )
    {
        std::string m = *it;

        if( mesh->hasMarker(m) )
        {
            auto r = Lh->dof()->markerToDof( m );
            size_type globalDof;
            auto map = Lh->dof()->mapGlobalProcessToGlobalCluster();
            if( r.first == r.second )
            {
                globalDof = SIZE_MAX;
            }
            else
            {
                auto localDof = r.first->second;
                globalDof = map[localDof];
            }
            int minGlobalDof;

            MPI_Allreduce( &globalDof, &minGlobalDof, 1, MPI_INT, MPI_MIN, Environment::worldComm());

            auto itToRemove = std::find(map.begin(), map.end(), minGlobalDof);
            if( itToRemove != map.end())
            {
                auto localDofToRemove = itToRemove - map.begin();
                auto indexeToRemove = std::find(boundaryIndexesToKeep.begin(),boundaryIndexesToKeep.end(),
                                                localDofToRemove + Nh->nLocalDofWithGhost());
                if( indexeToRemove != boundaryIndexesToKeep.end() )
                    boundaryIndexesToKeep.erase(indexeToRemove);

                if( Environment::isMasterRank() && ioption("offline.verbose") > 3)
                    std::cout << "#" << Environment::worldComm().globalRank()
                              << " remove index : " << localDofToRemove << std::endl;
            }
        }
        else
        {
            LOG(WARNING) << "Marker \"" << m << "\" not found !!" << std::endl;
            if( Environment::isMasterRank() )
                std::cout << "Marker \"" << m << "\" not found !!" << std::endl;
        }
    }
    // [remove]
    toc( "remove", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setMatrices()
{
    tic();
    // [fill]
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
    // [fill]

    // [submatrix]
    std::vector<size_type> rows( Nh->nLocalDof() );
    std::iota( rows.begin(), rows.end(), 0 );

    indexesToKeep = interiorIndexesToKeep;
    indexesToKeep.insert(indexesToKeep.end(),
                         boundaryIndexesToKeep.begin(),
                         boundaryIndexesToKeep.end() );
    std::sort(indexesToKeep.begin(), indexesToKeep.end() );

    C = cTilde->createSubMatrix(rows, indexesToKeep);
    C->close();
    // [submatrix]

    // [ptap]
    aHat = backend()->newMatrix(C->mapColPtr(), C->mapColPtr() );
    bHat = backend()->newMatrix(C->mapColPtr(), C->mapColPtr() );
    backend()->PtAP( matA, C, aHat );
    backend()->PtAP( matB, C, bHat );
    aHat->close();
    bHat->close();
    // [ptap]
    toc("matrices", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::solveEigen()
{
    tic();
    auto zhmodes = eigs( _matrixA=aHat, _matrixB=bHat );

    int i = 0;
    modes = eigenmodes_type(zhmodes.size(), std::make_tuple(0, Nh->element(), Sh->element()) );
    modes0 = eigen0_type(zhmodes.size(), Nh->element());
    for( auto const& pair: zhmodes )
    {
        if(Environment::isMasterRank() && ioption("offline.verbose") > 1)
            std::cout << i << " eigenvalue = " << boost::get<0>(pair.second) << std::endl;

        // zero eigenvalues discarded
        if( pair.first < tol )
            continue;

        std::get<0>(modes[i]) = boost::get<0>(pair.second);

        auto alphaHat = boost::get<2>(pair.second);

        // alpha = C*alphaHat
        auto tmpVec = backend()->newVector( Nh );
        C->multVector( alphaHat, tmpVec);
        tmpVec->close();
        std::get<1>(modes[i]) = *tmpVec;

        // g = g0 + grad(psi)
        // decomposition : keep only beta coefficients (psi)
        auto tmpVec2 = backend()->newVector( Lh );
        for(int i = 0; i < boundaryIndexesToKeep.size(); i++)
        {
            auto indexPtr = std::find(indexesToKeep.begin(), indexesToKeep.end(), boundaryIndexesToKeep[i]);
            auto index = indexPtr - indexesToKeep.begin();
            tmpVec2->set(boundaryIndexesToKeep[i] - Nh->nLocalDofWithGhost(), (*alphaHat)(index));
        }
        tmpVec2->close();
        auto tmp = Lh->element();
        tmp = *tmpVec2;
        std::get<2>(modes[i]) = vf::project(_space=Sh, _range=elements(mesh), _expr=idv(tmp));

        // decomposition : keep only alpha coefficients (g0)
        auto tmpVec3 = backend()->newVector( Nh );
        for(int i = 0; i < interiorIndexesToKeep.size(); i++)
        {
            auto indexPtr = std::find(indexesToKeep.begin(), indexesToKeep.end(), interiorIndexesToKeep[i]);
            auto index = indexPtr - indexesToKeep.begin();
            tmpVec3->set(interiorIndexesToKeep[i], (*alphaHat)(index));
        }
        tmpVec3->close();
        modes0[i] = *tmpVec3;

        i++;
    }

    LOG(INFO) << "number of converged eigenmodes : " << i;
    if( Environment::isMasterRank() && ioption("offline.verbose") > 1 )
        std::cout << "number of converged eigenmodes : " << i << std::endl;

    modes.resize(i);
    toc("solve", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::print()
{
    if( Environment::worldComm().isMasterRank() )
        std::cout << "Start printing" << std::endl;

    matA->printMatlab("a.m");
    matB->printMatlab("b.m");
    C->printMatlab("c.m");
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::save()
{
    tic();
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);

    for( int i = 0; i < modes.size(); i++)
    {
        if ( Environment::worldComm().isMasterRank() )
            s << std::get<0>(modes[i]) << std::endl;

        std::string path = (boost::format("mode-%1%")%i).str();
        std::get<1>(modes[i]).save(_path=path, _type=soption("eigen.format"));
        std::string pathP = (boost::format("psi-%1%")%i).str();
        std::get<2>(modes[i]).save(_path=pathP, _type=soption("eigen.format"));
    }

    if ( Environment::worldComm().isMasterRank() )
        s.close();
    toc( "save", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::load()
{
    tic();
    int nbMode = ioption("eigen.nb-mode");
    modes = eigenmodes_type(nbMode, make_tuple(0, Nh->element(), Sh->element()));

    std::fstream s;
    s.open ("lambda", std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found" << std::endl;
        exit(0);
    }

    for( int i = 0; i < nbMode && s.good(); i++ ){
        s >> std::get<0>(modes[i]);
        std::string path = (boost::format("mode-%1%")%i).str();
        std::get<1>(modes[i]).load(_path=path, _type=soption("eigen.format"));
        std::string pathP = (boost::format("psi-%1%")%i).str();
        std::get<2>(modes[i]).load(_path=pathP, _type=soption("eigen.format"));
    }

    s.close();
    toc("load", ioption("offline.verbose") > 2);
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::testEigs()
{
    for( int i = 0; i < modes.size(); i++)
    {
        auto gn = normL2(_range=boundaryfaces(mesh), _expr=trans(idv(std::get<1>(modes[i])))*N());
        auto divg = normL2(_range=elements(mesh), _expr=divv(std::get<1>(modes[i])));
        auto curlgn = normL2(_range=boundaryfaces(mesh), _expr=trans(curlv(std::get<1>(modes[i])))*N());
        auto g0xn = normL2(_range=boundaryfaces(mesh), _expr=cross(idv(modes0[i]),N()));
        auto g0n = normL2(_range=boundaryfaces(mesh), _expr=trans(idv(modes0[i]))*N());
        auto g0 = normL2(_range=boundaryfaces(mesh), _expr=idv(modes0[i]));
        auto e = normL2(_range=elements(mesh), _expr=idv(std::get<1>(modes[i]))-(idv(modes0[i])+trans(gradv(std::get<2>(modes[i])))));
        auto e2 = normL2(_range=elements(mesh), _expr=curlv(std::get<1>(modes[i]))-std::sqrt(std::get<0>(modes[i]))*idv(std::get<1>(modes[i])));
        auto curl2 = integrate(_range=elements(mesh), _expr=trans(curlv(std::get<1>(modes[i])))*curlv(std::get<1>(modes[i]))).evaluate()(0,0);
        auto norm = normL2(_range=elements(mesh), _expr=idv(std::get<1>(modes[i])));
        auto psi = integrate(_range=boundaryfaces(mesh), _expr=idv(std::get<2>(modes[i]))).evaluate()(0,0);
        auto psin = normL2(_range=boundaryfaces(mesh), _expr=gradv(std::get<2>(modes[i]))*N());

        if( Environment::isMasterRank())
        {
            std::cout << "i : " << i << std::endl
                      << "\t ||g.n||      = " << gn << std::endl
                      << "\t ||divg||     = " << divg << std::endl
                      << "\t ||curlg.n||  = " << curlgn << std::endl
                      << "\t ||g0xn||     = " << g0xn << std::endl
                      << "\t ||g0.n||     = " << g0n << std::endl
                      << "\t ||g0||       = " << g0 << std::endl
                      << "\t ||err||      = " << e << std::endl
                      << "\t ||curlg-lg|| = " << e2 << std::endl
                      << "\t (curl,curl)  = " << curl2 << std::endl
                      << "\t ||gi||       = " << norm << std::endl
                      << "\t ||Gpsi.n||   = " << psin << std::endl
                      << "\t int psi      = " << psi << std::endl;
        }
    }
}
