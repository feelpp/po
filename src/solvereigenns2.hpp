#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/solvereigen.hpp>

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

    double tol = 1.e-3;

    void setForms();
    void setInfo();
    void setDofsToRemove();
    void setMatrices();
    void solveEigen();
    void print();
    void save();
    void load();
    void logInfo();

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

    if( boption("solverns2.computeEigen"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute eigenmodes ----------\n";

        setForms();

        setInfo();
        setDofsToRemove();

        setMatrices();

        if( boption("solverns2.print") )
            print();

        toc( "matrices", ioption("solverns2.verbose") > 1);
        tic();

        solveEigen();
        toc( "eigs", ioption("solverns2.verbose") > 1);
        tic();

        save();
        toc( "save", ioption("solverns2.verbose") > 1);
    }
    else
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- load eigenmodes ----------\n";

        load();
        toc( "loadEigen", ioption("solverns2.verbose") > 1);
    }


    return modes;
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setForms()
{
    // [forms]
    auto u = Nh->element();
    auto v = Nh->element();

    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(u))*curl(v));
    matA = a.matrixPtr();
    matA->close();

    auto b = form2( _test=Nh, _trial=Nh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    matB = b.matrixPtr();
    matB->close();
    // [forms]
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setInfo()
{
    std::vector<bool> doneNh( Nh->nLocalDof(), false );
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
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setDofsToRemove()
{
    // [remove]
    auto markers = Environment::vm()["solverns2.markerList"].as<std::vector<std::string> >();
    auto s = markers.size();
    auto dofsToRemove = std::vector<int>();

    for( auto it = markers.begin(); it != markers.end(); ++it )
    {
        std::string m = *it;
        if( mesh->hasMarker(m) )
        {
            auto r = Lh->dof()->markerToDof( m );
            auto localDof = r.first->second;
            auto map = Lh->dof()->mapGlobalProcessToGlobalCluster();
            auto globalDof = map[localDof];
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

                std::cout << "#" << Environment::worldComm().globalRank() << " remove index : " << localDofToRemove << std::endl;
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
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::setMatrices()
{
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
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::solveEigen()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Eigs : nev = " << ioption(_name="solvereigen.nev")
                  << "\t ncv= " << ioption(_name="solvereigen.ncv") << std::endl;

    auto zhmodes = eigs( _matrixA=aHat,
                         _matrixB=bHat,
                         _solver=(EigenSolverType)EigenMap[soption("solvereigen.solver")],
                         _problem=(EigenProblemType)EigenMap[soption("solvereigen.problem")],
                         _transform=(SpectralTransformType)EigenMap[soption("solvereigen.transform")],
                         _spectrum=(PositionOfSpectrum)EigenMap[soption("solvereigen.spectrum")]
                         );

    int i = 0;
    modes = eigenmodes_type(zhmodes.size(), std::make_tuple(0, Nh->element(), Sh->element()) );
    for( auto const& pair: zhmodes )
    {
        if(Environment::isMasterRank())
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

        // decomposition : keep only beta coefficients
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

        i++;
    }
    LOG(INFO) << "number of converged eigenmodes : " << i;
    if( Environment::isMasterRank() )
        std::cout << "number of converged eigenmodes : " << i << std::endl;

    modes.resize(i);
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
    std::fstream s;
    if ( Environment::worldComm().isMasterRank() )
        s.open ("lambda", std::fstream::out);

    for( int i = 0; i < modes.size(); i++)
    {
        if ( Environment::worldComm().isMasterRank() )
            s << std::get<0>(modes[i]) << std::endl;

        std::string path = (boost::format("mode-%1%")%i).str();
        std::get<1>(modes[i]).save(_path=path);
        std::string pathP = (boost::format("psi-%1%")%i).str();
        std::get<2>(modes[i]).save(_path=pathP);
    }

    if ( Environment::worldComm().isMasterRank() )
        s.close();
}

template<typename T1, typename T2>
void
SolverEigenNS2<T1,T2>::load()
{
    int nbMode = ioption("solverns2.nbMode");
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
        std::get<1>(modes[i]).load(_path=path);
        std::string pathP = (boost::format("psi-%1%")%i).str();
        std::get<2>(modes[i]).load(_path=pathP);
    }

    s.close();
}
