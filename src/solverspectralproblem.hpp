#include <feel/feelalg/vectorblock.hpp>
#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpace1, typename FunctionSpace2>
class SolverSpectralProblem
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace1 vec_space_ptrtype;
    typedef typename vec_space_ptrtype::element_type vec_space_type;
    typedef typename vec_space_type::element_type vec_element_type;

    typedef FunctionSpace2 space_ptrtype;
    typedef typename space_ptrtype::element_type space_type;

    typedef SolverSpectralProblem<vec_space_ptrtype, space_ptrtype> solverspectralproblem_type;
    typedef typename boost::shared_ptr<solverspectralproblem_type> solverspectralproblem_ptrtype;

    typedef typename space_type::template sub_functionspace<0>::type space_edge_type;
    typedef boost::shared_ptr<space_edge_type> space_edge_ptrtype;
    typedef typename space_edge_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<1>::type space_vertex_type;
    typedef boost::shared_ptr<space_vertex_type> space_vertex_ptrtype;

    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    mesh_ptrtype mesh;
    vec_space_ptrtype Vh;
    space_ptrtype Xh;
    space_edge_ptrtype Nh;
    space_vertex_ptrtype Lh;

    vec_element_type a;

    sparse_matrix_ptrtype C;
    std::vector<size_type> indexesToKeep;

    double Re;

public:
    element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const space_ptrtype& Xh);
    void setA(const vec_element_type& a);
    void setEigen(const sparse_matrix_ptrtype _C, const std::vector<size_type> _indexesToKeep);
    element_type solve();
};

template<typename F, typename E>
typename SolverSpectralProblem<F,E>::solverspectralproblem_ptrtype
SolverSpectralProblem<F,E>::build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const space_ptrtype& Xh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F,E> );
    ap->mesh = mesh;
    ap->Vh = Vh;
    ap->Xh = Xh;
    ap->Nh = Xh->template functionSpace<0>();
    ap->Lh = Xh->template functionSpace<1>();
    return ap;
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::setA(const vec_element_type& a)
{
    this->a = a;
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::setEigen(const sparse_matrix_ptrtype _C, const std::vector<size_type> _indexesToKeep)
{
    this->C = _C;
    this->indexesToKeep = _indexesToKeep;
}


template<typename F, typename E>
typename SolverSpectralProblem<F,E>::element_type
SolverSpectralProblem<F,E>::solve()
{
    boost::mpi::timer t;
    u = Nh->element();

    auto P0 = Pch<0>(mesh);
    auto lambda  = P0->element();
    auto v = Nh->element();
    auto p = Lh->element();

    // Block A
    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(v))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();
    // Block B
    auto b = form2(_test=Lh, _trial=Nh);
    b = integrate( _range=elements(mesh), _expr=grad(p)*idt(v));
    auto matB = b.matrixPtr();
    matB->close();
    // Block Bt
    auto bt = form2(_test=Nh, _trial=Lh);
    bt = integrate( _range=elements(mesh), _expr=gradt(p)*id(v));
    auto matBt = bt.matrixPtr();
    matBt->close();
    // Block Ql
    auto ql = form2(_test=Lh, _trial=P0);
    ql = integrate( _range=elements(mesh), _expr=id(p)*idt(lambda));
    auto matQl = ql.matrixPtr();
    matQl->close();
    // Block Qlt
    auto qlt = form2(_test=P0, _trial=Lh);
    qlt = integrate( _range=elements(mesh), _expr=idt(p)*id(lambda));
    auto matQlt = qlt.matrixPtr();
    matQlt->close();

    // Matrix AA = ( A  Bt  0  )
    //             ( B  0   Ql )
    //             ( 0  Qlt 0  )
    BlocksBaseSparseMatrix<double> aaBlock(3,3);
    aaBlock(0,0) = matA;
    aaBlock(0,1) = matBt;
    aaBlock(1,0) = matB;
    aaBlock(1,2) = matQl;
    aaBlock(2,1) = matQlt;
    auto AA = backend()->newBlockMatrix(_block=aaBlock, _copy_values=true);
    AA->close();

    // Vector L = <alpha2,phi>_Gamma
    auto psi = Lh->element();
    auto alpha2 = expr(soption("solverns2.alpha2"));

    // vector in Nh x Lh
    auto ll = backend()->newVector( Xh );
    auto l = form1(_test=Xh, _vector=ll);
    l = integrate( _range=markedfaces(mesh,1), _expr=alpha2*id(psi));
    l += integrate( _range=markedfaces(mesh,2), _expr=-alpha2*id(psi));
    // vector in Nh(Omega\Gamma) x Lh(Gamma)
    auto L = ll->createSubVector(indexesToKeep);

    // Vector LL = ( L )   Zh
    //             ( 0 )   Lh
    //             ( 0 )   P0
    BlocksBaseVector<double> lBlock(3);
    lBlock(0) = L;
    lBlock(1) = backend()->newVector( _test=Lh );
    lBlock(2) = backend()->newVector( _test=P0 );
    auto LL = backend()->newBlockVector(_block=lBlock, _copy_values=true);

    // Matrix CC = ( C  0  )
    //             ( 0  Id )
    BlocksBaseSparseMatrix<double> ccBlock(3,3);
    ccBlock(0,0) = C;
    auto ccBlock11 = backend()->newMatrix( _test=Lh, _trial=Lh );
    auto d = backend()->newVector(Lh);
    d->setOnes();
    d->close();
    backend()->diag(d, ccBlock11);
    ccBlock(1,1) = ccBlock11;
    auto d2 = backend()->newMatrix( _test=P0, _trial=P0 );
    d2->set(0,0,1);
    ccBlock(2,2) = d2;
    auto CC = backend()->newBlockMatrix(_block=ccBlock, _copy_values=true);
    CC->close();

    // Matrix CC'*AA*CC
    auto aaHat = backend()->newMatrix(CC->mapColPtr(), CC->mapColPtr() );
    backend()->PtAP( AA, CC, aaHat );
    aaHat->close();

    if( Environment::numberOfProcessors() == 1 )
    {
        C->printMatlab("c.m");
        matA->printMatlab("a.m");
        matB->printMatlab("b.m");
        matBt->printMatlab("bt.m");
        AA->printMatlab("aa.m");
        CC->printMatlab("cc.m");
        L->printMatlab("l.m");
        LL->printMatlab("ll.m");
        aaHat->printMatlab("ah.m");
    }

    logTime(t, "matrices", ioption("solverns2.verbose") > 1);

    // solve aHat X = L (with precontditionner ?)
    auto s = LL;
    backend(_name="sp")->solve(_matrix=aaHat, _solution=s, _rhs=LL);

    logTime(t, "solve", ioption("solverns2.verbose") > 1);

    // vector C*s in Xh=Nh x Lh
    // pb tmpVec in Xh x P0
    BlocksBaseSparseMatrix<double> c2Block(2,3);
    c2Block(0,0) = C;
    c2Block(1,1) = ccBlock11;
    c2Block(0,2) = backend()->newMatrix( _test=P0, _trial=P0 );
    auto C2 = backend()->newBlockMatrix(_block=c2Block, _copy_values=true);
    C2->close();
    auto tmpVec = backend()->newVector( Xh );
    C2->multVector( s, tmpVec);
    tmpVec->close();
    auto U = Xh->element();
    U = *tmpVec;
    u = U.template element<0>();

    auto unG = form1(_test=Lh);
    unG = integrate(boundaryfaces(mesh), trans(idv(u))*N()*id(p));
    unG.close();
    auto unG1 = form1(_test=Lh);
    unG1 = integrate(markedfaces(mesh,1), trans(idv(u))*N()*id(p));
    unG1.close();
    auto unG2 = form1(_test=Lh);
    unG2 = integrate(markedfaces(mesh,2), trans(idv(u))*N()*id(p));
    unG2.close();
    auto unG3 = form1(_test=Lh);
    unG3 = integrate(markedfaces(mesh,3), trans(idv(u))*N()*id(p));
    unG3.close();

    if( Environment::numberOfProcessors() == 1 )
    {
        unG.vectorPtr()->printMatlab("unG.m");
        unG1.vectorPtr()->printMatlab("unG1.m");
        unG2.vectorPtr()->printMatlab("unG2.m");
        unG3.vectorPtr()->printMatlab("unG3.m");
    }

    auto du = normL2(elements(mesh), divv(u));
    auto intun = integrate(boundaryfaces(mesh),
                           trans(idv(u))*N()).evaluate()(0,0);
    auto nu = normL2(boundaryfaces(mesh), trans(idv(u))*N());
    auto cnu = normL2(boundaryfaces(mesh), trans(curlv(u))*N());
    if( Environment::isMasterRank() )
        std::cout << "||div u|| = " << du << std::endl
                  << "||u . n|| = " << nu << std::endl
                  << "int u.n   = " << intun <<  std::endl
                  << "||cu. n|| = " << cnu << std::endl;

    return u;
}
