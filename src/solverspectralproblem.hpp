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

    // Block matrices for AA
    auto v = Nh->element();
    auto p = Lh->element();
    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(v))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();
    auto b = form2(_test=Lh, _trial=Nh);
    b = integrate( _range=elements(mesh), _expr=divt(v)*id(p));
    auto matB = b.matrixPtr();
    matB->close();
    auto bt = form2(_test=Nh, _trial=Lh);
    bt = integrate( _range=elements(mesh), _expr=div(v)*idt(p));
    auto matBt = bt.matrixPtr();
    matBt->close();

    // Matrix AA = ( A  Bt )
    //             ( B  0  )
    BlocksBaseSparseMatrix<double> aaBlock(2,2);
    aaBlock(0,0) = matA;
    aaBlock(0,1) = matBt;
    aaBlock(1,0) = matB;
    aaBlock(1,1) = backend()->newMatrix( _test=Lh, _trial=Lh );
    auto AA = backend()->newBlockMatrix(_block=aaBlock, _copy_values=true);
    AA->close();

    // Vector L = <alpha2,phi>_Gamma
    auto psi = Lh->element();
    auto alpha2 = expr(soption("solverns2.alpha2"));

    // vector in Nh x Lh
    auto ll = backend()->newVector( Xh );
    auto l = form1(_test=Xh, _vector=ll);
    l = integrate( _range=boundaryfaces(mesh), _expr=alpha2*id(psi));
    // vector in Nh(Omega\Gamma) x Lh(Gamma)
    auto L = ll->createSubVector(indexesToKeep);

    // Vector LL = ( L )
    //             ( 0 )
    BlocksBaseVector<double> lBlock(2);
    lBlock(0) = L;
    lBlock(1) = backend()->newVector( _test=Lh );
    auto LL = backend()->newBlockVector(_block=lBlock, _copy_values=true);

    // Matrix CC = ( C  0  )
    //             ( 0  Id )
    BlocksBaseSparseMatrix<double> ccBlock(2,2);
    ccBlock(0,0) = C;//backend()->newMatrix( Nh->dof(), C->mapColPtr() );
    auto ccBlock11 = backend()->newMatrix( _test=Lh, _trial=Lh );
    auto d = backend()->newVector(Lh);
    d->setOnes();
    d->close();
    backend()->diag(d, ccBlock11);
    ccBlock(1,1) = ccBlock11;
    ccBlock.close();
    auto CC = backend()->newBlockMatrix(_block=ccBlock, _copy_values=true);
    CC->close();

    auto aaHat = backend()->newMatrix(CC->mapColPtr(), CC->mapColPtr() );
    backend()->PtAP( AA, CC, aaHat );
    aaHat->close();

    logTime(t, "matrices", ioption("solverns2.verbose") > 1);

    // solve aHat X = L (with precontditionner ?)
    auto s = LL;
    backend(_name="sp")->solve(_matrix=aaHat, _solution=s, _rhs=LL);
    auto tmpVec = backend()->newVector( Xh );
    CC->multVector( s, tmpVec);
    tmpVec->close();
    auto U = Xh->element();
    U = *tmpVec;
    u = U.template element<0>();
    return u;
}
