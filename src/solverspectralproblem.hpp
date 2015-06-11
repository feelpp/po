#include <Eigen/Dense>

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

    vec_element_type a;

    sparse_matrix_ptrtype C;
    std::vector<size_type> indexesToKeep;

    double Re;

public:
    element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const space_ptrtype& Xh);
    void setA(const vec_element_type& a);
    void setEigen(const sparse_matrix_ptrtype _C, const std::vector<size_type> _indexesToKeep);
    vec_element_type solve();
};

template<typename F, typename E>
typename SolverSpectralProblem<F,E>::solverspectralproblem_ptrtype
SolverSpectralProblem<F,E>::build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const space_ptrtype& Xh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F,E> );
    ap->mesh = mesh;
    ap->Vh = Vh;
    ap->Xh = Xh;
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
typename SolverSpectralProblem<F,E>::vec_element_type
SolverSpectralProblem<F,E>::solve()
{
    boost::mpi::timer t;

    auto U = Xh->element();
    auto psi = U.template element<1>();
    auto alpha2 = expr(soption("solverns2.alpha2"));

    auto ll = backend()->newVector( Xh );
    auto l = form1(_test=Xh, _vector=ll);
    l = integrate( _range=boundaryfaces(mesh), _expr=alpha2*id(psi));

    auto L = ll->createSubVector(indexesToKeep);

    auto v = Nh->element();
    auto a = form2( _test=Nh, _trial=Nh);
    a = integrate( _range=elements( mesh ), _expr=trans(curlt(v))*curl(v));
    auto matA = a.matrixPtr();
    matA->close();
    auto aHat = backend()->newMatrix(C->mapColPtr(), C->mapColPtr() );
    backend()->PtAP( matA, C, aHat );
    aHat->close();

    //solve aHat X = L (with precontditionner ?)
}
