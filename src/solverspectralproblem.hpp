#include <Eigen/Dense>

#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpace, typename EigenTuple>
class SolverSpectralProblem
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace vec_space_ptrtype;
    typedef typename vec_space_ptrtype::element_type vec_space_type;
    typedef typename vec_space_type::element_type vec_element_type;

    typedef EigenTuple eigentuple_type;
    typedef typename std::tuple_element<1, eigentuple_type>::type ned_element_type;
    typedef typename std::tuple_element<2, eigentuple_type>::type scalar_element_type;

    typedef SolverSpectralProblem<vec_space_ptrtype, eigentuple_type> solverspectralproblem_type;
    typedef typename boost::shared_ptr<solverspectralproblem_type> solverspectralproblem_ptrtype;

    typedef std::vector<eigentuple_type> eigenmodes_type;
    typedef std::vector<ned_element_type> eigenvec_type;
    typedef std::vector<scalar_element_type> psivec_type;

    mesh_ptrtype mesh;
    vec_space_ptrtype Vh;

    vec_element_type a;
    eigenvec_type g;
    psivec_type psi;
    VectorXd lambda;

    int M;
    MatrixXd j;
    VectorXd f;

    double Re;

    // [ri]
    Matrix<MatrixXd, Dynamic, 1 > Rijk;
    MatrixXd Riak;
    VectorXd Rfk;
    VectorXd Rpk;
    // [ri]

    void initRijk();
    void initRiak();
    void initRfk();
    void initRpk();


public:
    VectorXd c;
    vec_element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh);
    void setA(const vec_element_type& a);
    void setEigen(const eigenmodes_type& modes);
    void init();
    vec_element_type solve();
};

template<typename F, typename E>
typename SolverSpectralProblem<F,E>::solverspectralproblem_ptrtype
SolverSpectralProblem<F,E>::build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F,E> );
    ap->mesh = mesh;
    ap->Vh = Vh;
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
SolverSpectralProblem<F,E>::setEigen(const eigenmodes_type& modes)
{
    M = modes.size();
    lambda = VectorXd(M);

    g = eigenvec_type(M);
    psi = psivec_type(M);

    int i = 0;
    for(auto const& tuple : modes)
    {
        lambda(i) = std::get<0>(tuple);
        g[i] = std::get<1>(tuple);
        psi[i++] = std::get<2>(tuple);
    }
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::init()
{
    j = MatrixXd(M,M);
    f = VectorXd(M);

    double r = doption( "solverns2.radius" );
    double s = doption( "solverns2.speed" );
    double n = doption( "solverns2.nu" );
    Re = 2*r*s/n;
    u = Vh->element();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Initialization Spectral Problem -----" << std::endl
                  << "----- Re = " << Re << " -----" << std::endl;

    boost::mpi::timer t;

    // initRijk();
    // logTime(t, "Rijk", ioption("solverns2.verbose") > 1);
    initRiak();
    logTime(t, "Riak", ioption("solverns2.verbose") > 1);
    initRfk();
    logTime(t, "Rfk", ioption("solverns2.verbose") > 1);
    initRpk();
    logTime(t, "Rpk", ioption("solverns2.verbose") > 1);

    c = VectorXd::Ones(M);
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::initRijk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;

    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

    if( boption("solverns2.computeRijk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rijk", std::fstream::out);

        // Rijk = - Rikj && Rijj = 0
        for(int k = 0; k < M; k++)
        {
            for(int i = 0; i < M; i++)
            {
                for(int j = 0; j < k; j++)
                {
                    Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                                              _expr=inner( cross( curlv(g[i]),idv(g[j]) ), idv(g[k]) )
                                              ).evaluate()(0,0);
                    Rijk(j)(i,k) = -Rijk(k)(i,j);

                    if ( Environment::worldComm().isMasterRank() )
                    {
                        if( ioption("solverns2.verbose") > 2 )
                            std::cout << "Rijk(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
                        s << Rijk(k)(i,j) << std::endl;
                    }
                }
                Rijk(k)(i,k) = 0;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else
    {
        std::fstream s;
        s.open ("rijk", std::fstream::in);
        if( !s.is_open() )
        {
            std::cout << "Rijk not found\ntry to launch with --computeRijk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
        {
            for(int i = 0; i < M; i++)
            {
                for(int j = 0; j < k; j++)
                {
                    s >> Rijk(k)(i,j);
                    Rijk(j)(i,k) = -Rijk(k)(i,j);
                }
                Rijk(k)(i,k) = 0;
            }
        }
        s.close();
    }
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::initRiak()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Riak -----" << std::endl;

    Riak = MatrixXd(M,M);

    if( boption("solverns2.computeRiak") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("riak", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                Riak(k,i) = integrate( _range=elements( mesh ),
                                       _expr=inner( cross( idv(g[i]),idv(a) ), idv(g[k])) ).evaluate()(0,0);
                Riak(i,k) = -Riak(k,i);

                if ( Environment::worldComm().isMasterRank() )
                {
                    if( ioption("solverns2.verbose") > 2 )
                        std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
                    s << Riak(k,i) << std::endl;
                }
            }
            Riak(i,i) = 0;
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else
    {
        std::fstream s;
        s.open ("riak", std::fstream::in);
        if( !s.is_open() )
        {
            std::cout << "Riak not found\ntry to launch with --computeRiak=true" << std::endl;
            exit(0);
        }
        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                s >> Riak(k,i);
                Riak(i,k) = -Riak(k,i);
            }
            Riak(i,i) = 0;
        }
        s.close();
    }
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::initRfk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rfk -----" << std::endl;

    Rfk = VectorXd(M);

    if( boption("solverns2.computeRfk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rfk", std::fstream::out);

        auto ff = expr<3,1>(soption("solverns2.f"));
        for(int k = 0; k < M; k++)
        {
            Rfk(k) = integrate( _range=elements( mesh ),
                                _expr=trans(ff)*idv(g[k]) ).evaluate()(0,0);

            if ( Environment::worldComm().isMasterRank() )
            {
                if( ioption("solverns2.verbose") > 2)
                    std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
                s << Rfk(k) << std::endl;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else
    {
        std::fstream s;
        s.open ("rfk", std::fstream::in);
        if( !s.is_open() )
        {
            std::cout << "Rfk not found\ntry to launch with --computeRfk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
            s >> Rfk(k);
        s.close();
    }
}

template<typename F, typename E>
void
SolverSpectralProblem<F,E>::initRpk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rpk -----" << std::endl;

    Rpk = VectorXd(M);

    if( boption("solverns2.computeRpk") )
    {
        auto a2 = expr(soption("solverns2.alpha2"));

        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rpk", std::fstream::out);

        for(int k = 0; k < M; k++)
        {
            Rpk(k) = integrate( _range=markedfaces( mesh, 1 ),
                                _expr=-a2*idv(psi[k]) ).evaluate()(0,0);
            Rpk(k) += integrate( _range=markedfaces( mesh, 2 ),
                                 _expr=a2*idv(psi[k]) ).evaluate()(0,0);

            if ( Environment::worldComm().isMasterRank())
            {
                if( ioption("solverns2.verbose") > 2 )
                    std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
                s << Rpk(k) << std::endl;
            }
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else
    {
        std::fstream s;
        s.open ("rpk", std::fstream::in);
        if( !s.is_open() )
        {
            std::cout << "Rpk not found\ntry to launch with --computeRpk=true" << std::endl;
            exit(0);
        }
        for(int k = 0; k < M; k++)
            s >> Rpk(k);
        s.close();
    }
}

template<typename F, typename E>
typename SolverSpectralProblem<F,E>::vec_element_type
SolverSpectralProblem<F,E>::solve()
{
    boost::mpi::timer t;

    // [StokesA]
    MatrixXd A = MatrixXd(M,M);
    A = Riak;
    A += (lambda/Re).asDiagonal();
    // [StokesA]

    // [StokesB]
    VectorXd b = VectorXd(M);
    b = Rfk + Rpk;
    // [StokesB]

    // [StokesSolve]
    HouseholderQR<MatrixXd> qr(M,M);
    qr.compute(A);
    c = qr.solve(b);
    // [StokesSolve]

    logTime(t, "stokes", ioption("solverns2.verbose") > 1);

    for(int i=0; i<M; i++)
    {
        if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 2 )
            std::cout << "c(" << i << ") = " << c(i) << std::endl;
        u += vf::project( _space=Vh, _range=elements(mesh),
                          _expr = c(i)*idv(g[i]) );
    }

    return u;

    // // [NSInit]
    // VectorXd dc = VectorXd::Matrix(M);
    // double tol = 1.e-6;
    // // [NSInit]
    // // [NSSys1]
    // HouseholderQR<MatrixXd> qr(M,M);
    // // [NSSys1]

    // int i=0;
    // do{
    //     // [NSMatF]
    //     f = c.cwiseProduct(lambda)/Re + Riak*c - Rfk;
    //     for (int k = 0; k < M; k++)
    //         f(k) += c.transpose()*Rijk(k)*c;
    //     // [NSMatF]

    //     // [NSMatJ]
    //     for (int k = 0; k < M; k++)
    //         j.row(k) = c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
    //     j += Riak;
    //     j += lambda.asDiagonal();
    //     // [NSMatJ]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "j = " << j << std::endl << "f = " << f << std::endl;

    //     // [NSSys2]
    //     qr.compute(j);
    //     dc = qr.solve(-f);
    //     // [NSSys2]
    //     // [NSAdd]
    //     c += dc;
    //     // [NSAdd]

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << "iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

    //     if ( Environment::worldComm().isMasterRank() )
    //         std::cout << c << std::endl;

    //     i++;
    // } while(i < 10 && dc.norm() > tol);

    // if(i==10)
    //     std::cout << "Newton does not converge\n";

    // else{
    //     for( i = 0; i < M; i++){
    //         u += vf::project( _space=Vh, _range=elements(mesh),
    //                           _expr = c(i)*idv(g[i]) );
    //     }
    //     u += a;
    // }
}
