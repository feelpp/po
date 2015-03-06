#include <Eigen/Dense>

#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpaceType1, typename FunctionSpaceType2>
class SolverSpectralProblem
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpaceType1 ned_space_ptrtype;
    typedef typename ned_space_ptrtype::element_type ned_space_type;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef FunctionSpaceType2 vec_space_ptrtype;
    typedef typename vec_space_ptrtype::element_type vec_space_type;
    typedef typename vec_space_type::element_type vec_element_type;

    typedef SolverSpectralProblem<ned_space_ptrtype, vec_space_ptrtype> solverspectralproblem_type;
    typedef typename boost::shared_ptr<solverspectralproblem_type> solverspectralproblem_ptrtype;

    typedef std::pair<value_type, ned_element_type> eigenpair_type;
    typedef std::vector<eigenpair_type> eigenmodes_type;
    typedef std::vector<ned_element_type> eigenvec_type;

    mesh_ptrtype mesh;
    ned_space_ptrtype Nh;
    vec_space_ptrtype Vh;
    vec_element_type a0;
    eigenvec_type g;
    VectorXd lambda;

    int M;
    MatrixXd j;
    VectorXd f;

    double Re;

    // [ri]
    Matrix<MatrixXd, Dynamic, 1 > Rijk;
    MatrixXd Riak;
    MatrixXd Rik;
    VectorXd Rfk;
    // [ri]

    void initRijk();
    void initRiak();
    void initRfk();


public:
    VectorXd c;
    vec_element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const ned_space_ptrtype& Nh, const vec_space_ptrtype& Vh);
    void setA0(const vec_element_type& a0);
    void setEigen(const eigenmodes_type& modes);
    void init();
    vec_element_type solve();
};

template<typename T1, typename T2>
typename SolverSpectralProblem<T1,T2>::solverspectralproblem_ptrtype
SolverSpectralProblem<T1,T2>::build(const mesh_ptrtype& mesh, const ned_space_ptrtype& Nh, const vec_space_ptrtype& Vh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<T1,T2> );
    ap->mesh = mesh;
    ap->Nh = Nh;
    ap->Vh = Vh;
    return ap;
}

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::setA0(const vec_element_type& a0)
{
    this->a0 = a0;
}

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::setEigen(const eigenmodes_type& modes)
{
    M = modes.size();
    lambda = VectorXd(M);

    g = eigenvec_type(M);

    int i = 0;
    for(auto const& pair : modes)
    {
        lambda(i) = pair.first;
        g[i++] = pair.second;
    }
}

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::init()
{
    j = MatrixXd(M,M);
    f = VectorXd(M);

    double r = doption( _name="radius" );
    double s = doption( _name="speed" );
    double n = doption( _name="nu" );
    Re = 2*r*s/n;
    u = Vh->element();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Initialization Spectral Problem -----" << std::endl
                  << "----- Re = " << Re << " -----" << std::endl;

    boost::mpi::timer t;

    initRijk();
    logTime(t, "Rijk", FLAGS_v > 1);
    // initRiak();
    // logTime(t, "Riak", FLAGS_v > 1);
    initRfk();
    logTime(t, "Rfk", FLAGS_v > 1);

    c = VectorXd::Ones(M);
}

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::initRijk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;

    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

    if( boption("computeRijk") )
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

                    if ( Environment::worldComm().isMasterRank() && FLAGS_v > 2 )
                    {
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

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::initRiak() // to test i<->k mult ???
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Riak -----" << std::endl;

    Riak = MatrixXd(M,M);

    if( boption("computeRiak") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("riak", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                Riak(k,i) = integrate( _range=elements( mesh ),
                                       _expr=inner( cross( idv(g[i]),idv(a0) ), idv(g[k])) ).evaluate()(0,0);
                Riak(i,k) = -Riak(k,i);

                if ( Environment::worldComm().isMasterRank() && FLAGS_v > 2 )
                {
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

template<typename T1, typename T2>
void
SolverSpectralProblem<T1,T2>::initRfk()
{
    if ( Environment::worldComm().isMasterRank() )
        std::cout << "----- Rfk -----" << std::endl;

    Rfk = VectorXd(M);

    if( boption("computeRfk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rfk", std::fstream::out);

        auto ff = expr<3,1>(soption(_name="f"));
        for(int k = 0; k < M; k++)
        {
            Rfk(k) = integrate( _range=elements( mesh ),
                                _expr=trans(ff)*idv(g[k]) ).evaluate()(0,0);

            if ( Environment::worldComm().isMasterRank() && FLAGS_v > 2)
            {
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

template<typename T1, typename T2>
typename SolverSpectralProblem<T1,T2>::vec_element_type
SolverSpectralProblem<T1,T2>::solve()
{
    boost::mpi::timer t;

    // [StokesA]
    MatrixXd A = MatrixXd(M,M);
    A = Riak;
    A += (lambda/Re).asDiagonal();
    // [StokesA]

    // [StokesB]
    VectorXd b = VectorXd(M);
    b = Rfk;
    // [StokesB]

    // [StokesSolve]
    HouseholderQR<MatrixXd> qr(M,M);
    qr.compute(A);
    c = qr.solve(b);
    // [StokesSolve]

    logTime(t, "stokes", FLAGS_v > 1);

    for(int i=0; i<M; i++)
    {
        if ( Environment::worldComm().isMasterRank() && FLAGS_v > 2 )
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
