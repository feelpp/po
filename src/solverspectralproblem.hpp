#include <Eigen/Dense>

#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpace, typename FunctionSpace2, typename FunctionSpace3, typename EigenTuple>
class SolverSpectralProblem
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace vec_space_ptrtype;
    typedef typename vec_space_ptrtype::element_type vec_space_type;
    typedef typename vec_space_type::element_type vec_element_type;

    typedef FunctionSpace2 ned_space_ptrtype;
    typedef typename ned_space_ptrtype::element_type ned_space_type;
    typedef typename ned_space_type::element_type ned_element_type;

    typedef FunctionSpace3 scalar_space_ptrtype;
    typedef typename scalar_space_ptrtype::element_type scalar_space_type;
    typedef typename scalar_space_type::element_type scalar_element_type;

    typedef EigenTuple eigentuple_type;

    typedef SolverSpectralProblem<vec_space_ptrtype, ned_space_ptrtype, scalar_space_ptrtype, eigentuple_type> solverspectralproblem_type;
    typedef typename boost::shared_ptr<solverspectralproblem_type> solverspectralproblem_ptrtype;

    typedef std::vector<eigentuple_type> eigenmodes_type;
    typedef std::vector<ned_element_type> eigenvec_type;
    typedef std::vector<scalar_element_type> psivec_type;

    mesh_ptrtype mesh;
    vec_space_ptrtype Vh;
    ned_space_ptrtype Nh;
    scalar_space_ptrtype Sh;

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
    MatrixXd Raik;
    MatrixXd Riak;
    VectorXd Rfk;
    VectorXd Rpk;
    // [ri]

    void initRijk();
    void initRaik();
    void initRiak();
    void initRfk();
    void initRpk();


public:
    VectorXd c;
    vec_element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const ned_space_ptrtype& Nh, const scalar_space_ptrtype& Sh);
    void setA(const vec_element_type& a);
    void setEigen(const eigenmodes_type& modes);
    void init();
    vec_element_type solve();
};

template<typename F, typename F2, typename F3, typename E>
typename SolverSpectralProblem<F,F2,F3,E>::solverspectralproblem_ptrtype
SolverSpectralProblem<F,F2,F3,E>::build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const ned_space_ptrtype& Nh, const scalar_space_ptrtype& Sh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F,F2,F3,E> );
    ap->mesh = mesh;
    ap->Vh = Vh;
    ap->Nh = Nh;
    ap->Sh = Sh;
    return ap;
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::setA(const vec_element_type& a)
{
    this->a = a;
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::setEigen(const eigenmodes_type& modes)
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

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::init()
{
    j = MatrixXd(M,M);
    f = VectorXd(M);

    double r = doption( "solverns2.radius" );
    double s = doption( "solverns2.speed" );
    double n = doption( "solverns2.nu" );
    Re = 2*r*s/n;
    u = Vh->element();

    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Initialization Spectral Problem -----" << std::endl
                  << "----- Re = " << Re << " -----" << std::endl;

    if( ! boption("solverns2.stokes"))
    {
        tic();
        initRijk();
        toc( "Rijk", ioption("solverns2.verbose") > 1);
        tic();
        initRaik();
        toc( "Raik", ioption("solverns2.verbose") > 1);
        tic();
        initRiak();
        toc( "Riak", ioption("solverns2.verbose") > 1);
    }
    // tic();
    // initRfk();
    // toc( "Rfk", ioption("solverns2.verbose") > 1);
    tic();
    initRpk();
    toc( "Rpk", ioption("solverns2.verbose") > 1);

    c = VectorXd::Ones(M);
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::initRijk()
{
    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Rijk -----" << std::endl;

    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

    auto w = Nh->element();
    auto rijkForm = form2( _test=Nh, _trial=Nh );

    if( boption("solverns2.computeRijk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rijk", std::fstream::out);

        // Rijk = - Rikj && Rijj = 0
        for(int k = 0; k < M; k++)
        {
            rijkForm = integrate( elements(mesh),
                                  inner( cross( curl(w),idt(w) ), idv(g[k]) ));
            for(int i = 0; i < M; i++)
            {
                for(int j = 0; j < k; j++)
                {
                    //Rijk(k)(i,j) = integrate( _range=elements( mesh ),
                    //                          _expr=inner( cross( curlv(g[i]),idv(g[j]) ), idv(g[k]) )
                    //).evaluate()(0,0);
                    Rijk(k)(i,j) = rijkForm( g[i], g[j] );
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

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::initRaik()
{
    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Raik -----" << std::endl;

    Raik = MatrixXd(M,M);

    auto w = Nh->element();
    auto raikForm = form2( _trial=Nh, _test=Nh );
    raikForm = integrate( elements(mesh), inner( cross( curlv(a),id(w) ), idt(w)) );

    if( boption("solverns2.computeRaik") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("raik", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                //Raik(k,i) = integrate( _range=elements( mesh ),
                //                       _expr=inner( cross( curlv(a),idv(g[i]) ), idv(g[k])) ).evaluate()(0,0);
                Raik(k,i) = raikForm(g[i], g[k]);
                Raik(i,k) = -Raik(k,i);

                if ( Environment::worldComm().isMasterRank() )
                {
                    if( ioption("solverns2.verbose") > 2 )
                        std::cout << "Raik(" << k << "," << i << ") = " << Raik(k,i) << std::endl;
                    s << Raik(k,i) << std::endl;
                }
            }
            Raik(i,i) = 0;
        }
        if ( Environment::worldComm().isMasterRank() )
            s.close();
    }
    else
    {
        std::fstream s;
        s.open ("raik", std::fstream::in);
        if( !s.is_open() )
        {
            std::cout << "Raik not found\ntry to launch with --computeRaik=true" << std::endl;
            exit(0);
        }
        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                s >> Raik(k,i);
                Raik(i,k) = -Raik(k,i);
                if( ioption("solverns2.verbose") > 2 )
                    std::cout << "Raik(" << k << "," << i << ") = " << Raik(k,i) << std::endl;
            }
            Raik(i,i) = 0;
        }
        s.close();
    }
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::initRiak()
{
    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Riak -----" << std::endl;

    Riak = MatrixXd(M,M);

    auto w = Nh->element();
    auto riakForm = form2( _test=Nh, _trial=Nh );
    riakForm = integrate( elements(mesh), inner( cross( curl(w),idv(a) ), idt(w)) );
    if( boption("solverns2.computeRiak") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("riak", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < M; k++)
            {
                //Riak(k,i) = integrate( _range=elements( mesh ),
                //                       _expr=inner( cross( curlv(g[i]),idv(a) ), idv(g[k])) ).evaluate()(0,0);
                Riak(k,i) = riakForm( g[i], g[k] );

                if ( Environment::worldComm().isMasterRank() )
                {
                    if( ioption("solverns2.verbose") > 2 )
                        std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
                    s << Riak(k,i) << std::endl;
                }
            }
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
            for(int k = 0; k < M; k++)
            {
                s >> Riak(k,i);
                if( ioption("solverns2.verbose") > 2 )
                    std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
            }
        }
        s.close();
    }
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::initRfk()
{
    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Rfk -----" << std::endl;

    Rfk = VectorXd(M);

    auto ff = expr<3,1>(soption("solverns2.f"));
    auto w = Nh->element();
    auto rfkForm = form1( _test=Nh );
    rfkForm = integrate( elements(mesh), trans(ff)*id(w));

    if( boption("solverns2.computeRfk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rfk", std::fstream::out);

        for(int k = 0; k < M; k++)
        {
            //Rfk(k) = integrate( _range=elements( mesh ),
            //                    _expr=trans(ff)*idv(g[k]) ).evaluate()(0,0);
            Rfk(k) = rfkForm( g[k]);

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
        {
            s >> Rfk(k);
            if( ioption("solverns2.verbose") > 2)
                std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
        }

        s.close();
    }
}

template<typename F, typename F2, typename F3, typename E>
void
SolverSpectralProblem<F,F2,F3,E>::initRpk()
{
    if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
        std::cout << "----- Rpk -----" << std::endl;

    Rpk = VectorXd(M);

    auto a2 = expr(soption("solverns2.alpha2"));

    auto w = Sh->element();
    auto rpkForm = form1( _test=Sh );
    rpkForm = integrate( markedfaces(mesh, 1), -a2*id(w));
    rpkForm += integrate( markedfaces(mesh, 2), a2*id(w));

    if( boption("solverns2.computeRpk") )
    {
        std::fstream s;
        if ( Environment::worldComm().isMasterRank() )
            s.open ("rpk", std::fstream::out);

        for(int k = 0; k < M; k++)
        {
            // Rpk(k) = integrate( _range=markedfaces( mesh, 1 ),
            //                     _expr=-a2*idv(psi[k]) ).evaluate()(0,0);
            // Rpk(k) += integrate( _range=markedfaces( mesh, 2 ),
            //                      _expr=a2*idv(psi[k]) ).evaluate()(0,0);
            Rpk(k) = rpkForm(psi[k]);

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
        {
            s >> Rpk(k);
            if( ioption("solverns2.verbose") > 2 )
                std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
        }

        s.close();
    }
}

template<typename F, typename F2, typename F3, typename E>
typename SolverSpectralProblem<F,F2,F3,E>::vec_element_type
SolverSpectralProblem<F,F2,F3,E>::solve()
{

    if( boption("solverns2.stokes"))
    {
        tic();
        // [StokesA]
        MatrixXd A = MatrixXd(M,M);
        A = (lambda/Re).asDiagonal();
        // [StokesA]

        // [StokesB]
        VectorXd b = VectorXd(M);
        // b = Rfk;
        b = -Rpk;
        // [StokesB]

        // [StokesSolve]
        HouseholderQR<MatrixXd> qr(M,M);
        qr.compute(A);
        c = qr.solve(b);
        // [StokesSolve]

        toc( "stokes", ioption("solverns2.verbose") > 1);
    }
    else
    {
        tic();
        // [NSInit]
        VectorXd dc = VectorXd::Matrix(M);
        double tol = 1.e-8;
        // [NSInit]
        // [NSSys1]
        HouseholderQR<MatrixXd> qr(M,M);
        // [NSSys1]

        int i=0;
        do{
            // [NSMatF]
            f = c.cwiseProduct(lambda)/Re + Riak*c + Raik*c + Rpk;// - Rfk;
            for (int k = 0; k < M; k++)
                f(k) += c.transpose()*Rijk(k)*c;
            // [NSMatF]

            // [NSMatJ]
            for (int k = 0; k < M; k++)
                j.row(k) = c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
            j += Riak;
            j += Raik;
            j += lambda.asDiagonal();
            // [NSMatJ]

            // [NSSys2]
            qr.compute(j);
            dc = qr.solve(-f);
            // [NSSys2]
            // [NSAdd]
            c += dc;
            // [NSAdd]

            if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 1 )
                std::cout << "iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

            if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 3 )
                std::cout << c << std::endl;

            i++;
        } while(i < 10 && dc.norm() > tol);

        if(i==10)
            std::cout << "Newton does not converge\n";

        else{
            for( i = 0; i < M; i++){
                u += vf::project( _space=Vh, _range=elements(mesh),
                                  _expr = c(i)*idv(g[i]) );
            }
        }
        toc( "navier stokes", ioption("solverns2.verbose") > 1 );
    }

    tic();
    auto lhs = form2(_test=Vh,_trial=Vh);
    lhs = integrate(_range=elements(mesh), _expr=inner(idt(u),id(u)));
    auto rhs = form1(_test=Vh);
    for(int i=0; i<M; i++)
    {
        if ( Environment::worldComm().isMasterRank() && ioption("solverns2.verbose") > 2 )
            std::cout << "c(" << i << ") = " << c(i) << std::endl;
        rhs += integrate( _range=elements(mesh), _expr = c(i)*trans(idv(g[i]))*id(u) );
    }
    lhs.solve(_rhs=rhs, _solution=u);
    toc( "compute u", ioption("solverns2.verbose") > 1 );

    return u;
}
