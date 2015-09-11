#include <Eigen/Dense>
#include <feel/feeldiscr/expansion.hpp>
#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpace, typename FunctionSpace2, typename FunctionSpace3, typename FunctionSpace4, typename EigenTuple>
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

    typedef FunctionSpace4 rt_space_ptrtype;
    typedef typename rt_space_ptrtype::element_type rt_space_type;
    typedef typename rt_space_type::element_type rt_element_type;

    typedef EigenTuple eigentuple_type;

    typedef SolverSpectralProblem<vec_space_ptrtype, ned_space_ptrtype, scalar_space_ptrtype, rt_space_ptrtype, eigentuple_type> solverspectralproblem_type;
    typedef typename boost::shared_ptr<solverspectralproblem_type> solverspectralproblem_ptrtype;

    typedef std::vector<eigentuple_type> eigenmodes_type;
    typedef std::vector<ned_element_type> eigenvec_type;
    typedef std::vector<scalar_element_type> psivec_type;

    mesh_ptrtype mesh;
    vec_space_ptrtype Vh;
    ned_space_ptrtype Nh;
    scalar_space_ptrtype Sh;
    rt_space_ptrtype RTh;

    rt_element_type a;
    eigenvec_type g;
    psivec_type psi;
    VectorXd lambda;

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

    void initRaik();
    void initRiak();
    void initRfk( double t );
    void initRpk( double t );


public:
    int M;
    VectorXd c;
    VectorXd cNm1;
    ned_element_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const ned_space_ptrtype& Nh, const scalar_space_ptrtype& Sh, const rt_space_ptrtype& RThh);
    void setA(const rt_element_type& a);
    void setEigen();
    void setRijk();
    void init( double t );
    ned_element_type solve(double t);
};

template<typename F, typename F2, typename F3, typename F4, typename E>
typename SolverSpectralProblem<F,F2,F3,F4,E>::solverspectralproblem_ptrtype
SolverSpectralProblem<F,F2,F3,F4,E>::build(const mesh_ptrtype& mesh, const vec_space_ptrtype& Vh, const ned_space_ptrtype& Nh, const scalar_space_ptrtype& Sh, const rt_space_ptrtype& RTh)
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F,F2,F3,F4,E> );
    ap->mesh = mesh;
    ap->Vh = Vh;
    ap->Nh = Nh;
    ap->Sh = Sh;
    ap->RTh = RTh;

    double r = doption( "solverns2.radius" );
    double s = doption( "solverns2.speed" );
    double n = doption( "solverns2.nu" );
    ap->Re = 2*r*s/n;

    return ap;
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::setA(const rt_element_type& a)
{
    this->a = a;
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::setEigen()
{
    M = ioption("solverns2.nbMode");
    lambda = VectorXd(M);

    g = eigenvec_type(M);
    psi = psivec_type(M);

    std::fstream s;
    Feel::fs::path path(soption( _name="solverns2.path" ));
    auto basePathG = path.append("mode").string();
    auto basePathP = path.append("psi").string();
    auto basePathL = path.append("lambda").string();
    s.open (basePathL, std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found" << std::endl;
        exit(0);
    }

    for( int i = 0; i < M && s.good(); i++ ){
        s >> lambda(i);
        std::string pathG = (boost::format("%1%-%2%") % basePathG % i).str();
        g[i].load(_path=pathG);
        std::string pathP = (boost::format("%1%-%2%") % basePathP % i).str();
        psi[i].load(_path=pathP);
    }
    s.close();

    c = VectorXd::Ones(M);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::init( double t )
{
    if( ! boption("solverns2.stokes"))
    {
        initRaik();
        initRiak();
    }
    initRpk( t );
    initRfk( t );
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::setRijk()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rijk -----" << std::endl;

    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

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

    toc( "Rijk", ioption("solverns2.verbose") > 1);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::initRaik()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Raik -----" << std::endl;

    Raik = MatrixXd(M,M);

    auto w = Nh->element();
    auto raikForm = form2( _trial=Nh, _test=Nh );
    raikForm = integrate( elements(mesh), inner( cross( curlv(a),id(w) ), idt(w)) );

    if( boption("solverns2.computeRaik") )
    {
        std::fstream s;
        if ( Environment::isMasterRank() )
            s.open ("raik", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                Raik(k,i) = raikForm(g[i], g[k]);
                Raik(i,k) = -Raik(k,i);

                if ( Environment::isMasterRank() )
                {
                    if( ioption("solverns2.verbose") > 3 )
                        std::cout << "Raik(" << k << "," << i << ") = " << Raik(k,i) << std::endl;
                    s << Raik(k,i) << std::endl;
                }
            }
            Raik(i,i) = 0;
        }
        if ( Environment::isMasterRank() )
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
                if( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                    std::cout << "Raik(" << k << "," << i << ") = " << Raik(k,i) << std::endl;
            }
            Raik(i,i) = 0;
        }
        s.close();
    }
    toc( "Raik", ioption("solverns2.verbose") > 1);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::initRiak()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Riak -----" << std::endl;

    Riak = MatrixXd(M,M);

    auto w = Nh->element();
    auto riakForm = form2( _test=Nh, _trial=Nh );
    riakForm = integrate( elements(mesh), inner( cross( curl(w),idv(a) ), idt(w)) );

    if( boption("solverns2.computeRiak") )
    {
        std::fstream s;
        if ( Environment::isMasterRank() )
            s.open ("riak", std::fstream::out);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < M; k++)
            {
                Riak(k,i) = riakForm( g[i], g[k] );

                if ( Environment::isMasterRank() )
                {
                    if( ioption("solverns2.verbose") > 3 )
                        std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
                    s << Riak(k,i) << std::endl;
                }
            }
        }
        if ( Environment::isMasterRank() )
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
                if( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                    std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
            }
        }
        s.close();
    }
    toc( "Riak", ioption("solverns2.verbose") > 1);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::initRfk( double t )
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rfk -----" << std::endl;

    Rfk = VectorXd(M);

    auto ff = expr<3,1>(soption("solverns2.f"));
    ff.setParameterValues( {{"t",t}} );

    auto w = Nh->element();
    auto rfkForm = form1( _test=Nh );
    rfkForm = integrate( elements(mesh), trans(ff)*id(w));

    if( boption("solverns2.computeRfk") )
    {
        std::fstream s;
        if ( Environment::isMasterRank() )
            s.open ("rfk", std::fstream::out);

        for(int k = 0; k < M; k++)
        {
            Rfk(k) = rfkForm( g[k]);

            if ( Environment::isMasterRank() )
            {
                if( ioption("solverns2.verbose") > 3)
                    std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
                s << Rfk(k) << std::endl;
            }
        }
        if ( Environment::isMasterRank() )
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
            if( Environment::isMasterRank() && ioption("solverns2.verbose") > 3)
                std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
        }

        s.close();
    }
    toc( "Rfk", ioption("solverns2.verbose") > 1);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
void
SolverSpectralProblem<F,F2,F3,F4,E>::initRpk( double t )
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rpk -----" << std::endl;

    Rpk = VectorXd(M);

    auto a2 = expr(soption("solverns2.alpha2"));
    a2.setParameterValues( {{"t",t}} );

    auto w = Sh->element();
    auto rpkForm = form1( _test=Sh );
    rpkForm = integrate( markedfaces(mesh, 1), -a2*id(w));
    rpkForm += integrate( markedfaces(mesh, 2), a2*id(w));

    if( boption("solverns2.computeRpk") )
    {
        std::fstream s;
        if ( Environment::isMasterRank() )
            s.open ("rpk", std::fstream::out);

        for(int k = 0; k < M; k++)
        {
            Rpk(k) = rpkForm(psi[k]);

            if ( Environment::isMasterRank())
            {
                if( ioption("solverns2.verbose") > 3 )
                    std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
                s << Rpk(k) << std::endl;
            }
        }
        if ( Environment::isMasterRank() )
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
            if( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
        }

        s.close();
    }
    toc( "Rpk", ioption("solverns2.verbose") > 1);
}

template<typename F, typename F2, typename F3, typename F4, typename E>
typename SolverSpectralProblem<F,F2,F3,F4,E>::ned_element_type
SolverSpectralProblem<F,F2,F3,F4,E>::solve(double t)
{
    j = MatrixXd(M,M);
    f = VectorXd(M);

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
        cNm1 = c;
        c = VectorXd::Ones(M);

        // [NSInit]
        double dt = doption("solverns2.timeStep");

        VectorXd dc = VectorXd::Matrix(M);
        double tol = doption("solverns2.newton.tol");
        // [NSInit]
        // [NSSys1]
        HouseholderQR<MatrixXd> qr(M,M);
        // [NSSys1]

        int i=0;
        do{
            // [NSMatF]
            f = c.cwiseProduct(lambda)/Re + Riak*c + Raik*c + Rpk - Rfk;

            if( t > doption("solverns2.startTime") )
                f += c/dt - cNm1/dt;

            for (int k = 0; k < M; k++)
                f(k) += c.transpose()*Rijk(k)*c;
            // [NSMatF]

            // [NSMatJ]
            for (int k = 0; k < M; k++)
                j.row(k) = c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
            j += Riak;
            j += Raik;
            j += lambda.asDiagonal();
            if( t > doption("solverns2.startTime") )
                j += VectorXd::Constant(M, 1./dt).asDiagonal();
            // [NSMatJ]

            // [NSSys2]
            qr.compute(j);
            dc = qr.solve(-f);
            // [NSSys2]
            // [NSAdd]
            c += dc;
            // [NSAdd]

            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
                std::cout << "Newton iteration : " << i << " norm(dc) = " << dc.norm() << std::endl;

            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                std::cout << c << std::endl;

            i++;
        } while(i < ioption("solverns2.newton.maxIt") && dc.norm() > tol);

        if(i==ioption("solverns2.newton.maxIt"))
        {
            std::cout << "Newton does not converge\n";
        }
        toc( "navier stokes", ioption("solverns2.verbose") > 1 );
    }

    tic();
    u = Nh->element();
    u = expansion(g, c, c.size());
    toc ( "compute u", ioption("solverns2.verbose") > 1 );

    return u;
}
