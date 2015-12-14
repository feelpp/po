#include <Eigen/Dense>
#include <feel/feeldiscr/expansion.hpp>
#include <boost/mpi/timer.hpp>

using namespace Feel;
using namespace Eigen;

template<typename FunctionSpace, typename FunctionSpace2, typename AType>
class SolverSpectralProblem
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    using eigen_space_ptrtype = FunctionSpace;
    using eigen_type = typename eigen_space_ptrtype::element_type::element_type;
    using psi_space_ptrtype = FunctionSpace2;
    using psi_type = typename psi_space_ptrtype::element_type::element_type;
    using a_type = AType;

    using solverspectralproblem_type = SolverSpectralProblem<eigen_space_ptrtype,psi_space_ptrtype, a_type>;
    using solverspectralproblem_ptrtype = typename boost::shared_ptr<solverspectralproblem_type>;

    // using eigenmodes_type = std::vector<eigen_vec_type>;
    using eigenvec_type = std::vector<eigen_type>;
    using psivec_type = std::vector<psi_type>;

    mesh_ptrtype mesh;
    // vec_space_ptrtype Vh;
    eigen_space_ptrtype Nh;
    psi_space_ptrtype Sh;

    a_type a;
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
    int MMax;
    VectorXd c;
    VectorXd cNm1;
    eigen_type u;

    static solverspectralproblem_ptrtype build(const mesh_ptrtype& mesh, const eigen_space_ptrtype& Nh, const psi_space_ptrtype& Sh );
    void setA(const a_type& a);
    void setEigen();
    void setRijk();
    void init( double t );
    eigen_type solve(double t);
};

template<typename F1, typename F2, typename A1>
typename SolverSpectralProblem<F1,F2,A1>::solverspectralproblem_ptrtype
SolverSpectralProblem<F1,F2,A1>::build(const mesh_ptrtype& mesh, const eigen_space_ptrtype& Nh, const psi_space_ptrtype& Sh )
{
    solverspectralproblem_ptrtype ap( new SolverSpectralProblem<F1,F2,A1> );
    ap->mesh = mesh;
    ap->Nh = Nh;
    ap->Sh = Sh;

    ap->M = ioption("solverns2.nb-mode");

    double r = doption( "solverns2.radius" );
    double s = doption( "solverns2.speed" );
    double n = doption( "solverns2.nu" );
    ap->Re = 2*r*s/n;

    return ap;
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::setA(const a_type& a)
{
    this->a = a;
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::setEigen()
{
    tic();
    lambda = VectorXd(M);

    g = eigenvec_type(M, Nh->element());
    psi = psivec_type(M, Sh->element());

    std::fstream s;
    Feel::fs::path basePathG(soption( _name="solverns2.path" ));
    basePathG.append("mode");
    Feel::fs::path basePathP(soption( _name="solverns2.path" ));
    basePathP.append("psi");
    Feel::fs::path basePathL(soption( _name="solverns2.path" ));
    basePathL.append("lambda");

    s.open (basePathL.string(), std::fstream::in);
    if( !s.is_open() ){
        std::cout << "Eigen values not found" << std::endl;
        exit(0);
    }

    for( int i = 0; i < M && s.good(); i++ ){
        s >> lambda(i);
        std::string pathG = (boost::format("%1%-%2%") % basePathG.string() % i).str();
        g[i].load(_path=pathG, _type=soption("solverns2.format"));
        std::string pathP = (boost::format("%1%-%2%") % basePathP.string() % i).str();
        psi[i].load(_path=pathP, _type=soption("solverns2.format"));
    }
    s.close();

    c = VectorXd::Ones(M);
    toc( "eigen", ioption("solverns2.verbose") > 2);
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::init( double t )
{
    if( ! boption("solverns2.stokes"))
    {
        initRaik();
        initRiak();
    }
    initRpk( t );
    initRfk( t );
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::setRijk()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rijk -----" << std::endl;

    Feel::fs::path path(soption( _name="solverns2.path" ));
    auto filename = path.append("rijk").string();

    std::ifstream in(filename,std::ios::in | std::ios::binary);

    // load all the modes in Rijk
    in.read((char*) (&MMax),sizeof(typename decltype(Rijk)::Index));
    if( M > MMax)
    {
        LOG(WARNING) << "Number of modes " << M << " is greater than size of the saved coefficients " << MMax;
        if( Environment::isMasterRank() )
            std::cout << "WARNING : Number of modes " << M << " is greater than size of the saved coefficients " << MMax << std::endl;
    }
    Rijk = Matrix<MatrixXd, Dynamic, 1>(MMax,1);
    for(int k = 0; k < MMax; k++)
        Rijk(k) = MatrixXd(MMax,MMax);

    for( int k = 0; k < MMax; k++ )
        in.read( (char *) Rijk(k).data() , MMax*MMax*sizeof(typename decltype(Rijk)::Scalar::Scalar) );
    in.close();

    for( int k = 0; k< MMax; k++ )
        Rijk(k).conservativeResize(M,M);
    Rijk.conservativeResize(M,Eigen::NoChange);

    toc( "Rijk", ioption("solverns2.verbose") > 1);
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::initRaik()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Raik -----" << std::endl;

    auto w = Nh->element();
    auto raikForm = form2( _trial=Nh, _test=Nh );
    raikForm = integrate( elements(mesh), inner( cross( curlv(a),id(w) ), idt(w)) );

    if( boption("solverns2.compute-raik") )
    {
        Raik = MatrixXd(M,M);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < i; k++)
            {
                Raik(k,i) = raikForm(g[i], g[k]);
                Raik(i,k) = -Raik(k,i);

                if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                    std::cout << "Raik(" << k << "," << i << ") = " << Raik(k,i) << std::endl;
            }
            Raik(i,i) = 0;
        }

        if( Environment::isMasterRank() )
        {
            std::ofstream out("raik",std::ios::out | std::ios::binary | std::ios::trunc);
            out.write((char*) Raik.data(), M*M*sizeof(typename decltype(Raik)::Scalar) );
            out.close();
        }
    }
    else
    {
        Raik = MatrixXd(MMax,MMax);

        std::ifstream in("raik",std::ios::in | std::ios::binary);
        if( !in.is_open() )
        {
            std::cout << "Raik not found\ntry to launch with --solverns2.compute-raik=true" << std::endl;
            exit(0);
        }
        in.read( (char *) Raik.data() , MMax*MMax*sizeof(double) );
        in.close();

        Raik.conservativeResize(M,M);
    }
    toc( "Raik", ioption("solverns2.verbose") > 1);
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::initRiak()
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Riak -----" << std::endl;

    auto w = Nh->element();
    auto riakForm = form2( _test=Nh, _trial=Nh );
    riakForm = integrate( elements(mesh), inner( cross( curl(w),idv(a) ), idt(w)) );

    if( boption("solverns2.compute-riak") )
    {
        Riak = MatrixXd(M,M);

        for(int i = 0; i< M ; i++)
        {
            for(int k = 0; k < M; k++)
            {
                Riak(k,i) = riakForm( g[i], g[k] );

                if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                    std::cout << "Riak(" << k << "," << i << ") = " << Riak(k,i) << std::endl;
            }
        }
        if( Environment::isMasterRank() )
        {
            std::ofstream out("riak",std::ios::out | std::ios::binary | std::ios::trunc);
            out.write((char*) Riak.data(), M*M*sizeof(typename decltype(Riak)::Scalar) );
            out.close();
        }
    }
    else
    {
        Riak = MatrixXd(MMax,MMax);

        std::ifstream in("riak",std::ios::in | std::ios::binary);
        if( !in.is_open() )
        {
            std::cout << "Riak not found\ntry to launch with --solverns2.compute-riak=true" << std::endl;
            exit(0);
        }
        in.read( (char *) Riak.data() , MMax*MMax*sizeof(double) );
        in.close();

        Riak.conservativeResize(M,M);
    }
    toc( "Riak", ioption("solverns2.verbose") > 1);
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::initRfk( double t )
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rfk -----" << std::endl;

    auto ff = expr<3,1>(soption("solverns2.f"));
    ff.setParameterValues( {{"t",t}} );

    auto w = Nh->element();
    auto rfkForm = form1( _test=Nh );
    rfkForm = integrate( elements(mesh), trans(ff)*id(w));

    if( boption("solverns2.compute-rfk") )
    {
        Rfk = VectorXd(M);

        for(int k = 0; k < M; k++)
        {
            Rfk(k) = rfkForm( g[k]);

            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                std::cout << "Rfk(" << k << ") = " << Rfk(k) << std::endl;
        }
        if( Environment::isMasterRank() )
        {
            std::ofstream out("rfk",std::ios::out | std::ios::binary | std::ios::trunc);
            out.write((char*) Rfk.data(), M*sizeof(typename decltype(Rfk)::Scalar) );
            out.close();
        }
    }
    else
    {
        Rfk = VectorXd(MMax);

        std::ifstream in("rfk",std::ios::in | std::ios::binary);
        if( !in.is_open() )
        {
            std::cout << "Rfk not found\ntry to launch with --solverns2.compute-rfk=true" << std::endl;
            exit(0);
        }
        in.read( (char *) Rfk.data() , MMax*sizeof(double) );
        in.close();

        Rfk.conservativeResize(M);
    }
    toc( "Rfk", ioption("solverns2.verbose") > 1);
}

template<typename F1, typename F2, typename A1>
void
SolverSpectralProblem<F1,F2,A1>::initRpk( double t )
{
    tic();
    if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2 )
        std::cout << "----- Rpk -----" << std::endl;

    auto a2 = expr(soption("solverns2.alpha2"));
    a2.setParameterValues( {
            {"t",t},
            {"speed",doption("solverns2.speed")},
            {"radius",doption("solverns2.radius")}
        } );

    auto w = Sh->element();
    auto rpkForm = form1( _test=Sh );
    rpkForm = integrate( markedfaces(mesh, 1), -a2*id(w));
    rpkForm += integrate( markedfaces(mesh, 2), a2*id(w));

    if( boption("solverns2.compute-rpk") )
    {
        Rpk = VectorXd(M);

        for(int k = 0; k < M; k++)
        {
            Rpk(k) = rpkForm(psi[k]);

            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 3 )
                std::cout << "Rpk(" << k << ") = " << Rpk(k) << std::endl;
        }
        if( Environment::isMasterRank() )
        {
            std::ofstream out("rpk",std::ios::out | std::ios::binary | std::ios::trunc);
            out.write((char*) Rpk.data(), M*sizeof(typename decltype(Rpk)::Scalar) );
            out.close();
        }
    }
    else
    {
        Rpk = VectorXd(MMax);

        std::ifstream in("rpk",std::ios::in | std::ios::binary);
        if( !in.is_open() )
        {
            std::cout << "Rpk not found\ntry to launch with --solverns2.compute-rpk=true" << std::endl;
            exit(0);
        }
        in.read( (char *) Rpk.data() , MMax*sizeof(double) );
        in.close();

        Rpk.conservativeResize(M);
    }
    toc( "Rpk", ioption("solverns2.verbose") > 1);
}

template<typename F1, typename F2, typename A1>
typename SolverSpectralProblem<F1,F2,A1>::eigen_type
SolverSpectralProblem<F1,F2,A1>::solve(double t)
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
        double dt = doption("solverns2.time-step");

        VectorXd dc = VectorXd::Matrix(M);
        double tol = doption("solverns2.newton-tol");
        // [NSInit]
        // [NSSys1]
        HouseholderQR<MatrixXd> qr(M,M);
        // [NSSys1]

        int i=0;
        do{
            // [NSMatF]
            f = c.cwiseProduct(lambda/Re) + Riak*c + Raik*c + Rpk/Re - Rfk;

            if( t > doption("solverns2.start-time") )
                f += c/dt - cNm1/dt;

            for (int k = 0; k < M; k++)
                f(k) += c.transpose()*Rijk(k)*c;
            // [NSMatF]

            // [NSMatJ]
            for (int k = 0; k < M; k++)
                j.row(k) = c.transpose()*Rijk(k).transpose() + c.transpose()*Rijk(k);
            j += Riak;
            j += Raik;
            j += (lambda/Re).asDiagonal();
            if( t > doption("solverns2.start-time") )
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
        } while(i < ioption("solverns2.newton-max-it") && dc.norm() > tol);

        if(i==ioption("solverns2.newton-max-it"))
        {
            LOG(WARNING) << "Newton did not converge";
            if(Environment::isMasterRank())
                std::cout << "WARNING : Newton did not converge" << std::endl;
        }
        toc( "navier stokes", ioption("solverns2.verbose") > 1 );
    }

    tic();
    u = Nh->element();
    u = expansion(g, c, c.size());
    toc ( "compute u", ioption("solverns2.verbose") > 1 );

    return u;
}
