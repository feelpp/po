using namespace Feel;

template<typename FunctionSpaceType>
class SolverA
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpaceType ml_space_ptrtype;
    typedef typename ml_space_ptrtype::element_type ml_space_type;

    typedef typename ml_space_type::template sub_functionspace<0>::type rt_space_type;
    typedef boost::shared_ptr<rt_space_type> rt_space_ptrtype;
    typedef typename rt_space_type::element_type rt_element_type;

    typedef SolverA<ml_space_ptrtype> solvera_type;
    typedef typename boost::shared_ptr<solvera_type> solvera_ptrtype;

    using bilinear_form_type = typename Feel::meta::BilinearForm<ml_space_type, ml_space_type>::type;
    using linear_form_type = typename Feel::meta::LinearForm<ml_space_type>::type;

    mesh_ptrtype mesh;

    ml_space_ptrtype Mh;
    rt_space_ptrtype RTh;

    rt_element_type a;

    rt_element_type a0;
    bilinear_form_type a0Form2;
    bilinear_form_type a0Form2BC;
    linear_form_type a0Form1;

    rt_element_type a1;
    rt_element_type a2;

    void initA0();
    void computeA0( double t );
    void loadA0();
    void computeA1( double t );
    void loadA1();
    void computeA2( double t );
    void loadA2();

public:
    static solvera_ptrtype build(const mesh_ptrtype& mesh, const ml_space_ptrtype& Mh);
    rt_element_type solve( double t);
};

template<typename T>
typename SolverA<T>::solvera_ptrtype
SolverA<T>::build(const mesh_ptrtype& mesh, const ml_space_ptrtype& Mh)
{
    solvera_ptrtype solvera( new SolverA<T> );

    solvera->mesh = mesh;

    solvera->Mh = Mh;
    solvera->RTh = Mh->template functionSpace<0>();

    solvera->initA0();

    return solvera;
}

template<typename T>
typename SolverA<T>::rt_element_type
SolverA<T>::solve( double t )
{
    a = RTh->element();

    tic();
    if( boption("solverns2.compute-a0"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
            std::cout << " ---------- compute a0 ----------\n";
        computeA0( t );
    }
    else
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
            std::cout << " ---------- load a0 ----------\n";
        loadA0();
    }
    a += a0;
    toc( "a0", ioption("solverns2.verbose") > 2);

    if( boption("solverns2.need-a1") )
    {
        tic();
        if( boption("solverns2.compute-a1"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- compute a1 ----------\n";
            computeA1( t );
        }
        else
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- load a1 ----------\n";
            loadA1();
        }
        a += a1;
        toc( "a1", ioption("solverns2.verbose") > 2);
    }

    if( boption("solverns2.need-a2") )
    {
        tic();
        if( boption("solverns2.compute-a2"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- compute a2 ----------\n";
            computeA2( t );
        }
        else
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- load a2 ----------\n";
            loadA2();
        }
        a += a2;
        toc( "a2", ioption("solverns2.verbose") > 2);
    }

    return a;
}

template<typename T>
void
SolverA<T>::initA0()
{
    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.template element<0>() ;
    auto p = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto q = V.template element<1>() ;

    auto a0Form2 = form2( _test=Mh, _trial=Mh );
    auto a0Form2BC = form2( _test=Mh, _trial=Mh );
    auto a0Form1 = form1( _test=Mh );
    // a0Form2 = integrate( elements(mesh),
    //                      -trans(idt(u))*id(v)
    //                      -idt(p)*div(v)
    //                      -divt(u)*id(q)
    //                      -cst(0.5)*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q)))
    //                      -cst(0.5)*(divt(u)*div(v))
    //                      -cst(0.5)*(trans(curlt(u))*curl(v))
    //                      );
}

template<typename T>
void
SolverA<T>::computeA0( double t )
{
    // [option]
    auto g = expr(soption("solverns2.alpha0"));
    g.setParameterValues( {
            { "radius", doption( "solverns2.radius" ) },
            { "speed", doption( "solverns2.speed" ) },
            { "t", t} } );
    auto alpha0  = vec( cst(0.), cst(0.), g);
    // [option]

    // auto U = Mh->element();
    // auto u = U.template element<0>() ;
    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.template element<0>() ;
    auto p = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto q = V.template element<1>() ;

    // [bilinearA]
    auto a0Form2 = form2( _test=Mh, _trial=Mh );
    auto a0Form2BC = form2( _test=Mh, _trial=Mh );
    auto a0Form1 = form1( _test=Mh );
    a0Form2 = integrate( elements(mesh),
                         -trans(idt(u))*id(v)
                         -idt(p)*div(v)
                         -divt(u)*id(q)
                         -cst(0.5)*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q)))
                         -cst(0.5)*(divt(u)*div(v))
                         -cst(0.5)*(trans(curlt(u))*curl(v))
                         );
    a0Form1.zero();
    a0Form2BC = a0Form2;
    a0Form2BC += on( _range=boundaryfaces(mesh), _rhs=a0Form1, _element=u, _expr=alpha0);

    a0Form2BC.solveb( _backend=backend(_name="a0"), _rhs=a0Form1, _solution=U );

    a0 = RTh->element();
    a0 = u;

    std::string path = "a0";
    a0.save(_path=path);

    // auto diva = normL2(elements(mesh), divv(a0));
    // auto an = normL2(markedfaces(mesh,1), inner(idv(a0), N()) + g );

    // std::ofstream s;
    // if( Environment::isMasterRank() )
    // {
    //     s.open( "convergence_a.dat", std::ios::out | std::ios::app);
    //     s << "hsize" << "\t" << "nDof" << "\t" << "div" << "\t" << "an" << std::endl;
    //     s << doption("gmsh.hsize") << "\t" << RTh->nDof() << "\t" << diva << "\t" << an << std::endl;
    //     if( ioption("solverns2.verbose") > 2)
    //         std::cout << "diva = " << diva << std::endl
    //                   << "an = " << an << std::endl;
    // }
}

template<typename T>
void
SolverA<T>::loadA0()
{
    a0 = RTh->element();
    std::string path = "a0";
    a0.load(_path=path);
}

template<typename T>
void
SolverA<T>::computeA1( double t )
{
}

template<typename T>
void
SolverA<T>::loadA1()
{
    a1 = RTh->element();
    std::string path = "a1";
    a1.load(_path=path);
}

template<typename T>
void
SolverA<T>::computeA2( double t )
{
}

template<typename T>
void
SolverA<T>::loadA2()
{
    a2 = RTh->element();
    std::string path = "a2";
    a2.load(_path=path);
}
