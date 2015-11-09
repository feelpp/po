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

    void computeA0( );
    void loadA0();
    void computeA1( );
    void loadA1();
    void computeA2( );
    void loadA2();

public:
    static solvera_ptrtype build(const mesh_ptrtype& mesh, const ml_space_ptrtype& Mh);
    rt_element_type solve();
};

template<typename T>
typename SolverA<T>::solvera_ptrtype
SolverA<T>::build(const mesh_ptrtype& mesh, const ml_space_ptrtype& Mh)
{
    solvera_ptrtype solvera( new SolverA<T> );

    solvera->mesh = mesh;

    solvera->Mh = Mh;
    solvera->RTh = Mh->template functionSpace<0>();

    return solvera;
}

template<typename T>
typename SolverA<T>::rt_element_type
SolverA<T>::solve( )
{
    a = Mh->template functionSpace<0>()->element();

    tic();
    if( boption("solverns2.computeA0"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
            std::cout << " ---------- compute a0 ----------\n";
        computeA0( );
    }
    else
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
            std::cout << " ---------- load a0 ----------\n";
        loadA0();
    }
    a += a0;
    toc( "a0", ioption("solverns2.verbose") > 2);

    if( boption("solverns2.needA1") )
    {
        tic();
        if( boption("solverns2.computeA1"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- compute a1 ----------\n";
            computeA1( );
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

    if( boption("solverns2.needA2") )
    {
        tic();
        if( boption("solverns2.computeA2"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 2)
                std::cout << " ---------- compute a2 ----------\n";
            computeA2( );
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
SolverA<T>::computeA0( )
{
    // [option]
    auto g = expr(soption("solverns2.alpha0"));
    g.setParameterValues( {
            { "radius", doption( "solverns2.radius" ) },
            { "speed", doption( "solverns2.speed" ) } } );
    auto alpha0  = vec( cst(0.), cst(0.), g);
    // [option]

    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.template element<0>() ;
    auto p = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto q = V.template element<1>() ;

    auto a0Form2 = form2( _test=Mh, _trial=Mh );
    auto a0Form1 = form1( _test=Mh );
    a0Form2 = integrate( elements(mesh),
                         -trans(idt(u))*id(v)
                         -idt(p)*div(v)
                         -divt(u)*id(q)
                         -cst(0.5)*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q)))
                         -cst(0.5)*(divt(u)*div(v))
                         -cst(0.5)*(trans(curlt(u))*curl(v))
                         );

    a0Form2 += on( _range=boundaryfaces(mesh), _rhs=a0Form1, _element=u, _expr=alpha0);

    a0Form2.solveb( _backend=backend("a0"), _rhs=a0Form1, _solution=U );

    a0 = u.functionSpace()->element();
    a0 = u;

    std::string path = "a0";
    a0.save(_path=path);
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
SolverA<T>::computeA1( )
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
SolverA<T>::computeA2( )
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
