using namespace Feel;

template<typename FunctionSpaceType, typename FunctionSpaceType2>
class SolverA
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpaceType space_ptrtype;
    typedef typename space_ptrtype::element_type space_type;
    typedef typename space_type::element_type element_type;

    typedef FunctionSpaceType2 ml_space_ptrtype;
    typedef typename ml_space_ptrtype::element_type ml_space_type;

    typedef typename ml_space_type::template sub_functionspace<0>::type rt_space_type;
    typedef boost::shared_ptr<rt_space_type> rt_space_ptrtype;
    typedef typename rt_space_type::element_type rt_element_type;

    typedef SolverA<space_ptrtype, ml_space_ptrtype> solvera_type;
    typedef typename boost::shared_ptr<solvera_type> solvera_ptrtype;


    mesh_ptrtype mesh;
    space_ptrtype Vh;
    ml_space_ptrtype Mh;
    rt_space_ptrtype RTh;
    // scalar_space_ptrtype Sh;
    rt_element_type a;
    rt_element_type a0;
    rt_element_type a1;
    rt_element_type a2;

    void computeA0();
    void loadA0();
    void computeA1();
    void loadA1();
    void computeA2();
    void loadA2();

public:
    static solvera_ptrtype build(const mesh_ptrtype& mesh, const space_ptrtype& Vh, const ml_space_ptrtype& Mh);
    rt_element_type solve();

};

template<typename T, typename T2>
typename SolverA<T,T2>::solvera_ptrtype
SolverA<T,T2>::build(const mesh_ptrtype& mesh, const space_ptrtype& Vh, const ml_space_ptrtype& Mh)
{
    solvera_ptrtype solvera( new SolverA<T,T2> );
    solvera->mesh = mesh;
    solvera->Vh = Vh;
    solvera->Mh = Mh;
    solvera->RTh = Mh->template functionSpace<0>();
    return solvera;
}

template<typename T, typename T2>
typename SolverA<T,T2>::rt_element_type
SolverA<T,T2>::solve()
{
    a = RTh->element();

    tic();
    if( boption("solverns2.computeA0"))
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- compute a0 ----------\n";
        computeA0();
    }
    else
    {
        if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
            std::cout << " ---------- load a0 ----------\n";
        loadA0();
    }
    a += a0;
    toc( "a0", ioption("solverns2.verbose") > 1);

    if( boption("solverns2.needA1") )
    {
        tic();
        if( boption("solverns2.computeA1"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
                std::cout << " ---------- compute a1 ----------\n";
            computeA1();
        }
        else
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
                std::cout << " ---------- load a1 ----------\n";
            loadA1();
        }
        a += a1;
        toc( "a1", ioption("solverns2.verbose") > 1);
    }

    if( boption("solverns2.needA2") )
    {
        tic();
        if( boption("solverns2.computeA2"))
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
                std::cout << " ---------- compute a2 ----------\n";
            computeA2();
        }
        else
        {
            if ( Environment::isMasterRank() && ioption("solverns2.verbose") > 0)
                std::cout << " ---------- load a2 ----------\n";
            loadA2();
        }
        a += a2;
        toc( "a2", ioption("solverns2.verbose") > 1);
    }

    return a;
}

template<typename T, typename T2>
void
SolverA<T,T2>::computeA0()
{
    // [option]
    auto g_s = soption("solverns2.alpha0");
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", doption( "solverns2.radius" ) },
                { "speed", doption( "solverns2.speed" ) } } );
    auto alpha0  = vec( cst(0.), cst(0.), g);
    //auto f = div(alpha0);
    // [option]

    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.template element<0>() ;
    auto p = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto q = V.template element<1>() ;

    // [bilinearA]
    auto l = form1( _test=Mh );
    auto a = form2( _trial=Mh, _test=Mh );
    a = integrate( elements(mesh),
                   -trans(idt(u))*id(v)
                   -idt(p)*div(v)
                   -divt(u)*id(q)
                   -cst(0.5)*(trans(idt(u)) - gradt(p))*(-id(v)+trans(grad(q)))
                   -cst(0.5)*(divt(u)*div(v))
                   -cst(0.5)*(trans(curlt(u))*curl(v))
                   );
    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=alpha0);

    a.solve( _name="a0", _rhs=l, _solution=U );

    a0 = RTh->element();
    a0 = u;

    std::string path = "a0";
    a0.save(_path=path);

    auto diva = normL2(elements(mesh), divv(a0));
    auto an = normL2(markedfaces(mesh,1), inner(idv(a0), N()) + g );

    std::ofstream s;
    if( Environment::isMasterRank() )
    {
        s.open( "convergence_a.dat", std::ios::out | std::ios::app);
        s << "hsize" << "\t" << "nDof" << "\t" << "div" << "\t" << "an" << std::endl;
        s << doption("gmsh.hsize") << "\t" << RTh->nDof() << "\t" << diva << "\t" << an << std::endl;
        std::cout << "diva = " << diva << std::endl
                  << "an = " << an << std::endl;
    }
}

template<typename T, typename T2>
void
SolverA<T,T2>::loadA0()
{
    a0 = RTh->element();
    std::string path = "a0";
    a0.load(_path=path);
}

template<typename T, typename T2>
void
SolverA<T,T2>::computeA1()
{
}

template<typename T, typename T2>
void
SolverA<T,T2>::loadA1()
{
    a1 = RTh->element();
    std::string path = "a1";
    a1.load(_path=path);
}

template<typename T, typename T2>
void
SolverA<T,T2>::computeA2()
{
}

template<typename T, typename T2>
void
SolverA<T,T2>::loadA2()
{
    a2 = RTh->element();
    std::string path = "a2";
    a2.load(_path=path);
}
