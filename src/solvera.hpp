using namespace Feel;

template<typename FunctionSpaceType>
class SolverA
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpaceType space_ptrtype;
    typedef typename space_ptrtype::element_type space_type;
    typedef typename space_type::element_type element_type;

    typedef SolverA<space_ptrtype> solvera_type;
    typedef typename boost::shared_ptr<solvera_type> solvera_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar>, Lagrange<0, Scalar> > > ml_space_type;

    mesh_ptrtype mesh;
    space_ptrtype Vh;
    element_type a;
    element_type a0;
    element_type a1;
    element_type a2;

    void computeA0();
    void loadA0();
    void computeA1();
    void loadA1();
    void computeA2();
    void loadA2();

public:
    static solvera_ptrtype build(const mesh_ptrtype& mesh, const space_ptrtype& Vh);
    element_type solve();
};

template<typename T>
typename SolverA<T>::solvera_ptrtype
SolverA<T>::build(const mesh_ptrtype& mesh, const space_ptrtype& Vh)
{
    solvera_ptrtype solvera( new SolverA<T> );
    solvera->mesh = mesh;
    solvera->Vh = Vh;
    return solvera;
}

template<typename T>
typename SolverA<T>::element_type
SolverA<T>::solve()
{
    a = Vh->element();
    boost::mpi::timer t;

    if( boption("solverns2.needA0") )
    {
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
        logTime(t, "a0", ioption("solverns2.verbose") > 1);
    }

    if( boption("solverns2.needA1") )
    {
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
        logTime(t, "a1", ioption("solverns2.verbose") > 1);
    }

    if( boption("solverns2.needA2") )
    {
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
        logTime(t, "a2", ioption("solverns2.verbose") > 1);
    }

    return a;
}

template<typename T>
void
SolverA<T>::computeA0()
{
    // [option]
    auto g_s = soption("solverns2.alpha0");
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", doption( "solverns2.radius" ) },
                { "speed", doption( "solverns2.speed" ) } } );
    // [option]

    auto Mh = ml_space_type::New( mesh );

    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.template element<0>() ;
    auto lambda = U.template element<1>() ;
    auto v = V.template element<0>() ;
    auto nu = V.template element<1>() ;

    // [bilinearA]
    auto a = form2( _trial=Mh, _test=Mh );
    a = integrate( _range=elements(mesh),
                   _expr=inner(gradt(u),grad(v))
                   // [bilinearA]
                   // [bilinearB]
                   + id( v )*idt( lambda ) + idt( u )*id( nu )
    // [bilinearB]
                   );

    // [rhs]
    auto l = form1( _test=Mh );
    l = integrate( _range=markedfaces(mesh, 1), // inflow
                   _expr=-g*id(v) );
    l += integrate( _range=markedfaces(mesh, 2), // outflow
                    _expr=g*id(v) );
    // l += integrate( _range=markedfaces(mesh, 3), // wall
    //                 _expr=cst(0.)*id(v) );
    // [rhs]

    a.solve( _name="a0", _rhs=l, _solution=U );


    // [gradpsi0]
    auto w = Vh->element();
    a0 = Vh->element();
    auto b = form2( _trial=Vh, _test=Vh );
    b = integrate( _range=elements(mesh),
                   _expr=inner(id(w),idt(w)) );
    auto k = form1( _test=Vh );
    k = integrate( _range=elements(mesh),
                   _expr=inner(trans(gradv( U.template element<0>() )),id(w)) );
    // [gradpsi0]

    // a0 is the L2 projection of grad(psi0) over Vh
    b.solve( _name="grada0", _rhs=k, _solution=a0 );

    std::string path = "a0";
    a0.save(_path=path);
}

template<typename T>
void
SolverA<T>::loadA0()
{
    a0 = Vh->element();
    std::string path = "a0";
    a0.load(_path=path);
}

template<typename T>
void
SolverA<T>::computeA1()
{
}

template<typename T>
void
SolverA<T>::loadA1()
{
    a1 = Vh->element();
    std::string path = "a1";
    a1.load(_path=path);
}

template<typename T>
void
SolverA<T>::computeA2()
{
}

template<typename T>
void
SolverA<T>::loadA2()
{
    a2 = Vh->element();
    std::string path = "a2";
    a2.load(_path=path);
}

