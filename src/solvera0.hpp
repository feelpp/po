using namespace Feel;

template<typename FunctionSpaceType>
class SolverA0
{
    typedef double value_type;
    typedef Mesh<Simplex<3> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpaceType space_ptrtype;
    typedef typename space_ptrtype::element_type space_type;
    typedef typename space_type::element_type element_type;

    typedef SolverA0<space_ptrtype> solvera0_type;
    typedef typename boost::shared_ptr<solvera0_type> solvera0_ptrtype;

    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar>, Lagrange<0, Scalar> > > ml_space_type;

    mesh_ptrtype mesh;
    space_ptrtype Vh;
    element_type a0;

    void computeA0();
    void loadA0();

public:
    static solvera0_ptrtype build(const mesh_ptrtype& mesh, const space_ptrtype& Vh);
    element_type solve();
};

template<typename T>
typename SolverA0<T>::solvera0_ptrtype
SolverA0<T>::build(const mesh_ptrtype& mesh, const space_ptrtype& Vh)
{
    solvera0_ptrtype solvera0( new SolverA0<T> );
    solvera0->mesh = mesh;
    solvera0->Vh = Vh;
    return solvera0;
}

template<typename T>
typename SolverA0<T>::element_type
SolverA0<T>::solve()
{
    if( boption("computeA0"))
        computeA0();
    else
        loadA0();

    return a0;
}

template<typename T>
void
SolverA0<T>::computeA0()
{
    // [option]
    auto g_s = soption("alpha0");
    auto vars = Symbols{ "x", "y", "radius", "speed" };
    auto g_e = parse( g_s, vars );
    auto g = expr( g_e, vars );
    g.setParameterValues( {
            { "radius", doption( _name="radius" ) },
                { "speed", doption( _name="speed" ) } } );
    // [option]

    std::cout << g_s << std::endl;

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
    l += integrate( _range=markedfaces(mesh, 3), // wall
                    _expr=cst(0.)*id(v) );
    // [rhs]

    a.solve( _name="a0", _rhs=l, _solution=U );


    // [gradpsi0]
    auto w = Vh->element();
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
SolverA0<T>::loadA0()
{
    std::string path = "a0";
    a0.load(_path=path);
}
