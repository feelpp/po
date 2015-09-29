#include <Eigen/Dense>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Eigen;

template<typename EigenModes>
class InitCoeff
{
    using value_type = double;
    using mesh_type = Mesh<Simplex<3> >;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    using eigenmodes_type = EigenModes;
    using eigentuple_type = typename eigenmodes_type::value_type;
    using eigen_vec_type = typename std::tuple_element<1, eigentuple_type>::type;
    using eigen_space_ptrtype = typename eigen_vec_type::functionspace_ptrtype;
    using eigenvec_type = std::vector<eigen_vec_type>;

    using initcoeff_type = InitCoeff<eigenmodes_type>;
    using initcoeff_ptrtype = typename boost::shared_ptr<initcoeff_type>;

    mesh_ptrtype mesh;

    eigen_space_ptrtype Nh;
    eigenvec_type g;

    int M;
    Matrix<MatrixXd, Dynamic, 1 > Rijk;

    void write_binary_rijk(const char* filename);

public:
    static initcoeff_ptrtype build(const eigenmodes_type& modes);
    void initRijk();
};

template<typename E>
typename InitCoeff<E>::initcoeff_ptrtype
InitCoeff<E>::build(const eigenmodes_type& modes)
{
    initcoeff_ptrtype ap( new InitCoeff<E> );

    ap->M = modes.size();
    ap->g = eigenvec_type(ap->M);

    int i = 0;
    for(auto const& tuple : modes)
        ap->g[i++] = std::get<1>(tuple);

    ap->Nh = ap->g[0].functionSpace();
    ap->mesh = ap->Nh->mesh();

    return ap;
}

template<typename E>
void
InitCoeff<E>::write_binary_rijk(const char* filename)
{
    if( Environment::isMasterRank() )
    {
        std::ofstream out(filename,std::ios::out | std::ios::binary | std::ios::trunc);
        out.write((char*) (&M), sizeof(typename decltype(Rijk)::Index));
        for( int k = 0; k < M; k++ )
            out.write((char*) Rijk(k).data(), M*M*sizeof(typename decltype(Rijk)::Scalar::Scalar) );
        out.close();
    }
}

template<typename E>
void
InitCoeff<E>::initRijk()
{
    tic();
    if ( ioption("offline.verbose") > 2 && Environment::isMasterRank() )
        std::cout << "----- Rijk -----" << std::endl;

    tic();
    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

    auto w = Nh->element();
    auto rijkForm = form2( _test=Nh, _trial=Nh );

    // Rijk = - Rikj && Rijj = 0
    for (int k = 0; k < M; k++)
    {
        rijkForm = integrate( elements(mesh),
                              inner( cross( curl(w),idt(w) ), idv(g[k]) ));
        rijkForm.close();

        for(int i = 0; i < M; i++)
        {
            for(int j = 0; j < k; j++)
            {
                Rijk(k)(i,j) = rijkForm( g[i], g[j] );
                Rijk(j)(i,k) = -Rijk(k)(i,j);

                if ( ioption("offline.verbose") > 3 && Environment::isMasterRank() )
                    std::cout << "Rijk(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
            }
            Rijk(k)(i,k) = 0;
        }
    }
    toc("calcul rijk", ioption("offline.verbose") > 2);

    tic();
    write_binary_rijk( "rijk" );
    toc("write rijk", ioption("offline.verbose") > 1);

    toc( "Rijk", ioption("offline.verbose") > 1);
}
