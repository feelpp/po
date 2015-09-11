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

    return ap;
}

template<typename E>
void
InitCoeff<E>::initRijk()
{
    tic();
    if ( Environment::isMasterRank() && ioption("offline.verbose") > 2 )
        std::cout << "----- Rijk -----" << std::endl;

    Rijk = Matrix<MatrixXd, Dynamic, 1>(M,1);
    for(int k = 0; k < M; k++)
        Rijk(k) = MatrixXd(M,M);

    auto w = Nh->element();

    if( boption("coeff.computeRijk") )
    {
        std::fstream s;
        if ( Environment::isMasterRank() )
            s.open ("rijk", std::fstream::out);

        // Rijk = - Rikj && Rijj = 0
        for(int k = 0; k < M; k++)
        {
            auto rijkForm = form2( _test=Nh, _trial=Nh );
            rijkForm = integrate( elements(mesh),
                                  inner( cross( curl(w),idt(w) ), idv(g[k]) ));
            for(int i = 0; i < M; i++)
            {
                for(int j = 0; j < k; j++)
                {
                    Rijk(k)(i,j) = rijkForm( g[i], g[j] );
                    Rijk(j)(i,k) = -Rijk(k)(i,j);

                    if ( Environment::isMasterRank() )
                    {
                        if( ioption("offline.verbose") > 3 )
                            std::cout << "Rijk(" << k << "," << i << "," << j << ") = " << Rijk(k)(i,j) << std::endl;
                        s << Rijk(k)(i,j) << std::endl;
                    }
                }
                Rijk(k)(i,k) = 0;
            }
        }
        if ( Environment::isMasterRank() )
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
    toc( "Rijk", ioption("offline.verbose") > 1);

}
