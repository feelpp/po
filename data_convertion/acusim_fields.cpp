/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/functionspace.hpp>

//#include <feel/feelvf/vf.hpp>
//#include <feel/feeldiscr/pch.hpp>

using namespace Feel;

void loadPointsMeshMap( std::string const& filename, std::map<size_type,size_type> & pointsMap )
{
    std::ifstream __is( filename.c_str() );

    size_type counterPt = 0;

    node_type coords( 3 );
    int ptid = 0;

    while ( !__is.eof() )
    {
        __is >> ptid;
        if (__is.eof() )
            break;
        __is >> coords[0] >> coords[1] >> coords[2];

        pointsMap[ptid] = counterPt++;
    }

    __is.close();
}

template<typename ElementType>
void
loadAcusimFields( std::string const& filename, std::map<size_type,size_type> const& pointsMap, ElementType & field )
{
    if ( Environment::isMasterRank() )
        std::cout << "loadAcusimFields : " << filename << " start\n";
    std::ifstream __is( filename.c_str() );
    auto space = field.functionSpace();
    auto mesh = space->mesh();
    auto doftable = space->dof();
    uint16_type nComp = doftable->nComponents;
    std::vector< std::vector<double> > dataLoaded( pointsMap.size()/*mesh->numPoints()*/, std::vector<double>( nComp, 0.) ) ;
    size_type currentIndex = 0;
    std::vector<double> valueAtNode( nComp, 0.);


    while ( !__is.eof() )
    {
        for ( uint16_type c=0 ; c<nComp ; ++c)
        {
            __is >> valueAtNode[c];
        }
        dataLoaded[currentIndex++] = valueAtNode;
        if ( currentIndex == dataLoaded.size() )
            break;
    }
    __is.close();

    const size_type nLocDof = space->dof()->nLocalDof(true);
    for (auto const& elt : elements(mesh) )
    {
        const size_type eltId = elt.id();
        for ( int j=0 ; j<nLocDof ; ++j )
        {
            auto const& thepoint = elt.point(j);
            size_type ptId = thepoint.id();
            auto itFindPt = pointsMap.find(ptId );
            CHECK( itFindPt != pointsMap.end() ) << "pt not found";
            size_type indexDataLoaded = itFindPt->second;
            for ( uint16_type comp=0 ; comp<nComp ; ++comp )
            {
                size_type dofIndex = doftable->localToGlobal( eltId, j , comp ).index();
#if 0
                auto const& nodeAtNode = doftable->dofPoint( dofIndex ).template get<0>();
                CHECK( ( std::abs(nodeAtNode[0] - thepoint[0]) < 1e-9 ) &&
                       ( std::abs(nodeAtNode[1] - thepoint[1]) < 1e-9 ) &&
                       ( std::abs(nodeAtNode[2] - thepoint[2]) < 1e-9 ) ) << "points and dof must be same : " << nodeAtNode << " vs " << thepoint;
#endif
                CHECK( indexDataLoaded >= 0 && indexDataLoaded < dataLoaded.size() ) << "error " << indexDataLoaded << " with " << ptId;

                field.set( dofIndex, dataLoaded[indexDataLoaded][comp] );
            }
        }
    }
    if ( Environment::isMasterRank() )
        std::cout << "loadAcusimFields : " << filename << " finish\n";
}

int main( int argc, char** argv )
{
    po::options_description opts ( "Convertion of Acusim fields P1");
    opts.add_options()
        ( "input.mesh.filename", po::value<std::string>()->default_value(""), "input.mesh.filename" )
        ( "input.acusim.nodes", po::value<std::string>()->default_value(""), "acusim nodes filename" )
        ( "input.acusim.pressure", po::value<std::vector<std::string> >()->multitoken(), "(vector of string) acusim pressure" )
        ( "input.acusim.velocity", po::value<std::vector<std::string> >()->multitoken(), "(vector of string) acusim velocity" )
        ( "output.directory", po::value<std::string>()->default_value("DATA_PO_CONVERTED"), "output directory" )
        ( "do-export.feel-format", po::value<bool>()->default_value(true), "do-export.feel-format" )
        ( "do-export.visu-format", po::value<bool>()->default_value(false), "do-export.visu-format" )
        ;

    Environment env( _argc=argc, _argv=argv,
                     _desc=opts,
                     _about=about( _name="convertion_acusim_fields" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    typedef Mesh<Simplex<3>> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    size_type updateComponentsMesh = MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES;
    //size_type updateComponentsMesh = MESH_UPDATE_FACES_MINIMAL;
    //size_type updateComponentsMesh = MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES;

    std::string meshFilename = soption("input.mesh.filename");
    CHECK( fs::exists( meshFilename ) ) << "mesh file does not exists : " << meshFilename;
    auto mesh = loadMesh(_mesh=new mesh_type,
                         _filename=meshFilename,
                         _update=updateComponentsMesh );

    if ( Environment::isMasterRank() )
        std::cout << "loadMesh done" << std::endl;


    std::vector<std::string> filenameFieldsPressure;
    if ( Environment::vm().count("input.acusim.pressure"))
        filenameFieldsPressure = Environment::vm()["input.acusim.pressure"].as<std::vector<std::string> >();
    std::vector<std::string> filenameFieldsVelocity;
    if ( Environment::vm().count("input.acusim.velocity"))
        filenameFieldsVelocity = Environment::vm()["input.acusim.velocity"].as<std::vector<std::string> >();

    for ( std::string const& filename : filenameFieldsPressure )
        CHECK( fs::exists( filename ) ) << "pressure file does not exists : " << filename;
    for ( std::string const& filename : filenameFieldsVelocity )
        CHECK( fs::exists( filename ) ) << "velocity file does not exists : " << filename;


    int nFieldPressure = filenameFieldsPressure.size();
    int nFieldVelocity = filenameFieldsVelocity.size();

    if ( nFieldVelocity == 0 && nFieldPressure == 0 )
    {
        if ( Environment::isMasterRank() )
            std::cout << "nothing to export : exit\n";
        return 0;
    }

    std::string filenameNodes = soption(_name="input.acusim.nodes");
    CHECK( fs::exists( filenameNodes ) ) << "acusim nodes filename does not exists : " << filenameNodes;
    std::map<size_type,size_type> pointsMap;
    loadPointsMeshMap( filenameNodes, pointsMap );


    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Vectorial,Continuous> > > space_velocity_type;
    typedef boost::shared_ptr<space_velocity_type> space_velocity_ptrtype;
    typedef typename space_velocity_type::component_functionspace_type space_pressure_type;
    typedef boost::shared_ptr<space_pressure_type> space_pressure_ptrtype;

    space_velocity_ptrtype VhVelocity;
    space_pressure_ptrtype VhPressure;

    if ( nFieldVelocity > 0 )
    {
        VhVelocity = space_velocity_type::New( _mesh=mesh );
        if ( nFieldPressure > 0 )
            VhPressure = VhVelocity->compSpace();
    }
    else if ( nFieldPressure > 0 )
    {
        VhPressure = space_pressure_type::New( _mesh=mesh );
    }


    std::vector<typename space_pressure_type::element_ptrtype> fieldsPressure( nFieldPressure );
    for (int k=0;k<nFieldPressure;++k)
    {
        fieldsPressure[k] = VhPressure->elementPtr();
        loadAcusimFields( filenameFieldsPressure[k], pointsMap, *fieldsPressure[k] );
    }

    std::vector<typename space_velocity_type::element_ptrtype> fieldsVelocity( nFieldVelocity );
    for (int k=0;k<nFieldVelocity;++k)
    {
        fieldsVelocity[k] = VhVelocity->elementPtr();
        loadAcusimFields( filenameFieldsVelocity[k], pointsMap, *fieldsVelocity[k] );
    }
    std::string outputDir = soption(_name="output.directory");
    if ( outputDir.empty() )
        outputDir = fs::current_path().string();
    else if ( fs::path(outputDir).is_relative() )
        outputDir = (fs::path(Feel::Environment::rootRepository())/fs::path(outputDir)).string();

    if ( Environment::worldComm().isMasterRank() )
        if ( !fs::exists( outputDir ) )
            fs::create_directories( outputDir );
    Environment::worldComm().globalComm().barrier();

    if ( boption(_name="do-export.feel-format" ) )
    {
        for ( int k=0;k<nFieldPressure;++k )
        {
            std::string nameFieldExport = fs::path(filenameFieldsPressure[k]).stem().string();
            std::string pathFieldExport = (fs::path(outputDir)/fs::path("feel-format")/fs::path(nameFieldExport)).string();
            if ( Environment::isMasterRank() )
                std::cout << "save pressure field in feel-format : " << pathFieldExport << "\n";
            fieldsPressure[k]->save( _path=pathFieldExport,_type="hdf5" );
        }
        for ( int k=0;k<nFieldVelocity;++k )
        {
            std::string nameFieldExport = fs::path(filenameFieldsVelocity[k]).stem().string();
            std::string pathFieldExport = (fs::path(outputDir)/fs::path("feel-format")/fs::path(nameFieldExport)).string();
            if ( Environment::isMasterRank() )
                std::cout << "save velocity field in feel-format : " << pathFieldExport << "\n";
            fieldsVelocity[k]->save( _path=pathFieldExport,_type="hdf5" );
        }
    }

    if ( boption(_name="do-export.visu-format") )
    {
        std::string nameExporter = fs::path(meshFilename).stem().string();
        std::string pathExporter = (fs::path(outputDir)/fs::path("visu-format")).string();
        if ( Environment::isMasterRank() )
            std::cout << "save accusim fields in visu-format : " << pathExporter << "\n";
        auto e = exporter( _mesh=mesh,_name=nameExporter,_path=pathExporter );
        for (int k=0;k<nFieldPressure;++k)
        {
            std::string nameExport = fs::path(filenameFieldsPressure[k]).stem().string();
            e->add( "pressure_"+nameExport, *fieldsPressure[k] );
        }
        for (int k=0;k<nFieldVelocity;++k)
        {
            std::string nameExport = fs::path(filenameFieldsVelocity[k]).stem().string();
            e->add( "velocity_"+nameExport, *fieldsVelocity[k] );
        }
        e->save();
    }

    return 0;
}
