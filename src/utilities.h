#include <boost/mpi/timer.hpp>

using namespace Feel;

void
logTime(boost::mpi::timer& t, std::string s, bool verbose)
{
    LOG(INFO) << "[timer] " << s << " = " << t.elapsed() << std::endl;
    if ( verbose && Environment::worldComm().isMasterRank() )
        std::cout << "[timer] " << s << " = " << t.elapsed() << std::endl;
    t.restart();
}

template<typename MeshPtrType >
std::ostream&
logMesh(std::ostream& out, const MeshPtrType& mesh)
{
    if(Environment::isMasterRank())
    {
        out << "[mesh]  - mesh entities" << std::endl
            << "[mesh]  number of elements : " << mesh->numGlobalElements() << std::endl
            << "[mesh]  number of faces : " << mesh->numGlobalFaces() << std::endl
            << "[mesh]  number of edges : " << mesh->numGlobalEdges() << std::endl
            << "[mesh]  number of points : " << mesh->numGlobalPoints() << std::endl
            << "[mesh]  number of vertices : " << mesh->numGlobalVertices() << std::endl
            << "[mesh]  - mesh sizes" << std::endl
            << "[mesh]  h max : " << mesh->hMax() << std::endl
            << "[mesh]  h min : " << mesh->hMin() << std::endl
            << "[mesh]  h avg : " << mesh->hAverage() << std::endl
            << "[mesh]  measure : " << mesh->measure() << std::endl;

        for( auto marker: mesh->markerNames() )
            out << "[mesh]  - Marker " << marker.first << std::endl;
    }

    return out;
}
