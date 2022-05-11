#include "mesh_def.h"

typedef typename PolyMesh::ScalarType ScalarType;
typedef vcg::Point3<ScalarType> CoordType; 

ScalarType EdgeDev(PolyMesh& mesh)
{
    std::set<std::pair<size_t, size_t>> edgeSet;
    for (size_t i = 0; i < mesh.fn; i++)
    {
        if (mesh.face[i].VN() != 4)
        {
            std::cout << "In Edge: mesh not pure quad\n";
        }
        for (size_t j = 0; j < 4; j++)
        {
            size_t indexV0 = vcg::tri::Index(mesh, mesh.face[i].V0(j));
            size_t indexV1 = vcg::tri::Index(mesh, mesh.face[i].V1(j));
            
            edgeSet.insert(std::pair<size_t, size_t>(std::min(indexV0, indexV1), std::max(indexV0, indexV1)));
        }
    }

    ScalarType avgEdge = 0;
    for (auto edge: edgeSet)
    {
        CoordType vec = mesh.vert[edge.first].P() - mesh.vert[edge.second].P();
        avgEdge += vec.Norm();
    }
    avgEdge /= edgeSet.size();
    
    ScalarType avg = 0;
    for (auto edge: edgeSet)
    {
        CoordType vec = mesh.vert[edge.first].P() - mesh.vert[edge.second].P();
        ScalarType length = vec.Norm();
        avg += abs(length - avgEdge);
    }

    return avg / (edgeSet.size() * avgEdge);
}

ScalarType AngleDev(PolyMesh& mesh)
{
    ScalarType avg = 0;
    size_t count = 0;

    for (size_t i = 0; i < mesh.fn; i++)
    {
        if (mesh.face[i].VN() != 4)
        {
            std::cout << "In Angle: mesh not pure quad\n";
            continue;
        }
        CoordType edge[4];
        for (size_t j = 0; j < 4; j++)
        {
            edge[j] = (mesh.face[i].P1(j)-mesh.face[i].P0(j)).normalized();
        }
        for (size_t j = 0; j < 4; j++)
        {
            ScalarType cosvalue = -edge[j] * edge[(j+1)%4];
            if (cosvalue < -1) cosvalue = -1;
            if (cosvalue > 1) cosvalue = 1;
            
            ScalarType angle = acos(cosvalue) * 180 / M_PI;
            ScalarType bias = abs(angle - 90);

            avg += bias;
            count++;
        }
    }

    return avg / count;
}

ScalarType IrrNum(PolyMesh& mesh)
{
    std::vector<int> ncount(mesh.vert.size(), 0);
    for (size_t i = 0; i < mesh.fn; i++)
    {
        if (mesh.face[i].VN() != 4)
            std::cout << "In Irr: mesh not pure quad\n";
        for (size_t j = 0; j < 4; j++)
        {
            size_t indexV = vcg::tri::Index(mesh, mesh.face[i].V(j));
            ncount[indexV]++;
        }
    }
    int count = 0;
    for (size_t i = 0; i < ncount.size(); i++)
    {
        if (!mesh.vert[i].IsB() && ncount[i] != 4)
            count++;
    }
    return count;
}

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Must input a mesh's path as arg...\n";
        return 0;
    }

    std::string pathM = std::string(argv[1]);
    std::cout << "Loading: " << pathM.c_str() << std::endl;

    PolyMesh quadMesh;
    int mask=vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
	int err = vcg::tri::io::ImporterOBJ<PolyMesh>::Open(quadMesh, pathM.c_str(), mask);
	if (err)
	{
		std::cerr << "Import Error: " << vcg::tri::io::ImporterOBJ<PolyMesh>::ErrorMsg(err) << std::endl;
		return err;
	}

    quadMesh.UpdateDataStructures();

    if (quadMesh.fn != quadMesh.face.size())
        std::cout << "Warning: face number not match\n";
    if (quadMesh.vn != quadMesh.vert.size())
        std::cout << "Warning: vertex number not match\n";

    std::cout << "Num of face: " << quadMesh.fn << std::endl;
    std::cout << "Num of vertices: " << quadMesh.vn << std::endl;
    std::cout << "NUm of irregular: " << IrrNum(quadMesh) << std::endl;
    std::cout << "Angle deviation: " << AngleDev(quadMesh) << std::endl;
    std::cout << "Edge deviation: " << EdgeDev(quadMesh) << std::endl;

    return 0;   
}