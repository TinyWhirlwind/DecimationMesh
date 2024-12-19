#include <iostream>
#include <DecimationMesh.h>

int main()
{
    for (float i = 0.1f; i < 0.2f; i = i + 0.1f)
    {
        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "D:/data/Decimation/toothHole.stl");
        DecimationMesh dm(mesh);
        dm.setT(i);
        dm.applyDecimation();
        MyMesh output_mesh = dm.getMesh();
        std::string output_file = "D:/data/Decimation/output_" + std::to_string(i) + ".stl";
        OpenMesh::IO::write_mesh(output_mesh, output_file);
    }
}
