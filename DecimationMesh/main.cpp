#include <iostream>
#include <DecimationMesh.h>

int main()
{
    for (float i = 0.5f; i < 0.6f; i = i + 0.1f)
    {
        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "Model/tooth.stl");
        DecimationMesh dm(mesh);
        dm.setT(i);
        dm.applyDecimation();
        MyMesh output_mesh = dm.getMesh();
        std::string output_file = "Model/output_" + std::to_string(i) + ".obj";
        OpenMesh::IO::write_mesh(output_mesh, output_file);
    }
}
