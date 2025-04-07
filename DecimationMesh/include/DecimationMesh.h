#ifndef DECIMATION_H
#define DECIMATION_H

#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
enum DecimationType
{
	EDGE_LENGTH = 0,
	CURVATURE,
	POINT_NUMBER
};
struct DecimationPara
{
	float paramT;
	bool paramLocal;
	bool paramPerBoundary;
};
class DecimationMesh
{
public:
	DecimationMesh();
	DecimationMesh(MyMesh mesh);
	~DecimationMesh();

public:
	void setMesh(MyMesh mesh);
	void setT(float t = 0);
	void setLocal(bool isLocal);
	//void setPreserveBoundarie(bool isBoundary = true);
	//void setPreserveFeature(bool isFeature = false);
	void applyDecimation();
	MyMesh getMesh();
private:
	struct CostEdge;
	class PImpl;
	std::shared_ptr<PImpl> impl_;
};
#endif