﻿#include "DecimationMesh.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <queue>
typedef std::map<MyMesh::EdgeHandle, std::set<MyMesh::FaceHandle>> PlaneV;

struct DecimationMesh::CostEdge
{
	MyMesh::EdgeHandle edge;
	MyMesh::VertexHandle v0;
	MyMesh::VertexHandle v1;
	double cost;
	bool operator>(const CostEdge& other) const
	{
		return cost > other.cost;
	}
};
class DecimationMesh::PImpl
{
public:
	PImpl(MyMesh mesh)
	{
		mesh_ = mesh;
		origin_mesh_ = mesh;
		v_bar_.clear();
		avg_edge_length_ = 0.0;
	}
	PImpl()
	{
		avg_edge_length_ = 0.0;
		v_bar_.clear();
	}
	~PImpl() {}
public:
	MyMesh mesh_;
	MyMesh origin_mesh_;
	MyMesh pairMesh_;
	double avg_edge_length_;
	float t_;
	int vert_number_limit_ = 0;
	std::map<MyMesh::EdgeHandle, MyMesh::Point> v_bar_;

public:
	std::set<MyMesh::FaceHandle> getPlaneSet(MyMesh::EdgeHandle e_it);
	Eigen::MatrixXd getPlaneFunction(MyMesh::FaceHandle face);//获取法向归一化的平面方程
	void applyContraction();
	void applyContraction2();
	void getVerEdgeCost(CostEdge& ct);
	Eigen::MatrixXd pseudoInverse(Eigen::MatrixXd mat);
	bool isLegalCollapseEdge(CostEdge& ct);
};

DecimationMesh::DecimationMesh()
{
	impl_.reset(new PImpl());
}
DecimationMesh::DecimationMesh(MyMesh mesh)
{
	impl_.reset(new PImpl(mesh));
	impl_->mesh_.request_edge_status();
	impl_->mesh_.request_face_status();
	impl_->mesh_.request_vertex_status();
	impl_->mesh_.request_face_normals();
	for (auto e_it = impl_->origin_mesh_.edges_begin(); e_it != impl_->origin_mesh_.edges_end(); ++e_it)
	{
		float edge_length = impl_->origin_mesh_.calc_edge_length(*e_it);
		impl_->avg_edge_length_ += edge_length;
	}
	impl_->avg_edge_length_ /= impl_->origin_mesh_.n_edges();
}
DecimationMesh::~DecimationMesh()
{
	impl_->mesh_.release_edge_status();
	impl_->mesh_.release_face_status();
	impl_->mesh_.release_vertex_status();
	impl_->mesh_.release_face_normals();
};

void DecimationMesh::setMesh(MyMesh mesh)
{
	impl_->mesh_ = mesh;
	impl_->origin_mesh_ = mesh;
	for (auto e_it = impl_->origin_mesh_.edges_begin(); e_it != impl_->origin_mesh_.edges_end(); ++e_it)
	{
		float edge_length = impl_->origin_mesh_.calc_edge_length(*e_it);
		impl_->avg_edge_length_ += edge_length;
	}
	impl_->avg_edge_length_ /= impl_->origin_mesh_.n_edges();
	impl_->mesh_.request_edge_status();
	impl_->mesh_.request_face_status();
	impl_->mesh_.request_vertex_status();
	impl_->mesh_.request_face_normals();
}

void DecimationMesh::setT(float t)
{
	impl_->t_ = t;
	impl_->vert_number_limit_ = impl_->mesh_.n_vertices() * impl_->t_;
}

void DecimationMesh::setLocal(bool isLocal)
{

}

void DecimationMesh::applyDecimation()
{
	impl_->applyContraction();
}

MyMesh DecimationMesh::getMesh()
{
	return impl_->mesh_;
}

std::set<MyMesh::FaceHandle> DecimationMesh::PImpl::getPlaneSet(MyMesh::EdgeHandle e_it)
{
	std::set<MyMesh::FaceHandle> plane_set;
	MyMesh::HalfedgeHandle half_edge = mesh_.halfedge_handle(e_it, 0);
	MyMesh::VertexHandle v0 = mesh_.from_vertex_handle(half_edge);
	MyMesh::VertexHandle v1 = mesh_.to_vertex_handle(half_edge);
	for (auto vf_it = mesh_.vf_begin(v0); vf_it != mesh_.vf_end(v0); ++vf_it)
	{
		plane_set.insert(*vf_it);
	}
	for (auto vf_it = mesh_.vf_begin(v1); vf_it != mesh_.vf_end(v1); ++vf_it)
	{
		plane_set.insert(*vf_it);
	}
	return plane_set;
}

Eigen::MatrixXd DecimationMesh::PImpl::getPlaneFunction(MyMesh::FaceHandle face)
{
	double a, b, c = 0.0;
	double x0, y0, z0 = 0.0;
	MyMesh::Normal n = mesh_.calc_face_normal(face);
	double length = n.norm();
	/*if (length > 1e-6)
	{
		n /= length;
	}
	else
	{
		n = MyMesh::Normal(0.0, 0.0, 1.0);
	}*/
	MyMesh::FaceVertexIter fv_it = mesh_.fv_iter(face);
	MyMesh::Point point = mesh_.point(*fv_it);
	Eigen::MatrixXd pVector(4, 1);
	a = n[0];
	b = n[1];
	c = n[2];
	x0 = point[0];
	y0 = point[1];
	z0 = point[2];
	double d = -(a * x0 + b * y0 + c * z0);
	pVector(0, 0) = a;
	pVector(1, 0) = b;
	pVector(2, 0) = c;
	pVector(3, 0) = d;
	//std::cout << pVector(0, 0) <<"--" << pVector(1, 0)<<"--" << pVector(2, 0) << "--" << pVector(3, 0)<< std::endl;
	return pVector;
}

void DecimationMesh::PImpl::applyContraction()
{
	std::priority_queue<CostEdge, std::vector<CostEdge>, std::greater<CostEdge> > small_heap;
	for (auto e_it = mesh_.edges_begin(); e_it != mesh_.edges_end(); ++e_it)
	{
		MyMesh::HalfedgeHandle h0 = mesh_.halfedge_handle(*e_it, 0);
		if (!mesh_.is_valid_handle(h0))
		{
			h0 = mesh_.halfedge_handle(*e_it, 1);
		}
		if (!mesh_.is_valid_handle(h0))
		{
			continue;
		}
		CostEdge ct;
		ct.edge = *e_it;
		ct.v0 = mesh_.from_vertex_handle(h0);;
		ct.v1 = mesh_.to_vertex_handle(h0);
		ct.cost = FLT_MAX;
		getVerEdgeCost(ct);
		small_heap.push(ct);
	}

	OpenMesh::IO::write_mesh(pairMesh_, "D:/data/Decimation/pairMesh_.obj");
	int vert_num = mesh_.n_vertices();
	while (!small_heap.empty() && vert_num > vert_number_limit_)
	{
		auto ct = small_heap.top();
		small_heap.pop();
		if (!ct.edge.is_valid() || !ct.v0.is_valid() || !ct.v1.is_valid())
		{
			continue;
		}
		if (!isLegalCollapseEdge(ct))
		{
			continue;
		}
		//mesh_.request_vertex_status();
		//mesh_.request_edge_status();
		//mesh_.request_face_status();
		MyMesh::Point new_v = v_bar_[ct.edge];
		MyMesh::HalfedgeHandle heh0 = mesh_.halfedge_handle(ct.edge, 0);
		if (!mesh_.is_collapse_ok(heh0))
			continue;
		MyMesh::VertexHandle v_keep = mesh_.to_vertex_handle(heh0);
		mesh_.collapse(heh0);
		mesh_.set_point(v_keep, new_v);
		vert_num--;
		for (auto ve_it = mesh_.ve_iter(v_keep); ve_it.is_valid(); ++ve_it)
		{
			MyMesh::HalfedgeHandle h0 = mesh_.halfedge_handle(*ve_it, 0);
			if (!mesh_.is_valid_handle(h0))
			{
				h0 = mesh_.halfedge_handle(*ve_it, 1);
			}
			if (!mesh_.is_valid_handle(h0))
			{
				continue;
			}
			CostEdge ct;
			ct.edge = *ve_it;
			ct.v0 = mesh_.from_vertex_handle(h0);;
			ct.v1 = mesh_.to_vertex_handle(h0);
			ct.cost = FLT_MAX;
			getVerEdgeCost(ct);
			small_heap.push(ct);
		}
	}
	mesh_.garbage_collection();
	OpenMesh::IO::write_mesh(mesh_, "D:/data/Decimation/output.obj");
}


void DecimationMesh::PImpl::getVerEdgeCost(CostEdge& ct)
{
	if (mesh_.is_boundary(ct.edge) || mesh_.is_boundary(ct.v0) || mesh_.is_boundary(ct.v1))
		ct.cost = FLT_MAX;
	//收缩
	Eigen::MatrixXd QMat(4, 4);
	Eigen::MatrixXd zeroV(4, 1);
	Eigen::MatrixXd vBar(4, 1);

	float cur_cost = 0.0;
	QMat.setZero();
	zeroV.setZero();
	zeroV(3, 0) = 1;

	std::set<MyMesh::FaceHandle> plane_set = getPlaneSet(ct.edge);
	for (auto f_it : plane_set)
	{
		Eigen::MatrixXd pVector(4, 1);
		pVector = getPlaneFunction(f_it);
		QMat += (pVector * pVector.transpose());
	}

	vBar = pseudoInverse(QMat) * zeroV;
	if (std::abs(vBar(3, 0)) > 1e-6)
	{
		vBar /= vBar(3, 0);
	}
	else
	{
		MyMesh::Point v0 = mesh_.point(ct.v0);
		MyMesh::Point v1 = mesh_.point(ct.v1);
		MyMesh::Point mid_point = (v0 + v1) * 0.5;
		vBar(0, 0) = mid_point[0];
		vBar(1, 0) = mid_point[1];
		vBar(2, 0) = mid_point[2];
		vBar(3, 0) = 1;
	}
	v_bar_[ct.edge] = MyMesh::Point{ vBar(0, 0), vBar(1, 0), vBar(2, 0) };

	MyMesh::Point v0 = mesh_.point(ct.v0);
	MyMesh::Point v1 = mesh_.point(ct.v1);
	//std::cout <<"V0" << "X:" << v0[0] << "Y:" << v0[1] << "Z:" << v0[2] << std::endl;
	//std::cout <<"V1" << "X:" << v1[0] << "Y:" << v1[1] << "Z:" << v1[2] << std::endl;
	std::cout << "VBar" << "X:" << vBar(0, 0) << "Y:" << vBar(1, 0) << "Z:" << vBar(2, 0) << "T:" << vBar(3, 0) << std::endl;

	int vertex_num = pairMesh_.n_vertices();
	MyMesh::VertexHandle vh0 = pairMesh_.add_vertex(v0);
	MyMesh::VertexHandle vh1 = pairMesh_.add_vertex(v1);
	MyMesh::VertexHandle vh2 = pairMesh_.add_vertex(v_bar_[ct.edge]);
	pairMesh_.add_face(vh0, vh1, vh2);

	ct.cost = (vBar.transpose() * QMat * vBar)(0, 0);
}

Eigen::MatrixXd DecimationMesh::PImpl::pseudoInverse(Eigen::MatrixXd mat)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::VectorXd singularValuesInv = svd.singularValues();
	for (int i = 0; i < singularValuesInv.size(); ++i) {
		if (singularValuesInv(i) > 1e-10) {
			singularValuesInv(i) = 1.0 / singularValuesInv(i);
		}
		else {
			singularValuesInv(i) = 0.0; // 忽略过小的奇异值
		}
	}
	return svd.matrixV() * singularValuesInv.asDiagonal() * svd.matrixU().transpose();
}

bool DecimationMesh::PImpl::isLegalCollapseEdge(CostEdge& ct)
{
	if (!ct.edge.is_valid())
	{
		return false;
	}
	//判断是否合法塌缩边
	std::vector<MyMesh::VertexHandle> all_Vertex;
	std::set<MyMesh::VertexHandle> ring_Vertex;
	for (auto vv_it = mesh_.vv_iter(ct.v0); vv_it.is_valid(); vv_it++)
	{
		ring_Vertex.insert(*vv_it);
		all_Vertex.push_back(*vv_it);
	}
	for (auto vv_it = mesh_.vv_iter(ct.v1); vv_it.is_valid(); vv_it++)
	{
		ring_Vertex.insert(*vv_it);
		all_Vertex.push_back(*vv_it);
	}

	int num = all_Vertex.size() - ring_Vertex.size();
	if (num == 2)
	{
		return true;
	}
	else
	{
		return false;
	}
}