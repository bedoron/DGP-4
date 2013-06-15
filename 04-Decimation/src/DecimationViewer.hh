//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS QualityViewer
//
//=============================================================================


#ifndef DECIMATIONVIEWERWIDGET_HH
#define DECIMATIONVIEWERWIDGET_HH


//== INCLUDES =================================================================


#include "MeshViewer.hh"
#include "QuadricT.hh"
#include <set>
#include <vector>
using std::vector;
//== CLASS DEFINITION =========================================================


class DecimationViewer : public MeshViewer
{

protected:

//	typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
	typedef vector<Mesh::VertexHandle> VList;
	typedef vector<Mesh::EdgeHandle> EList;
	int percentage_;
public:

	/// default constructor
	DecimationViewer(const char* _title, int _width, int _height);
	~DecimationViewer() {};
	virtual void keyboard(int key, int x, int y);

	//mesh decimation actions
	void  init();
	bool  is_collapse_legal(Mesh::HalfedgeHandle _hh);
	bool  is_collapse_legal2(Mesh::HalfedgeHandle _hh);

	float priority(Mesh::HalfedgeHandle _heh);
	float priority2(Mesh::HalfedgeHandle _heh);

//	bool  is_collapse_legal(Mesh::EdgeHandle _eh);
//	float priority(Mesh::EdgeHandle _eh);

	void  decimate(unsigned int _n_vertices);
	void  blowup(unsigned int _n_vertices);

	void  enqueue_vertex(Mesh::VertexHandle vh);
	void  enqueue_vertex2(Mesh::VertexHandle vh);
//	void enqueue_edge(Mesh::EdgeHandle eh);

	//bool  is_collapse_legal_special(Mesh::HalfedgeHandle _hh);
	//float priority(Mesh::HalfedgeHandle _heh);
	//void  decimate(unsigned int _n_vertices);
	//void  enqueue_vertex(Mesh::VertexHandle vh);

	// access quadric of vertex _vh
	Quadricd& quadric(Mesh::VertexHandle _vh)  { return mesh_.property(vquadric, _vh); }
//	Quadricd quadric(Mesh::EdgeHandle _eh);

	// access priority of vertex _vh
	float& priority(Mesh::VertexHandle _vh)				{ return mesh_.property(vprio, _vh); }

	// access target halfedge of vertex _vh
	Mesh::HalfedgeHandle& target(Mesh::VertexHandle _vh) { return mesh_.property(vtarget, _vh); }
//	Vec3f& target(Mesh::EdgeHandle _eh) { return mesh_.property(starget, _eh); }

	void collapse_edge_and_move_to(Mesh::EdgeHandle eh, Mesh::Point newPosition);
	//Vec3f calculateOptimalPoint(Mesh::EdgeHandle eh);
	Vec3f calculateOptimalPoint(Quadricd q0, Quadricd q1);

	//inline void edgeToPoints(Mesh::EdgeHandle _eh, Mesh::VertexHandle &p0, Mesh::VertexHandle &p1) {
	//		Mesh::HalfedgeHandle _heh = mesh_.halfedge_handle(_eh,0);
	//		p0 = mesh_.from_vertex_handle(_heh);
	//		p1 = mesh_.to_vertex_handle(_heh);
	//};

	//void getOneRing(Mesh::EdgeHandle _eh, VList ringers); // not used
	//void getOneRing(Mesh::EdgeHandle _eh, EList ringers);

protected:

//	Mesh                       mesh_; // Y U DO DIS ?!
	std::vector<unsigned int>  indices_;

	OpenMesh::VPropHandleT<Quadricd>                vquadric;
	OpenMesh::VPropHandleT<float>                   vprio;
	OpenMesh::VPropHandleT<Mesh::HalfedgeHandle>   vtarget;

	OpenMesh::EPropHandleT<Vec3f>					starget;

	struct QueueVertex{
		Mesh::VertexHandle v;
		float prio;
	};

	// compare functor for priority queue
	struct VertexCmp
	{

		bool operator()(QueueVertex _v0, QueueVertex _v1) const
		{
			// std::set needs UNIQUE keys -> handle equal priorities
			return (( _v0.prio ==  _v1.prio) ? 
				(_v0.v.idx() < _v1.v.idx()) :
			(  _v0.prio <   _v1.prio));
		}
	};


	std::set<QueueVertex, VertexCmp>  queue;


	struct QEdge {
		Mesh::EdgeHandle eh;
		float prio;

		bool operator<(const QEdge &rhs) const {
			if(prio == rhs.prio)
				return eh.idx() < rhs.eh.idx();
			return prio < rhs.prio;
		}
	};

	std::set<QEdge> qEdge;
};


//=============================================================================
#endif // QUALVIEWERWIDGET_HH defined
//=============================================================================

