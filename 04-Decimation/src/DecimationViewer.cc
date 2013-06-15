//=============================================================================
//                                                
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Applied Geometry Group and 
//							Computer Graphics Laboratory, ETH Zurich
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
//  CLASS ReconViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>

#include <vector>
#include <float.h>
#include <stack>

#include <windows.h>
#include "DecimationViewer.hh"

using std::stack;

#undef max
#undef min
#include <gmm.h>

typedef gmm::dense_matrix<double>  gmmMatrix;
typedef std::vector<double>        gmmVector;

static 	void solve_linear_system( gmmMatrix& _M, 
								 gmmVector& _b, 
								 gmmVector& _x ) {
	unsigned int N = _b.size();
	_x.resize(N);
	std::vector< size_t >  ipvt(N);
	gmm::lu_factor( _M, ipvt );
	//try {
		gmm::lu_solve( _M, ipvt, _x, _b );
	//} catch(...) {
	//	cout << "Something is fishy :(\n";
	//}
}

//== IMPLEMENTATION ========================================================== 

DecimationViewer::
	DecimationViewer(const char* _title, int _width, int _height) : 
	MeshViewer(_title, _width, _height) {

	MeshViewer::init();
	mesh_.add_property(vprio);
	mesh_.add_property(vquadric);
	mesh_.add_property(vtarget);
	mesh_.add_property(starget);
	mesh_.request_edge_status();
	mesh_.request_vertex_status();
	mesh_.request_face_status();
	mesh_.request_halfedge_status();

	init();
};
//-----------------------------------------------------------------------------




//-----------------------------------------------------------------------------


void DecimationViewer::keyboard(int key, int x, int y) 
{
	switch (key)
	{
	case 'o':
		{
			OPENFILENAME ofn={0};
			char szFileName[MAX_PATH]={0};
			ofn.lStructSize=sizeof(OPENFILENAME);
			ofn.Flags=OFN_ALLOWMULTISELECT|OFN_EXPLORER;
			ofn.lpstrFilter="All Files (*.*)\0*.*\0";
			ofn.lpstrFile=szFileName;
			ofn.nMaxFile=MAX_PATH;
			if(GetOpenFileName(&ofn))
			{
				mesh_.clear();
				MeshViewer::open_mesh(szFileName);
				init();
			}
		}
		break;
	case 'z':  //up percentage by 5%
		percentage_=percentage_+5;
		if (percentage_>95) percentage_=95;
		std::cout<<"Percentage is %"<<percentage_<<"\n";
		break;

	case 'x':  //up percentage by 5%
		percentage_=percentage_-5;
		if (percentage_<5) percentage_=5;
		std::cout<<"Percentage is %"<<percentage_<<"\n";
		break;

	case 'd': //decimate
		// compute normals & quadrics
		init();

		// decimate
		decimate((percentage_/100.0)*mesh_.n_vertices());
		std::cout << "#vertices: " << mesh_.n_vertices() << std::endl;

		break;
	case 'b': // blowup!
		
		blowup((percentage_/100.0)*mesh_.n_vertices());
		std::cout << "#vertices: " << mesh_.n_vertices() << std::endl;
		break;
	default:
		GlutExaminer::keyboard(key, x, y);
		break;
	}
}


//=============================================================================
//DECIMATION IMPLEMENTATION FUNCTIONS
//==============================================================================
void DecimationViewer::init()
{
	// compute face normals
	mesh_.update_face_normals();


	Mesh::VertexIter  v_it, v_end = mesh_.vertices_end();
	Mesh::Point       n;
	double              a, b, c, d;

	for (v_it=mesh_.vertices_begin(); v_it != v_end; ++v_it)
	{
		priority(v_it) = -1.0;
		quadric(v_it).clear();


		// Exercise 4.1 --------------------------------------
		// INSERT CODE:
		// calc vertex quadrics from incident triangles
		// ---------------------------------------------------
		
		Quadricd qd(0,0,0,0);
		qd.clear(); // Just in case...

		for(Mesh::VFIter vf_iter = mesh_.vf_iter(v_it); vf_iter; ++vf_iter) {
			Vec3f fnormal = mesh_.normal(vf_iter.handle());
			float d_i = -(mesh_.point(v_it)|fnormal);
			qd += Quadricd(fnormal[0],fnormal[1],fnormal[2], d_i);
		}
		quadric(v_it) = qd;
	}
}


//-----------------------------------------------------------------------------


bool DecimationViewer::is_collapse_legal(Mesh::HalfedgeHandle _hh)
{
	// collect vertices
	Mesh::VertexHandle v0, v1;
	v0 = mesh_.from_vertex_handle(_hh);
	v1 = mesh_.to_vertex_handle(_hh);


	// collect faces
	Mesh::FaceHandle fl = mesh_.face_handle(_hh);
	Mesh::FaceHandle fr = mesh_.face_handle(mesh_.opposite_halfedge_handle(_hh));


	// backup point positions
	Mesh::Point p0 = mesh_.point(v0);
	Mesh::Point p1 = mesh_.point(v1);


	// topological test
	if (!mesh_.is_collapse_ok(_hh))
		return false;


	// test boundary stuff
	if (mesh_.is_boundary(v0) && !mesh_.is_boundary(v1))
		return false;


	// Exercise 4.2 -----------------------------------------------
	// INSERT CODE:
	// test normal flipping:
	//   if normal vector of a (non-degenerate) triangle changes by 
	//   more than pi/4 degrees, return false.
	// ------------------------------------------------------------
	for(Mesh::VFIter taround = mesh_.vf_iter(v0); taround; ++taround) {
		if(!taround.handle().is_valid()) continue;
		if(taround.handle() == fl || taround.handle() == fr) continue;
		Mesh::FaceVertexIter fv_iter = mesh_.fv_iter(taround);
		vector<Mesh::Point> pts(2); 
		for(; fv_iter; ++fv_iter) {
			if(fv_iter.handle() == v0) continue;
			pts.push_back(mesh_.point(fv_iter));
		}
		Vec3f origin_normal = ((pts[0] - p0)%(pts[1] - p0)).normalize();
		Vec3f target_normal = ((pts[0] - p1)%(pts[1] - p1)).normalize();
		float projection = abs(origin_normal|target_normal);
		if(projection<(0.70710678118)) // 1/sqrt(2)
			return false; 
	}
	return true;
}

bool DecimationViewer::is_collapse_legal2(Mesh::HalfedgeHandle _hh)
{
	// collect vertices
	Mesh::VertexHandle v0, v1;
	v0 = mesh_.from_vertex_handle(_hh);
	v1 = mesh_.to_vertex_handle(_hh);


	// collect faces
	Mesh::FaceHandle fl = mesh_.face_handle(_hh);
	Mesh::FaceHandle fr = mesh_.face_handle(mesh_.opposite_halfedge_handle(_hh));


	// backup point positions
	Mesh::Point p0 = mesh_.point(v0);
	Mesh::Point p1 = mesh_.point(v1);


	// topological test
	if (!mesh_.is_collapse_ok(_hh))
		return false;


	// test boundary stuff
	Vec3f newPoint = calculateOptimalPoint(quadric(v0),quadric(v1));

	if (mesh_.is_boundary(v0) && !mesh_.is_boundary(v1))
		return false;

	for(Mesh::VertexHandle curr = v0 ; curr != v1; curr = v1) { // Cycle points

		for(Mesh::VFIter taround = mesh_.vf_iter(curr); taround; ++taround) {
			if(!taround.handle().is_valid()) continue;
			if(taround.handle() == fl || taround.handle() == fr) continue;
			Mesh::FaceVertexIter fv_iter = mesh_.fv_iter(taround);
			vector<Mesh::Point> pts(2); 
			for(; fv_iter; ++fv_iter) {
				if(mesh_.point(fv_iter) == newPoint) continue;
				pts.push_back(mesh_.point(fv_iter));
			}
			Mesh::Point curr_point = mesh_.point(curr);
			Vec3f origin_normal = ((pts[0] - curr_point)%(pts[1] - curr_point)).normalize();
			Vec3f target_normal = ((pts[0] - newPoint)%(pts[1] - newPoint)).normalize();
			float projection = abs(origin_normal|target_normal);
			if(projection<(0.70710678118)) // 1/sqrt(2)
				return false; 
		}
	}
	return true;

}


//
//void DecimationViewer::getOneRing(Mesh::EdgeHandle _eh, EList ringers) {
//	VertexHandle p0, p1;
//	edgeToPoints(_eh, p0, p1);
//
//	for(VertexHandle vh = p0; vh != p1; vh = p1) {
//		for(Mesh::VertexEdgeIter oneRing = mesh_.ve_iter(vh); oneRing; ++oneRing) {
//			if(oneRing.handle() == _eh) 
//				continue;
//			ringers.push_back(oneRing);
//		}
//	}
//}
//
//void DecimationViewer::getOneRing(Mesh::EdgeHandle _eh, VList ringers) {
//	VertexHandle p0, p1;
//	edgeToPoints(_eh, p0, p1);
//
//	for(VertexHandle vh = p0; vh != p1; vh = p1) {
//		for(Mesh::VVIter oneRing = mesh_.vv_iter(vh); oneRing ; ++oneRing) {
//			if(mesh_.vertex_handle(oneRing) == p0 || mesh_.vertex_handle(oneRing) == p1)
//				continue;
//			ringers.push_back(mesh_.vertex_handle(oneRing));
//		}
//	}
//}
// TODO _ BLA
//bool DecimationViewer::is_collapse_legal(Mesh::EdgeHandle _eh) {
//	VertexHandle p0, p1;
//	Vec3f dest = target(_eh);
//	edgeToPoints(_eh, p0, p1);
//
//	
//	if (!mesh_.is_collapse_ok(mesh_.halfedge_handle(_eh,0)))
//		return false;
//
//	if (!mesh_.is_collapse_ok(mesh_.halfedge_handle(_eh,1)))
//		return false;
//
//	// test boundary stuff
//	if (mesh_.is_boundary(p0) && !mesh_.is_boundary(p1))
//		return false;
//
//	// Collect "1Ring"
//	std::vector<Mesh::Point> rings;
//	for(VertexHandle vh = p0; vh != p1; vh = p1) {
//		for(Mesh::VVIter oneRing = mesh_.vv_iter(vh); oneRing ; ++oneRing) {
//			if(mesh_.vertex_handle(oneRing) == p0 || mesh_.vertex_handle(oneRing) == p1)
//				continue;
//			rings.push_back(mesh_.point(oneRing));
//		}
//	}
//
//	std::vector< std::vector<Vec3f> > edges;
//	std::vector<Mesh::Point>::iterator vit;
//	int valence = rings.size();
//
//	for(vit = rings.begin(); vit != rings.end(); ++vit) {
//		edges[0].push_back(*vit - mesh_.point(p0).normalize());
//		edges[1].push_back(*vit - mesh_.point(p1).normalize());
//		edges[2].push_back((*vit - dest).normalize());
//	}
//
//	for(int i=0; i < valence; ++i) {
//		Vec3f plane_normal_0 = (edges[0][i]%edges[0][(i+1)%valence]).normalize(); 
//		Vec3f plane_normal_1 = (edges[1][i]%edges[1][(i+1)%valence]).normalize(); 
//		Vec3f target_normal =  (edges[2][i]%edges[2][(i+1)%valence]).normalize(); 
//
//		float projection_0 = abs(plane_normal_0|target_normal);
//		float projection_1 = abs(plane_normal_1|target_normal);
//
//		if( (projection_0 < (0.70710678118) ) || (projection_1 < (0.70710678118)))
//			return false;
//	}
//}

//-----------------------------------------------------------------------------
//Quadricd DecimationViewer::quadric(Mesh::EdgeHandle _eh) {
//	Mesh::VertexHandle p0;
//	Mesh::VertexHandle p1;
//
//	edgeToPoints(_eh, p0, p1);
//
//	Quadricd errors(0,0,0);
//	errors += quadric(p0);
//	errors += quadric(p1);
//
//	return errors;
//}

float DecimationViewer::priority(Mesh::HalfedgeHandle _heh)
{
	// Exercise 4.3 ----------------------------------------------
	// INSERT CODE:
	// return priority: the smaller the better
	// use quadrics to estimate approximation error
	// -----------------------------------------------------------
	Mesh::VertexHandle v0, v1;
	v0 = mesh_.from_vertex_handle(_heh);
	v1 = mesh_.to_vertex_handle(_heh);

	Mesh::Point p0 = mesh_.point(v0);
	Mesh::Point p1 = mesh_.point(v1);

	return quadric(v0)(p1)+quadric(v1)(p1);
}

float DecimationViewer::priority2(Mesh::HalfedgeHandle _heh)
{
	Mesh::VertexHandle v0, v1;
	v0 = mesh_.from_vertex_handle(_heh);
	v1 = mesh_.to_vertex_handle(_heh);

	Mesh::Point p1 = calculateOptimalPoint(quadric(v0), quadric(v1));

	return quadric(v0)(p1)+quadric(v1)(p1);
}
//float DecimationViewer::priority(Mesh::EdgeHandle _eh) {
//	Vec3f new_point = target(_eh);
//	Quadricd errors = quadric(_eh);
//	return errors(new_point);
//}

Vec3f DecimationViewer::calculateOptimalPoint(Quadricd q0, Quadricd q1) {
	Quadricd errors;
	errors.clear();
	errors += q0;
	errors += q1;

	gmmMatrix qmat(4,4);
	qmat(0,0) = errors.a; qmat(0,1) = errors.b; qmat(0,2) = errors.c; qmat(0,3) = errors.d;
	qmat(1,0) = errors.b; qmat(1,1) = errors.e; qmat(1,2) = errors.f; qmat(1,3) = errors.g;
	qmat(2,0) = errors.c; qmat(2,1) = errors.f; qmat(2,2) = errors.h; qmat(2,3) = errors.i;
	qmat(3,0) = 0; qmat(3,1) = 0; qmat(3,2) = 0; qmat(3,3) = 1;

	gmmVector constraints(4);
	constraints[0] = 0;
	constraints[1] = 0;
	constraints[2] = 0;
	constraints[3] = 1;

	gmmVector solution;

	solve_linear_system(qmat, constraints, solution);
	Vec3f new_point(solution[0]/solution[3], solution[1]/solution[3], solution[2]/solution[3]);
	return new_point;
}

void DecimationViewer::enqueue_vertex(Mesh::VertexHandle _vh)
{
	float                   prio, min_prio(FLT_MAX);
	Mesh::HalfedgeHandle  min_hh;

	// find best out-going halfedge
	for (Mesh::VOHIter vh_it(mesh_, _vh); vh_it; ++vh_it)
	{
		if (is_collapse_legal(vh_it))
		{
			prio = priority(vh_it);	

			if (prio >= -1.0 && prio < min_prio)
			{
				min_prio = prio;
				min_hh   = vh_it.handle();
			}
		}
	}


	// update queue
	QueueVertex qv;
	qv.v=_vh; qv.prio=priority(_vh);
	if (priority(_vh) != -1.0) 
	{
		queue.erase(qv);
		priority(_vh) = -1.0;
	}

	if (min_hh.is_valid()) 
	{
		priority(_vh) = min_prio;
		target(_vh)   = min_hh;
		qv.prio=min_prio;
		queue.insert(qv);
	}
}

void DecimationViewer::enqueue_vertex2(Mesh::VertexHandle _vh)
{
	float                   prio, min_prio(FLT_MAX);
	Mesh::HalfedgeHandle  min_hh;

	// find best out-going halfedge
	for (Mesh::VOHIter vh_it(mesh_, _vh); vh_it; ++vh_it)
	{
		if (is_collapse_legal2(vh_it))
		{
			prio = priority2(vh_it);	
			if (prio >= -1.0 && prio < min_prio)
			{
				min_prio = prio;
				min_hh   = vh_it.handle();
			}
		}
	}


	// update queue
	QueueVertex qv;
	qv.v=_vh; qv.prio=priority(_vh);
	if (priority(_vh) != -1.0) 
	{
		queue.erase(qv);
		priority(_vh) = -1.0;
	}

	if (min_hh.is_valid()) 
	{
		priority(_vh) = min_prio;
		target(_vh)   = min_hh;
		qv.prio=min_prio;
		queue.insert(qv);
	}
}
//void DecimationViewer::enqueue_edge(Mesh::EdgeHandle eh) {
//	float prio, min_prio(FLT_MAX);
//	QEdge qe;
//	qe.eh = eh;
//	qe.prio = priority(eh);
//	qEdge.insert(qe);
//
//	target(eh) = calculateOptimalPoint(eh);
//}


//-----------------------------------------------------------------------------
void  DecimationViewer::blowup(unsigned int _n_vertices) {
	unsigned int nv(mesh_.n_vertices());

	Mesh::HalfedgeHandle hh;
	Mesh::VertexHandle   to, from;
	Mesh::VVIter         vv_it;

	std::vector<Mesh::VertexHandle>            one_ring;
	std::vector<Mesh::VertexHandle>::iterator  or_it, or_end;



	// build priority queue
	Mesh::VertexIter  v_it  = mesh_.vertices_begin(), 
	v_end = mesh_.vertices_end();

	queue.clear();
	for (; v_it!=v_end; ++v_it)
		enqueue_vertex2(v_it.handle());
	
	cout << "Queue contains : " << queue.size() << " elements \n";

	std::cout << "BLOWUP!";
	while (nv > _n_vertices && !queue.empty())
	{
		// Exercise 4.3 ----------------------------------------------
		// INSERT CODE:
		// Decimate using priority queue:
		//   1) take 1st element of queue
		//   2) collapse this halfedge
		//   3) update queue
		// -----------------------------------------------------------
		std::set<QueueVertex, VertexCmp>::iterator iter_minimal = queue.begin();
		QueueVertex qV_minimal = (*iter_minimal);
		queue.erase(iter_minimal); // bye bye

		Mesh::HalfedgeHandle victim_hh = target(qV_minimal.v);

		if(is_collapse_legal2(victim_hh)) {
			--nv;
			// Collect the 1Ring
			std::stack<Mesh::VertexHandle> oneRingers;
			Mesh::VertexHandle to = mesh_.to_vertex_handle(victim_hh);

			for(Mesh::VertexHandle lvh = qV_minimal.v; lvh != to; lvh = to) {
				for(Mesh::VVIter oneRing = mesh_.vv_iter(lvh); oneRing; ++oneRing) {
					if((to == oneRing)||(qV_minimal.v == oneRing)) continue;
					if(!oneRing.handle().is_valid()) continue;
					oneRingers.push(oneRing);
				}
			}

			mesh_.collapse(victim_hh);
			mesh_.set_point(to, calculateOptimalPoint(quadric(to), quadric(qV_minimal.v)));
			quadric(to) += quadric(qV_minimal.v); // The article said so
			
			while(!oneRingers.empty()) {
				enqueue_vertex2(oneRingers.top()); // Update all priorities around
				oneRingers.pop();
			}

			enqueue_vertex2(to);

		}
	}
//	cout << "nv > _n_vertices: " << (nv > _n_vertices) << " !queue.empty() : " << !queue.empty() << "\n";

	// clean up
	queue.clear();

	// now, delete the items marked to be deleted
	mesh_.garbage_collection();

	// re-compute face & vertex normals
	mesh_.update_normals();


	// re-update face indices for faster rendering
	update_face_indices();
}

void DecimationViewer::decimate(unsigned int _n_vertices)
{
	unsigned int nv(mesh_.n_vertices());

	Mesh::HalfedgeHandle hh;
	Mesh::VertexHandle   to, from;
	Mesh::VVIter         vv_it;

	std::vector<Mesh::VertexHandle>            one_ring;
	std::vector<Mesh::VertexHandle>::iterator  or_it, or_end;



	// build priority queue
	Mesh::VertexIter  v_it  = mesh_.vertices_begin(), 
		v_end = mesh_.vertices_end();

	queue.clear();
	for (; v_it!=v_end; ++v_it)
		enqueue_vertex(v_it.handle());


	std::cout << "DECIMATION!";
	while (nv > _n_vertices && !queue.empty())
	{
		// Exercise 4.3 ----------------------------------------------
		// INSERT CODE:
		// Decimate using priority queue:
		//   1) take 1st element of queue
		//   2) collapse this halfedge
		//   3) update queue
		// -----------------------------------------------------------
		std::set<QueueVertex, VertexCmp>::iterator iter_minimal = queue.begin();
		QueueVertex qV_minimal = (*iter_minimal);
		queue.erase(iter_minimal); // bye bye

		Mesh::HalfedgeHandle victim_hh = target(qV_minimal.v);

		if(is_collapse_legal(victim_hh)) {
			--nv;
			// Collect the 1Ring
			std::stack<Mesh::VertexHandle> oneRingers;
			Mesh::VertexHandle to = mesh_.to_vertex_handle(victim_hh);

			for(Mesh::VertexHandle lvh = qV_minimal.v; lvh != to; lvh = to) {
				for(Mesh::VVIter oneRing = mesh_.vv_iter(lvh); oneRing; ++oneRing) {
					if((to == oneRing)||(qV_minimal.v == oneRing)) continue;
					if(!oneRing.handle().is_valid()) continue;
					oneRingers.push(oneRing);
				}
			}

			quadric(to) += quadric(qV_minimal.v); // The article said so
			mesh_.collapse(victim_hh);

			while(!oneRingers.empty()) {
				enqueue_vertex(oneRingers.top()); // Update all priorities around
				oneRingers.pop();
			}

			enqueue_vertex(to);

		}
	}



	// clean up
	queue.clear();

	// now, delete the items marked to be deleted
	mesh_.garbage_collection();

	// re-compute face & vertex normals
	mesh_.update_normals();


	// re-update face indices for faster rendering
	update_face_indices();
}