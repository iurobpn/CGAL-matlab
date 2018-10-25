#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h> 
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/HalfedgeDS_vertex_max_base_with_id.h>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef typename HalfedgeDS::Vertex   Vertex;
typedef typename Vertex::Point Point;
typedef Polyhedron::Vertex_iterator        Vertex_iterator;
typedef Polyhedron::Facet_iterator        Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;
typedef std::vector<int> Facet;
typedef std::vector<Point> VertexList;
typedef std::vector<Facet> FacetList;

struct MeshPolyhedron {
	VertexList vertices;
	FacetList facets;

	MeshPolyhedron(){}
	MeshPolyhedron(VertexList vs, FacetList fs) {
		this->vertices = vs;
		this->facets = fs;
	}
};
/* typedef std::list<Polyhedron> PolyhedronList; */

template <class HDS>
class BuildPolyhedron : public CGAL::Modifier_base<HDS> {
private:
	VertexList vertices;
	FacetList facets;
public:
	BuildPolyhedron(VertexList vertices, FacetList facets) 
	{

		this->vertices.insert(std::end(this->vertices),
				std::begin(vertices), std::end(vertices));
		this->facets.insert(std::end(this->facets),
				std::begin(facets), std::end(facets));
	}

    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
	int num_vertices = vertices.size();
	int num_facets = facets.size();
	// std::cout << "num vertices: " << num_vertices << ", num facets: " << num_facets << std::endl;
        B.begin_surface( num_vertices, num_facets);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
	for (auto v : vertices)
		B.add_vertex( v );

	int k = 0;
	for (auto facet : facets) {
		k++;
		bool isfacet = B.test_facet(facet.begin(),facet.end());
		if (!isfacet) {
			std::cout << "obs, facet not well defined " << k << "\n.";
		}
		B.begin_facet();
		for (auto point : facet) {
			B.add_vertex_to_facet( point );
			// std::cout << point << std::endl;
		}
		B.end_facet();
	}
        B.end_surface();
    }
};

Polyhedron createPolyhedron(double vs[][3], int num_vertices, int fs[][4], int num_facets);
std::list<Polyhedron> decompose(Nef_polyhedron N);
std::list<Polyhedron> decompose(Polyhedron P);
/* MeshPolyhedron poly2mesh(Polyhedron P); */
void poly2mesh(Polyhedron P, MeshPolyhedron &MP);
int find(VertexList vs, Point p);

int main() {
	// obstacle info
	int num_vertices = 8;
	int num_facets = 6;
	// double v[8][3] = {80,0,20,80,20,20,80,20,0,80,0,0,120,0,20,120,20,20,120,20,0,120,0,0};
	double v[8][3] = {60,0,20,60,20,20,60,20,0,60,0,0,90,0,20,90,20,20,90,20,0,90,0,0};
	int f[6][4] = {1,2,3,4,8,7,6,5,3,2,6,7,4,3,7,8,1,4,8,5,2,1,5,6};
	Polyhedron Pobs = createPolyhedron(v,num_vertices,f,num_facets);
	/* decompose(P); */
	// simplified velodyne range of vision
	num_vertices = 8;
	num_facets = 6;
	double vs[8][3] = {10,-90,10,10,10,110,10,110,10,10,10,-90,110,-90,10,110,10,110,110,110,10,110,10,-90};
	int fs[6][4] = {1,2,3,4,8,7,6,5,1,4,8,5,4,3,7,8,6,7,3,2,5,6,2,1};
	Polyhedron P = createPolyhedron(vs,num_vertices,fs,num_facets);
	/* decompose(P); */
	Nef_polyhedron Nobs(Pobs);
	Nef_polyhedron N(P);
	Nef_polyhedron Nf = N-Nobs;
	std::list<Polyhedron> convex_parts = decompose(Nf);
	Polyhedron P1 = convex_parts.front();

	MeshPolyhedron MP;
	poly2mesh(P1,MP);
	
	std::cout << "polyhedro 1 all facets\n";
	for (auto f : MP.facets) {
		std::cout << "nova face\n";
		for (auto i : f) {
			/* std::cout << i << " ";	 */
			std::cout << MP.vertices[i] << std::endl;	
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void poly2mesh(Polyhedron P, MeshPolyhedron &MP)
{
	std::size_t i = 0; 
	for(Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
		  MP.vertices.push_back(v->point());

	for ( Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f) {
		i++;
		HF_circulator hinit = f->facet_begin();
		HF_circulator h = hinit;
		Facet facet;
		do {
			Point p = h->vertex()->point();
			int index = find(MP.vertices,p);
			facet.push_back(index);
			h++;
		} while(h!=hinit);
		MP.facets.push_back(facet);
	}

	/* MeshPolyhedron Pout(vertices,facets); */
}

int find(VertexList vs, Point p) 
{
	int index = 0;
	for(auto v : vs) {
		if (p[0]==v[0] && p[1]==v[1] && p[2]==v[2])
			return index;
		index++;
	}
	return -1;
}

Polyhedron createPolyhedron(double vs[][3], int num_vertices, int fs[][4], int num_facets)
{
	VertexList vertices;
	for (int i = 0; i < num_vertices; ++i)
		vertices.push_back(Point(vs[i][0], vs[i][1], vs[i][2]));
	FacetList facets;
	for (int i = 0; i < num_facets; ++i) {
		facets.push_back(Facet());
		for (int j = 0; j < 4; ++j)
			facets[i].push_back(fs[i][j]-1);
	}

	Polyhedron P;
	BuildPolyhedron<HalfedgeDS> poly(vertices,facets);
	P.delegate(poly);

	return P;
}

std::list<Polyhedron> decompose(Polyhedron P)
{
	Nef_polyhedron N(P);
	return decompose(N);
}

std::list<Polyhedron> decompose(Nef_polyhedron N)
{
	/* CGAL_assertion( P.is_triangle( P.halfedges_begin())); */

	CGAL::convex_decomposition_3(N);
	std::list<Polyhedron> convex_parts;

	// the first volume is the outer volume, which is 
	// ignored in the decomposition
	Volume_const_iterator ci = ++N.volumes_begin();
	for( ; ci != N.volumes_end(); ++ci) {
		if(ci->mark()) {
			Polyhedron P;
			N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
			convex_parts.push_back(P);
		}
	}
	std::cout << "decomposition into " << convex_parts.size() << " convex parts " << std::endl;

	return convex_parts;
}

