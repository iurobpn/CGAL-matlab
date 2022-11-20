/** MyMEXFunction
 * c = MyMEXFunction(a,b);
 * Adds offset argument a to each element of double array b and
 * returns the modified array c.
 */

#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h> 
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <list>
#include <vector>
#include <thread>
#include "mex.hpp"
#include "mexAdapter.hpp"


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
};

template <class HDS>
class BuildPolyhedron : public CGAL::Modifier_base<HDS> {
private:
public:
	VertexList vertices;
	FacetList facets;
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
			std::cout << "obs, facet not well defined " << k << "\n";

			std::cout << "facet size: " << facet.size() << "\n";
			for (auto i : facet) {
				std::cout << i << " ";
			}
			std::cout << std::endl;
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

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
	void operator()(ArgumentList outputs, ArgumentList inputs) {
		checkArguments(outputs, inputs);
		// TypedArray<double> inx = std::move(inputs[0]);
		// const size_t numRows = inx.getDimensions()[0];
		// double* x = new double[inx.getNumberOfElements()];
		// memcpy(x,&*inx.begin(),sizeof(float)*inx.getNumberOfElements());

		// Polygon_List partition_polys = 	optimal_convex_partition_2(numRows,x);
		// int N = 2;partition_polys.size();

		double *verticesRaw;
		int Nin = inputs.size();

		std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
		ArrayFactory f;
		/* Check for proper number of arguments. */
		size_t num_poly_in = Nin/2.0; //number of input polygons

		Nef_polyhedron Nf(Nef_polyhedron::EMPTY); 
		for (int i = 0; i < Nin; i = i + 2) {
			size_t numRows = inputs[i].getDimensions()[0];
			size_t numCols = inputs[i].getDimensions()[1];
			if (inputs[i].getType() != ArrayType::DOUBLE ||
					inputs[i].getType() == ArrayType::COMPLEX_DOUBLE ||
					!(numRows > 4 && numCols == 3)) {
				matlabPtr->feval(u"error",
						0,
						std::vector<Array>({ f.createScalar("It must be at least a 4x3 double matrix") }));
			}
			TypedArray<double> Vi = std::move(inputs[i]);
			verticesRaw = new double[Vi.getNumberOfElements()];
			memcpy(verticesRaw,&*Vi.begin(),sizeof(float)*Vi.getNumberOfElements());

			CellArray facets = inputs[i+1];

			Polyhedron P = createPolyhedron(verticesRaw,numRows,facets);
			Nef_polyhedron N(P);
			if (i==0) {
				Nf = Nf + N;
			} else {
				Nf = Nf - N;
			}
		}

		// for (int i = 0; i < N; ++i) {
		// 	Polygon_2 polyi = partition_polys[i];
		// 	int M = polyi.size();
		// 	TypedArray<double> outi =  f.createArray<double>({ (long unsigned int)M, 2});
        //
		// 	int j = 0;
		// 	for (auto &elem : outi) {
		// 		int k = j%M;
		// 		Point_2 p = polyi.vertex(k);
		// 		if (j < M)
		// 			elem = (double) p.x();
		// 		else
		// 			elem = (double) p.y();
		// 		j++;
		// 	}
		// 	out[i] = outi;
		// }

		// std::cout << "num poly in: " << num_poly_in << std::endl;

		std::list<Polyhedron> convex_parts = decompose(Nf);
		size_t num_out = convex_parts.size();
		CellArray out = f.createCellArray({1,num_out});

		int cell_index = 0;
		int i = 0;
		for (auto poly : convex_parts) {
			MeshPolyhedron MP;
			poly2mesh(poly,MP);
			size_t num_vertices = poly.size_of_vertices();
			int k = 0;
			TypedArray<double> Vi = f.createArray<double>({num_vertices,3});
			for (auto v : MP.vertices) {
				Vi[k] = CGAL::to_double(v[0]);
				Vi[k + num_vertices] = CGAL::to_double(v[1]);
				Vi[k + num_vertices*2] = CGAL::to_double(v[2]);
				k++;
			}

			size_t num_facets = MP.facets.size();
			CellArray facetsCell = f.createCellArray({1,num_facets});
			size_t num;
			int facet_index = 0;
			for (auto facet : MP.facets) {
				num = facet.size();
				TypedArray<double> F = f.createArray<double>({num,1});
				k=0;
				for (auto i : facet) {
					F[k] = (double)i+1; //add 1 before returnin to matlab(indices start at 1)
					k++;
				}
				facetsCell[facet_index] = F;
				facet_index++;
			}

			num = 2;
			CellArray polyCell = f.createCellArray({1,num});
			polyCell[0] = Vi;
			polyCell[1] = facetsCell;
			out[cell_index++] = polyCell;
		}
		outputs[0] = out;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs) {
		// Get pointer to engine
		std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

		// Get array factory
		ArrayFactory factory;

		if (inputs.size() < 2)
		{
			matlabPtr->feval(u"error",
					0,
					std::vector<Array>({ factory.createScalar("Error: it must have 1 one input") }));
		}

		// Check number of outputs
		if (outputs.size() > 1) 
		{
			matlabPtr->feval(u"error",
					0,
					std::vector<Array>({ factory.createScalar("Too many output arguments") }));
		}
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

	/**
	 * verticesRaw:  
	 *
	 */
	Polyhedron createPolyhedron(double *verticesRaw, int num_vertices, const CellArray facets)
	{
		size_t num_facets = facets.getNumberOfElements();
		// Assign pointers to each input and output.
		double *facetRaw;

		FacetList facets_lst;
		for (int i = 0; i < num_facets; ++i) {
			std::vector<int> facet;
			TypedArray<double> Fi = facets[i];

			for (auto elem : Fi)
				facet.push_back(elem-1); //compensate for matlab stating with 1 on its indexes
			facets_lst.push_back(facet);
		}


		VertexList vertices;
		for (int i = 0; i < num_vertices; ++i)
			vertices.push_back(Point(verticesRaw[i], verticesRaw[num_vertices+i], verticesRaw[2*num_vertices+i]));

		Polyhedron P;
		BuildPolyhedron<HalfedgeDS> poly(vertices,facets_lst);
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
		// std::cout << "decomposition into " << convex_parts.size() << " convex parts " << std::endl;

		return convex_parts;
	}
};

