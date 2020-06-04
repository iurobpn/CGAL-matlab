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
#include "mex.h"
#include "mex_main.h"
#include "matrix.h"

class mystream : public std::streambuf {
protected:
	virtual std::streamsize xsputn(const char *s, std::streamsize n) { mexPrintf("%.*s", n, s); return n; }
	virtual int overflow(int c=EOF) { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
};

class scoped_redirect_cout {
public:
	scoped_redirect_cout() { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
	~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
private:
	mystream mout;
	std::streambuf *old_buf;
};
scoped_redirect_cout mycout_redirect;

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
	/* MeshPolyhedron(VertexList vs, FacetList fs) { */
	/* 	this->vertices = vs; */
	/* 	this->facets = fs; */
	/* } */
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

Polyhedron createPolyhedron(double *verticesRaw, int num_vertices, const mxArray *prhs);
std::list<Polyhedron> decompose(Nef_polyhedron N);
std::list<Polyhedron> decompose(Polyhedron P);
void poly2mesh(Polyhedron P, MeshPolyhedron &MP);
int find(VertexList vs, Point p);


void __mexFunction__( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	double *verticesRaw, *y;
	size_t mrowsV, ncolsV, mrowsi, ncolsi;

	/* Check for proper number of arguments. */
	mwSize num_poly_in = nrhs/2; //number of input polygons
	// std::cout << "num poly in: " << num_poly_in << std::endl;
	if(nrhs<2) {
		mexErrMsgIdAndTxt( "MATLAB:convex_decomposition_3:invalidNumInputs",
	    				"At least two inputs required.");
	} else if(nlhs > 1) {
		mexErrMsgIdAndTxt( "MATLAB:convex_decomposition_3:maxlhs",
	    				"Too many output arguments.");
	}
	Nef_polyhedron Nf(Nef_polyhedron::EMPTY); 
	for (int i = 0; i < nrhs; i=i+2) {
		/* The input must be a noncomplex 4x3 double matrix.*/
		mrowsV = mxGetM(prhs[i]);
		ncolsV = mxGetN(prhs[i]);

		// mexPrintf("rows:%d, cols:%d\n",mrows,ncols);
		if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || !(mrowsV>4 && ncolsV==3) ) {
			mexErrMsgIdAndTxt( "MATLAB:convex_decomposition_3:inputNotBigEnough",
					    "Input must be at least 4x3 double matrix.");
		}

		verticesRaw = mxGetPr(prhs[i]);

		Polyhedron P = createPolyhedron(verticesRaw,mrowsV,prhs[i+1]);
		Nef_polyhedron N(P);
		if (i==0) {
			Nf = Nf + N;
		} else {
			Nf = Nf - N;
		}

	}
	std::list<Polyhedron> convex_parts = decompose(Nf);
	mwSize num_polyhedrons = convex_parts.size();

	/* plhs[0] = mxCreateDoubleMatrix( 0, 0, mxREAL ); */
	/* return; */

	plhs[0] = mxCreateCellArray(1, &num_polyhedrons);
	/* y = mxGetPr(plhs[0]); */
	int cell_index = 0;
	int i = 0;
	// for (auto poly = convex_parts.begin(); poly != convex_parts.end(); ++poly) {
	// 	MeshPolyhedron MP;
	// 	poly2mesh(*poly,MP);
	// 	std::cout << "Size of polyhedron " << i << ": " << MP.vertices.size() << std::endl;	
	// 	i++;
	// }
	for (auto poly : convex_parts) {
		MeshPolyhedron MP;
		poly2mesh(poly,MP);
		int num_vertices = poly.size_of_vertices();
		mxArray *V = mxCreateDoubleMatrix( num_vertices, 3, mxREAL);
		double *Vptr = mxGetPr(V);
		/* std::cout << "size of vertices: " << num_vertices << std::endl; */
		int k = 0;
		for (auto v : MP.vertices) {
			/* std::cout << v << " \n"; */
			Vptr[k] = CGAL::to_double(v[0]);
			Vptr[k + num_vertices] = CGAL::to_double(v[1]);
			Vptr[k + num_vertices*2] = CGAL::to_double(v[2]);
			k++;
		}

		mwSize num_facets = MP.facets.size();
		/* std::cout << "size of vertices: " << num_vertices << std::endl; */
		mxArray *facetsCell = mxCreateCellArray(1, &num_facets);
		int facet_index = 0;
		for (auto facet : MP.facets) {
			int num = facet.size();
			mxArray *F = mxCreateDoubleMatrix( num, 1, mxREAL);
			double *Fptr = mxGetPr(F);
			k=0;
			for (auto i : facet) {
				Fptr[k] = (double)i+1; //add 1 before returnin to matlab(indices start at 1)
				k++;
			}
			mxSetCell(facetsCell,facet_index,F);
			facet_index++;
		}

		mwSize num = 2;
		mxArray *polyCell = mxCreateCellArray(1, &num);
		mxSetCell(polyCell,0,V);
		mxSetCell(polyCell,1,facetsCell);
		mxSetCell(plhs[0],cell_index,polyCell);
		cell_index++;
	}

	return;
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
Polyhedron createPolyhedron(double *verticesRaw, int num_vertices, const mxArray *mxFacetsRaw)
{
	int num_facets = mxGetNumberOfElements(mxFacetsRaw);
	// Assign pointers to each input and output.
	const mxArray *cell_element_ptr;
	double *facetRaw;
	
	FacetList facets;
	for (int i = 0; i < num_facets; ++i) {
		std::vector<int> facet;
		cell_element_ptr = mxGetCell(mxFacetsRaw,i);
		mwSize mrowsi = mxGetM(cell_element_ptr);
		mwSize ncolsi = mxGetN(cell_element_ptr);
		double *facetRaw = mxGetPr(cell_element_ptr);
		
		for (int k = 0; k < mrowsi; ++k) {
			int index = facetRaw[k];
			facet.push_back(index-1); //compensate for matlab stating with 1 on its indexes
		}
		facets.push_back(facet);
	}


	VertexList vertices;
	for (int i = 0; i < num_vertices; ++i)
		vertices.push_back(Point(verticesRaw[i], verticesRaw[num_vertices+i], verticesRaw[2*num_vertices+i]));
	
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
	// std::cout << "decomposition into " << convex_parts.size() << " convex parts " << std::endl;

	return convex_parts;
}

