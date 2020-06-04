#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <cassert>
#include <list>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include "mex.h"
#include "mex_main.h"
#include "matrix.h"
  // mxArray *mxCreateCellArray(mwSize ndim, const mwSize *dims);

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef CGAL::Is_convex_2<Traits>                           Is_convex_2;
typedef Traits::Polygon_2                                   Polygon_2;
typedef Traits::Point_2                                     Point_2;
// typedef Polygon_2::Vertex_const_iterator                    Vertex_iterator;
typedef Polygon_2::Vertex_iterator 			    VertexIterator;
typedef std::list<Polygon_2>                                Polygon_list;
typedef CGAL::Partition_is_valid_traits_2<Traits, Is_convex_2>
                                                            Validity_traits;
typedef CGAL::Creator_uniform_2<int, Point_2>               Creator;
typedef CGAL::Random_points_in_square_2<Point_2, Creator>   Point_generator;

void print_polygon(const Polygon_2& p);
Polygon_list optimal_convex_partition_2(const size_t &rows, double *x);
Polygon_list greene_approx_convex_partition(const size_t &rows, double *x);

Polygon_2 make_polygon(double *x,const int &M);

void __mexFunction__( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	double *x,*y;
	size_t mrows,ncols;

	/* Check for proper number of arguments. */
	if(nrhs!=1) {
		mexErrMsgIdAndTxt( "MATLAB:timesfour:invalidNumInputs",
	    				"One input required.");
	} else if(nlhs>1) {
		mexErrMsgIdAndTxt( "MATLAB:timesfour:maxlhs",
	    				"Too many output arguments.");
	}

	/* The input must be a noncomplex 3x2 double matrix.*/
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	// mexPrintf("rows:%d, cols:%d\n",mrows,ncols);
	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(mrows>2 && ncols==2) ) {
		mexErrMsgIdAndTxt( "MATLAB:optimal_convex_partition_2:inputNotBigEnough",
				    "Input must be at least 3x2 double matrix.");
	}


	// Assign pointers to each input and output.
	x = mxGetPr(prhs[0]);

	/* Call the timesfour subroutine. */
	Polygon_list pl = greene_approx_convex_partition(mrows,x);
	// mxArray *mxCreateCellArray(mwSize ndim, const mwSize *dims);
	/* Create matrix for the return argument. */
	mwSize numpoly[1];
	numpoly[0] = pl.size();
	// mexPrintf("numpoly: %d\n",numpoly[0]);

	plhs[0] = mxCreateCellArray(1, numpoly);
	y = mxGetPr(plhs[0]);

	Polygon_list::iterator it;
	int index = 0;
	for(it = pl.begin(); it != pl.end(); it++) {
		Polygon_2 p = *it;
		/*  set the output pointer to the output matrix */
		int n = p.size();
		// mexPrintf("polygon size: %d\n",n);
		mxArray *c = mxCreateDoubleMatrix( n, 2, mxREAL);
		double *cptr = mxGetPr(c);

		int i = 0;
		for(VertexIterator vi = p.vertices_begin(); vi != p.vertices_end();
			++vi) {
			cptr[i] = vi->x();
			cptr[i+n] = vi->y();
			i++;
		}
		mxSetCell(plhs[0],index,c);
		index++;
	}

	return;
}

Polygon_list optimal_convex_partition_2(const size_t &rows, double *x)
{
	Polygon_2             polygon;
	Polygon_list          partition_polys;
	Traits                partition_traits;
	Validity_traits       validity_traits;

	polygon = make_polygon(x,rows);
	// print_polygon(polygon);
	CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
				    polygon.vertices_end(),
				    std::back_inserter(partition_polys),
				    partition_traits);
	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
			     polygon.vertices_end(),
			     partition_polys.begin(),
			     partition_polys.end(),
			     validity_traits));


	return partition_polys;
}

Polygon_list greene_approx_convex_partition(const size_t &rows, double *x)
{
	Polygon_2    polygon;
	Polygon_list partition_polys;
	Traits       partition_traits;

	polygon = make_polygon(x,rows);
	CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
											polygon.vertices_end(),
											std::back_inserter(partition_polys),
											partition_traits);

	assert(CGAL::convex_partition_is_valid_2(polygon.vertices_begin(),
											polygon.vertices_end(),
											partition_polys.begin(),
											partition_polys.end(),
											partition_traits));

	return partition_polys;
}

Polygon_2 make_polygon(double *x, const int &M)
{
	Polygon_2 polygon;

	for(int i=0; i<M; i++)
		polygon.push_back(Point_2(x[i],x[i+M]));

	return polygon;
}



void print_polygon(const Polygon_2& p)
{
	std::cout << "P: ";
	for(VertexIterator vi = p.vertices_begin(); vi != p.vertices_end();
			++vi) {
		mexPrintf("(%f,%f)\n", vi->x(), vi->y());
	}
	std::cout << "\n";
}

