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

void make_polygon(Polygon_2& polygon);

void __mexFunction__( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
}

int main_in(int argc, char *argv[])
{
   Polygon_2             polygon;
   Polygon_list          partition_polys;
   Traits                partition_traits;
   Validity_traits       validity_traits;
/*
   CGAL::random_polygon_2(50, std::back_inserter(polygon),
                          Point_generator(100));
*/

	make_polygon(polygon);
	std::chrono::high_resolution_clock::time_point to = std::chrono::
		high_resolution_clock::now();
	CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
				    polygon.vertices_end(),
				    std::back_inserter(partition_polys),
				    partition_traits);
	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
			     polygon.vertices_end(),
			     partition_polys.begin(),
			     partition_polys.end(),
			     validity_traits));
   	

	// std::string outfilename("result.log");
	// if (argc > 1)
	// 	outfilename = argv[1];

	// unsigned int nthreads = std::thread::hardware_concurrency();
	// std::ofstream outfile (outfilename);

	// if (!outfile.is_open()) {
	// 	printf("Unable to open file");
	// 	exit(1);
	// }

	// outfile << "Case #" << i+1 << ": " << n << std::endl;
	// outfile.close();
	


	std::chrono::high_resolution_clock::time_point tf = std::chrono::
		high_resolution_clock::now();
	std::chrono::duration<double,std::ratio<1l,1000000l>> delta_t = 
		std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1l,1000000l>>>(tf - to);
	// print_polygon(polygon);
	std::list<Polygon_2>::iterator it;
	for(it = partition_polys.begin(); it != partition_polys.end(); it++) {
	    // std::cout << *it << ' ';
            //
	// for(int i=0; i< partition_polys.size(); i++) {
		print_polygon(*it);
	}

	std::cout << "Elapsed: " << delta_t.count() << " us.\n";


   return 0;
}

void make_polygon(Polygon_2& polygon)
{
   polygon.push_back(Point_2(2.2699525, -4.4550328));
   polygon.push_back(Point_2(2.7527637,	-2.0));
   polygon.push_back(Point_2(4.4550328,	-2.2699525));
   polygon.push_back(Point_2(5.0, 0.0));
   polygon.push_back(Point_2(4.7552829,	1.5450850));
   polygon.push_back(Point_2(2.5000000,	1.2738137));
   polygon.push_back(Point_2(2.0381019,	4.0));
   polygon.push_back(Point_2(-1.7484555e-07, 4.0));
}



void print_polygon(const Polygon_2& p)
{
	std::cout << "P: ";
	for(VertexIterator vi = p.vertices_begin(); vi != p.vertices_end();
			++vi) { 
		std::cout << "(" << vi->x() << "," << vi->y() << ") ";
	}
	std::cout << "\n";
}

