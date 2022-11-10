/** MyMEXFunction
 * c = MyMEXFunction(a,b);
 * Adds offset argument a to each element of double array b and
 * returns the modified array c.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/partition_2.h>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#include "mex.hpp"
#include "mexAdapter.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef Traits::Polygon_2                                   Polygon_2;
typedef Traits::Point_2                                     Point_2;
typedef CGAL::Is_convex_2<Traits>                           Is_convex_2;
typedef std::vector<Polygon_2>                                Polygon_List;
typedef CGAL::Partition_is_valid_traits_2<Traits, Is_convex_2>
                                                            Validity_traits;
typedef Polygon_2::Vertex_iterator 			    VertexIterator;

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
		checkArguments(outputs, inputs);
        TypedArray<double> inx = std::move(inputs[0]);
		const size_t numRows = inx.getDimensions()[0];
		double* x = new double[inx.getNumberOfElements()];
		memcpy(x,&*inx.begin(),sizeof(float)*inx.getNumberOfElements());

		Polygon_List partition_polys = 	optimal_convex_partition_2(numRows,x);
		int N = partition_polys.size();
		ArrayFactory f;
		CellArray out = f.createCellArray({1,N});
			

		for (int i = 0; i < N; ++i) {
			Polygon_2 polyi = partition_polys[i];
			int M = polyi.size();
			TypedArray<double> outi =  f.createArray<double>({ M, 2});

			int j = 0;
			for (auto &elem : outi) {
				int k = j%M;
				Point_2 p = polyi.vertex(k);
				if (j < M)
					elem = (double) p.x();
				else
					elem = (double) p.y();
				j++;
			}
			out[i] = outi;
		}
		outputs[0] = out;
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check offset argument: First input must be scalar double
        if (inputs.size() != 1)
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
		// Check offset argument: First input must be scalar double
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be double") }));
        }
		const size_t numRows = inputs[0].getDimensions()[0];
        const size_t numColumns = inputs[0].getDimensions()[1];
        if (numRows > 2 && numColumns == 2)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be at lesat 3x2 double matrix") }));
        }

    }

	Polygon_List optimal_convex_partition_2(const size_t &rows, double *x)
	{
		Polygon_2             polygon;
		Polygon_List          partition_polys;
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

	Polygon_2 make_polygon(double *x, const int &M)
	{
		Polygon_2 polygon;

		for(int i=0; i<M; i++)
			polygon.push_back(Point_2(x[i],x[i+M]));

		return polygon;
	}

	// void print_polygon(const Polygon_2& p)
	// {
	// 	std::cout << "P: ";
	// 	for(VertexIterator vi = p.vertices_begin(); vi != p.vertices_end();
	// 			++vi) { 
	// 		printf("(%f,%f)\n", vi->x(), vi->y());
	// 	}
	// 	std::cout << "\n";
	// }
};

