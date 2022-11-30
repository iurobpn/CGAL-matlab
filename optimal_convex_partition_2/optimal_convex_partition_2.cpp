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
        TypedArray<double> vertex = std::move(inputs[0]);
		const size_t numRows = vertex.getDimensions()[0];
		// double* x = new double[vertex.getNumberOfElements()];
		// memcpy(x,&*vertex.begin(),sizeof(float)*vertex.getNumberOfElements());

		Polygon_List partition_polys = 	optimal_convex_partition_2(numRows,vertex);
		int N = partition_polys.size();
		ArrayFactory f;
		CellArray out = f.createCellArray({1,(long unsigned int)N});
			

		for (int i = 0; i < N; ++i) {
			Polygon_2 polyi = partition_polys[i];
			int M = polyi.size();
			TypedArray<double> outi =  f.createArray<double>({ (long unsigned int)M, 2});

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
        if (numRows < 3 || numColumns != 2)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be at least 3x2 double matrix") }));
        }
    }

	Polygon_List optimal_convex_partition_2(const size_t &rows, const TypedArray<double> &vertex)
	{
		Polygon_2             polygon;
		Polygon_List          partition_polys;
		Traits                partition_traits;
		Validity_traits       validity_traits;

		polygon = make_polygon(vertex,rows);
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

	Polygon_2 make_polygon(const TypedArray<double> &vertex, const int &M)
	{
		Polygon_2 polygon;

		for(int i=0; i<M; i++)
			polygon.push_back(Point_2(vertex[i][0], vertex[i][1]));

		return polygon;
	}
};

