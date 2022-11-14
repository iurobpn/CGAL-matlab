/* MyMEXFunction
 * c = MyMEXFunction(a,b);
 * Adds offset argument a to each element of double array b and
 * returns the modified array c.
*/

#include <cassert>
#include <list>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
		checkArguments(outputs, inputs);
        TypedArray<double> inx = std::move(inputs[0]);
		// const size_t numRows = inx.getDimensions()[0];
		// float* x = new float[inx.getNumberofElements()];
		// memcpy(x,&*inx.begin,sizeof(float)*inx.getNumberofElements());

        outputs[0] = inx;
		// Poligon_List Ps = 	optimal_convex_partition_2(numRows,x);
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

	// Polygon_list optimal_convex_partition_2(const size_t &rows, double *x)
	// {
	// 	Polygon_2             polygon;
	// 	Polygon_list          partition_polys;
	// 	Traits                partition_traits;
	// 	Validity_traits       validity_traits;
    //
	// 	polygon = make_polygon(x,rows);
	// 	// print_polygon(polygon);
	// 	CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
	// 					polygon.vertices_end(),
	// 					std::back_inserter(partition_polys),
	// 					partition_traits);
	// 	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
	// 				 polygon.vertices_end(),
	// 				 partition_polys.begin(),
	// 				 partition_polys.end(),
	// 				 validity_traits));
	// 	
    //
	// 	return partition_polys;
	// }

	// Polygon_2 make_polygon(double *x, const int &M)
	// {
	// 	Polygon_2 polygon;
    //
	// 	for(int i=0; i<M; i++)
	// 		polygon.push_back(Point_2(x[i],x[i+M]));
    //
	// 	return polygon;
	// }



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

