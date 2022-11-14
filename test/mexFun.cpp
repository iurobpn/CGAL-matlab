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
		std::ostringstream stream;
		stream << "dims: " << numRows << "," << numColumns << "\n";
		matlabPtr->feval(u"fprintf", 0,
				std::vector<Array>({ factory.createScalar(stream.str()) }));
        if (numRows < 3 || numColumns != 2)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be at least 3x2 double matrix") }));
        }

    }

	void display(std::shared_ptr<matlab::engine::MATLABEngine> &matlabPtr, std::ostringstream& stream) {

        ArrayFactory factory;
		// Pass stream content to MATLAB fprintf function
		matlabPtr->feval(u"fprintf", 0,
				std::vector<Array>({ factory.createScalar(stream.str()) }));
		// Clear stream buffer
		stream.str("");
	}
};

