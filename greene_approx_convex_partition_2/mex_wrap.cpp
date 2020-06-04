#include "mex.h"
#include "mex_main.h"
#include <time.h>
#include <thread>

/* static void at_exit(); */

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    __mexFunction__(nlhs, plhs, nrhs, prhs);

}

