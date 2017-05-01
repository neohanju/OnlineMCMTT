/*******************************************************************************
* Piotr's Image&Video Toolbox      Version 3.00
* Copyright 2012 Piotr Dollar.  [pdollar-at-caltech.edu]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the Simplified BSD License [see external/bsd.txt]
*******************************************************************************/
#ifndef _WRAPPERS_HPP_
#define _WRAPPERS_HPP_

#include "mxArray.h"

inline void wrError(const char *errormsg) { throw errormsg; };

// platform independent aligned memory allocation (see also alFree)
void* alMalloc(size_t size, int alignment);

// platform independent alignned memory de-allocation (see also alMalloc)
void alFree(void* aligned);

#endif
