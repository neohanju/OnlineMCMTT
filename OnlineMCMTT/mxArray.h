#pragma once

#define NOMINMAX

#include <iostream>
#include <vector>

using namespace std;

// Type definitions
typedef enum
{
	mxUNKNOWN_CLASS = 0,
	mxCELL_CLASS,
	mxSTRUCT_CLASS,
	mxLOGICAL_CLASS,
	mxCHAR_CLASS,
	mxVOID_CLASS,
	mxDOUBLE_CLASS,
	mxSINGLE_CLASS,
	mxINT8_CLASS,
	mxUINT8_CLASS,
	mxINT16_CLASS,
	mxUINT16_CLASS,
	mxINT32_CLASS,
	mxUINT32_CLASS,
	mxINT64_CLASS,
	mxUINT64_CLASS,
	mxFUNCTION_CLASS,
	mxOPAQUE_CLASS,
	mxOBJECT_CLASS, /* keep the last real item in the list */
#if defined(_LP64) || defined(_WIN64)
	mxINDEX_CLASS = mxUINT64_CLASS,
#else
	mxINDEX_CLASS = mxUINT32_CLASS,
#endif
	/* TEMPORARY AND NASTY HACK UNTIL mxSPARSE_CLASS IS COMPLETELY ELIMINATED */
	mxSPARSE_CLASS = mxVOID_CLASS /* OBSOLETE! DO NOT USE */
}
mxClassID;


typedef enum
{
	mxREAL,
	mxCOMPLEX
}
mxComplexity;



//mxArray structure
class mxArray {
public:
	bool bDataAlloc;
	void* data;

	bool bIsMats;
	mxArray* mats[256];
	
	char fieldName[256][1024];
	int nField;
	
	size_t dim[256];
	size_t nDim;
	mxClassID classID;
	mxComplexity complexity;

};


// mxArray get functions
mxArray* mxGetField(mxArray* pa, size_t i, char* fieldname);
size_t *mxGetDimensions(const mxArray *pa);
size_t mxGetNumberOfDimensions(const mxArray *pa);
void *mxGetData(const mxArray *pa);
double mxGetScalar(const mxArray *pa);
void *mxGetPr(const mxArray *pa);
mxClassID mxGetClassID(const mxArray *pa);
int mxGetString(const mxArray *pa, char *buf, size_t buflen);
size_t mxGetNumberOfElements(const mxArray *pa);
size_t mxGetM(const mxArray *pa);
size_t mxGetN(const mxArray *pa);


//mxArray set functions
int mxSetDimensions(mxArray *pa, const size_t *pdims, size_t ndims);
void mxSetData(mxArray *pa, void  *newdata);


//mxArray memory functions
mxArray *mxCreateNumericArray(size_t ndim, const size_t *dims, mxClassID classid, mxComplexity flag);
mxArray *mxCreateString(char *str);
mxArray *mxCreateDoubleScalar(double value);
mxArray *mxCreateNumericMatrix(size_t m, size_t n, mxClassID classid, mxComplexity flag);
void mxDestroyArray(mxArray *pa);
void *mxMalloc(size_t	n);
void *mxCalloc(size_t	n, size_t size);


//MATFile functions
mxArray* jwmatOpen(char *filename);
mxArray* jwmatOpen(FILE* fp);

void mexErrMsgTxt(const char	*error_msg);

void mxFree(void *pa);
//MATFile functions
mxArray* jwmatOpen(char *filename);

bool mxIsDouble(const mxArray *pa);


