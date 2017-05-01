#include "mxArray.h"

mxArray *mxCreateNumericArray(size_t ndim, const size_t *dims, mxClassID classid, mxComplexity flag) {

	if (ndim > 256)
		return NULL;

	mxArray* output;
	output = (mxArray*)malloc(sizeof(mxArray));
	memset(output, 0, sizeof(mxArray));

	output->nDim = ndim;
	memcpy(output->dim, dims, sizeof(size_t)*ndim);
	output->classID = classid;
	output->complexity = flag;

	size_t nSize = 1;
	for (int ii = 0; ii < ndim; ii++)
		nSize *= dims[ii];

	/*if (ndim > 0) {
		switch (classid) {
			case mxCHAR_CLASS:			output->data = (char*)malloc(nSize*sizeof(char));				break;
			case mxINT8_CLASS:			output->data = (char*)malloc(nSize*sizeof(char));				break;
			case mxUINT8_CLASS:			output->data = (unsigned char*)malloc(nSize*sizeof(char));		break;
			case mxINT16_CLASS:			output->data = (short*)malloc(nSize*sizeof(short));				break;
			case mxUINT16_CLASS:		output->data = (unsigned short*)malloc(nSize*sizeof(short));		break;
			case mxINT32_CLASS:			output->data = (int*)malloc(nSize*sizeof(int));					break;
			case mxUINT32_CLASS:		output->data = (unsigned int*)malloc(nSize*sizeof(int));		break;
			case mxINT64_CLASS:			output->data = (long long*)malloc(nSize*sizeof(long long));			break;
			case mxUINT64_CLASS:		output->data = (unsigned long long*)malloc(nSize*sizeof(long long));	break;
			case mxSINGLE_CLASS:		output->data = (float*)malloc(nSize*sizeof(float));				break;
			case mxDOUBLE_CLASS:		output->data = (double*)malloc(nSize*sizeof(double));				break;
		}
		output->bDataAlloc = 1;
	}*/

	if (ndim > 0) {
		switch (classid) {
			case mxCHAR_CLASS:			output->data = (char*)calloc(nSize, sizeof(char));				break;
			case mxINT8_CLASS:			output->data = (char*)calloc(nSize, sizeof(char));				break;
			case mxUINT8_CLASS:			output->data = (unsigned char*)calloc(nSize, sizeof(unsigned char));		break;
			case mxINT16_CLASS:			output->data = (short*)calloc(nSize, sizeof(short));				break;
			case mxUINT16_CLASS:		output->data = (unsigned short*)calloc(nSize, sizeof(short));		break;
			case mxINT32_CLASS:			output->data = (int*)calloc(nSize, sizeof(int));					break;
			case mxUINT32_CLASS:		output->data = (unsigned int*)calloc(nSize, sizeof(int));		break;
			case mxINT64_CLASS:			output->data = (long long*)calloc(nSize, sizeof(long long));			break;
			case mxUINT64_CLASS:		output->data = (unsigned long long*)calloc(nSize, sizeof(long long));	break;
			case mxSINGLE_CLASS:		output->data = (float*)calloc(nSize, sizeof(float));				break;
			case mxDOUBLE_CLASS:		output->data = (double*)calloc(nSize, sizeof(double));				break;
		}
		output->bDataAlloc = 1;
	}
	else
		output->bDataAlloc = 0;

	output->bIsMats = 0;
	
	return output;

}





mxArray *mxCreateString(char *str){

	size_t length = 0;
	char* curr = str;
	while (*(curr++) != NULL)
		length++;
	length++;
	length = 1024;
	mxArray* output = mxCreateNumericArray(1, (const size_t *)&length, mxINT8_CLASS, mxREAL);
	memset(output->data, 0, sizeof(char)*1024);
	memcpy(output->data, str, length);
	 
	return output;

}





mxArray *mxCreateDoubleScalar(double value){

	size_t length[1]; length[0] = 1;

	mxArray* output = mxCreateNumericArray(1, length, mxDOUBLE_CLASS, mxREAL);

	memcpy(output->data, &value, sizeof(double));

	return output;

}





mxArray *mxCreateNumericMatrix(size_t m, size_t n, mxClassID classid, mxComplexity flag){

	size_t length[2]; length[0] = m; length[1] = n;

	mxArray* output = mxCreateNumericArray(2, length, classid, mxREAL);

	return output;

}





void mxDestroyArray(mxArray *pa){
	if (pa == NULL) { return;}
	if (pa->bIsMats) {
		for (int ii = 0; ii < pa->nField; ii++) {
				mxDestroyArray(pa->mats[ii]);
		}
	}
	if (pa->bDataAlloc)
		free(pa->data);
	free(pa);
	pa = NULL;
}





void *mxMalloc(size_t	n){

	return (void*)malloc(n);

}





void *mxCalloc(size_t	n, size_t	size){

	return (void*)calloc(n,size);

}








mxArray* mxGetField(mxArray* pa, size_t i, char* fieldname) {

	if (pa->bIsMats) {

		for (int ii = 0; ii < pa->nField; ii++) {
			if (!strcmp(fieldname, pa->fieldName[ii])) {

				mxArray* output;
				output = pa->mats[ii];
				
				return output;

			}
		}

	}

	return NULL;

}



size_t *mxGetDimensions(const mxArray *pa) {

	return (size_t*)(pa->dim);

}




size_t mxGetNumberOfDimensions(const mxArray *pa) {

	return (size_t)(pa->nDim);

}




void *mxGetData(const mxArray *pa) {

	return (void*)(pa->data);

}




double mxGetScalar(const mxArray *pa) {

	double* temp = (double*)pa->data;

	return temp[0];

}




void *mxGetPr(const mxArray *pa) {

	return (void*)pa->data;

}




mxClassID mxGetClassID(const mxArray *pa) {

	return pa->classID;

}




int mxGetString(const mxArray *pa, char *buf, size_t buflen) {

	memcpy(buf, (char*)pa->data, buflen);

	return 0;

}




size_t mxGetNumberOfElements(const mxArray *pa) {

	size_t length = 1;
	for (int ii = 0; ii < pa->nDim; ii++)
		length *= pa->dim[ii];

	return (size_t)length;

}




size_t mxGetM(const mxArray *pa) {

	return pa->dim[0];

}




size_t mxGetN(const mxArray *pa) {

	return pa->dim[1];

}




int mxSetDimensions(mxArray *pa, const size_t *pdims, size_t ndims) {
	if(pa == NULL)
	{
		pa = mxCreateNumericArray(ndims, pdims, mxSINGLE_CLASS, mxREAL);
	}
	else {
	memcpy(pa->dim, pdims, sizeof(size_t)*ndims);
	pa->nDim = ndims;
	}
	return 0;

}





void mxSetData(mxArray *pa, void  *newdata) {
	if(pa->bDataAlloc)
		if(pa->data != newdata)
			free(pa->data);
	pa->data = newdata;
	pa->bDataAlloc = 1;
}




void mexErrMsgTxt(const char	*error_msg) {
	printf(error_msg);
}


// 수정할 때 아래 void jwmatOpen도 동일하게 수정할 것.
mxArray* jwmatOpen(char *filename) {
	
	FILE* fp;
	fopen_s(&fp, filename, "r+");

	unsigned int nVars, sz;

	fscanf_s(fp, "%d\n", &nVars);
	fscanf_s(fp, "%d\n", &sz);

	mxArray* output;
	output = (mxArray*)malloc(sizeof(mxArray));
	memset(output, 0, sizeof(mxArray));

	output->bIsMats = 1;
	
	char currInput[1024];
	int currVar = 0;
	int currPos = 0;
	int nDims;
	int dim[2];
	int typ, cls;
	char cc;
	int ii = 0;
	int kk = 0;
	while (1) {
		ii = 0;
		while (1) {
			fscanf_s(fp, "%c", &cc,1);
			//printf("%c %d ", cc, ii);
			if (cc == '\n') break;
			currInput[ii++] = cc;
		}		
		currInput[ii] = '\0';
		unsigned int length = ii+1;
		
		if (!strcmp(currInput, "end -1\0")) break;

		memcpy(output->fieldName[kk], currInput, length);

		fscanf_s(fp, "%d %d %d\n", &typ, &cls, &nDims);
		sz = 1;
		if (typ == 0)		//structure
		{
			mxArray* tt;
			tt = jwmatOpen(fp);
			output->mats[kk++] = tt;
		}
		else if (typ == 1)	//string
		{
			for (int ii = 0; ii < nDims; ii++) {
				fscanf_s(fp, "%d\n", &(dim[ii]));
			}
			mxArray* tt;
			tt = mxCreateString("");
			
			output->mats[kk++] = tt;

			for (int ii = 0; ii < dim[0] * dim[1]; ii++) {
				fscanf_s(fp, "%c\n", (char*)(tt->data) + ii, 1);
				//printf("%lf\n", *((double*)(output->data) + (currPos)));
				currPos++;
			}
			memset((char*)tt->data + dim[0]*dim[1], 0, sizeof(char)*(1024-dim[0]*dim[1]));
		}
		else if (typ == 2)	//numeric matrix
		{
			for (int ii = 0; ii < nDims; ii++) {
				fscanf_s(fp, "%d\n", &(dim[ii]));
			}
			mxArray* tt;
			if (cls == 6)
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxDOUBLE_CLASS, mxREAL);
			else if (cls == 7)
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxSINGLE_CLASS, mxREAL);
			else
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxUINT32_CLASS, mxREAL);

			output->mats[kk++] = tt;

			double reading;
			for (int ii = 0; ii < dim[0] * dim[1]; ii++) {
				fscanf_s(fp, "%lf\n", &reading);
				//printf("%lf\n", *((double*)(output->data) + (currPos)));
				if (cls == 6)
					memcpy(((double*)(tt->data) + ii), &reading, sizeof(double));
				else if (cls == 7) {
					float aa = (float)reading;
					memcpy(((float*)(tt->data) + ii), &aa, sizeof(float));
				}
				else {
					unsigned long aa = (unsigned int)reading;
					memcpy(((unsigned int*)(tt->data) + ii), &aa, sizeof(unsigned int));
				}

				currPos++;
			}
		}
		
	}

	fclose(fp);
	
	output->nField = kk;
	return output;

}

mxArray* jwmatOpen(FILE* fp) {
	
	unsigned int nVars, sz;

	fscanf_s(fp, "%d\n", &nVars);
	fscanf_s(fp, "%d\n", &sz);

	mxArray* output;
	output = (mxArray*)malloc(sizeof(mxArray));
	memset(output, 0, sizeof(mxArray));

	output->bIsMats = 1;

	char currInput[1024];
	int currVar = 0;
	int currPos = 0;
	int nDims;
	int dim[2];
	int typ, cls;
	char cc;
	int ii = 0;
	int kk = 0;
	while (1) {
		ii = 0;
		while (1) {
			fscanf_s(fp, "%c", &cc, 1);
			if (cc == '\n') break;
			currInput[ii++] = cc;
		}
		currInput[ii] = '\0';
		unsigned int length = ii + 1;

		if (!strcmp(currInput, "end -1\0")) break;

		memcpy(output->fieldName[kk], currInput, length);

		fscanf_s(fp, "%d %d %d\n", &typ, &cls, &nDims);
		sz = 1;
		if (typ == 0)		//structure
		{
			mxArray* tt;
			tt = jwmatOpen(fp);
			output->mats[kk++] = tt;
		}
		else if (typ == 1)	//string
		{
			for (int ii = 0; ii < nDims; ii++) {
				fscanf_s(fp, "%d\n", &(dim[ii]));
			}
			mxArray* tt;
			tt = mxCreateString("");
			output->mats[kk++] = tt;

			for (int ii = 0; ii < dim[0]*dim[1]; ii++) {
				fscanf_s(fp, "%c\n", (char*)(tt->data) + ii, 1);				
				//printf("%lf\n", *((double*)(output->data) + (currPos)));
				currPos++;
			}
			memset((char*)tt->data + dim[0]*dim[1], 0, sizeof(char)*(1024-dim[0]*dim[1]));
		}
		else if (typ == 2)	//numeric matrix
		{
			for (int ii = 0; ii < nDims; ii++) {
				fscanf_s(fp, "%d\n", &(dim[ii]));
			}
			mxArray* tt;
			if (cls == 6)
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxDOUBLE_CLASS, mxREAL);
			else if (cls == 7)
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxSINGLE_CLASS, mxREAL);
			else
				tt = mxCreateNumericMatrix(dim[0], dim[1], mxUINT32_CLASS, mxREAL);

			output->mats[kk++] = tt;

			double reading;
			for (int ii = 0; ii < dim[0] * dim[1]; ii++) {
				fscanf_s(fp, "%lf\n", &reading);
				//printf("%lf\n", *((double*)(output->data) + (currPos)));
				if (cls == 6)
					memcpy(((double*)(tt->data) + ii), &reading, sizeof(double));
				else if (cls == 7) {
					float aa = (float)reading;
					memcpy(((float*)(tt->data) + ii), &aa, sizeof(float));
				}
				else {
					unsigned long aa = (unsigned int)reading;
					memcpy(((unsigned int*)(tt->data) + ii), &aa, sizeof(unsigned int));
				}

				currPos++;
			}
		}



	}
	
	output->nField = kk;
	return output;

}

void mxFree(void *pa) {
	free(pa);
}

bool mxIsDouble(const mxArray *pa) {
	if(pa->classID == mxDOUBLE_CLASS) {
		return 1;
	}
	return 0;
}

