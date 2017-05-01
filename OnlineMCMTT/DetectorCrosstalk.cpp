#pragma once

#include <math.h>
#include <time.h>
#include <string.h>

#include "xmmintrin.h."
#include "sse.hpp"

#include "haanju_math.hpp"
#include "DetectorCrosstalk.hpp"

#define PSN_2D_MAX_HEIGHT      (2300)
#define PSN_2D_MIN_HEIGHT      (1400)
#define PSN_2D_MAX_HEAD_HEIGHT (800)
#define PSN_2D_MIN_HEAD_HEIGHT (330)

typedef unsigned int uint32;

// 반올림함수
#define PI 3.14159265f

inline void getChild( float *chns1, uint32 *cids, uint32 *fids, float *thrs, uint32 offset, uint32 &k0, uint32 &k )
{
  float ftr = chns1[cids[fids[k]]];
  k = (ftr<thrs[k]) ? 1 : 2;
  k0=k+=k0*2; k+=offset;
}

inline bool SortDetectionScoreDescend(hj::CDetection &val1, hj::CDetection &val2)
{
	return val1.score > val2.score? true : false;
}

inline void VECTOR_MAX(std::vector<double> &outputVector, const std::vector<double> *vector1, const std::vector<double> *vector2)
{
	assert(vector1->size() == vector2->size());
	outputVector.reserve(vector1->size());
	for (size_t elementIdx = 0; elementIdx < vector1->size(); elementIdx++)
	{
		double maxElement =  (*vector1)[elementIdx] > (*vector2)[elementIdx] ? (*vector1)[elementIdx] : (*vector2)[elementIdx];
		outputVector.push_back(maxElement);
	}
}
inline int MOD (int a, int b)
{
   int ret = a % b;
   if (ret < 0)
	   ret+=b;
   return ret;
}

inline double Rounding( double x, int digit )
{
    return (floor( (x) * pow( float(10), digit) + 0.5f) / pow( float(10), digit ));
}

CDetectorCrosstalk::CDetectorCrosstalk(void)
{
	stCrosstalkParams_.mxArray_opts  = NULL;
	stCrosstalkParams_.mxArray_clf   = NULL;
	stCrosstalkParams_.mxArray_info1 = NULL;
	stCrosstalkParams_.mxArray_info2 = NULL;
	stCrosstalkParams_.mxArray_info3 = NULL;
}

CDetectorCrosstalk::~CDetectorCrosstalk(void)
{
}

bool 
CDetectorCrosstalk::Initialize(hj::stParamDetect2D _params)
{
	////////////////////////////// detector 읽기(opts, clf, info 따로) ///////////////////////////////

	char file_opts[300];
	char file_clf[300];
	char file_info1[300];
	char file_info2[300];
	char file_info3[300];

	stParams_ = _params;	

	sprintf_s(file_opts,  "%s/opts.jwmat",  stParams_.strDetectorDir.c_str());
	sprintf_s(file_clf,   "%s/clf.jwmat",   stParams_.strDetectorDir.c_str());
	sprintf_s(file_info1, "%s/info1.jwmat", stParams_.strDetectorDir.c_str());
	sprintf_s(file_info2, "%s/info2.jwmat", stParams_.strDetectorDir.c_str());
	sprintf_s(file_info3, "%s/info3.jwmat", stParams_.strDetectorDir.c_str());

	stCrosstalkParams_.mxArray_clf = jwmatOpen(file_clf);
	stCrosstalkParams_.mxArray_opts = jwmatOpen(file_opts);
	stCrosstalkParams_.mxArray_info1 = jwmatOpen(file_info1);
	stCrosstalkParams_.mxArray_info2 = jwmatOpen(file_info2);
	stCrosstalkParams_.mxArray_info3 = jwmatOpen(file_info3);

	stCrosstalkParams_.pr_pPyramid     = mxGetField(stCrosstalkParams_.mxArray_opts, 0, "pPyramid");
	stCrosstalkParams_.pr_pChns        = mxGetField(stCrosstalkParams_.pr_pPyramid, 0, "pChns");
	stCrosstalkParams_.pr_shrink       = mxGetField(stCrosstalkParams_.pr_pChns, 0, "shrink");
	stCrosstalkParams_.pr_pColor       = mxGetField(stCrosstalkParams_.pr_pChns, 0, "pColor");
	stCrosstalkParams_.pr_colorSpace   = mxGetField(stCrosstalkParams_.pr_pColor,    0, "colorSpace");
	stCrosstalkParams_.pr_nPerOct      = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "nPerOct");
	stCrosstalkParams_.pr_nOctUp       = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "nOctUp");
	stCrosstalkParams_.pr_minDs        = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "minDs");
	stCrosstalkParams_.pr_nApprox      = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "nApprox");
	stCrosstalkParams_.pr_smooth       = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "smooth");
	stCrosstalkParams_.pr_enabled      = mxGetField(stCrosstalkParams_.pr_pColor,    0, "enabled");
	stCrosstalkParams_.pr_smooth_color = mxGetField(stCrosstalkParams_.pr_pColor,    0, "smooth");
	stCrosstalkParams_.pr_lambdas      = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "lambdas");
	stCrosstalkParams_.pr_pad          = mxGetField(stCrosstalkParams_.pr_pPyramid,  0, "pad");
	stCrosstalkParams_.pr_padWith0     = mxGetField(stCrosstalkParams_.mxArray_info1,0,"padWith"); // unsigned char
	stCrosstalkParams_.pr_padWith1     = mxGetField(stCrosstalkParams_.mxArray_info2,1,"padWith");
	stCrosstalkParams_.pr_padWith2     = mxGetField(stCrosstalkParams_.mxArray_info3,2,"padWith");

	ptCalibrationInfo_ = _params.pCalibrationInfo;

	return true;
}

bool
CDetectorCrosstalk::Finalize()
{	
	mxDestroyArray(stCrosstalkParams_.mxArray_clf);
	mxDestroyArray(stCrosstalkParams_.mxArray_info1);
	mxDestroyArray(stCrosstalkParams_.mxArray_info2);
	mxDestroyArray(stCrosstalkParams_.mxArray_info3);
	mxDestroyArray(stCrosstalkParams_.mxArray_opts);
	return true;
}

vector<hj::CDetection>
CDetectorCrosstalk::Detect(cv::Mat _targetImage, int _frameIndex)
{
	//clear previous detection
	m_DetectionResultVectorOneCam.clear();
	m_RawDetectionResultVector.clear();

	if (1.0 != stParams_.dImageRescale)
	{
		cv::resize(_targetImage, _targetImage, 
			cv::Size((int)((double)_targetImage.cols * stParams_.dImageRescale), 
			         (int)((double)_targetImage.rows * stParams_.dImageRescale)));
	}

	///////////////////////////////// Mat -> mxArray 타입변환 /////////////////////////////////////
	mxArray* mxArray_matCurrentFrame_k = this->CvMat2MxArray(&_targetImage);
	// create and set output array
	const size_t *dims = mxGetDimensions(mxArray_matCurrentFrame_k);
	m_ImgSize.height = (int) dims[0];
	m_ImgSize.width = (int) dims[1];
		
	///////////////////////////////// RGB -> LUV //////////////////////////////////////////////////
	mxArray* mxArray_testimg_LUV = this->CrosstalkRGBConvert(mxArray_matCurrentFrame_k);
		
	char * mxArray_colorSpace_pr = (char *)mxGetData(stCrosstalkParams_.pr_colorSpace);
	mxArray_colorSpace_pr = "orig";
		
		
	this->CrosstalkGetScales();
	m_pyr_chn.resize((int)m_nScales);
	for (int scale_idx = 0; scale_idx < m_nScales; scale_idx++) {
		m_pyr_chn[scale_idx].resize(3);
		m_pyr_chn[scale_idx][0] = NULL;
		m_pyr_chn[scale_idx][1] = NULL;
		m_pyr_chn[scale_idx][2] = NULL;
	}
	m_pyr_chn_concat.resize((int)m_nScales, NULL);
		
	this->CrosstalkComputeImagePyramid(m_pyr_chn, mxArray_testimg_LUV); // deleaker memory leak			
	this->CrosstalkComputeImagePyramidApp(m_pyr_chn); 		
	this->CrosstalkChannelSmoothETC(m_pyr_chn);
	this->CrosstalkComputeImPad(m_pyr_chn);
	this->CrosstalkConcat(m_pyr_chn_concat, m_pyr_chn);	

	m_DetectionResultVectorOneCam = this->CrosstalkComputeDetect(m_pyr_chn_concat);

	m_RawDetectionResultVector = m_DetectionResultVectorOneCam;
	m_DetectionResultVectorOneCam = bbNms_maxg(m_DetectionResultVectorOneCam, stParams_.dNMSOverlapRatio_1);

	int count = 0; // detection threshold validation
	int OriginalResultSize = (int)m_DetectionResultVectorOneCam.size() ;
	for(int idx = 0 ; idx < OriginalResultSize; ++idx) {
		if(m_DetectionResultVectorOneCam[count].score < stParams_.nDetectionThreshold) {
			m_DetectionResultVectorOneCam.erase(m_DetectionResultVectorOneCam.begin() + count);
			continue;
		}
		count++;
	}

	////// image resize restoration 
	//OriginalResultSize = (int)m_DetectionResultVectorOneCam.size() ;
	//for(int idx = 0 ; idx < OriginalResultSize; ++idx) {
	//	m_DetectionResultVectorOneCam[idx].box.x = m_DetectionResultVectorOneCam[idx].box.x/resize_ratio;
	//	m_DetectionResultVectorOneCam[idx].box.y = m_DetectionResultVectorOneCam[idx].box.y/resize_ratio;
	//	m_DetectionResultVectorOneCam[idx].box.w = m_DetectionResultVectorOneCam[idx].box.w/resize_ratio;
	//	m_DetectionResultVectorOneCam[idx].box.h = m_DetectionResultVectorOneCam[idx].box.h/resize_ratio;
	//}

	// detection size validation
	if (NULL != ptCalibrationInfo_)
	{
		if (stParams_.nDetectionThreshold == 0 || stParams_.nDetectionThreshold == 2) {
			hj::Point2D headBottomCenter(0.0, 0.0), bottomCenter(0.0, 0.0), topCenter(0.0, 0.0), boxCenter(0.0, 0.0);
			hj::Point3D WorldBoxcenter, WorldBottomCenter;
			double curHeight = 0.0;
			count = 0;
			int OriginalResultSize2 = (int)m_DetectionResultVectorOneCam.size();
			for (size_t detectionIdx = 0; detectionIdx < OriginalResultSize2; detectionIdx++)
			{
				headBottomCenter = m_DetectionResultVectorOneCam[count].box.bottomCenter();
				bottomCenter = m_DetectionResultVectorOneCam[count].box.bottomCenter();
				topCenter = bottomCenter;
				topCenter.y -= m_DetectionResultVectorOneCam[count].box.h;

				boxCenter = topCenter;
				boxCenter.y += m_DetectionResultVectorOneCam[count].box.h / 2;
				WorldBoxcenter.z = 1550;
				ptCalibrationInfo_->cCamModel.imageToWorld(boxCenter.x, boxCenter.y, WorldBoxcenter.z, WorldBoxcenter.x, WorldBoxcenter.y);
				WorldBottomCenter = WorldBoxcenter;
				WorldBottomCenter.z = 0;
				ptCalibrationInfo_->cCamModel.worldToImage(WorldBottomCenter.x, WorldBottomCenter.y, WorldBottomCenter.z, bottomCenter.x, bottomCenter.y);

				m_DetectionResultVectorOneCam[count].box.x = m_DetectionResultVectorOneCam[count].box.x - m_DetectionResultVectorOneCam[count].box.w / 2;
				m_DetectionResultVectorOneCam[count].box.w = 2 * m_DetectionResultVectorOneCam[count].box.w;
				m_DetectionResultVectorOneCam[count].box.h = abs(bottomCenter.y - topCenter.y);

				curHeight = hj::EstimateBoxHeight(m_DetectionResultVectorOneCam[count].box, *ptCalibrationInfo_);

				if (PSN_2D_MAX_HEAD_HEIGHT < curHeight || PSN_2D_MIN_HEAD_HEIGHT > curHeight || m_DetectionResultVectorOneCam[count].box.h > 0) { // || m_DetectionResultVectorOneCam[count].box.h >120
					m_DetectionResultVectorOneCam.erase(m_DetectionResultVectorOneCam.begin() + count);
					continue;
				}
				count++;
			}
		}

		if (stParams_.nDetectionThreshold == 1) {
			hj::Point2D headBottomCenter(0.0, 0.0), bottomCenter(0.0, 0.0), topCenter(0.0, 0.0), boxCenter(0.0, 0.0);
			hj::Point3D WorldTopcenter, WorldBottomCenter;
			double curHeight = 0.0;
			count = 0;
			int OriginalResultSize2 = (int)m_DetectionResultVectorOneCam.size();
			for (size_t detectionIdx = 0; detectionIdx < OriginalResultSize2; detectionIdx++)
			{
				curHeight = hj::EstimateBoxHeight(m_DetectionResultVectorOneCam[count].box, *ptCalibrationInfo_);
				if (PSN_2D_MAX_HEIGHT < curHeight || PSN_2D_MIN_HEIGHT > curHeight) { // || m_DetectionResultVectorOneCam[count].box.h >120
					m_DetectionResultVectorOneCam.erase(m_DetectionResultVectorOneCam.begin() + count);
					continue;
				}
				count++;
			}
		}
	}
	
	m_RawDetectionResultVector = m_DetectionResultVectorOneCam;
	m_DetectionResultVectorOneCam = bbNms_maxg(m_DetectionResultVectorOneCam, stParams_.dNMSOverlapRatio_2);

	// detection size shrink
	double shrinkRatio = 0.9;
	double shrinkWidth, shrinkHeight;
	for(int idx = 0 ; idx < (int)m_DetectionResultVectorOneCam.size(); ++idx) {
		shrinkWidth = ((1.0-shrinkRatio)/2.0)*m_DetectionResultVectorOneCam[idx].box.w;
		shrinkHeight = ((1.0-shrinkRatio)/2.0)*m_DetectionResultVectorOneCam[idx].box.h;
		m_DetectionResultVectorOneCam[idx].box.x += shrinkWidth;
		m_DetectionResultVectorOneCam[idx].box.y += shrinkHeight;
		m_DetectionResultVectorOneCam[idx].box.w *= shrinkRatio;
		m_DetectionResultVectorOneCam[idx].box.h *= shrinkRatio;

		// recover size
		m_DetectionResultVectorOneCam[idx].box *= stParams_.dImageRescaleRecover;
	}

	mxDestroyArray(mxArray_matCurrentFrame_k);
	mxDestroyArray(mxArray_testimg_LUV);
	this->CrosstalkClearMxarray();
	int NumDetection = (int)m_DetectionResultVectorOneCam.size();

	return this->m_DetectionResultVectorOneCam;
}

vector<hj::CDetection>
CDetectorCrosstalk::CrosstalkGetRawDetection(void) {
	return this->m_RawDetectionResultVector;
}

void
CDetectorCrosstalk::CrosstalkClearMxarray(void) {
	for (int scale_idx = 0; scale_idx < m_pyr_chn.size(); scale_idx++) {
		for(int j = 0; j < 3; j ++) {			
			mxDestroyArray(m_pyr_chn[scale_idx][j]);
		}
	}
	m_pyr_chn.clear();
	for (int scale_idx = 0; scale_idx < (int) m_pyr_chn_concat.size(); scale_idx++) {
		mxDestroyArray(m_pyr_chn_concat[scale_idx]);
	}
	m_pyr_chn_concat.clear();
	m_scales.clear();
	m_scaleshw.clear();
	m_isR.clear();
	m_isA.clear();
	m_isJ.clear();
	m_isN.clear();
}


void 
CDetectorCrosstalk::CrosstalkConcat(vector<mxArray *> &pyr_chn_concat,vector<vector<mxArray*>> &pyr_chn) {
	for(int i = 0 ; i < pyr_chn.size();i++) {
		float *pyr_chn0 = (float *)mxGetData(pyr_chn[i][0]);
		size_t *dims = (size_t *)mxGetDimensions(pyr_chn[i][0]);
		dims[2] = 10;
		pyr_chn_concat[i] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
		float *chn_concat_pr = (float *)mxGetData(pyr_chn_concat[i]);
		float *pyr_chn1 = (float *)mxGetData(pyr_chn[i][1]);
		float *pyr_chn2 = (float *)mxGetData(pyr_chn[i][2]);
		memcpy(chn_concat_pr,pyr_chn0,dims[0]*dims[1]*3*sizeof(float));
		memcpy(chn_concat_pr+dims[0]*dims[1]*3,pyr_chn1,dims[0]*dims[1]*sizeof(float));
		memcpy(chn_concat_pr+dims[0]*dims[1]*4,pyr_chn2,dims[0]*dims[1]*6*sizeof(float));
	}
	//memory clearing
	for(int i = 0 ; i < pyr_chn.size();i++) {
		for(int j = 0 ; j < pyr_chn[i].size(); j++) {
			mxDestroyArray(pyr_chn[i][j]);
		}
		pyr_chn[i].clear();
	}
	pyr_chn.clear();
}

vector<hj::CDetection> 
CDetectorCrosstalk::CrosstalkComputeDetect(vector<mxArray *> &pyr_chn_concat) {
	vector<mxArray*> vec_bbs((int)m_nScales); // for output
	vector<hj::CDetection> st_bbs;
	int nDs = 1; // Fix nDs to 1
	mxArray *mxArray_modelDsPad = mxGetField(stCrosstalkParams_.mxArray_opts,0,"modelDsPad");
	double *modelDsPad = (double*) mxGetData(mxArray_modelDsPad);
	mxArray *mxArray_modelDs = mxGetField(stCrosstalkParams_.mxArray_opts,0,"modelDs");
	double *modelDs = (double*) mxGetData(mxArray_modelDs);
	double *pad = (double*) mxGetData(stCrosstalkParams_.pr_pad);
	
	//int modified_scale = (int)m_nScales;
	int modified_scale = (int)m_nScales-4;
	m_scaleshw.erase(m_scaleshw.end()-1);
	m_scaleshw.erase(m_scaleshw.end()-1);
	m_scaleshw.erase(m_scaleshw.end()-1);
	m_scaleshw.erase(m_scaleshw.end()-1);

	for(int ii=0;ii < modified_scale;++ii) {
		for(int jj=1;jj<=nDs;jj++) {
			double curHeightRatio = m_scaleshw[ii].height;
			mxArray *bb = CrosstalkAcfDetect(pyr_chn_concat[ii], ii);
			double shift[2];
			shift[0] = (modelDsPad[0]-modelDs[0])/2-pad[0];
			shift[1] = (modelDsPad[1]-modelDs[1])/2-pad[1];

			const size_t *bb_size = mxGetDimensions(bb);
			double *temp_bb = (double*) mxGetData(bb);
			// 1st column
			for(int row=0; row<bb_size[0]; row++) {
				temp_bb[row] = (temp_bb[row]+shift[1])/m_scaleshw[ii].width;
			}
			// 2nd column
			for(int row=0; row<bb_size[0]; row++) {
				temp_bb[bb_size[0]+row] = (temp_bb[bb_size[0]+row]+shift[0])/m_scaleshw[ii].height;
			}
			// 3rd column
			for(int row=0; row<bb_size[0]; row++) {
				temp_bb[bb_size[0]*2+row] = modelDs[1]/m_scales[ii];
			}
			// 4th column
			for(int row=0; row<bb_size[0]; row++) {
				temp_bb[bb_size[0]*3+row] = modelDs[0]/m_scales[ii];
			}
						
			//vec_bbs.push_back(bb);
			vec_bbs[ii] = bb;

			for(int row=0; row<bb_size[0]; row++) {
				hj::CDetection st_temp;
				st_temp.box.x = temp_bb[row];
				st_temp.box.y = temp_bb[row + bb_size[0]];
				st_temp.box.w = temp_bb[row + bb_size[0]*2];
				st_temp.box.h = temp_bb[row + bb_size[0]*3];
				st_temp.score = (float)temp_bb[row + bb_size[0]*4];

				st_bbs.push_back(st_temp);
			}
		}
	}

	for(int i=0; i< vec_bbs.size(); i++) {
		mxDestroyArray(vec_bbs[i]);
	}

	return st_bbs;
}
// J = rgbConvertMex(I,flag,single); see rgbConvert.m for usage details
// pr[0] = I, pr[1] = flag, pr[2] =single 
mxArray*
CDetectorCrosstalk::CrosstalkRGBConvert(mxArray* mxArray_InputImage, int flag, bool single) {
	int nDims, n, d; void *I; void *J;
	mxClassID idIn, idOut;

	const size_t *dims = mxGetDimensions(mxArray_InputImage);
	n=((int) dims[0])*((int) dims[1]);
	nDims = (int)mxGetNumberOfDimensions(mxArray_InputImage);
	d = 1; 
	for( int i=2; i<nDims; i++ ) 
		d*=(int) dims[i];
	// extract input arguments
	I = mxGetPr(mxArray_InputImage);
	idIn = mxGetClassID(mxArray_InputImage);

	idOut = single ? mxSINGLE_CLASS : mxDOUBLE_CLASS;
	mxArray* mxArray_OutputImage = mxCreateNumericMatrix(0,0,idOut,mxREAL); // create pointer
	J = mxGetData(mxArray_OutputImage);

	// call rgbConvert() based on type of input and output array
	if(!((d==1 && flag==0) || flag==1 || (d/3)*3==d))
	mexErrMsgTxt("I must have third dimension d==1 or (d/3)*3==d.");
	if( idIn == mxSINGLE_CLASS && !single )
		J = (void*) rgbConvert( (float*) I, n, d, flag, 1.0f );
	else if( idIn == mxSINGLE_CLASS && single )
		J = (void*) rgbConvert( (float*) I, n, d, flag, 1.0f );
	else if( idIn == mxDOUBLE_CLASS && !single )
		J = (void*) rgbConvert( (double*) I, n, d, flag, 1.0 );
	else if( idIn == mxDOUBLE_CLASS && single )
		J = (void*) rgbConvert( (double*) I, n, d, flag, 1.0 );
	else if( idIn == mxUINT8_CLASS && !single )
		J = (void*) rgbConvert( (unsigned char*) I, n, d, flag, 1.0/255);
	else if( idIn == mxUINT8_CLASS && single )
		J = (void*) rgbConvert( (unsigned char*) I, n, d, flag, 1.0f/255);
	else
	mexErrMsgTxt("Unsupported image type.");

	mxSetData(mxArray_OutputImage,J);  // allocate memory
	mxSetDimensions(mxArray_OutputImage,(const size_t*) dims,3); // set dimension
	return mxArray_OutputImage;
}

// Output : m_scales, m_scaleshw
void
CDetectorCrosstalk::CrosstalkGetScales(void)
{
	//주의 : 정확한타입으로 캐스팅 해주지 않으면 엉뚱한 값이 넘어간다.
	double shrink = mxGetScalar(stCrosstalkParams_.pr_shrink);
	double nPerOct = mxGetScalar(stCrosstalkParams_.pr_nPerOct);
	double nOctUp = mxGetScalar(stCrosstalkParams_.pr_nOctUp);
	double nApprox = mxGetScalar(stCrosstalkParams_.pr_nApprox);
	double *minDs = (double *)mxGetData(stCrosstalkParams_.pr_minDs);
	m_nScales = floor(nPerOct*(nOctUp+log(MIN(m_ImgSize.height/minDs[0],m_ImgSize.width/minDs[1]))/log(2.0))+1);
	double d0, d1;
	for(double i = 0;i < m_nScales ;i++) {
		m_scales.push_back(pow(2.0,-(i/nPerOct)+nOctUp));
	}
	if(m_ImgSize.height < m_ImgSize.width) {
		d0 = m_ImgSize.height;
		d1 = m_ImgSize.width;
	}
	else {
		d0 = m_ImgSize.height;
		d1 = m_ImgSize.width;
	}
	for(int i = 0;i < (int)m_nScales ;i++) {
		double s = m_scales[i];
		double s0 = (Rounding(d0*s/shrink,0)*shrink-0.25*shrink)/d0;
		double s1 = (Rounding(d0*s/shrink,0)*shrink+0.25*shrink)/d0;
		vector<double> es0, es1, ss;
		double minValue = DBL_MAX;
		int minIdx = 0;
		double constJ = 0.0;
		for(int j = 0; j <= 100;j++) {
			if(100 == j){ constJ -= 1.0E-32; }
			double ss_temp = constJ*(s1-s0)+s0;
			double es0_temp = d0*ss_temp;
			double es1_temp = d1*ss_temp;

			es0_temp = abs(es0_temp-Rounding(es0_temp/shrink,0)*shrink);
			es1_temp = abs(es1_temp-Rounding(es1_temp/shrink,0)*shrink);

			ss.push_back(ss_temp);
			es0.push_back(es0_temp);
			es1.push_back(es1_temp);
			constJ += 0.01;

			double curMaxValue = es0_temp > es1_temp ? es0_temp : es1_temp;
			if(curMaxValue > minValue){	continue; }
			minValue = curMaxValue;
			minIdx = j;			
		}
		m_scales[i] = ss[minIdx];
	}
	vector<bool> kp((int)m_nScales, true);
	for(int i = 0;i < (int)m_nScales-1 ;i++) {
		kp[i] = (m_scales[i] == m_scales[i+1]) ? false : true;
	}
	int count = 0;
	for(int i = 0;i < (int)m_nScales ;i++) {
		if(kp[count] == false) {
			m_scales.erase(m_scales.begin() + count);
		}
		count++;
	}
	count = 0;
	for(int i = 0; i < (int)m_scales.size(); i++)
	{
		Size2d m_scaleshw_temp;
		m_scaleshw_temp.height = Rounding(m_ImgSize.height*m_scales[count]/shrink,0)*shrink/m_ImgSize.height;
		m_scaleshw_temp.width = Rounding(m_ImgSize.width*m_scales[count]/shrink,0)*shrink/m_ImgSize.width;
		m_scaleshw.push_back(m_scaleshw_temp);
		count ++;
	}
	m_nScales = (double) m_scales.size();
	
	for(int i = 0 ; i < (int)m_nScales ; i = i + (int)(nApprox) + 1) {
		m_isR.push_back(i);
	}
	for(int i = 0 ; i < (int)m_nScales ; i = i++) {
		m_isA.push_back(i);
		m_isN.push_back(i);
	}
	for(int i = 0 ; i < (int)m_isR.size() ; i++) {
		int aa = m_isR[(int)m_isR.size()-1-i];
		m_isA.erase(m_isA.begin()+aa);
	}
	for(int i = 0 ; i < (int)m_isR.size() + 1; i = i++) {
		if(i == 0) m_isJ.push_back(-1);
		else if(i == (int)m_isR.size())  m_isJ.push_back((int)m_nScales-1);
		else m_isJ.push_back((m_isR[i-1]+m_isR[i])/2);
	}
	for(int i = 0 ; i < (int)m_isR.size(); i = i++) {
		for(int j = m_isJ[i]+1 ; j <= m_isJ[i+1];j++) {
				m_isN[j] = m_isR[i];
		}
	}
	//this->m_scales.push_back();
}

void
CDetectorCrosstalk::CrosstalkComputeImagePyramid(vector<vector<mxArray*>> &pyr_chns, mxArray*& InputImage) {
	double shrink = mxGetScalar(stCrosstalkParams_.pr_shrink);
	double nApprox = mxGetScalar(stCrosstalkParams_.pr_nApprox);
	double nPerOct = mxGetScalar(stCrosstalkParams_.pr_nPerOct);
	for(int i = 0 ; i < (int)m_isR.size(); i++) {
	//for(int i = 0 ; i < 1; i++) {
		bool flagI1Delete = true;
		int j = m_isR[i];
		double s = m_scales[j];
		cv::Size2d sz1(Rounding(m_ImgSize.width*s/shrink,0)*shrink, Rounding(m_ImgSize.height*s/shrink,0)*shrink);
		mxArray* I1 = NULL;
		if(sz1.height == m_ImgSize.height && sz1.width == m_ImgSize.width)
		{
			I1 = InputImage;
			flagI1Delete = false;
		}
		else
			I1 = CrosstalkIMResample(InputImage,(int) sz1.height,(int) sz1.width,1.0);				
		if(s == 0.5 && (nApprox > 0 || nPerOct == 1))
		{
			mxDestroyArray(InputImage);
			InputImage = I1;
			flagI1Delete = false;
		}
		mxArray* temp[3];
		chnsCompute(temp, I1);
		pyr_chns[j][0] = temp[0];
		pyr_chns[j][1] = temp[1];
		pyr_chns[j][2] = temp[2];
		if(flagI1Delete)
			mxDestroyArray(I1);
	}
}
void 
CDetectorCrosstalk::CrosstalkComputeImagePyramidApp(vector<vector<mxArray*>> &pyr_chns) {
	double shrink = mxGetScalar(stCrosstalkParams_.pr_shrink);
	double *lambdas = (double *)mxGetData(stCrosstalkParams_.pr_lambdas);
	for(int i = 0 ; i < (int)m_isA.size(); i++) {
		int j = m_isA[i];
		int iR = m_isN[j];
		cv::Size2d sz1(Rounding(m_ImgSize.width*m_scales[j]/shrink,0), Rounding(m_ImgSize.height*m_scales[j]/shrink,0));
		int nTypes = 3;
		for(int k = 0 ; k < nTypes; k++) {
			float ratio = (float) pow(m_scales[j]/m_scales[iR],-lambdas[k]);
			float * prr = (float *)mxGetData(pyr_chns[iR][k]); 
			pyr_chns[j][k] = CrosstalkIMResample(pyr_chns[iR][k],(int)sz1.height,(int)sz1.width,ratio); // allocate memory
		}
	}
}
void 
CDetectorCrosstalk::CrosstalkChannelSmoothETC(vector<vector<mxArray*>> &pyr_chns) {
	for(int i = 0 ; i < pyr_chns.size();i++) {
		for(int j = 0 ; j < 3;j++) {
			mxArray* smoothedResult = this->ConvTri(pyr_chns[i][j],1);
			mxDestroyArray(pyr_chns[i][j]);
			pyr_chns[i][j] = smoothedResult;
		}
	}
}
void 
CDetectorCrosstalk::CrosstalkComputeImPad(vector<vector<mxArray*>> &pyr_chns) {
	double *pad_pr = (double *)mxGetData(stCrosstalkParams_.pr_pad);
	double shrink = mxGetScalar(stCrosstalkParams_.pr_shrink);
	pad_pr[0] = pad_pr[0]/shrink;
	pad_pr[1] = pad_pr[1]/shrink;
	vector<mxArray *> vec_padwith;
	vec_padwith.push_back(stCrosstalkParams_.pr_padWith0);
	vec_padwith.push_back(stCrosstalkParams_.pr_padWith1);
	vec_padwith.push_back(stCrosstalkParams_.pr_padWith2);
	for(int i = 0 ; i < pyr_chns.size();i++) {
		for(int j = 0 ; j < 3;j++) {
			mxArray* ImPadResult = this->CrosstalkImPadMex(pyr_chns[i][j], stCrosstalkParams_.pr_pad, vec_padwith[j]);
			mxDestroyArray(pyr_chns[i][j]);
			pyr_chns[i][j] = ImPadResult;
		}
	}
	pad_pr[0] = pad_pr[0]*shrink; // 다시 되돌리기
	pad_pr[1] = pad_pr[1]*shrink;
}
mxArray* 
CDetectorCrosstalk::CrosstalkIMResample(mxArray* mxArray_InputImage, int sz1_1, int sz1_2, float ratio){
	// sz1_1 height, sz1_2 width
	int ms[3], n, m, nCh;
	size_t nDims;
	const size_t *ns;
	void *A, *B; mxClassID id; double nrm;

	// Error checking on arguments
	//if( nrhs!=4) mexErrMsgTxt("Four inputs expected.");
	//if( nlhs>1 ) mexErrMsgTxt("One output expected.");

	nDims = mxGetNumberOfDimensions(mxArray_InputImage); id=mxGetClassID(mxArray_InputImage);
	ns = (const size_t *) mxGetDimensions(mxArray_InputImage); nCh=((int)nDims == 2) ? 1 : (int)ns[2];
	if( (nDims!=2 && nDims!=3) ||
	(id!=mxSINGLE_CLASS && id!=mxDOUBLE_CLASS && id!=mxUINT8_CLASS) )
	mexErrMsgTxt("A should be 2D or 3D single, double or uint8 array.");
	ms[0]=sz1_1; ms[1]=sz1_2; ms[2]=nCh;
	if( ms[0]<=0 || ms[1]<=0 ) mexErrMsgTxt("downsampling factor too small.");
	nrm=ratio;

	const size_t ms_[3] = {(size_t)ms[0], (size_t)ms[1], (size_t)ms[2]};
	// create output array
	mxArray* mxArray_OutputImage = mxCreateNumericArray(3, (const size_t*) ms_, id, mxREAL);
	n=((int)ns[0])* ((int)ns[1])*nCh; m=ms[0]*ms[1]*nCh;

	// perform resampling (w appropriate type)
	A=mxGetData(mxArray_InputImage); 
	B=mxGetData(mxArray_OutputImage);
	if( id==mxDOUBLE_CLASS ) {
		resample((double*)A, (double*)B, (int)ns[0], ms[0], (int)ns[1], ms[1], nCh, nrm);
	} else if( id==mxSINGLE_CLASS ) {
		resample((float*)A, (float*)B, (int)ns[0], ms[0], (int)ns[1], ms[1], nCh, float(nrm));
	} else if( id==mxUINT8_CLASS ) {
		uchar *A1 = (uchar*) mxMalloc(n*sizeof(uchar));
		uchar *B1 = (uchar*) mxCalloc(m,sizeof(uchar));
		for(int i=0; i<n; i++) A1[i]=(uchar) ((uchar*)A)[i];
		resample(A1, B1, (int)ns[0], ms[0], (int)ns[1], ms[1], nCh, uchar(nrm));
		for(int i=0; i<m; i++) ((uchar*)B)[i]=(uchar) (B1[i]+.5);
	} else {
	mexErrMsgTxt("Unsupported type.");
	}
	//mxArray *mxArray_OutputImage = mxCreateNumericMatrix(0, 0, id, mxREAL);
	mxSetData(mxArray_OutputImage,B);
	mxSetDimensions(mxArray_OutputImage, ms_, 3);
	//mxDestroyArray(mxArray_InputImage);
	return mxArray_OutputImage;
}  
void CDetectorCrosstalk::chnsCompute(mxArray* chns[], mxArray* InputImage) {
	int cr[2]; size_t sz2[2];

	double shrink = mxGetScalar(stCrosstalkParams_.pr_shrink);
	int smooth = (int)mxGetScalar(stCrosstalkParams_.pr_smooth);
	int enabled = (int)mxGetScalar(stCrosstalkParams_.pr_enabled);
	
	/*mxSetData(pChns,stCrosstalkParams_.pr_pChns);
	mxAddField(info,"name"); mxAddField(info,"pChns"); mxAddField(info,"nChns"); mxAddField(info,"padWith");
	mxAddField(chns,"pChns");  mxAddField(chns,"nTypes"); mxAddField(chns,"data"); mxAddField(chns,"info");*/
	const size_t *sz1 = mxGetDimensions(InputImage);
	cr[0] = MOD((int)sz1[0],(int)shrink);
	cr[1] = MOD((int)sz1[1],(int)shrink);
	if(cr[0] !=0 || cr[1] !=0) {
		sz2[0] = sz1[0] - (size_t)cr[0];
		sz2[1] = sz1[1] - (size_t)cr[1];
	}
	else {
		sz2[0] = sz1[0];
		sz2[1] = sz1[1];
	}
	mxArray *AfterResize = this->MxArrayResize(InputImage,(int)sz2[0],(int)sz2[1]); // allocate memory
	sz2[0] = (size_t)(sz2[0]/shrink); sz2[1] = (size_t)(sz2[1]/shrink);

	mxArray* AfterResize_Convert = CrosstalkRGBConvert(AfterResize,1,1); // flag 1 = orig
	mxDestroyArray(AfterResize);

	mxArray *AfterResize_ConvTri = this->ConvTri(AfterResize_Convert,(float) smooth);
	mxDestroyArray(AfterResize_Convert);// added by BMS

	mxArray* AfterResize_gradientMag[2];	
	this->CorsstalkgradientMag(AfterResize_gradientMag,AfterResize_ConvTri); // 0: M, 1 : O
	mxArray* AfterResize_H = CrosstalkgradientHist(AfterResize_gradientMag[0],AfterResize_gradientMag[1],4, 6, 0, 0, 0.2f, 0);
	mxDestroyArray(AfterResize_gradientMag[1]);
	
	if(enabled) {
		this->CrosstalkAddChn(chns, AfterResize_ConvTri, sz2, 0);
		this->CrosstalkAddChn(chns, AfterResize_gradientMag[0], sz2, 1);
		this->CrosstalkAddChn(chns, AfterResize_H, sz2, 2);
	}
	return;
}

void 
CDetectorCrosstalk::CrosstalkAddChn(mxArray **chns, mxArray *AfterResize_Convert, size_t *sz2, int idx)
{
	const size_t *sz1 = mxGetDimensions(AfterResize_Convert);
	if(sz1[0] != sz2[0] || sz1[1] != sz2[1]) {
		mxArray *AfterResize_Convert_chns = CrosstalkIMResample(AfterResize_Convert,(int)sz2[0],(int)sz2[1]);
		chns[idx] = AfterResize_Convert_chns;
		mxDestroyArray(AfterResize_Convert);
	}
	else chns[idx] = AfterResize_Convert;
}

mxArray* 
CDetectorCrosstalk::ConvTri(mxArray* I, float r, int s, int nomex)
{
	mxArray* output = NULL;
	if(nomex==0)
	{
		if( r>0 && r<=1 && s<=2 )		
			output = convConst("convTri1", I, 12/r/(r+2)-2, s);		
		else
			output = convConst("convTri", I, r, s);	
	}
	else
	{
		printf("error at ConvTri");
	}	
	return output;
}

mxArray* CDetectorCrosstalk::convConst(char* type, mxArray* I, float p, int s)
{
  int nDims, d, m, r; float *A, *B;
  const size_t* ns;
  size_t ms[3];
  mxClassID id;
  mxArray* output;

  // error checking on arguments
  //if(nrhs!=4) mexErrMsgTxt("Four inputs required.");
  //if(nlhs > 1) mexErrMsgTxt("One output expected.");
  nDims = (int) mxGetNumberOfDimensions(I);
  id = mxGetClassID(I);
  ns = (size_t*) mxGetDimensions(I);
  d = (nDims == 3) ? (int)ns[2] : 1;
  m = ((int)ns[0] < (int)ns[1]) ? (int)ns[0] : (int)ns[1];
  if( (nDims!=2 && nDims!=3) || id!=mxSINGLE_CLASS || m<4 )
    mexErrMsgTxt("A must be a 4x4 or bigger 2D or 3D float array.");

  // extract inputs
  //if(mxGetString(prhs[0],type,1024))
  //  mexErrMsgTxt("Failed to get type.");
  A = (float*) mxGetData(I);
  r = (int) p;
  //p = (float) mxGetScalar(prhs[2]);
  //r = (int) mxGetScalar(prhs[2]);
  //s = (int) mxGetScalar(prhs[3]);A
  if( s<1 ) mexErrMsgTxt("Invalid sampling value s");
  if( r<0 ) mexErrMsgTxt("Invalid radius r");

  // create output array (w/o initializing to 0)
  ms[0]=ns[0]/s; ms[1]=ns[1]/s; ms[2]=d;
  output = mxCreateNumericArray(3, ms, mxSINGLE_CLASS, mxREAL);
  B = (float*) mxGetData(output);
  //B = (float*) mxMalloc(ms[0]*ms[1]*d*sizeof(float));
  //output = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
  //mxSetData(output, B); mxSetDimensions(output,(mwSize*)ms,nDims);

  // perform appropriate type of convolution
  if(!strcmp(type,"convBox")) {
    if(r>=m/2) mexErrMsgTxt("mask larger than image (r too large)");
    convBox( A, B, (int)ns[0], (int)ns[1], d, r, s );
  } else if(!strcmp(type,"convTri")) {
    if(r>=m/2) mexErrMsgTxt("mask larger than image (r too large)");
    convTri( A, B, (int)ns[0], (int)ns[1], d, r, s );
  } else if(!strcmp(type,"conv11")) {
    if( s>2 ) mexErrMsgTxt("conv11 can sample by at most s=2");
    conv11( A, B, (int)ns[0], (int)ns[1], d, r, s );
  } else if(!strcmp(type,"convTri1")) {
    if( s>2 ) mexErrMsgTxt("convTri1 can sample by at most s=2");
    convTri1( A, B, (int)ns[0], (int)ns[1], d, p, s );
  } else if(!strcmp(type,"convMax")) {
    if( s>1 ) mexErrMsgTxt("convMax cannot sample");
    convMax( A, B, (int)ns[0], (int)ns[1], d, r );
  } else {
    mexErrMsgTxt("Invalid type.");
  }

  return output;
}

void CDetectorCrosstalk::CorsstalkgradientMag(mxArray* pl[], mxArray* I, int channel, int normRad, float NormConst, int full)
{
	int nl,nr; 
	const size_t *sz1; 
	size_t sz2[2];

	sz1 = mxGetDimensions(I);
	sz2[0] = sz1[0];
	sz2[1] = sz1[1];
	nl = 2;
	if(nl==2)
	{
		nr = 4;				 
		mxArray* type = mxCreateString("gradientMag");	
		mxArray* _channel = mxCreateDoubleScalar(channel);
		mxArray* _full = mxCreateDoubleScalar(full);
		const mxArray* pr[4] = {type, I, _channel, _full};
		gradientMex( nl, pl, nr, pr);
		mxDestroyArray(_channel);
		mxDestroyArray(_full);
		mxDestroyArray(type);
	}
	else if(nl==1)	
	{
		//TODO
	}
	mxSetDimensions(pl[0],sz2,2);

	float * pr_ = (float *)mxGetData(pl[0]);
	//normalize
	if(normRad==0){
		return;
	}

	mxArray*S = ConvTri(pl[0], (float)normRad);	
	{
		nr = 4;				 
		mxArray* type = mxCreateString("gradientMagNorm");	
		mxArray* M = pl[0];
		mxArray* _NormConst = mxCreateDoubleScalar(NormConst);
		const mxArray* pr[4] = {type, M, S, _NormConst};
		gradientMex( NULL, NULL, nr, pr);
		mxDestroyArray(_NormConst);
		mxDestroyArray(type);
	}		
	mxDestroyArray(S);
	return;
	
}

void CDetectorCrosstalk::gradientMex( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  int f; char action[1024]; f=mxGetString(pr[0],action,1024); nr--; pr++;
  if(f) mexErrMsgTxt("Failed to get action.");
  else if(!strcmp(action,"gradient2")) mGrad2(nl,pl,nr,pr);
  else if(!strcmp(action,"gradientMag")) mGradMag(nl,pl,nr,pr);
  else if(!strcmp(action,"gradientMagNorm")) mGradMagNorm(nl,pl,nr,pr);
  else if(!strcmp(action,"gradientHist")) mGradHist(nl,pl,nr,pr);
  else mexErrMsgTxt("Invalid action.");
    
}

mxArray* CDetectorCrosstalk::CrosstalkgradientHist(mxArray* M,mxArray* O, int binSize, int nOrients, int softBin, int useHog, float clipHog, int full)
{     
   int nl = 1;
   mxArray* pl[1];   
   mxArray* type = mxCreateString("gradientHist");      
   mxArray*_binSize = mxCreateDoubleScalar(binSize);
   mxArray*_nOrients = mxCreateDoubleScalar(nOrients);
   mxArray*_softBin = mxCreateDoubleScalar(softBin);
   mxArray*_useHog = mxCreateDoubleScalar(useHog);
   mxArray*_clipHog = mxCreateDoubleScalar(clipHog);
   mxArray*_full = mxCreateDoubleScalar(full);

   int nr = 9;
   mxArray* pr[9] = {type, M, O, _binSize,_nOrients,_softBin,_useHog,_clipHog,_full};
   gradientMex( nl, pl, nr, (const mxArray**) pr);
   mxDestroyArray(type);
   mxDestroyArray(_binSize);
   mxDestroyArray(_nOrients);
   mxDestroyArray(_softBin);
   mxDestroyArray(_useHog);
   mxDestroyArray(_clipHog);
   mxDestroyArray(_full);
   //free(type); free(_binSize); free(_nOrients); free(_softBin); free(_useHog); free(_clipHog); free(_full);
   return pl[0];
}

// B = imPadMex(A,pad,type); see imPad.m for usage details
mxArray * 
CDetectorCrosstalk::CrosstalkImPadMex(mxArray *chn, mxArray* pad, mxArray* padWidth) {
  int ms[3], nCh, nDims, pt, pb, pl, pr, flag, k; double *p;
  void *A, *B; mxClassID id; double val=0;
  const size_t* ns;
  //// Error checking on arguments
  //if( nrhs!=3 ) mexErrMsgTxt("Three inputs expected.");
  //if( nlhs>1 ) mexErrMsgTxt("One output expected.");
  nDims=(int)mxGetNumberOfDimensions(chn); id=mxGetClassID(chn);
  ns = mxGetDimensions(chn); nCh=(nDims==2) ? 1 : (int)ns[2];
  if( (nDims!=2 && nDims!=3) ||
    (id!=mxSINGLE_CLASS && id!=mxDOUBLE_CLASS && id!=mxUINT8_CLASS) )
    mexErrMsgTxt("A should be 2D or 3D single, double or uint8 array.");
  if( !mxIsDouble(pad) ) mexErrMsgTxt("Input pad must be a double array.");

  // extract padding amounts
  k = (int) mxGetNumberOfElements(pad);
  p = (double*) mxGetData(pad);
  if(k==1) { pt=pb=pl=pr=int(p[0]); }
  else if (k==2) { pt=pb=int(p[0]); pl=pr=int(p[1]); }
  else if (k==4) { pt=int(p[0]); pb=int(p[1]); pl=int(p[2]); pr=int(p[3]); }
  else mexErrMsgTxt( "Input pad must have 1, 2, or 4 values.");

 // // figure out padding type (flag and val)
 // if( !mxGetString(padWidth,type,1024) ) {
 //   if(!strcmp(type,"replicate")) flag=1;
 //   else if(!strcmp(type,"symmetric")) flag=2;
 //   else if(!strcmp(type,"circular")) flag=3;
	//else {flag=0; val=(double)mxGetScalar(padWidth);}
 //   //else mexErrMsgTxt("Invalid pad value.");
 // } else {
 //   flag=0; val=(double)mxGetScalar(padWidth);
 // }
    // figure out padding type (flag and val)
  if(!strcmp((char*)padWidth->data,"replicate")) flag=1;
  else if(!strcmp((char*)padWidth->data,"symmetric")) flag=2;
  else if(!strcmp((char*)padWidth->data,"circular")) flag=3;
  else {flag=0; val=(double)mxGetScalar(padWidth);}
    //else mexErrMsgTxt("Invalid pad value.");
  if( ns[0]==0 || ns[1]==0 ) flag=0;

  // create output array
  ms[0]=(int)ns[0]+pt+pb; ms[1]=(int)ns[1]+pl+pr; ms[2]=nCh;
  if( ms[0]<0 || (int)ns[0]<=-pt || (int)ns[0]<=-pb ) ms[0]=0;
  if( ms[1]<0 || (int)ns[1]<=-pl || (int)ns[1]<=-pr ) ms[1]=0;
  size_t ms_[3]; ms_[0] = (size_t)ms[0]; ms_[1] = (size_t)ms[1]; ms_[2] = (size_t)ms[2];
  mxArray *chn_pad = mxCreateNumericArray(3, (const size_t*) ms_, id, mxREAL);
  if( ms[0]==0 || ms[1]==0 ) return chn_pad;

  // pad array
  A=mxGetData(chn); 
  B=mxGetData(chn_pad);
  if( id==mxDOUBLE_CLASS ) {
    imPad( (double*)A,(double*)B,(int)ns[0],(int)ns[1],nCh,pt,pb,pl,pr,flag,val );
  } else if( id==mxSINGLE_CLASS ) {
    imPad( (float*)A,(float*)B,(int)ns[0],(int)ns[1],nCh,pt,pb,pl,pr,flag,float(val) );
  } else if( id==mxUINT8_CLASS ) {
    imPad( (uchar*)A,(uchar*)B,(int)ns[0],(int)ns[1],nCh,pt,pb,pl,pr,flag,uchar(val) );
  } else {
    mexErrMsgTxt("Unsupported image type.");
  }
  return chn_pad;
}

// 코 //
mxArray*
CDetectorCrosstalk::CrosstalkAcfDetect(mxArray * chn_concat, int curStep)
{
	double *shrink_ = (double*) mxGetData(stCrosstalkParams_.pr_shrink);
	mxArray *mxArray_modelDsPad = mxGetField(stCrosstalkParams_.mxArray_opts,0,"modelDsPad");
	double *modelDsPad = (double*) mxGetData(mxArray_modelDsPad);
	mxArray *mxArray_stride = mxGetField(stCrosstalkParams_.mxArray_opts,0,"stride");
	double *stride_ = (double*) mxGetData(mxArray_stride);
	mxArray *mxArray_cascThr = mxGetField(stCrosstalkParams_.mxArray_opts,0,"cascThr");
	double *cascThr_ = (double*) mxGetData(mxArray_cascThr);

	// get inputs
	float *chns = (float*) mxGetData(chn_concat);
	mxArray *trees = (mxArray*) stCrosstalkParams_.mxArray_clf;
	const int shrink = (int) shrink_[0];
	const int modelHt = (int)  modelDsPad[0];
	const int modelWd = (int) modelDsPad[1];
	const int stride = (int) stride_[0];
	const float cascThr = (float) cascThr_[0];

	// extract relevant fields from trees
	float *thrs = (float*) mxGetData(mxGetField(trees,0,"thrs"));
	float *hs = (float*) mxGetData(mxGetField(trees,0,"hs"));
	uint32 *fids = (uint32*) mxGetData(mxGetField(trees,0,"fids"));
	uint32 *child = (uint32*) mxGetData(mxGetField(trees,0,"child"));
	int* treeDepth_ = mxGetField(trees,0,"treeDepth")==NULL ? 0 :
		(int*) mxGetData(mxGetField(trees,0,"treeDepth"));
	const int treeDepth = treeDepth_[0];

	// get dimensions and constants
	// const mwSize *chnsSize = mxGetDimensions(prhs[0]);
	const size_t *chnsSize = mxGetDimensions(chn_concat);
	const int height = (int) chnsSize[0];
	int start_r;
	int end_r;
	if(curStep < 4)
	{
		start_r = 0;
		end_r = height/2;
	}
	else if(curStep < 10)
	{
		start_r = 0;
		end_r = (int)(height/1.5);
	}
	else
	{
		start_r = height/4;
		end_r = height;
	}
	const int width = (int) chnsSize[1];
	// const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
	const int nChns = mxGetNumberOfDimensions(chn_concat)<=2 ? 1 : (int) chnsSize[2];
	const size_t *fidsSize = mxGetDimensions(mxGetField(trees,0,"fids"));
	const int nTreeNodes = (int) fidsSize[0];
	const int nTrees = (int) fidsSize[1];
	//const int height1 = (int) ceil(float(height*shrink-modelHt+1)/stride
	const int height1 = (int) ceil(float(end_r*shrink-modelHt+1)/stride);
	const int width1 = (int) ceil(float(width*shrink-modelWd+1)/stride);

	// construct cids array
	int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
	uint32 *cids = new uint32[nFtrs]; int m=0;
	for( int z=0; z<nChns; z++ )
	for( int c=0; c<modelWd/shrink; c++ )
		for( int r=0; r<modelHt/shrink; r++ )
		cids[m++] = z*width*height + c*height + r;

	// apply classifier to each patch
	vector<int> rs, cs; vector<float> hs1;
	for( int c=0; c<width1; c++ ) for( int r=start_r; r<height1; r++ ) {
	float h=0, *chns1=chns+(r*stride/shrink) + (c*stride/shrink)*height;
	if( treeDepth==1 ) {
		// specialized case for treeDepth==1
		for( int t = 0; t < nTrees; t++ ) {
		uint32 offset=t*nTreeNodes, k=offset, k0=0;
		getChild(chns1,cids,fids,thrs,offset,k0,k);
		h += hs[k]; if( h<=cascThr ) break;
		}
	} else if( treeDepth==2 ) {
		// specialized case for treeDepth==2
		for( int t = 0; t < nTrees; t++ ) {
		uint32 offset=t*nTreeNodes, k=offset, k0=0;
		getChild(chns1,cids,fids,thrs,offset,k0,k);
		getChild(chns1,cids,fids,thrs,offset,k0,k);
		h += hs[k]; if( h<=cascThr ) break;
		}
	} else if( treeDepth>2) {
		// specialized case for treeDepth>2
		for( int t = 0; t < nTrees; t++ ) {
		uint32 offset=t*nTreeNodes, k=offset, k0=0;
		for( int i=0; i<treeDepth; i++ )
			getChild(chns1,cids,fids,thrs,offset,k0,k);
		h += hs[k]; if( h<=cascThr ) break;
		}
	} else {
		// general case (variable tree depth)
		for( int t = 0; t < nTrees; t++ ) {
		uint32 offset=t*nTreeNodes, k=offset, k0=k;
		while( child[k] ) {
			float ftr = chns1[cids[fids[k]]];
			k = (ftr<thrs[k]) ? 1 : 0;
			k0 = k = child[k0]-k+offset;
		}
		h += hs[k]; if( h<=cascThr ) break;
		}
	}
	if(h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }
	}
	delete [] cids; m=(int)cs.size();

	// convert to bbs
	//plhs[0] = mxCreateNumericMatrix(m,5,mxDOUBLE_CLASS,mxREAL);
	mxArray* bb = mxCreateNumericMatrix(m,5,mxDOUBLE_CLASS,mxREAL);
	//double *bbs = (double*) mxGetData(plhs[0]);
	double *bbs = (double*) mxGetData(bb);
	for( int i=0; i<m; i++ ) {
	bbs[i+0*m]=cs[i]*stride; bbs[i+2*m]=modelWd;
	bbs[i+1*m]=rs[i]*stride; bbs[i+3*m]=modelHt;
	bbs[i+4*m]=hs1[i];
	}
	return bb;
}

// 코 //
vector<hj::CDetection>
CDetectorCrosstalk::bbNms_maxg(vector<hj::CDetection> st_bbs, double overlap)
{
	double thr = -1000;
	double overDnm = 1;
	double greedy = 1;
	
	sort(st_bbs.begin(), st_bbs.end(), SortDetectionScoreDescend);
	bool* kp = (bool*)malloc(sizeof(bool)*st_bbs.size());
	double* as = (double*)malloc(sizeof(double)*st_bbs.size());
	for(int i=0;i<st_bbs.size();i++)
	{
		kp[i] = true;
		as[i] = st_bbs[i].box.w * st_bbs[i].box.h;
	}

	for(int i=0;i<st_bbs.size();i++)
	{
		if(greedy && (kp[i]==false)) continue;
		for(int j=i+1;j<st_bbs.size();j++)
		{
			if(kp[j]==false) continue;
			double iw = min(st_bbs[i].box.x+st_bbs[i].box.w,st_bbs[j].box.x+st_bbs[j].box.w)-max(st_bbs[i].box.x,st_bbs[j].box.x);
			if(iw<=0) continue;
			double ih = min(st_bbs[i].box.y+st_bbs[i].box.h,st_bbs[j].box.y+st_bbs[j].box.h)-max(st_bbs[i].box.y,st_bbs[j].box.y);
			if(ih<=0) continue;
			double o = iw*ih;
			double u;
			if(overDnm) u = as[i]+as[j]-o;
			else u = min(as[i],as[j]);
			o = o/u;
			if(o>overlap) kp[j]=false;
		}
	}

	vector<hj::CDetection> st_bbs_valid;
	for(int i=0;i<st_bbs.size();i++)
	{
		if(kp[i]==true) st_bbs_valid.push_back(st_bbs[i]);
	}

	free(kp);
	free(as);
	return st_bbs_valid;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//convConst
///////////////////////////////////////////////////////////////////////////////////////////////
//convolve one column of I by a 2rx1 ones filter
void 
CDetectorCrosstalk::convBoxY( float *I, float *O, int h, int r, int s ) {
  float t; int j, p=r+1, q=2*h-(r+1), h0=r+1, h1=h-r, h2=h;
  t=0; for(j=0; j<=r; j++) t+=I[j]; t=2*t-I[r]; j=0;
  if( s==1 ) {
    for(; j<h0; j++) O[j]=t-=I[r-j]-I[r+j];
    for(; j<h1; j++) O[j]=t-=I[j-p]-I[r+j];
    for(; j<h2; j++) O[j]=t-=I[j-p]-I[q-j];
  } else {
    int k=(s-1)/2; h2=(h/s)*s; if(h0>h2) h0=h2; if(h1>h2) h1=h2;
    for(; j<h0; j++) { t-=I[r-j]-I[r+j]; k++; if(k==s) { k=0; *O++=t; } }
    for(; j<h1; j++) { t-=I[j-p]-I[r+j]; k++; if(k==s) { k=0; *O++=t; } }
    for(; j<h2; j++) { t-=I[j-p]-I[q-j]; k++; if(k==s) { k=0; *O++=t; } }
  }
}

// convolve I by a 2r+1 x 2r+1 ones filter (uses SSE)
void 
CDetectorCrosstalk::convBox( float *I, float *O, int h, int w, int d, int r, int s ) {
  float nrm = 1.0f/((2*r+1)*(2*r+1)); int i, j, k=(s-1)/2, h0, h1, w0;
  if(h%4==0) h0=h1=h; else { h0=h-(h%4); h1=h0+4; } w0=(w/s)*s;
  float *T=(float*) alMalloc(h1*sizeof(float),16);
  while(d-- > 0) {
    // initialize T
    memset( T, 0, h1*sizeof(float) );
    for(i=0; i<=r; i++) for(j=0; j<h0; j+=4) INC(T[j],LDu(I[j+i*h]));
    for(j=0; j<h0; j+=4) STR(T[j],MUL(nrm,SUB(MUL(2,LD(T[j])),LDu(I[j+r*h]))));
    for(i=0; i<=r; i++) for(j=h0; j<h; j++ ) T[j]+=I[j+i*h];
    for(j=h0; j<h; j++ ) T[j]=nrm*(2*T[j]-I[j+r*h]);
    // prepare and convolve each column in turn
    k++; if(k==s) { k=0; convBoxY(T,O,h,r,s); O+=h/s; }
    for( i=1; i<w0; i++ ) {
      float *Il=I+(i-1-r)*h; if(i<=r) Il=I+(r-i)*h;
      float *Ir=I+(i+r)*h; if(i>=w-r) Ir=I+(2*w-r-i-1)*h;
      for(j=0; j<h0; j+=4) DEC(T[j],MUL(nrm,SUB(LDu(Il[j]),LDu(Ir[j]))));
      for(j=h0; j<h; j++ ) T[j]-=nrm*(Il[j]-Ir[j]);
      k++; if(k==s) { k=0; convBoxY(T,O,h,r,s); O+=h/s; }
    }
    I+=w*h;
  }
  alFree(T);
}

// convolve one column of I by a [1; 1] filter (uses SSE)
void
CDetectorCrosstalk::conv11Y( float *I, float *O, int h, int side, int s ) {
  #define C4(m,o) ADD(LDu(I[m*j-1+o]),LDu(I[m*j+o]))
  int j=0, k=((~((size_t) O) + 1) & 15)/4;
  const int d = (side % 4 >= 2) ? 1 : 0, h2=(h-d)/2;
  if( s==2 ) {
    for( ; j<k; j++ ) O[j]=I[2*j+d]+I[2*j+d+1];
    for( ; j<h2-4; j+=4 ) STR(O[j],_mm_shuffle_ps(C4(2,d+1),C4(2,d+5),136));
    for( ; j<h2; j++ ) O[j]=I[2*j+d]+I[2*j+d+1];
    if(d==1 && h%2==0) O[j]=2*I[2*j+d];
  } else {
    if(d==0) { O[0]=2*I[0]; j++; if(k==0) k=4; }
    for( ; j<k; j++ ) O[j]=I[j-1+d]+I[j+d];
    for( ; j<h-4-d; j+=4 ) STR(O[j],C4(1,d) );
    for( ; j<h-d; j++ ) O[j]=I[j-1+d]+I[j+d];
    if(d==1) { O[j]=2*I[j]; j++; }
  }
  #undef C4
}

// convolve I by a [1 1; 1 1] filter (uses SSE)
void CDetectorCrosstalk::conv11( float *I, float *O, int h, int w, int d, int side, int s ) {
  const float nrm = 0.25f; int i, j;
  float *I0, *I1, *T = (float*) alMalloc(h*sizeof(float),16);
  for( int d0=0; d0<d; d0++ ) for( i=s/2; i<w; i+=s ) {
    I0=I1=I+i*h+d0*h*w; if(side%2) { if(i<w-1) I1+=h; } else { if(i) I0-=h; }
    for( j=0; j<h-4; j+=4 ) STR( T[j], MUL(nrm,ADD(LDu(I0[j]),LDu(I1[j]))) );
    for( ; j<h; j++ ) T[j]=nrm*(I0[j]+I1[j]);
    conv11Y(T,O,h,side,s); O+=h/s;
  }
  alFree(T);
}

// convolve one column of I by a 2rx1 triangle filter
void CDetectorCrosstalk::convTriY( float *I, float *O, int h, int r, int s ) {
  r++; float t, u; int j, r0=r-1, r1=r+1, r2=2*h-r, h0=r+1, h1=h-r+1, h2=h;
  u=t=I[0]; for( j=1; j<r; j++ ) u+=t+=I[j]; u=2*u-t; t=0;
  if( s==1 ) {
    O[0]=u; j=1;
    for(; j<h0; j++) O[j] = u += t += I[r-j]  + I[r0+j] - 2*I[j-1];
    for(; j<h1; j++) O[j] = u += t += I[j-r1] + I[r0+j] - 2*I[j-1];
    for(; j<h2; j++) O[j] = u += t += I[j-r1] + I[r2-j] - 2*I[j-1];
  } else {
    int k=(s-1)/2; h2=(h/s)*s; if(h0>h2) h0=h2; if(h1>h2) h1=h2;
    if(++k==s) { k=0; *O++=u; } j=1;
    for(;j<h0;j++) { u+=t+=I[r-j] +I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h1;j++) { u+=t+=I[j-r1]+I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h2;j++) { u+=t+=I[j-r1]+I[r2-j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
  }
}

// convolve I by a 2rx1 triangle filter (uses SSE)
void CDetectorCrosstalk::convTri( float *I, float *O, int h, int w, int d, int r, int s ) {
  r++; float nrm = 1.0f/(r*r*r*r); int i, j, k=(s-1)/2, h0, h1, w0;
  if(h%4==0) h0=h1=h; else { h0=h-(h%4); h1=h0+4; } w0=(w/s)*s;
  float *T=(float*) alMalloc(2*h1*sizeof(float),16), *U=T+h1;
  while(d-- > 0) {
    // initialize T and U
    for(j=0; j<h0; j+=4) STR(U[j], STR(T[j], LDu(I[j])));
    for(i=1; i<r; i++) for(j=0; j<h0; j+=4) INC(U[j],INC(T[j],LDu(I[j+i*h])));
    for(j=0; j<h0; j+=4) STR(U[j],MUL(nrm,(SUB(MUL(2,LD(U[j])),LD(T[j])))));
    for(j=0; j<h0; j+=4) STR(T[j],0);
    for(j=h0; j<h; j++ ) U[j]=T[j]=I[j];
    for(i=1; i<r; i++) for(j=h0; j<h; j++ ) U[j]+=T[j]+=I[j+i*h];
    for(j=h0; j<h; j++ ) { U[j] = nrm * (2*U[j]-T[j]); T[j]=0; }
    // prepare and convolve each column in turn
    k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    for( i=1; i<w0; i++ ) {
      float *Il=I+(i-1-r)*h; if(i<=r) Il=I+(r-i)*h; float *Im=I+(i-1)*h;
      float *Ir=I+(i-1+r)*h; if(i>w-r) Ir=I+(2*w-r-i)*h;
      for( j=0; j<h0; j+=4 ) {
        INC(T[j],ADD(LDu(Il[j]),LDu(Ir[j]),MUL(-2,LDu(Im[j]))));
        INC(U[j],MUL(nrm,LD(T[j])));
      }
      for( j=h0; j<h; j++ ) U[j]+=nrm*(T[j]+=Il[j]+Ir[j]-2*Im[j]);
      k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    }
    I+=w*h;
  }
  alFree(T);
}

// convolve one column of I by a [1 p 1] filter (uses SSE)
void CDetectorCrosstalk::convTri1Y( float *I, float *O, int h, float p, int s ) {
  #define C4(m,o) ADD(ADD(LDu(I[m*j-1+o]),MUL(p,LDu(I[m*j+o]))),LDu(I[m*j+1+o]))
  int j=0, k=((~((size_t) O) + 1) & 15)/4, h2=(h-1)/2;
  if( s==2 ) {
    for( ; j<k; j++ ) O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    for( ; j<h2-4; j+=4 ) STR(O[j],_mm_shuffle_ps(C4(2,1),C4(2,5),136));
    for( ; j<h2; j++ ) O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    if( h%2==0 ) O[j]=I[2*j]+(1+p)*I[2*j+1];
  } else {
    O[j]=(1+p)*I[j]+I[j+1]; j++; if(k==0) k=(h<=4) ? h-1 : 4;
    for( ; j<k; j++ ) O[j]=I[j-1]+p*I[j]+I[j+1];
    for( ; j<h-4; j+=4 ) STR(O[j],C4(1,0));
    for( ; j<h-1; j++ ) O[j]=I[j-1]+p*I[j]+I[j+1];
    O[j]=I[j-1]+(1+p)*I[j];
  }
  #undef C4
}

// convolve I by a [1 p 1] filter (uses SSE)
void CDetectorCrosstalk::convTri1( float *I, float *O, int h, int w, int d, float p, int s ) {
  const float nrm = 1.0f/((p+2)*(p+2)); int i, j, h0=h-(h%4);
  float *Il, *Im, *Ir, *T=(float*) alMalloc(h*sizeof(float),16);
  for( int d0=0; d0<d; d0++ ) for( i=s/2; i<w; i+=s ) {
    Il=Im=Ir=I+i*h+d0*h*w; if(i>0) Il-=h; if(i<w-1) Ir+=h;
    for( j=0; j<h0; j+=4 )
      STR(T[j],MUL(nrm,ADD(ADD(LDu(Il[j]),MUL(p,LDu(Im[j]))),LDu(Ir[j]))));
    for( j=h0; j<h; j++ ) T[j]=nrm*(Il[j]+p*Im[j]+Ir[j]);
    convTri1Y(T,O,h,p,s); O+=h/s;
  }
  alFree(T);
}

// convolve one column of I by a 2rx1 max filter
void CDetectorCrosstalk::convMaxY( float *I, float *O, float *T, int h, int r ) {
  int y, y0, y1, yi, m=2*r+1;
  #define max1(a,b) a>b ? a : b;
  #define maxk(y0,y1) { O[y]=I[y0]; \
    for( yi=y0+1; yi<=y1; yi++ ) { if(I[yi]>O[y]) O[y]=I[yi]; }}
  for( y=0; y<r; y++ ) { y1=y+r; if(y1>h-1) y1=h-1; maxk(0,y1); }
  for( ; y<=h-m-r; y+=m ) {
    T[m-1] = I[y+r];
    for( yi=1; yi<m; yi++ ) T[m-1-yi] = max1( T[m-1-yi+1], I[y+r-yi] );
    for( yi=1; yi<m; yi++ ) T[m-1+yi] = max1( T[m-1+yi-1], I[y+r+yi] );
    for( yi=0; yi<m; yi++ ) O[y+yi] = max1( T[yi], T[yi+m-1] );
  }
  for( ; y<h-r; y++ ) { maxk(y-r,y+r); }
  for( ; y<h; y++ ) { y0=y-r; if(y0<0) y0=0; maxk(y0,h-1); }
  #undef maxk
  #undef max1
}

// convolve I by a 2rx1 max filter
void CDetectorCrosstalk::convMax( float *I, float *O, int h, int w, int d, int r ) {
  if( r>w-1 ) r=w-1; if( r>h-1 ) r=h-1; int m=2*r+1;
  float *T=(float*) alMalloc(m*2*sizeof(float),16);
  for( int d0=0; d0<d; d0++ ) for( int x=0; x<w; x++ ) {
    float *Oc=O+d0*h*w+h*x, *Ic=I+d0*h*w+h*x;
    convMaxY(Ic,Oc,T,h,r);
  }
  alFree(T);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//gradientMag
///////////////////////////////////////////////////////////////////////////////////////////////
// compute x and y gradients for just one column (uses sse)
void CDetectorCrosstalk::grad1( float *I, float *Gx, float *Gy, int h, int w, int x ) {
  int y, y1; float *Ip, *In, r; __m128 *_Ip, *_In, *_G, _r;
  // compute column of Gx
  Ip=I-h; In=I+h; r=.5f;
  if(x==0) { r=1; Ip+=h; } else if(x==w-1) { r=1; In-=h; }
  if( h<4 || h%4>0 || (size_t(I)&15) || (size_t(Gx)&15) ) {
    for( y=0; y<h; y++ ) *Gx++=(*In++-*Ip++)*r;
  } else {
    _G=(__m128*) Gx; _Ip=(__m128*) Ip; _In=(__m128*) In; _r = SET(r);
    for(y=0; y<h; y+=4) *_G++=MUL(SUB(*_In++,*_Ip++),_r);
  }
  // compute column of Gy
  #define GRADY(r) *Gy++=(*In++-*Ip++)*r;
  Ip=I; In=Ip+1;
  // GRADY(1); Ip--; for(y=1; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  y1=((~((size_t) Gy) + 1) & 15)/4; if(y1==0) y1=4; if(y1>h-1) y1=h-1;
  GRADY(1); Ip--; for(y=1; y<y1; y++) GRADY(.5f);
  _r = SET(.5f); _G=(__m128*) Gy;
  for(; y+4<h-1; y+=4, Ip+=4, In+=4, Gy+=4)
    *_G++=MUL(SUB(LDu(*In),LDu(*Ip)),_r);
  for(; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  #undef GRADY
}

// compute x and y gradients at each location (uses sse)
void CDetectorCrosstalk::grad2( float *I, float *Gx, float *Gy, int h, int w, int d ) {
  int o, x, c, a=w*h; for(c=0; c<d; c++) for(x=0; x<w; x++) {
    o=c*a+x*h; grad1( I+o, Gx+o, Gy+o, h, w, x );
  }
}

// build lookup table a[] s.t. a[x*n]~=acos(x) for x in [-1,1]
float* CDetectorCrosstalk::acosTable() {
  const int n=10000, b=10; int i;
  static float a[n*2+b*2]; static bool init=false;
  float *a1=a+n+b; if( init ) return a1;
  for( i=-n-b; i<-n; i++ )   a1[i]=PI;
  for( i=-n; i<n; i++ )      a1[i]=float(acos(i/float(n)));
  for( i=n; i<n+b; i++ )     a1[i]=0;
  for( i=-n-b; i<n/10; i++ ) if( a1[i] > PI-1e-6f ) a1[i]=PI-1e-6f;
  init=true; return a1;
}

// compute gradient magnitude and orientation at each location (uses sse)
void CDetectorCrosstalk::gradMag( float *I, float *M, float *O, int h, int w, int d, bool full ) {
  int x, y, y1, c, h4, s; float *Gx, *Gy, *M2; __m128 *_Gx, *_Gy, *_M2, _m;
  float *acost = acosTable(), acMult=10000.0f;
  // allocate memory for storing one column of output (padded so h4%4==0)
  h4=(h%4==0) ? h : h-(h%4)+4; s=d*h4*sizeof(float);
  M2=(float*) alMalloc(s,16); _M2=(__m128*) M2;
  Gx=(float*) alMalloc(s,16); _Gx=(__m128*) Gx;
  Gy=(float*) alMalloc(s,16); _Gy=(__m128*) Gy;
  // compute gradient magnitude and orientation for each column
  for( x=0; x<w; x++ ) {
    // compute gradients (Gx, Gy) with maximum squared magnitude (M2)
    for(c=0; c<d; c++) {
      grad1( I+x*h+c*w*h, Gx+c*h4, Gy+c*h4, h, w, x );
      for( y=0; y<h4/4; y++ ) {
        y1=h4/4*c+y;
        _M2[y1]=ADD(MUL(_Gx[y1],_Gx[y1]),MUL(_Gy[y1],_Gy[y1]));
        if( c==0 ) continue; _m = CMPGT( _M2[y1], _M2[y] );
        _M2[y] = OR( AND(_m,_M2[y1]), ANDNOT(_m,_M2[y]) );
        _Gx[y] = OR( AND(_m,_Gx[y1]), ANDNOT(_m,_Gx[y]) );
        _Gy[y] = OR( AND(_m,_Gy[y1]), ANDNOT(_m,_Gy[y]) );
      }
    }
    // compute gradient mangitude (M) and normalize Gx
    for( y=0; y<h4/4; y++ ) {      
	  _m = MIN_( RCPSQRT(_M2[y]), SET(1e10f) );		
      _M2[y] = RCP(_m);
      if(O) _Gx[y] = MUL( MUL(_Gx[y],_m), SET(acMult) );
      if(O) _Gx[y] = XOR( _Gx[y], AND(_Gy[y], SET(-0.f)) );
    };
    memcpy( M+x*h, M2, h*sizeof(float) );
    // compute and store gradient orientation (O) via table lookup
    if( O!=0 ) for( y=0; y<h; y++ ) O[x*h+y] = acost[(int)Gx[y]];
    if( O!=0 && full ) {
      y1=((~size_t(O+x*h)+1)&15)/4; y=0;
      for( ; y<y1; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
      for( ; y<h-4; y+=4 ) STRu( O[y+x*h],
        ADD( LDu(O[y+x*h]), AND(CMPLT(LDu(Gy[y]),SET(0.f)),SET(PI)) ) );
      for( ; y<h; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
    }
  }
  alFree(Gx); alFree(Gy); alFree(M2);
}

// normalize gradient magnitude at each location (uses sse)
void CDetectorCrosstalk::gradMagNorm( float *M, float *S, int h, int w, float norm ) {
  __m128 *_M, *_S, _norm; int i=0, n=h*w, n4=n/4;
  _S = (__m128*) S; _M = (__m128*) M; _norm = SET(norm);
  bool sse = !(size_t(M)&15) && !(size_t(S)&15);
  if(sse) { for(; i<n4; i++) *_M++=MUL(*_M,RCP(ADD(*_S++,_norm))); i*=4; }
  for(; i<n; i++) M[i] /= (S[i] + norm);
}

// helper for gradHist, quantize O and M into O0, O1 and M0, M1 (uses sse)
void CDetectorCrosstalk::gradQuantize( float *O, float *M, int *O0, int *O1, float *M0, float *M1,
  int nb, int n, float norm, int nOrients, bool full, bool interpolate )
{
  // assumes all *OUTPUT* matrices are 4-uchar aligned
  int i, o0, o1; float o, od, m;
  __m128i _o0, _o1, *_O0, *_O1; __m128 _o, _od, _m, *_M0, *_M1;
  // define useful constants
  const float oMult=(float)nOrients/(full?2*PI:PI); const int oMax=nOrients*nb;
  const __m128 _norm=SET(norm), _oMult=SET(oMult), _nbf=SET((float)nb);
  const __m128i _oMax=SET(oMax), _nb=SET(nb);
  // perform the majority of the work with sse
  _O0=(__m128i*) O0; _O1=(__m128i*) O1; _M0=(__m128*) M0; _M1=(__m128*) M1;
  if( interpolate ) for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(_o); _od=SUB(_o,CVT(_o0));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    _o1=ADD(_o0,_nb); _o1=AND(CMPGT(_oMax,_o1),_o1); *_O1++=_o1;
    _m=MUL(LDu(M[i]),_norm); *_M1=MUL(_od,_m); *_M0++=SUB(_m,*_M1); _M1++;
  } else for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(ADD(_o,SET(.5f)));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    *_M0++=MUL(LDu(M[i]),_norm); *_M1++=SET(0.f); *_O1++=SET(0);
  }
  // compute trailing locations without sse
  if( interpolate ) for( i; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) o; od=o-o0;
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    o1=o0+nb; if(o1==oMax) o1=0; O1[i]=o1;
    m=M[i]*norm; M1[i]=od*m; M0[i]=m-M1[i];
  } else for( i; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) (o+.5f);
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    M0[i]=M[i]*norm; M1[i]=0; O1[i]=0;
  }
}

// compute nOrients gradient histograms per bin x bin block of pixels
void CDetectorCrosstalk::gradHist( float *M, float *O, float *H, int h, int w,
  int bin, int nOrients, int softBin, bool full )
{
  const int hb=h/bin, wb=w/bin, h0=hb*bin, w0=wb*bin, nb=wb*hb;
  const float s=(float)bin, sInv=1/s, sInv2=1/s/s;
  float *H0, *H1, *M0, *M1; int x, y; int *O0, *O1;
  O0=(int*)alMalloc(h*sizeof(int),16); M0=(float*) alMalloc(h*sizeof(float),16);
  O1=(int*)alMalloc(h*sizeof(int),16); M1=(float*) alMalloc(h*sizeof(float),16);
  // main loop
  for( x=0; x<w0; x++ ) {
    // compute target orientation bins for entire column - very fast
    gradQuantize(O+x*h,M+x*h,O0,O1,M0,M1,nb,h0,sInv2,nOrients,full,softBin>=0);

    if( softBin<0 && softBin%2==0 ) {
      // no interpolation w.r.t. either orienation or spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) for(y=0; y<h0;) { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else if( softBin%2==0 || bin==1 ) {
      // interpolate w.r.t. orientation only, not spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; H1[O1[y]]+=M1[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) for(y=0; y<h0;) { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else {
      // interpolate using trilinear interpolation
      float ms[4], xyd, xb, yb, xd, yd, init; __m128 _m, _m0, _m1;
      bool hasLf, hasRt; int xb0, yb0;
      if( x==0 ) { init=(0+.5f)*sInv-0.5f; xb=init; }
      hasLf = xb>=0; xb0 = hasLf?(int)xb:-1; hasRt = xb0 < wb-1;
      xd=xb-xb0; xb+=sInv; yb=init; y=0;
      // macros for code conciseness
      #define GHinit yd=yb-yb0; yb+=sInv; H0=H+xb0*hb+yb0; xyd=xd*yd; \
        ms[0]=1-xd-yd+xyd; ms[1]=yd-xyd; ms[2]=xd-xyd; ms[3]=xyd;
      #define GH(H,ma,mb) H1=H; STRu(*H1,ADD(LDu(*H1),MUL(ma,mb)));
      // leading rows, no top bin
      for( ; y<bin/2; y++ ) {
        yb0=-1; GHinit;
        if(hasLf) { H0[O0[y]+1]+=ms[1]*M0[y]; H0[O1[y]+1]+=ms[1]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb+1]+=ms[3]*M0[y]; H0[O1[y]+hb+1]+=ms[3]*M1[y]; }
      }
      // main rows, has top and bottom bins, use SSE for minor speedup
      if( softBin<0 ) for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit; _m0=SET(M0[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]); GH(H0+O0[y],_m,_m0); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]); GH(H0+O0[y]+hb,_m,_m0); }
      } else for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit;
        _m0=SET(M0[y]); _m1=SET(M1[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]);
          GH(H0+O0[y],_m,_m0); GH(H0+O1[y],_m,_m1); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]);
          GH(H0+O0[y]+hb,_m,_m0); GH(H0+O1[y]+hb,_m,_m1); }
      }
      // final rows, no bottom bin
      for( ; y<h0; y++ ) {
        yb0 = (int) yb; GHinit;
        if(hasLf) { H0[O0[y]]+=ms[0]*M0[y]; H0[O1[y]]+=ms[0]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb]+=ms[2]*M0[y]; H0[O1[y]+hb]+=ms[2]*M1[y]; }
      }
      #undef GHinit
      #undef GH
    }
  }
  alFree(O0); alFree(O1); alFree(M0); alFree(M1);
  // normalize boundary bins which only get 7/8 of weight of interior bins
  if( softBin%2!=0 ) for( int o=0; o<nOrients; o++ ) {
    x=0; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=0; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    x=wb-1; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=hb-1; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
  }
}

/******************************************************************************/

// HOG helper: compute 2x2 block normalization values (padded by 1 pixel)
float* CDetectorCrosstalk::hogNormMatrix( float *H, int nOrients, int hb, int wb, int bin ) {
  float *N, *N1, *n; int o, x, y, dx, dy, hb1=hb+1, wb1=wb+1;
  float eps = 1e-4f/4/bin/bin/bin/bin; // precise backward equality
  N = (float*)calloc(hb1*wb1,sizeof(float)); N1=N+hb1+1;
  for( o=0; o<nOrients; o++ ) for( x=0; x<wb; x++ ) for( y=0; y<hb; y++ )
    N1[x*hb1+y] += H[o*wb*hb+x*hb+y]*H[o*wb*hb+x*hb+y];
  for( x=0; x<wb-1; x++ ) for( y=0; y<hb-1; y++ ) {
    n=N1+x*hb1+y; *n=1/float(sqrt(n[0]+n[1]+n[hb1]+n[hb1+1]+eps)); }
  x=0;     dx= 1; dy= 1; y=0;                  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=0;     dx= 1; dy= 0; for(y=0; y<hb1; y++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=0;     dx= 1; dy=-1; y=hb1-1;              N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy= 1; y=0;                  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy= 0; for( y=0; y<hb1; y++) N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy=-1; y=hb1-1;              N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  y=0;     dx= 0; dy= 1; for(x=0; x<wb1; x++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  y=hb1-1; dx= 0; dy=-1; for(x=0; x<wb1; x++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  return N;
}

// HOG helper: compute HOG or FHOG channels
void CDetectorCrosstalk::hogChannels( float *H, const float *R, const float *N,
  int hb, int wb, int nOrients, float clip, int type )
{
  #define GETT(blk) t=R1[y]*N1[y-(blk)]; if(t>clip) t=clip; c++;
  const float r=.2357f; int o, x, y, c; float t;
  const int nb=wb*hb, nbo=nOrients*nb, hb1=hb+1;
  for( o=0; o<nOrients; o++ ) for( x=0; x<wb; x++ ) {
    const float *R1=R+o*nb+x*hb, *N1=N+x*hb1+hb1+1;
    float *H1 = (type<=1) ? (H+o*nb+x*hb) : (H+x*hb);
    if( type==0) for( y=0; y<hb; y++ ) {
      // store each orientation and normalization (nOrients*4 channels)
      c=-1; GETT(0); H1[c*nbo+y]=t; GETT(1); H1[c*nbo+y]=t;
      GETT(hb1); H1[c*nbo+y]=t; GETT(hb1+1); H1[c*nbo+y]=t;
    } else if( type==1 ) for( y=0; y<hb; y++ ) {
      // sum across all normalizations (nOrients channels)
      c=-1; GETT(0); H1[y]+=t*.5f; GETT(1); H1[y]+=t*.5f;
      GETT(hb1); H1[y]+=t*.5f; GETT(hb1+1); H1[y]+=t*.5f;
    } else if( type==2 ) for( y=0; y<hb; y++ ) {
      // sum across all orientations (4 channels)
      c=-1; GETT(0); H1[c*nb+y]+=t*r; GETT(1); H1[c*nb+y]+=t*r;
      GETT(hb1); H1[c*nb+y]+=t*r; GETT(hb1+1); H1[c*nb+y]+=t*r;
    }
  }
  #undef GETT
}

// compute HOG features
void CDetectorCrosstalk::hog( float *M, float *O, float *H, int h, int w, int binSize,
  int nOrients, int softBin, bool full, float clip )
{
  float *N, *R; const int hb=h/binSize, wb=w/binSize, nb=hb*wb;
  // compute unnormalized gradient histograms
  R = (float*)calloc(wb*hb*nOrients, sizeof(float));
  gradHist( M, O, R, h, w, binSize, nOrients, softBin, full );
  // compute block normalization values
  N = hogNormMatrix( R, nOrients, hb, wb, binSize );
  // perform four normalizations per spatial block
  hogChannels( H, R, N, hb, wb, nOrients, clip, 0 );
  free(N); free(R);
}

// compute FHOG features
void CDetectorCrosstalk::fhog( float *M, float *O, float *H, int h, int w, int binSize,
  int nOrients, int softBin, float clip )
{
  const int hb=h/binSize, wb=w/binSize, nb=hb*wb, nbo=nb*nOrients;
  float *N, *R1, *R2; int o, x;
  // compute unnormalized constrast sensitive histograms
  R1 = (float*)calloc(wb*hb*nOrients*2, sizeof(float));
  gradHist( M, O, R1, h, w, binSize, nOrients*2, softBin, true );
  // compute unnormalized contrast insensitive histograms
  R2 = (float*)calloc(wb*hb*nOrients, sizeof(float));
  for( o=0; o<nOrients; o++ ) for( x=0; x<nb; x++ )
    R2[o*nb+x] = R1[o*nb+x]+R1[(o+nOrients)*nb+x];
  // compute block normalization values
  N = hogNormMatrix( R2, nOrients, hb, wb, binSize );
  // normalized histograms and texture channels
  hogChannels( H+nbo*0, R1, N, hb, wb, nOrients*2, clip, 1 );
  hogChannels( H+nbo*2, R2, N, hb, wb, nOrients*1, clip, 1 );
  hogChannels( H+nbo*3, R1, N, hb, wb, nOrients*2, clip, 2 );
  free(N); mxFree(R1); free(R2);
}

// Create [hxwxd] mxArray array, initialize to 0 if c=true
mxArray* CDetectorCrosstalk::mxCreateMatrix3( int h, int w, int d, mxClassID id, bool c, void **I ){
  const size_t dims[3]={(size_t)h, (size_t)w, (size_t)d}, n=h*w*d; int b; mxArray* M;
  if( id==mxINT32_CLASS ) b=sizeof(int);
  else if( id==mxDOUBLE_CLASS ) b=sizeof(double);
  else if( id==mxSINGLE_CLASS ) b=sizeof(float);
  else mexErrMsgTxt("Unknown mxClassID.");
  *I = c ? mxCalloc(n,b) : mxMalloc(n*b);
  M = mxCreateNumericMatrix(0,0,id,mxREAL);
  mxSetData(M,*I); mxSetDimensions(M,dims,3); return M;
}

// Check inputs and outputs to mex, retrieve first input I
void CDetectorCrosstalk::checkArgs( int nl, mxArray *pl[], int nr, const mxArray *pr[], int nl0,
  int nl1, int nr0, int nr1, int *h, int *w, int *d, mxClassID id, void **I )
{
  const size_t *dims; int nDims;
  if( nl<nl0 || nl>nl1 ) mexErrMsgTxt("Incorrect number of outputs.");
  if( nr<nr0 || nr>nr1 ) mexErrMsgTxt("Incorrect number of inputs.");
  nDims = (int)mxGetNumberOfDimensions(pr[0]); 
  dims = mxGetDimensions(pr[0]);
  *h=(int)dims[0]; *w=(int)dims[1]; *d=(nDims==2) ? 1 : (int)dims[2]; 
  *I = mxGetPr(pr[0]);
  if( nDims!=2 && nDims!=3 ) mexErrMsgTxt("I must be a 2D or 3D array.");
  if( mxGetClassID(pr[0])!=id ) mexErrMsgTxt("I has incorrect type.");
}

// [Gx,Gy] = grad2(I) - see gradient2.m
void CDetectorCrosstalk::mGrad2( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
  int h, w, d; float *I, *Gx, *Gy;
  checkArgs(nl,pl,nr,pr,1,2,1,1,&h,&w,&d,mxSINGLE_CLASS,(void**)&I);
  if(h<2 || w<2) mexErrMsgTxt("I must be at least 2x2.");
  pl[0]= mxCreateMatrix3( h, w, d, mxSINGLE_CLASS, 0, (void**) &Gx );
  pl[1]= mxCreateMatrix3( h, w, d, mxSINGLE_CLASS, 0, (void**) &Gy );
  grad2( I, Gx, Gy, h, w, d );
}

// [M,O] = gradMag( I, channel, full ) - see gradientMag.m
void CDetectorCrosstalk::mGradMag( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
  int h, w, d, c, full; float *I, *M, *O=0;
  checkArgs(nl,pl,nr,pr,1,2,3,3,&h,&w,&d,mxSINGLE_CLASS,(void**)&I);
  if(h<2 || w<2) mexErrMsgTxt("I must be at least 2x2.");
  c = (int) mxGetScalar(pr[1]); full = (int) mxGetScalar(pr[2]);
  if( c>0 && c<=d ) { I += h*w*(c-1); d=1; }
  pl[0] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&M);
  if(nl>=2) pl[1] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&O);
  gradMag(I, M, O, h, w, d, full>0 );
}

// gradMagNorm( M, S, norm ) - operates on M - see gradientMag.m
void CDetectorCrosstalk::mGradMagNorm( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
  int h, w, d; float *M, *S, norm;
  checkArgs(nl,pl,nr,pr,0,0,3,3,&h,&w,&d,mxSINGLE_CLASS,(void**)&M);
  if( mxGetM(pr[1])!=h || mxGetN(pr[1])!=w || d!=1 ||
    mxGetClassID(pr[1])!=mxSINGLE_CLASS ) mexErrMsgTxt("M or S is bad.");
  S = (float*) mxGetPr(pr[1]); norm = (float) mxGetScalar(pr[2]);
  gradMagNorm(M,S,h,w,norm);
}

// H=gradHist(M,O,[...]) - see gradientHist.m
void CDetectorCrosstalk::mGradHist( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {
  int h, w, d, hb, wb, nChns, binSize, nOrients, softBin, useHog;
  bool full; float *M, *O, *H, clipHog;
  checkArgs(nl,pl,nr,pr,1,3,2,8,&h,&w,&d,mxSINGLE_CLASS,(void**)&M);
  O = (float*) mxGetPr(pr[1]);
  if( mxGetM(pr[1])!=h || mxGetN(pr[1])!=w || d!=1 ||
    mxGetClassID(pr[1])!=mxSINGLE_CLASS ) mexErrMsgTxt("M or O is bad.");
  binSize  = (nr>=3) ? (int)   mxGetScalar(pr[2])    : 8;
  nOrients = (nr>=4) ? (int)   mxGetScalar(pr[3])    : 9;
  softBin  = (nr>=5) ? (int)   mxGetScalar(pr[4])    : 1;
  useHog   = (nr>=6) ? (int)   mxGetScalar(pr[5])    : 0;
  clipHog  = (nr>=7) ? (float) mxGetScalar(pr[6])    : 0.2f;
  full     = (nr>=8) ? (bool) (mxGetScalar(pr[7])>0) : false;
  hb = h/binSize; wb = w/binSize;
  nChns = useHog== 0 ? nOrients : (useHog==1 ? nOrients*4 : nOrients*3+5);
  pl[0] = mxCreateMatrix3(hb,wb,nChns,mxSINGLE_CLASS,1,(void**)&H);
  if( nOrients==0 ) return;
  if( useHog==0 ) {
    gradHist( M, O, H, h, w, binSize, nOrients, softBin, full );
  } else if(useHog==1) {
    hog( M, O, H, h, w, binSize, nOrients, softBin, full, clipHog );
  } else {
    fhog( M, O, H, h, w, binSize, nOrients, softBin, clipHog );
  }
}

// Crosstalk UTIL FUNCTIONs
template<class T> void 
localimshow(mxArray* I, T type, const char* windowName)
{
	// error가 나는 두가지 경우 1. 인풋이미지의 포인터가 잘못돼었을때, 2. dimension이 잘 들어오지 않았을 경우
	mxClassID I_id = mxGetClassID(I);
	T* ptr = (T*)mxGetData(I);

	const size_t* ns = mxGetDimensions(I);
	int w = (int) ns[1];
	int h = (int) ns[0];
	int nChannel = (int) ns[2];

	Mat matI;
	if(nChannel==1)
	{
		matI = Mat::zeros(h,w, CV_8UC1);
	}
	else if(nChannel==3)
	{
		matI = Mat::zeros(h,w, CV_8UC3);      
	}
	else
	{
		return;
	}
	int pixelOffset = 0;
	for (int ch = 0; ch < nChannel; ch++)
	{            		
		for (int j = 0; j < w; j++)
		{ 
			for (int i = 0; i < h; i++, pixelOffset++)
			{
				matI.row(i).col(j).data[nChannel - 1 - ch] = (uchar)ptr[pixelOffset];	
			}
		}
	}
	cv::imshow(windowName,matI);
	waitKey(0);
}

void CDetectorCrosstalk::imshow(mxArray* I, const char* windowName)
{
	// error가 나는 두가지 경우 1. 인풋이미지의 포인터가 잘못돼었을때, 2. dimension이 잘 들어오지 않았을 경우
	mxClassID I_id = mxGetClassID(I);

	if( I_id == mxSINGLE_CLASS)
	{
		// float
		float type = 0.0;
		localimshow(I, type, windowName);
	}
	else if( I_id == mxDOUBLE_CLASS)
	{
		// double
		double type = 0.0;
		localimshow(I, type, windowName);
	}	
	else if( I_id == mxUINT8_CLASS)
	{
		// unsigned char
		unsigned char type = 0;
		localimshow(I, type, windowName);
	}	
	else
	{
		mexErrMsgTxt("Unsupported image type.");
	}
}

template<class T> cv::Mat 
localMxArray2CvMat(mxArray* I, T type)
{
	// error가 나는 두가지 경우 1. 인풋이미지의 포인터가 잘못돼었을때, 2. dimension이 잘 들어오지 않았을 경우
	mxClassID I_id = mxGetClassID(I);
	T* ptr = (T*)mxGetData(I);

	const size_t* ns = mxGetDimensions(I);
	int w = (int) ns[1];
	int h = (int) ns[0];
	int nChannel = (int) ns[2];

	Mat matI;
	if(nChannel==1)
	{
		if(I_id == mxSINGLE_CLASS)
			matI = Mat::zeros(h,w, CV_32FC1);
		else if(I_id == mxDOUBLE_CLASS)
			matI = Mat::zeros(h,w, CV_64FC1);
		else if(I_id == mxUINT8_CLASS)
			matI = Mat::zeros(h,w, CV_8UC1);
		else
			mexErrMsgTxt("Unsupported image type.");
	}
	else if(nChannel==3)
	{
		if(I_id == mxSINGLE_CLASS)
			matI = Mat::zeros(h,w, CV_64FC3);
		else if(I_id == mxDOUBLE_CLASS)
			matI = Mat::zeros(h,w, CV_64FC3);
		else if(I_id == mxUINT8_CLASS)
			matI = Mat::zeros(h,w, CV_8UC3);
		else
			mexErrMsgTxt("Unsupported image type.");   
	}
	int pixelOffset = 0;
	for (int ch = 0; ch < nChannel; ch++)
	{            		
		for (int j = 0; j < w; j++)
		{ 
			for (int i = 0; i < h; i++, pixelOffset++)
			{
				matI.row(i).col(j).data[nChannel - 1 - ch] = (uchar)ptr[pixelOffset];	// 이게 지금 float인데 UINT8로 선언되어있음.....
			}
		}
	}
	return matI;
}

cv::Mat 
CDetectorCrosstalk::MxArray2CvMat(mxArray* I)
{
	// error가 나는 두가지 경우 1. 인풋이미지의 포인터가 잘못돼었을때, 2. dimension이 잘 들어오지 않았을 경우
	mxClassID I_id = mxGetClassID(I);
	cv::Mat MatI;
	if( I_id == mxSINGLE_CLASS)
	{
		// float
		float type = 0.0f;
		MatI = localMxArray2CvMat(I, type);
	}
	else if( I_id == mxDOUBLE_CLASS)
	{
		// double
		double type = 0.0;
		MatI = localMxArray2CvMat(I, type);
	}	
	else if( I_id == mxUINT8_CLASS)
	{
		// unsigned char
		unsigned char type = 0;
		MatI = localMxArray2CvMat(I, type);
	}	
	else
	{
		mexErrMsgTxt("Unsupported image type.");
	}
	return MatI;
}

mxArray* CDetectorCrosstalk::CvMat2MxArray(cv::Mat* I, int flag) {
	mxArray *mxArray_testimg= NULL;
	uchar* mxArray_testimg_pr = NULL;
	int h = (*I).rows;
	int w = (*I).cols;
	int imgSrc_Channels = (*I).channels();
	if(imgSrc_Channels == 1) { // gray image
		size_t dims[2] = {(size_t)h, (size_t)w};
		mxArray_testimg = mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL);
		mxArray_testimg_pr = (uchar *)mxGetData(mxArray_testimg);
		//mxArray_testimg->data = I->data;
		memcpy(mxArray_testimg->data,I->data,h*w);
	}
	else if (imgSrc_Channels == 3) {// 3-channel image
		size_t dims[3] = {(size_t)h, (size_t)w, (size_t)imgSrc_Channels};
		mxArray_testimg = mxCreateNumericArray(imgSrc_Channels, dims, mxUINT8_CLASS, mxREAL);

		vector<Mat> channels(3);
		split(*I, channels);
		transpose(channels[0], channels[0]);
		transpose(channels[1], channels[1]);
		transpose(channels[2], channels[2]);

		mxArray_testimg_pr = (uchar *)mxGetData(mxArray_testimg);
		memcpy((char *)(mxArray_testimg->data),channels[2].data,h*w);
		memcpy((char *)(mxArray_testimg->data)+h*w,channels[1].data,h*w);
		memcpy((char *)(mxArray_testimg->data)+h*w*2,channels[0].data,h*w);
	}
	mxSetData(mxArray_testimg, mxArray_testimg_pr);
	return mxArray_testimg;
}

mxArray* CDetectorCrosstalk::MxArrayResize(mxArray *BeforeResize,int crop_height, int crop_width)
{
	const size_t *sz1 = mxGetDimensions(BeforeResize);
	mxArray* AfterResize = mxCreateNumericArray(sz1[2],sz1,mxSINGLE_CLASS,mxREAL);
	float *BeforeResizeData = (float *)mxGetData(BeforeResize);
	float *AfterResizeData = (float *)mxGetData(AfterResize);
	int pixelOffset = 0;
	for(int ch = 0; ch < (int)sz1[2];ch++) {
		for(int col = 0; col < crop_width;col++) {
			for(int row = 0; row < crop_height;row++,pixelOffset++) {
				//AfterResizeData[pixelOffset] = BeforeResizeData[pixelOffset];
				AfterResizeData[pixelOffset] = BeforeResizeData[ch*sz1[0]*sz1[1] + col*sz1[0] + row];
			}
		}
	}
	//mxSetData(AfterResize,AfterResizeData);
	return AfterResize;
}

 