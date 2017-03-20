/**************************************************************************
* Title        : CSCTracker
* Author       : Haanju Yoo
* Initial Date : 2014.03.01 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  :
*	Single camera multiple target tracker
**************************************************************************/

#pragma once

#include <vector>
#include <list>
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/video/tracking.hpp"
#include "types.hpp"

namespace hj
{

/////////////////////////////////////////////////////////////////////////
// INPUT DETECTION
/////////////////////////////////////////////////////////////////////////
class CDetectedObject
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CDetectedObject();
	~CDetectedObject();

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	CDetection   detection;
	Point3D      location;
	double       height;
	bool         bMatchedWithTracker;
	bool         bCoveredByOtherDetection;

	/* back tracking related */
	std::vector<std::vector<cv::Point2f>> vecvecTrackedFeatures; // current -> past order
	std::vector<Rect> boxes; // current -> past order

	/* appearance related */
	cv::Mat patchGray;
	cv::Mat patchRGB;
	
};

/////////////////////////////////////////////////////////////////////////
// SINGLE TARGET TRACKER
/////////////////////////////////////////////////////////////////////////
class CTracker2D
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CTracker2D();
	~CTracker2D();

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	unsigned int timeStart;
	unsigned int timeEnd;
	unsigned int timeLastUpdate;
	unsigned int duration;
	unsigned int numStatic;
	double confidence;
	std::deque<hj::Rect> boxes;
	std::deque<hj::Rect> heads;
	std::vector<cv::Point2f> featurePoints;
	std::vector<cv::Point2f> trackedPoints;
	Point3D lastPosition;
	hj::Rect estimatedBox;
	double height;
};


/////////////////////////////////////////////////////////////////////////
// MULTI-TARGET TRACKER
/////////////////////////////////////////////////////////////////////////
class CSCMTTracker
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CSCMTTracker();
	~CSCMTTracker();

	void Initialize(int _nCamID, stParamTrack2D &_stParams);
	void Finalize(void);
	CTrack2DResult& Track(
		std::vector<CDetection> vecInputDetections, 
		cv::Mat curFrame, 
		int frameIdx);

private:
	/* MAIN OPERATIONS */	
	void Track2D_GenerateDetectedObjects(
		const std::vector<CDetection> &vecDetections, 
		std::vector<CDetectedObject> &vecDetectedObjects);
	void Track2D_BackwardTracking(std::vector<CDetectedObject> &vecDetectedObjects);
	void Track2D_ForwardTracking(std::deque<CTracker2D*> &queueTrackers);
	void Track2D_MatchingAndUpdating(
		const std::vector<CDetectedObject> &vecDetectedObjects, 
		std::deque<CTracker2D*> &queueTrackers);
	void Track2D_ResultPackaging();

	/* TRACKING RELATED */
	bool FeatureExtraction(
		const hj::Rect inputBox, 
		const cv::Mat inputImage, 
		std::vector<cv::Point2f> &vecFeaturePoints);
	bool FeatureTracking(
		const hj::Rect inputBox, 
		const cv::Mat inputImage, 
		const cv::Mat targetImage, 
		std::vector<cv::Point2f> &vecInputFeatures, 
		std::vector<cv::Point2f> &vecOutputFeatures, 
		std::vector<int> &vecFeatureInlierIndex, 
		hj::Rect &trackingResult);
	std::vector<cv::Point2f> FindInlierFeatures(
		std::vector<cv::Point2f> *vecInputFeatures, 
		std::vector<cv::Point2f> *vecOutputFeatures, 
		std::vector<unsigned char> *vecPointStatus);
	Rect LocalSearchKLT(
		Rect preBox, 
		std::vector<cv::Point2f> &preFeatures, 
		std::vector<cv::Point2f> &curFeatures, 
		std::vector<int> &inlierFeatureIndex);	
	static double BoxMatchingCost(Rect &box1, Rect &box2);
	static double GetTrackingConfidence(Rect &box, std::vector<cv::Point2f> &vecTrackedFeatures);	

	/* ETC */
	void ResultWithTracker(CTracker2D *curTracker, CObject2DInfo &outObjectInfo);

	/* FOR DEBUGGING */
	void VisualizeResult();

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	bool           bInit_;
	stParamTrack2D stParam_;
	unsigned int   nCamID_;
	unsigned int   nCurrentFrameIdx_;

	/* calibration related */
	CCalibrationInfo *pCalibrationInfo_;
	unsigned int     nInputWidth_;
	unsigned int     nInputHeight_;
	double           dDefaultBottomZ_;

	/* input related */
	std::vector<CDetectedObject> vecDetectedObjects_;
	CMatFIFOBuffer cImageBuffer_;
	cv::Size sizeBufferImage_;
	cv::Mat  matGrayImage_;
	cv::Mat  matResizedGrayImage_;	
	
	/* traker related */
	unsigned int            nNewTrackerID_;
	std::list<CTracker2D>   listCTracker2D_;
	std::deque<CTracker2D*> queueActiveTracker2D_;

	/* matching related */
	std::vector<float> arrMatchingCost_;

	/* feature tracking related */
	cv::Ptr<cv::AgastFeatureDetector> featureDetector_;
	cv::Mat matFeatureExtractionMask_;

	/* result related */
	CTrack2DResult trackingResult_;

	/* visualization related */
	bool        bVisualizeResult_;
	cv::Mat     matTrackingResult_;
	std::string strVisWindowName_;
	std::vector<cv::Scalar> vecColors_;

	// TEMPORAL
	bool bFirstDraw_;
};

}

//()()
//('')HAANJU.YOO
