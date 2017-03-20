#include <limits>
#include <assert.h>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "HungarianMethod.h"
#include "haanju_math.hpp"
#include "haanju_visualize.hpp"
#include "SCMTTracker.h"

// TEMPORAL
#include "haanju_3D.hpp"


namespace hj
{

/////////////////////////////////////////////////////////////////////////
// CDetectedObject MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////

CDetectedObject::CDetectedObject()
	: id(0)
	, location(0.0, 0.0, 0.0)
	, height(0.0)
	, bMatchedWithTracker(false)
	, bCoveredByOtherDetection(false)
{
}


CDetectedObject::~CDetectedObject()
{
	for (size_t vecIdx = 0; vecIdx < vecvecTrackedFeatures.size(); vecIdx++)
	{
		vecvecTrackedFeatures[vecIdx].clear();
	}
	vecvecTrackedFeatures.clear();
	boxes.clear();

	if (!patchGray.empty()) { patchGray.release(); }
	if (!patchRGB.empty())  { patchRGB.release(); }
}

/////////////////////////////////////////////////////////////////////////
// CTracker2D MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////

CTracker2D::CTracker2D()
	: id(0)
	, timeStart(0)
	, timeEnd(0)
	, timeLastUpdate(0)
	, duration(0)
	, numStatic(0)
	, confidence(0.0)
	, lastPosition(0.0, 0.0, 0.0)
	, estimatedBox(0.0, 0.0, 0.0, 0.0)
	, height(0.0)
{
}


CTracker2D::~CTracker2D()
{
	boxes.clear();
	heads.clear();
	featurePoints.clear();
	trackedPoints.clear();
}


/////////////////////////////////////////////////////////////////////////
// CSCMTTracker MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////
CSCMTTracker::CSCMTTracker()
	: bInit_(false)
	, bVisualizeResult_(false)
	, strVisWindowName_("")
{
}


CSCMTTracker::~CSCMTTracker()
{
	Finalize();
}


/************************************************************************
 Method Name: Initialize
 Description:
	- Initialize the multiple target tracker with parameters
 Input Arguments:
	- nCamID:      index of the camera
	- stParams:    parameters for tracking algorithm
	- ptCalibInfo: calibration information
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Initialize(int _nCamID, stParamTrack2D &_stParams)
{
	// check duplicated initialization
	if (bInit_) { Finalize(); }

	stParam_ = _stParams;
	stParam_.dImageRescaleRecover = 1.0 / stParam_.dImageRescale;

	nCamID_ = _nCamID;
	nCurrentFrameIdx_ = 0;
	trackingResult_.camID = nCamID_;
	trackingResult_.frameIdx = 0;
	trackingResult_.object2DInfos.clear();

	// calibration related
	pCalibrationInfo_ = _stParams.pCalibrationInfo;	
	nInputWidth_      = pCalibrationInfo_->cCamModel.width();
	nInputHeight_     = pCalibrationInfo_->cCamModel.height();
	dDefaultBottomZ_  = 0.0;

	// detection related
	vecDetectedObjects_.clear();

	// tracker related
	nNewTrackerID_ = 0;
	listCTracker2D_.clear();
	queueActiveTracker2D_.clear();

	// input related
	sizeBufferImage_ = cv::Size(
		(int)((double)nInputWidth_  * stParam_.dImageRescale),
		(int)((double)nInputHeight_ * stParam_.dImageRescale));
	matGrayImage_ = cv::Mat(nInputHeight_, nInputWidth_, CV_8UC1);
	cImageBuffer_.set(stParam_.nBackTrackingLength);	

	// feature tracking related		
	featureDetector_          = cv::AgastFeatureDetector::create();
	matFeatureExtractionMask_ = cv::Mat(sizeBufferImage_, CV_8UC1, cv::Scalar(0));

	// visualization related
	bVisualizeResult_ = stParam_.bVisualize;
	strVisWindowName_ = "2D Tracking with cam " + std::to_string(nCamID_);
	if (bVisualizeResult_)
	{
		vecColors_ = hj::GenerateColors(400);
	}

	// TEMPORAL
	bFirstDraw_ = true;

	// initialization flag
	bInit_ = true;
}


/************************************************************************
 Method Name: Finalize
 Description:
	- terminate the class with memory clean up
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Finalize(void)
{	
	if (!bInit_) { return; }

	/* detection related */
	this->vecDetectedObjects_.clear();

	/* tracker related */
	for (std::list<CTracker2D>::iterator trackerIter = listCTracker2D_.begin();
		trackerIter != listCTracker2D_.end();
		trackerIter++)
	{
		(*trackerIter).boxes.clear();
		(*trackerIter).featurePoints.clear();
		(*trackerIter).trackedPoints.clear();
	}
	listCTracker2D_.clear();
	queueActiveTracker2D_.clear();

	/* input related */
	cImageBuffer_.clear();
	if (!matGrayImage_.empty()) { matGrayImage_.release(); }
	if (!matTrackingResult_.empty()) { matTrackingResult_.release(); }

	/* matching related */
	arrMatchingCost_.clear();

	/* result related */
	trackingResult_.object2DInfos.clear();

	/* visualize related */
	if (bVisualizeResult_) { cv::destroyWindow(strVisWindowName_); }

	/* initialization flag */
	bInit_ = false;
}


/************************************************************************
 Method Name: Track
 Description:
	- Run the tracking algorithm on the current input frame
 Input Arguments:
	- vecInputDetections: input detections of the current frame
	- curFrame: current input frame image
	- frameIdx: current frame index
 Return Values:
	- CTrack2DResult: tracking result of the current frame
************************************************************************/
CTrack2DResult& CSCMTTracker::Track(std::vector<CDetection> vecInputDetections, cv::Mat curFrame, int frameIdx)
{
	assert(bInit_ && curFrame.rows == matGrayImage_.rows && curFrame.cols == matGrayImage_.cols);	
	
	nCurrentFrameIdx_ = frameIdx;

	/* buffering */
	cv::cvtColor(curFrame, matGrayImage_, CV_BGR2GRAY);
	cImageBuffer_.insert_resize(matGrayImage_, sizeBufferImage_);
	if (!matTrackingResult_.empty()) { matTrackingResult_.release(); }
	matTrackingResult_ = curFrame.clone();

	/* input pre-processing */
	Track2D_GenerateDetectedObjects(vecInputDetections, vecDetectedObjects_);

	/* bi-directional tracking */
	Track2D_BackwardTracking(vecDetectedObjects_);
	Track2D_ForwardTracking(queueActiveTracker2D_);	
	Track2D_MatchingAndUpdating(vecDetectedObjects_, queueActiveTracker2D_);

	/* result packaging */
	Track2D_ResultPackaging();

	/* visualize */
	if (bVisualizeResult_) { VisualizeResult(); }

	return this->trackingResult_;
}


/************************************************************************
 Method Name: Track2D_GenerateDetectedObjects
 Description:
	- Generate set of 'CDetectedObject' with input detection.
 Input Arguments:
	- vecDetections     : input detections
	- vecDetectedObjects: (output) pre-processed detections
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Track2D_GenerateDetectedObjects(
	const std::vector<CDetection> &vecDetections,
	std::vector<CDetectedObject> &vecDetectedObjects)
{
	/* reset 'vecDetectedObjecs' for usage of the current frame */
	size_t detectionID = 0;
	vecDetectedObjects_.clear();
	vecDetectedObjects_.reserve(vecDetections.size());

	double  estimatedHeight = 0.0;
	Point3D estimatedLocation;

	for (size_t detectionIdx = 0; detectionIdx < vecDetections.size(); detectionIdx++)
	{
		/* validation with box height (or scale) */
		estimatedHeight = hj::EstimateBoxHeight(
			vecDetections[detectionIdx].box, 
			*pCalibrationInfo_,
			dDefaultBottomZ_, 
			&estimatedLocation);
		if (stParam_.dDetectionMaxHeight < estimatedHeight || stParam_.dDetectionMinHeight > estimatedHeight)
		{ 
			continue; 
		}

		/* generate detection information */
		CDetectedObject curDetection;
		curDetection.id             = (unsigned int)detectionID++;
		curDetection.detection      = vecDetections[detectionIdx];
		curDetection.detection.box *= stParam_.dImageRescale;
		curDetection.location       = estimatedLocation;
		curDetection.height         = estimatedHeight;
		curDetection.bMatchedWithTracker      = false;
		curDetection.bCoveredByOtherDetection = false;
		curDetection.vecvecTrackedFeatures.reserve(stParam_.nBackTrackingLength);
		curDetection.boxes.reserve(stParam_.nBackTrackingLength);
		curDetection.boxes.push_back(curDetection.detection.box);

		vecDetectedObjects.push_back(curDetection);
	}
}


/************************************************************************
 Method Name: Track2D_BackwardTracking
 Description:
	- Track the input detections in a backward direction with input frame
	  buffers.
 Input Arguments: 
	- vecDetectedObjects: (in/out) pre-processed detections with backward tracking results
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Track2D_BackwardTracking(std::vector<CDetectedObject> &vecDetectedObjects)
{
	//---------------------------------------------------
	// BACKWARD FEATURE TRACKING
	//---------------------------------------------------
	for (size_t dIdx = 0; dIdx < vecDetectedObjects.size(); dIdx++)
	{			
		hj::Rect rectRescaledDetectionBox = vecDetectedObjects[dIdx].detection.box.scale(stParam_.dImageRescale);
		hj::Rect rectEstimatedBox;

		std::vector<cv::Point2f> vecInputFeatures, vecTrackedFeatures;
		std::vector<int>         vecInlireIndices;

		// feature extraction
		FeatureExtraction(vecDetectedObjects[dIdx].boxes.front(), *cImageBuffer_.rbegin(), vecInputFeatures);
		vecDetectedObjects[dIdx].vecvecTrackedFeatures.push_back(vecInputFeatures);

		/* feature tracking */
		for (hj::CMatFIFOBuffer::reverse_iterator bufferIter = cImageBuffer_.rbegin();
			bufferIter != cImageBuffer_.rend() - 1;
			bufferIter++)
		{
			if (!FeatureTracking(rectRescaledDetectionBox,
				                 *bufferIter,
				                 *(bufferIter + 1),
				                 vecInputFeatures,
				                 vecTrackedFeatures,
				                 vecInlireIndices,
				                 rectEstimatedBox))
			{
				break;
			}

			vecDetectedObjects[dIdx].boxes.push_back(rectEstimatedBox.scale(stParam_.dImageRescaleRecover));

			/* save inliers */
			vecInputFeatures.clear();
			vecInputFeatures.reserve(vecInlireIndices.size());
			for (size_t indexIdx = 0; indexIdx < vecInlireIndices.size(); indexIdx++)
			{
				vecInputFeatures.push_back(vecTrackedFeatures[vecInlireIndices[indexIdx]]);
			}
			vecDetectedObjects[dIdx].vecvecTrackedFeatures.push_back(vecInputFeatures);
		}
	}

	//---------------------------------------------------
	// CHECK OVERLAP
	//---------------------------------------------------
	for (int detect1Idx = 0; detect1Idx < vecDetectedObjects.size(); detect1Idx++)
	{
		if (vecDetectedObjects[detect1Idx].bCoveredByOtherDetection) { continue; }
		for (int detect2Idx = detect1Idx + 1; detect2Idx < vecDetectedObjects.size(); detect2Idx++)
		{
			if (vecDetectedObjects[detect1Idx].detection.box.overlap(vecDetectedObjects[detect2Idx].detection.box))
			{
				vecDetectedObjects[detect1Idx].bCoveredByOtherDetection = true;
				break;
			}
		}
	}
}


/************************************************************************
 Method Name: Track2D_ForwardTracking
 Description:
	- Estimate trackers positions at the current frame. Estimated positions
	  are inserted at the end of box array of each tracker.
 Input Arguments:
	- queueTrackers: (in/out) trackers of which we want to estimate the location.
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Track2D_ForwardTracking(std::deque<CTracker2D*> &queueTrackers)
{
	if (2 > cImageBuffer_.num_elements())
	{ 
		// there is no frame to track
		return; 
	}

	//---------------------------------------------------
	// FORWARD FEATURE TRACKING
	//---------------------------------------------------
	for (size_t trackIdx = 0; trackIdx < queueTrackers.size(); trackIdx++)
	{
		std::vector<cv::Point2f> vecTrackedFeatures;
		std::vector<int>         vecInlierIndices;
		hj::Rect estimatedBox;
		if (!FeatureTracking(queueTrackers[trackIdx]->boxes.back(),
			                 *(cImageBuffer_.rbegin() + 1),
			                 *cImageBuffer_.rbegin(),
			                 queueTrackers[trackIdx]->featurePoints,
			                 vecTrackedFeatures,
			                 vecInlierIndices,
			                 estimatedBox))
		{
			continue;
		}

		queueTrackers[trackIdx]->trackedPoints.clear();
		if (stParam_.nMinNumFeatures > vecInlierIndices.size()) { continue; }
		for (size_t featureIdx = 0; featureIdx < vecInlierIndices.size(); featureIdx++)
		{
			queueTrackers[trackIdx]->trackedPoints.push_back(vecTrackedFeatures[vecInlierIndices[featureIdx]]);
		}		
		queueTrackers[trackIdx]->estimatedBox = estimatedBox;		
	}
}


/************************************************************************
 Method Name: Track2D_MatchingAndUpdating
 Description:
	- Do matching between detections and trackers. And then, update trackers
 Input Arguments:
	- vecDetectedObjects: (pre-processed) detections
	- queueTrackers: trackers
 Return Values:
	- none
************************************************************************/
void CSCMTTracker::Track2D_MatchingAndUpdating(
	const std::vector<CDetectedObject> &vecDetectedObjects,
	std::deque<CTracker2D*> &queueTrackers)
{
	/////////////////////////////////////////////////////////////////////////////
	// CALCULATE MATCHING COSTS
	/////////////////////////////////////////////////////////////////////////////

	// matching cost matrix: default score is an infinite
	arrMatchingCost_.clear();
	arrMatchingCost_.resize(vecDetectedObjects.size() * queueTrackers.size(), std::numeric_limits<float>::infinity());
	//std::vector<float> arrMatchingCost_(vecDetectedObjects.size() * queueTrackers.size(), std::numeric_limits<float>::infinity());

	// to determine occlusion
	std::vector<std::deque<int>> featuresInDetectionBox(vecDetectedObjects.size());

	//---------------------------------------------------
	// COST WITH BI-DIRECTIONAL TRACKING
	//---------------------------------------------------
	for (size_t trackIdx = 0; trackIdx < queueTrackers.size(); trackIdx++)
	{
		CTracker2D *curTracker = queueTrackers[trackIdx];


		for (size_t detectIdx = 0, costPos = trackIdx; 
			detectIdx < vecDetectedObjects.size(); 
			detectIdx++, costPos += queueTrackers.size())
		{
			// validate with backward tracking result
			if (!curTracker->estimatedBox.overlap(vecDetectedObjects[detectIdx].detection.box)) { continue; }

			// count feature points inside the detection box
			for (int featureIdx = 0; featureIdx < curTracker->trackedPoints.size(); featureIdx++)
			{
				if (!vecDetectedObjects[detectIdx].detection.box.contain(curTracker->trackedPoints[featureIdx])) { continue; }
				featuresInDetectionBox[detectIdx].push_back((int)trackIdx);
			}		

			// determine the possible longest comparison interval
			size_t lengthForCompare = std::min((size_t)stParam_.nBackTrackingLength, 
				std::min(curTracker->boxes.size(), vecDetectedObjects[detectIdx].boxes.size()));			
			
			// croll tracker boxes for comparison (reverse ordering)
			size_t numBoxCopy = lengthForCompare - 1;
			std::vector<hj::Rect> vecTrackerBoxes;
			vecTrackerBoxes.reserve(lengthForCompare);
			vecTrackerBoxes.push_back(curTracker->estimatedBox);
			if (0 < numBoxCopy)
			{
				vecTrackerBoxes.insert(vecTrackerBoxes.begin() + 1,
					curTracker->boxes.rbegin(), 
					curTracker->boxes.rbegin() + lengthForCompare - 1);
			}

			double boxCost = 0.0;
			Rect detectionBox, trackerBox;
			for (size_t boxIdx = 0; boxIdx < lengthForCompare; boxIdx++)
			{
				detectionBox = vecDetectedObjects[detectIdx].boxes[boxIdx];
				trackerBox   = vecTrackerBoxes[boxIdx];

				if (!detectionBox.overlap(trackerBox) // rejection criterion
					|| stParam_.dMaxBoxDistance < detectionBox.distance(trackerBox)
					|| stParam_.dMinBoxOverlapRatio > detectionBox.overlappedArea(trackerBox) / std::min(detectionBox.area(), trackerBox.area())
					|| stParam_.dMaxBoxCenterDiffRatio * std::max(detectionBox.w, trackerBox.w) < (detectionBox.center() - trackerBox.center()).norm_L2())
				{
					boxCost = std::numeric_limits<double>::infinity();
					break;
				}
				boxCost += BoxMatchingCost(trackerBox, detectionBox);
			}
			if (std::numeric_limits<double>::infinity() == boxCost) { continue; }
			boxCost /= (double)lengthForCompare;

			arrMatchingCost_[costPos] = (float)boxCost;
		}
	}

	//---------------------------------------------------
	// OCCLUSION HANDLING
	//---------------------------------------------------	
	// If a detection box contains feature points from more than one tracker, we examine whether there exists a
	// dominant tracker or not. If there is no dominant tracker, that means the ownership of the detection is 
	// not clear, we set the scores between that detection and trackers to infinite. This yields termination of
	// all related trackers.

	int numFeatureFromMajorTracker = 0,
		numFeatureFromCurrentTracker = 0,
		majorTrackerIdx = 0,
		currentTrackerIdx = 0;
	for (size_t detectIdx = 0, costPos = 0; detectIdx < vecDetectedObjects.size(); detectIdx++, costPos += queueTrackers.size())
	{
		if (0 == featuresInDetectionBox[detectIdx].size()) { continue; }

		// find dominant tracker of the detection
		numFeatureFromMajorTracker = numFeatureFromCurrentTracker = 0;
		majorTrackerIdx = featuresInDetectionBox[detectIdx].front();
		currentTrackerIdx = featuresInDetectionBox[detectIdx].front();
		for (int featureIdx = 0; featureIdx < featuresInDetectionBox[detectIdx].size(); featureIdx++)
		{
			if (currentTrackerIdx == featuresInDetectionBox[detectIdx][featureIdx])
			{
				// we assume that the same tracker indices in 'featuresInDetectionBox' are grouped together
				numFeatureFromCurrentTracker++;
				continue;
			}

			if (numFeatureFromCurrentTracker > numFeatureFromMajorTracker)
			{
				majorTrackerIdx = currentTrackerIdx;
				numFeatureFromMajorTracker = numFeatureFromCurrentTracker;
			}
			currentTrackerIdx = featuresInDetectionBox[detectIdx][featureIdx];
			numFeatureFromCurrentTracker = 0;
		}

		// case 1: only one tracker has its features points in the detecion box
		if (featuresInDetectionBox[detectIdx].front() == currentTrackerIdx)
		{
			continue;
		}

		// case 2: there is a domninant tracker among related trackers
		if (numFeatureFromMajorTracker > featuresInDetectionBox[detectIdx].size() * stParam_.dMinOpticalFlowMajorityRatio)
		{
			continue;
		}

		// case 3: more than one trackers are related and there is no dominant tracker
		for (size_t infCostPos = costPos; infCostPos < costPos + queueTrackers.size(); infCostPos++)
		{
			arrMatchingCost_[infCostPos] = std::numeric_limits<float>::infinity();
		}
	}


	//---------------------------------------------------
	// INFINITE HANDLING
	//---------------------------------------------------
	// To ensure a proper operation of our Hungarian implementation, we convert infinite to the finite value
	// that is little bit (=100.0f) greater than the maximum finite cost in the original cost function.
	float maxCost = -1000.0f;
	for (int costIdx = 0; costIdx < arrMatchingCost_.size(); costIdx++)
	{
		if (!_finitef(arrMatchingCost_[costIdx])) { continue; }
		if (maxCost < arrMatchingCost_[costIdx]) { maxCost = arrMatchingCost_[costIdx]; }
	}
	maxCost = maxCost + 100.0f;
	for (int costIdx = 0; costIdx < arrMatchingCost_.size(); costIdx++)
	{
		if (_finitef(arrMatchingCost_[costIdx])) { continue; }
		arrMatchingCost_[costIdx] = maxCost;
	}

	/////////////////////////////////////////////////////////////////////////////
	// MATCHING
	/////////////////////////////////////////////////////////////////////////////
	trackingResult_.object2DInfos.clear();
	size_t numDetection = this->vecDetectedObjects_.size();
	CHungarianMethod cHungarianMatcher;
	cHungarianMatcher.Initialize(arrMatchingCost_, (unsigned int)numDetection, (unsigned int)this->queueActiveTracker2D_.size());
	stMatchInfo *curMatchInfo = cHungarianMatcher.Match();
	for (size_t matchIdx = 0; matchIdx < curMatchInfo->rows.size(); matchIdx++)
	{
		if (maxCost == curMatchInfo->matchCosts[matchIdx]) { continue; }
		CDetectedObject *curDetection = &this->vecDetectedObjects_[curMatchInfo->rows[matchIdx]];
		CTracker2D *curTracker = this->queueActiveTracker2D_[curMatchInfo->cols[matchIdx]];

		//---------------------------------------------------
		// MATCHING VALIDATION
		//---------------------------------------------------
		// distance in 3D space
		if ((curDetection->location - curTracker->lastPosition).norm_L2() > stParam_.dMaxDetectionDistance) { continue; }
		// height difference
		if (std::abs(curDetection->height - curTracker->height) > stParam_.dMaxHeightDifference) { continue; }
		// max tracklet length
		if (curTracker->duration >= (unsigned int)stParam_.nMaxTrackletLength) { continue; }
		double curConfidence = 1.0;

		//---------------------------------------------------
		// TRACKER UPDATE
		//---------------------------------------------------		
		curTracker->timeEnd        = this->nCurrentFrameIdx_;
		curTracker->timeLastUpdate = this->nCurrentFrameIdx_;
		curTracker->duration       = curTracker->timeEnd - curTracker->timeStart + 1;
		curTracker->numStatic      = 0;		
		curTracker->confidence     = curConfidence;
		curTracker->lastPosition   = curDetection->location;
		curTracker->height         = curDetection->height;
		curTracker->boxes.push_back(curDetection->detection.box);

		curDetection->bMatchedWithTracker = true;
		//queueNewActiveTrackers.push_back(curTracker);

		////---------------------------------------------------
		//// RESULT PACKAGING
		////---------------------------------------------------
		//CObject2DInfo objectInfo;
		//ResultWithTracker(curTracker, objectInfo);
		//trackingResult_.object2DInfos.push_back(objectInfo);

		// update features with detection (after result packaging)
		curTracker->featurePoints = curDetection->vecvecTrackedFeatures.front();
		curTracker->trackedPoints.clear();
	}
	cHungarianMatcher.Finalize();

	/////////////////////////////////////////////////////////////////////////////
	// TRACKER GENERATION
	/////////////////////////////////////////////////////////////////////////////
	for (std::vector<CDetectedObject>::iterator detectionIter = this->vecDetectedObjects_.begin();
		detectionIter != this->vecDetectedObjects_.end();
		detectionIter++)
	{
		if ((*detectionIter).bMatchedWithTracker) { continue; }

		CTracker2D newTracker;
		newTracker.id = this->nNewTrackerID_++;
		newTracker.timeStart = this->nCurrentFrameIdx_;
		newTracker.timeEnd = this->nCurrentFrameIdx_;
		newTracker.timeLastUpdate = this->nCurrentFrameIdx_;
		newTracker.duration = 1;
		newTracker.numStatic = 0;
		newTracker.boxes.push_back((*detectionIter).detection.box);
		newTracker.featurePoints = (*detectionIter).vecvecTrackedFeatures.front();
		newTracker.trackedPoints.clear();
		newTracker.confidence = 1.0;
		newTracker.lastPosition = (*detectionIter).location;
		newTracker.height = (*detectionIter).height;

		// generate tracklet instance
		this->listCTracker2D_.push_back(newTracker);
		//queueNewActiveTrackers.push_back(&this->listCTracker2D_.back());

		////---------------------------------------------------
		//// RESULT PACKAGING
		////---------------------------------------------------
		//CObject2DInfo objectInfo;
		//ResultWithTracker(&newTracker, objectInfo);
		//trackingResult_.object2DInfos.push_back(objectInfo);
	}

	/////////////////////////////////////////////////////////////////////////////
	// TRACKER TERMINATION
	/////////////////////////////////////////////////////////////////////////////
	std::deque<CTracker2D*> queueNewActiveTrackers;
	for (std::list<CTracker2D>::iterator trackerIter = listCTracker2D_.begin();
		trackerIter != listCTracker2D_.end();
		/*trackerIter++*/)
	{
		if ((*trackerIter).timeLastUpdate == nCurrentFrameIdx_)
		{
			// keep
			queueNewActiveTrackers.push_back(&(*trackerIter));
			trackerIter++;
			continue;
		}
		// termination
		(*trackerIter).boxes.clear();
		(*trackerIter).featurePoints.clear();
		(*trackerIter).trackedPoints.clear();
		trackerIter = this->listCTracker2D_.erase(trackerIter);
	}
	queueActiveTracker2D_ = queueNewActiveTrackers;

	//// matching cost
	//int cost_pos = 0;
	//if (!this->trackingResult_.matMatchingCost.empty()) { trackingResult_.matMatchingCost.release(); }
	//trackingResult_.matMatchingCost = 
	//	cv::Mat((int)vecDetectedObjects_.size(), (int)trackingResult_.vecTrackerRects.size(), CV_32F);
	//for (int detectionIdx = 0; detectionIdx < vecDetectedObjects_.size(); detectionIdx++)
	//{
	//	for (int trackIdx = 0; trackIdx < trackingResult_.vecTrackerRects.size(); trackIdx++)
	//	{
	//		this->trackingResult_.matMatchingCost.at<float>(detectionIdx, trackIdx) = arrMatchingCost_[cost_pos];
	//		cost_pos++;
	//	}
	//}
}


/************************************************************************
 Method Name: ResultWithTracker
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
void CSCMTTracker::Track2D_ResultPackaging()
{
	time_t timePackaging = clock();
	CObject2DInfo objectInfo;

	trackingResult_.camID     = nCamID_;
	trackingResult_.frameIdx  = nCurrentFrameIdx_;
	trackingResult_.timeStamp = (unsigned int)timePackaging;
	trackingResult_.object2DInfos.clear();	
	for (size_t trackerIdx = 0; trackerIdx < queueActiveTracker2D_.size(); trackerIdx++)
	{
		ResultWithTracker(queueActiveTracker2D_[trackerIdx], objectInfo);
		trackingResult_.object2DInfos.push_back(objectInfo);
	}

	int cost_pos = 0;
	if (!this->trackingResult_.matMatchingCost.empty()) { trackingResult_.matMatchingCost.release(); }
	trackingResult_.matMatchingCost =
		cv::Mat((int)vecDetectedObjects_.size(), (int)trackingResult_.vecTrackerRects.size(), CV_32F);
	for (int detectionIdx = 0; detectionIdx < vecDetectedObjects_.size(); detectionIdx++)
	{
		for (int trackIdx = 0; trackIdx < trackingResult_.vecTrackerRects.size(); trackIdx++)
		{
			this->trackingResult_.matMatchingCost.at<float>(detectionIdx, trackIdx) = arrMatchingCost_[cost_pos];
			cost_pos++;
		}
	}
}


/************************************************************************
 Method Name: FeatureExtraction
 Description:
	- Tracks the input feature points of the input frame in the target frame.
 Input Arguments:
	- boundingBox      : bounding box of the target at the input frame
	- inputImage       : image containing the input feature points.
	- targetImage      : target image of feature point tracking.
	- vecInputFeatures : input feature points.
	- vecOutputFeatures: points of tracking result. Actually it is an output.
	- vecFeatureInlierIndex: index of features that are inliers of estimated motion.
	- trackingResult   : (output) estimated box at the target frame.
 Return Values:
	- bool: true = success / false = failure
************************************************************************/
bool CSCMTTracker::FeatureExtraction(
	const hj::Rect inputBox, 
	const cv::Mat  inputImage, 
	std::vector<cv::Point2f> &vecFeaturePoints)
{
	vecFeaturePoints.clear();

	//---------------------------------------------------
	// EXTRACT FEATURE POINTS
	//---------------------------------------------------
	std::vector<cv::KeyPoint> newKeypoints;
	cv::Rect rectROI = inputBox.cropWithSize(inputImage.cols, inputImage.rows).cv();
	matFeatureExtractionMask_(rectROI) = cv::Scalar(255); // masking with the bounding box
	featureDetector_->detect(inputImage, newKeypoints, matFeatureExtractionMask_);
	matFeatureExtractionMask_(rectROI) = cv::Scalar(0);   // restore the mask image

	if (stParam_.nMinNumFeatures > newKeypoints.size())
	{
		// it is impossible to track the target because there are insufficient number of feature points
		return false;
	}

	//---------------------------------------------------
	// EXTRACT SELECTION
	//---------------------------------------------------	
	std::random_shuffle(newKeypoints.begin(), newKeypoints.end());
	for (int pointIdx = 0; pointIdx < std::min((int)newKeypoints.size(), stParam_.nMaxNumFeatures); pointIdx++)
	{
		vecFeaturePoints.push_back(newKeypoints[pointIdx].pt);
	}

	return true;
}


/************************************************************************
 Method Name: FeatureTracking
 Description:
	- Tracks the input feature points of the input frame in the target frame.
 Input Arguments:
	- boundingBox      : bounding box of the target at the input frame
	- inputImage       : image containing the input feature points.
	- targetImage      : target image of feature point tracking.
	- vecInputFeatures : input feature points.
	- vecOutputFeatures: points of tracking result. Actually it is an output.
	- vecFeatureInlierIndex: index of features that are inliers of estimated motion.
	- trackingResult   : (output) estimated box at the target frame.
 Return Values:
	- bool: true = success / false = failure
************************************************************************/
bool CSCMTTracker::FeatureTracking(
	const hj::Rect inputBox,
	const cv::Mat  inputImage, 
	const cv::Mat  targetImage,
	std::vector<cv::Point2f> &vecInputFeatures, 
	std::vector<cv::Point2f> &vecOutputFeatures,
	std::vector<int>         &vecFeatureInlierIndex,
	hj::Rect &trackingResult)
{
	if (0 == vecInputFeatures.size()) { return false; }

	trackingResult = hj::Rect(0.0, 0.0, 0.0, 0.0);

	//---------------------------------------------------
	// CONVERT TO GRAY SCALE IMAGES
	//---------------------------------------------------
	cv::Mat currImage, nextImage;
	
	if (1 == inputImage.channels()) { currImage = inputImage; }
	else { cv::cvtColor(inputImage, currImage, CV_BGR2GRAY); }
	
	if (1 == targetImage.channels()) { nextImage = targetImage; }
	else { cv::cvtColor(targetImage, nextImage, CV_BGR2GRAY); }

	//---------------------------------------------------
	// EXTRACT FEATURE POINTS
	//---------------------------------------------------
	if (vecInputFeatures.empty())
	{
		if (!FeatureExtraction(inputBox, currImage, vecInputFeatures))
		{
			return false;
		}
	}

	//---------------------------------------------------
	// FEATURE TRACKING
	//---------------------------------------------------
	std::vector<uchar> vecFeatureStatus;
	vecOutputFeatures.clear();

	cv::Mat vecErrors;
	cv::calcOpticalFlowPyrLK(
		currImage,
		nextImage,
		vecInputFeatures,
		vecOutputFeatures,
		vecFeatureStatus,
		cv::noArray(),
		//vecErrors,
		cv::Size((int)(inputBox.w * stParam_.dFeatureTrackWindowSizeRatio),
		         (int)(inputBox.w * stParam_.dFeatureTrackWindowSizeRatio)));

	//---------------------------------------------------
	// BOX ESTIMATION
	//---------------------------------------------------	
	hj::Rect newRect = LocalSearchKLT(inputBox, vecInputFeatures, vecOutputFeatures, vecFeatureInlierIndex);
	if (stParam_.nMinNumFeatures > vecFeatureInlierIndex.size())
	{ 
		// tracking failure because of insufficient number of tracked feature points
		return false; 
	}
	else
	{
		trackingResult = newRect;
	}

	return true;
}


/************************************************************************
 Method Name: FindInlierFeatures
 Description:
	- Find inlier feature points
 Input Arguments:
	- vecInputFeatures:
	- vecOutputFeatures:
	- vecPointStatus:
 Return Values:
	- portion of inlier
************************************************************************/
std::vector<cv::Point2f> CSCMTTracker::FindInlierFeatures(
	std::vector<cv::Point2f> *vecInputFeatures,
	std::vector<cv::Point2f> *vecOutputFeatures,
	std::vector<unsigned char> *vecPointStatus)
{
	size_t numTrackedFeatures = 0;
	// find center of disparity
	cv::Point2f disparityCenter(0.0f, 0.0f);
	std::vector<cv::Point2f> vecDisparity;
	std::vector<size_t> vecInlierIndex;
	for (size_t pointIdx = 0; pointIdx < vecPointStatus->size(); pointIdx++)
	{
		if (!(*vecPointStatus)[pointIdx]) { continue; }
		vecDisparity.push_back((*vecOutputFeatures)[pointIdx] - (*vecInputFeatures)[pointIdx]);
		disparityCenter += vecDisparity.back();
		vecInlierIndex.push_back(pointIdx);
		numTrackedFeatures++;
	}
	disparityCenter = (1 / (float)numTrackedFeatures) * disparityCenter;

	// find distribution of disparity norm
	float norm;
	float normAverage = 0.0f;
	float normSqauredAverage = 0.0f;
	float normStd = 0.0;
	std::vector<float> vecNorm;
	for (size_t pointIdx = 0; pointIdx < vecDisparity.size(); pointIdx++)
	{
		norm = (float)cv::norm(vecDisparity[pointIdx] - disparityCenter);
		vecNorm.push_back(norm);
		normAverage += norm;
		normSqauredAverage += norm * norm;
	}
	normAverage /= (float)numTrackedFeatures;
	normSqauredAverage /= (float)numTrackedFeatures;
	normStd = sqrtf(((float)numTrackedFeatures / ((float)numTrackedFeatures - 1)) * (normSqauredAverage - (normAverage * normAverage)));

	std::vector<cv::Point2f> vecInlierFeatures;
	for (size_t pointIdx = 0; pointIdx < vecNorm.size(); pointIdx++)
	{
		if (abs(vecNorm[pointIdx] - normAverage) > 1 * normStd) { continue; }
		vecInlierFeatures.push_back((*vecOutputFeatures)[vecInlierIndex[pointIdx]]);
	}

	return vecInlierFeatures;
}


/************************************************************************
 Method Name: LocalSearchKLT
 Description:
	- estimate current box location with feature tracking result
 Input Arguments:
	- preFeatures  : feature positions at the previous frame
	- curFeatures  : feature positions at the current frame
	- featureStatus: (output) indicates inlier features
 Return Values:
	- Rect: estimated box
************************************************************************/
#define PSN_LOCAL_SEARCH_PORTION_INLIER (false)
#define PSN_LOCAL_SEARCH_MINIMUM_MOVEMENT (0.1)
#define PSN_LOCAL_SEARCH_NEIGHBOR_WINDOW_SIZE_RATIO (0.2)
Rect CSCMTTracker::LocalSearchKLT(
	Rect preBox,
	std::vector<cv::Point2f> &preFeatures,
	std::vector<cv::Point2f> &curFeatures,
	std::vector<int> &inlierFeatureIndex)
{
	size_t numFeatures = preFeatures.size();
	size_t numMovingFeatures = 0;
	inlierFeatureIndex.clear();
	inlierFeatureIndex.reserve(numFeatures);

	// find disparity of moving features
	std::vector<Point2D> vecMovingVector;
	std::vector<int> vecMovingFeatuerIdx;
	std::vector<double> vecDx;
	std::vector<double> vecDy;
	vecMovingVector.reserve(numFeatures);
	vecMovingFeatuerIdx.reserve(numFeatures);
	vecDx.reserve(numFeatures);
	vecDy.reserve(numFeatures);
	Point2D movingVector;
	double disparity = 0.0;
	for (int featureIdx = 0; featureIdx < (int)numFeatures; featureIdx++)
	{
		movingVector = curFeatures[featureIdx] - preFeatures[featureIdx];
		disparity = movingVector.norm_L2();
		if (disparity < PSN_LOCAL_SEARCH_MINIMUM_MOVEMENT * stParam_.dImageRescale) { continue; }

		vecMovingVector.push_back(movingVector);
		vecMovingFeatuerIdx.push_back(featureIdx);
		vecDx.push_back(movingVector.x);
		vecDy.push_back(movingVector.y);

		numMovingFeatures++;
	}

	// check static movement
	if (numMovingFeatures < numFeatures * 0.5)
	{ 
		for (int featureIdx = 0; featureIdx < (int)numFeatures; featureIdx++)
		{
			if (preBox.contain(curFeatures[featureIdx]))
			{
				inlierFeatureIndex.push_back(featureIdx);
			}
		}
		return preBox; 
	}

	std::sort(vecDx.begin(), vecDx.end());
	std::sort(vecDy.begin(), vecDy.end());

	// estimate major disparity
	double windowSize = preBox.w * PSN_LOCAL_SEARCH_NEIGHBOR_WINDOW_SIZE_RATIO * stParam_.dImageRescale;
	size_t maxNeighborX = 0, maxNeighborY = 0;
	Point2D estimatedDisparity;
	for (size_t disparityIdx = 0; disparityIdx < numMovingFeatures; disparityIdx++)
	{
		size_t numNeighborX = 0;
		size_t numNeighborY = 0;
		// find neighbors in each axis
		for (size_t compIdx = 0; compIdx < numMovingFeatures; compIdx++)
		{
			if (std::abs(vecDx[disparityIdx] - vecDx[compIdx]) < windowSize) { numNeighborX++; } // X		
			if (std::abs(vecDy[disparityIdx] - vecDy[compIdx]) < windowSize) { numNeighborY++; } // Y
		}
		// disparity in X axis
		if (maxNeighborX < numNeighborX)
		{
			estimatedDisparity.x = vecDx[disparityIdx];
			maxNeighborX = numNeighborX;
		}
		// disparity in Y axis
		if (maxNeighborY < numNeighborY)
		{
			estimatedDisparity.y = vecDy[disparityIdx];
			maxNeighborY = numNeighborY;
		}
	}

	// find inliers
	for (int vectorIdx = 0; vectorIdx < (int)numMovingFeatures; vectorIdx++)
	{
		if ((vecMovingVector[vectorIdx] - estimatedDisparity).norm_L2() < windowSize)
		{
			inlierFeatureIndex.push_back(vecMovingFeatuerIdx[vectorIdx]);
		}
	}

	// estimate box
	Rect estimatedBox = preBox;
	estimatedBox.x += estimatedDisparity.x;
	estimatedBox.y += estimatedDisparity.y;

	return estimatedBox;
}


/************************************************************************
 Method Name: BoxDistance
 Description:
	- Calculate the distance between two boxes
 Input Arguments:
	- box1: the first box
	- box2: the second box
 Return Values:
	- double: distance between two boxes
************************************************************************/
double CSCMTTracker::BoxMatchingCost(Rect &box1, Rect &box2)
{
	double nominator = (box1.center() - box2.center()).norm_L2();
	double denominator = (box1.w + box2.w) / 2.0;
	double boxDistance = (nominator * nominator) / (denominator * denominator);

	return boxDistance;

	//double probability = std::exp(-boxDistance);
	//double cost = -std::numeric_limits<double>::infinity();
	//if(1.0 > probability)
	//{
	//	cost = boxDistance + std::log(1.0 - probability); 
	//}
	//return cost;
}


/************************************************************************
 Method Name: GetTrackingConfidence
 Description:
	- Calculate the confidence of tracking by the number of features lay in the box
 Input Arguments:
	- box: target position
	- vecTrackedFeatures: tracked features
 Return Values:
	- tracking confidence
************************************************************************/
double CSCMTTracker::GetTrackingConfidence(Rect &box, std::vector<cv::Point2f> &vecTrackedFeatures)
{
	double numFeaturesInBox = 0.0;
	for (std::vector<cv::Point2f>::iterator featureIter = vecTrackedFeatures.begin();
		featureIter != vecTrackedFeatures.end();
		featureIter++)
	{
		if (box.contain(*featureIter))
		{
			numFeaturesInBox++;
		}
	}

	return numFeaturesInBox / (double)vecTrackedFeatures.size();
}


/************************************************************************
 Method Name: ResultWithTracker
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
void CSCMTTracker::ResultWithTracker(CTracker2D *curTracker, CObject2DInfo &outObjectInfo)
{
	Rect curBox = curTracker->boxes.back();
	curBox.x *= (float)stParam_.dImageRescaleRecover;
	curBox.y *= (float)stParam_.dImageRescaleRecover;
	curBox.w *= (float)stParam_.dImageRescaleRecover;
	curBox.h *= (float)stParam_.dImageRescaleRecover;
	outObjectInfo.prevFeatures = curTracker->featurePoints;
	outObjectInfo.currFeatures = curTracker->trackedPoints;
	outObjectInfo.id = curTracker->id;
	outObjectInfo.box = curBox;
	outObjectInfo.score = 0;
	for (int pointIdx = 0; pointIdx < outObjectInfo.prevFeatures.size(); pointIdx++)
	{
		outObjectInfo.prevFeatures[pointIdx].x *= (float)stParam_.dImageRescaleRecover;
		outObjectInfo.prevFeatures[pointIdx].y *= (float)stParam_.dImageRescaleRecover;
		if (pointIdx >= outObjectInfo.currFeatures.size()) { continue; }
		outObjectInfo.currFeatures[pointIdx].x *= (float)stParam_.dImageRescaleRecover;
		outObjectInfo.currFeatures[pointIdx].y *= (float)stParam_.dImageRescaleRecover;
	}
}


/************************************************************************
 Method Name: VisualizeResult
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
void CSCMTTracker::VisualizeResult()
{
	/* detections */
	for (int detIdx = 0; detIdx < vecDetectedObjects_.size(); detIdx++)
	{
		cv::rectangle(
			matTrackingResult_, 
			vecDetectedObjects_[detIdx].detection.box.cv(),
			cv::Scalar(255, 255, 255), 
			1);
	}

	/* tracklets */
	for (int objIdx = 0; objIdx < trackingResult_.object2DInfos.size(); objIdx++)
	{
		CObject2DInfo *curObject = &trackingResult_.object2DInfos[objIdx];

		// feature points
		for (int pointIdx = 0; pointIdx < curObject->prevFeatures.size(); pointIdx++)
		{
			if (pointIdx < curObject->currFeatures.size())
			{
				cv::circle(
					matTrackingResult_,
					curObject->prevFeatures[pointIdx],
					1, cv::Scalar(0, 255, 0), 1);
				cv::line(
					matTrackingResult_,
					curObject->prevFeatures[pointIdx],
					curObject->currFeatures[pointIdx],
					cv::Scalar(255, 255, 255), 1);
				cv::circle(
					matTrackingResult_,
					curObject->currFeatures[pointIdx],
					1, cv::Scalar(0, 255, 0), 1);
			}
			else
			{
				cv::circle(
					matTrackingResult_,
					curObject->prevFeatures[pointIdx],
					1, cv::Scalar(0, 0, 255), 1);
			}
		}

		// tracklet box
		hj::DrawBoxWithID(matTrackingResult_, curObject->box, curObject->id, 0, 0, &vecColors_);
	}

	cv::namedWindow(strVisWindowName_);
	cv::moveWindow(strVisWindowName_, 10 + nCamID_*640, 10);
	cv::imshow(strVisWindowName_, matTrackingResult_);
	cv::waitKey(1);
	matTrackingResult_.release();
}

}

//()()
//('')HAANJU.YOO

