/*************************************************************************************
NOTES:
- There are some inaccurate calculations in 3D point reconstruction and vision prior
at branching stage for conveniency for computing.
**************************************************************************************/
#include <fstream>
#include <iostream>
#include <string>
#include <omp.h>
#include "opencv2/highgui/highgui.hpp"
#include "HungarianMethod.h"
#include "SGSmoother.h"
#include "haanju_math.hpp"
#include "haanju_misc.hpp"
#include "haanju_visualize.hpp"
#include "Associator3D.h"

#define NOMINMAX
#include <process.h>
#include <windows.h>


namespace hj
{

/////////////////////////////////////////////////////////////////////////
// PARAMETERS (unit: mm)
/////////////////////////////////////////////////////////////////////////

#define PRINT_TRACKS (0) // 0: no print / 1: print
#define MAX_NUM_SOLVING_THREADS (2000)

// smoothing related
#define MIN_SMOOTHING_LENGTH (SGS_DEFAULT_SPAN / 2)
std::vector<Qset> precomputedQsets;

/////////////////////////////////////////////////////////////////////////
// THREAD RELATED
/////////////////////////////////////////////////////////////////////////
struct stSolvingThreadParams
{
	int                threadID;
	HypothesisSet      *outBranchHypotheses;
	stParamAssociate3D *params;
	TrackSet           *tracks;
	TrackSet           *initialSolutionTracks;
};
static stSolvingThreadParams gVecSolvingThreadParams[MAX_NUM_SOLVING_THREADS];

unsigned int __stdcall AssociationSolvingWork(void *data)
{
	stSolvingThreadParams *pParams = (stSolvingThreadParams*)data;

	hj::printf_debug("  >> solving thread is started, id: %d\n", pParams->threadID);

	hj::CAssociator3D::Hypothesis_BranchHypotheses(
		*pParams->outBranchHypotheses,
		pParams->params,
		pParams->tracks,
		pParams->initialSolutionTracks);

	hj::printf_debug("  >> solving thread is terminated, id: %d\n", pParams->threadID);

	return 0;
}

/////////////////////////////////////////////////////////////////////////
// LOCAL FUNCTIONS
/////////////////////////////////////////////////////////////////////////

// compartors for sorting
bool hjTrackIDAscendComparator(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	return track1->id < track2->id;
}

bool hjTrackAscendentCostComparator(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	return track1->costTotal < track2->costTotal;
}

bool hjTrackGPDescendComparator(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	return track1->GTProb > track2->GTProb;
}

bool hjTrackGPDescendAndDurationAscendComparator(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	// higher GTProb, shorter duration!
	if (track1->GTProb > track2->GTProb) { return true; }
	else if (track1->GTProb == track2->GTProb)
	{
		if (track1->duration < track2->duration) { return true; }
	}
	return false;
}

bool hjTrackGTPandLLDescend(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	// higher GTProb, lower cost!
	if (track1->GTProb > track2->GTProb) { return true; }
	else if (track1->GTProb == track2->GTProb)
	{
		if (track1->loglikelihood > track2->loglikelihood) { return true; }
	}
	return false;
}

bool hjTrackGPandLengthComparator(
	const CTrack3D *track1, 
	const CTrack3D *track2)
{
	if (track1->tree->bConfirmed && track2->tree->bConfirmed)
	{
		if (track1->GTProb > track2->GTProb) { return true; }
		else if (track1->GTProb < track2->GTProb) { return false; }
		if (track1->duration > track2->duration) { return true; }
		return false;
	}
	else if (track1->tree->bConfirmed && track2->tree->bConfirmed) { return false; }
	return true;
}

bool hjSolutionLogLikelihoodDescendComparator(
	const stGlobalHypothesis &solution1, 
	const stGlobalHypothesis &solution2)
{
	return solution1.logLikelihood > solution2.logLikelihood;
}


/////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////


/************************************************************************
 Method Name: CAssociator3D
 Description:
	- class constructor
 Input Arguments:
	- none
 Return Values:
	- class instance
************************************************************************/
CAssociator3D::CAssociator3D()
	: bInit_(false)	
{
}


/************************************************************************
 Method Name: ~CAssociator3D
 Description:
	- class destructor
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
CAssociator3D::~CAssociator3D()
{
}


/************************************************************************
 Method Name: Initialize
 Description:
	- initialization routine for 3D association
 Input Arguments:
	- datasetPath: the string for a dataset path
	- stCalibInfo: calibration information
 Return Values:
	- none
************************************************************************/
void CAssociator3D::Initialize(stParamAssociate3D &_stParam)
{
	if (bInit_) { return; }
	stParam_ = _stParam;
	nNumCameras_ = stParam_.stViewInfo.nNumCameras;
	bSnapshotReaded_ = false;

	nCurrentFrameIdx_       = 0;
	nNumFramesForProc_      = 0;
	nCountForPenalty_       = 0;
	dCurrentProcessingTime_ = 0.0;
	dCurrentSolvingTime_    = 0.0;

	// image buffer related
	vecMatCurrentFrames_.resize(nNumCameras_);

	// 2D tracklet related
	vecTracklet2DSet_.resize(nNumCameras_);
	nNumTotalActive2DTracklet_ = 0;

	// 3D track related
	nNewTrackID_ = 0;
	nNewTreeID_  = 0;
	nLastPrintedDeferredResultFrameIdx_ = 0;
	nLastPrintedInstantResultFrameIdx_  = 0;
	bReceiveNewMeasurement_ = false;
	bInitiationPenaltyFree_ = true;

	queueNewSeedTracks_.clear();
	queueActiveTrack_.clear();
	queuePausedTrack_.clear();
	queueTracksInWindow_.clear();
	queueTracksInBestSolution_.clear();

	// optimization related
	queuePrevGlobalHypotheses_.clear();
	queueCurrGlobalHypotheses_.clear();

	// result related
	nNewTargetID_ = 0;
	queuePairTreeIDToTargetID_.clear();

	// visualize related
	bVisualizeResult_ = stParam_.bVisualize;
	strVisWindowName_ = "Association Result";
	if (bVisualizeResult_)
	{
		vecColors_ = hj::GenerateColors(400);
	}
	
	// smoothing related
	for (int windowSize = 1; windowSize <= SGS_DEFAULT_SPAN; windowSize++)
	{
		precomputedQsets.push_back(CSGSmoother::CalculateQ(windowSize));
	}

	// calibration related
	vecMatProjectionSensitivity_.resize(nNumCameras_);
	vecMatDistanceFromBoundary_.resize(nNumCameras_);
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		vecMatDistanceFromBoundary_[camIdx] = 
			stParam_.vecPCalibrationInfo[camIdx]->matDistanceFromBoundary;
		vecMatProjectionSensitivity_[camIdx] = 
			stParam_.vecPCalibrationInfo[camIdx]->matProjectionSensitivity;
	}

	// file path
	//strLogFileName_ = std::string(RESULT_SAVE_PATH) + "/" + std::string(strParameter) + "/" + strTime_ + "_log.txt";
	//strTrackLogFileName_ = std::string(RESULT_SAVE_PATH) + "/" + "hj_tracks_" + strTime_ + "txt";

	// processing time
	queueProcessingTime_.clear();

	// init flag
	bInit_ = true;
}

/************************************************************************
 Method Name: Finalize
 Description:
	- termination routine for 3D association
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CAssociator3D::Finalize(void)
{
	if (!bInit_) { return; }
#ifdef SAVE_SNAPSHOT_
	hj::CreateDirectoryForWindows(std::string(SNAPSHOT_PATH));
	this->SaveSnapshot(SNAPSHOT_PATH);
#endif
	/////////////////////////////////////////////////////////////////////////////
	// PRINT RESULT
	/////////////////////////////////////////////////////////////////////////////

	// tracks in solutions
	std::deque<CTrack3D*> queueResultTracks;
	for (std::deque<CTrack3D*>::iterator trackIter = queueTracksInBestSolution_.begin();
		trackIter != queueTracksInBestSolution_.end();
		trackIter++)
	{
		queueResultTracks.push_back(*trackIter);
	}

	// print results
	//this->PrintResult(strInstantResultFileName_.c_str(), &queueTrackingResult_);	

//	/////////////////////////////////////////////////////////////////////////////
//	// EVALUATION
//	/////////////////////////////////////////////////////////////////////////////
//#ifndef LOAD_SNAPSHOT_	
//	for (int evalIdx = 0; evalIdx < vecEvaluator_->size(); evalIdx++)
//	{
//		// fill the result
//		for (unsigned int timeIdx = (unsigned int)std::max(0, (int)nCurrentFrameIdx_ - (*vecEvaluator_)[evalIdx].second + 1);
//			timeIdx <= nCurrentFrameIdx_;
//			timeIdx++)
//		{
//			(*vecEvaluator_)[evalIdx].first.SetResult(queueTracksInBestSolution_, timeIdx);
//		}
//	}
//#endif

	///////////////////////////////////////////////////////////////////////////////
	//// LOGGING
	///////////////////////////////////////////////////////////////////////////////	
	//std::string strLog = "";
	//for (int fIdx = 0; fIdx < queueProcessingTime_.size(); fIdx++)
	//{
	//	strLog += std::to_string(queueProcessingTime_[fIdx]) + "\n";
	//}
	//hj::printLog(strLogFileName_.c_str(), strLog);

	/////////////////////////////////////////////////////////////////////////////
	// FINALIZE
	/////////////////////////////////////////////////////////////////////////////
	// clean-up camera model
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		vecMatCurrentFrames_[camIdx].release();
		stParam_.vecPCalibrationInfo[camIdx] = NULL;		
		vecMatProjectionSensitivity_[camIdx].release();
		vecMatDistanceFromBoundary_[camIdx].release();
	}
	vecMatCurrentFrames_.clear();
	vecMatProjectionSensitivity_.clear();
	vecMatDistanceFromBoundary_.clear();

	dCurrentProcessingTime_ = 0.0;
	nCurrentFrameIdx_ = 0;

	// 2D tracklet related
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		vecTracklet2DSet_[camIdx].activeTracklets.clear();
		for (std::list<CTracklet2D>::iterator trackletIter = vecTracklet2DSet_[camIdx].tracklets.begin();
			trackletIter != vecTracklet2DSet_[camIdx].tracklets.end();
			trackletIter++)
		{
			trackletIter->rects.clear();
		}
		vecTracklet2DSet_[camIdx].tracklets.clear();
	}
	vecTracklet2DSet_.clear();

	// 3D track related
	nNewTrackID_ = 0;
	nNewTreeID_ = 0;

	/* visualize related */
	if (bVisualizeResult_) { cv::destroyWindow(strVisWindowName_); }
}


/************************************************************************
 Method Name: Run
 Description:
	- main process for 3D association
 Input Arguments:
	-
	-
 Return Values:
	- CTrack3DResult: 3D tracking result from the processing
************************************************************************/
CTrack3DResult CAssociator3D::Run(
	std::vector<CTrack2DResult> &curTrack2DResult,
	std::vector<cv::Mat> curFrames,
	int frameIdx)
{
	assert(bInit_);

#ifdef LOAD_SNAPSHOT_
	if (!bSnapshotReaded_) { bSnapshotReaded_ = this->LoadSnapshot(SNAPSHOT_PATH); }
	if (bSnapshotReaded_ && frameIdx <= (int)nCurrentFrameIdx_)
	{
		CTrack3DResult curResult;
		if (frameIdx < queueTrackingResult_.size()) { curResult = queueTrackingResult_[frameIdx]; }
		return curResult;
	}
#endif
	clock_t timer_start;
	clock_t timer_end;
	timer_start = clock();
	double processingTime;

	/////////////////////////////////////////////////////////////////////////////
	// PRE-PROCESSING
	/////////////////////////////////////////////////////////////////////////////
	// get frames
	nCurrentFrameIdx_ = frameIdx;
	nNumFramesForProc_++;
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++) 
	{ 
		vecMatCurrentFrames_[camIdx] = curFrames[camIdx];
	}
	// enterance/exit penalty
	if (bInitiationPenaltyFree_ && (int)nCountForPenalty_++ > stParam_.nEnterPenaltyFreeLength) 
	{ 
		bInitiationPenaltyFree_ = false; 
	}

	/////////////////////////////////////////////////////////////////////////////
	// 2D TRACKLET
	/////////////////////////////////////////////////////////////////////////////
	Tracklet2D_UpdateTracklets(curTrack2DResult, nCurrentFrameIdx_);

	/////////////////////////////////////////////////////////////////////////////
	// 3D TRACK
	/////////////////////////////////////////////////////////////////////////////	
	Track3D_Management(queueNewSeedTracks_, nCurrentFrameIdx_);

	/////////////////////////////////////////////////////////////////////////////
	// GLOBAL HYPOTHESES
	/////////////////////////////////////////////////////////////////////////////	
	// solve MHT
	Hypothesis_UpdateHypotheses(queuePrevGlobalHypotheses_, &queueNewSeedTracks_);
	Hypothesis_Formation(queueCurrGlobalHypotheses_, &queuePrevGlobalHypotheses_);
	// post-pruning
	Hypothesis_PruningNScanBack(
		nCurrentFrameIdx_, 
		stParam_.nProcWindowSize, 
		&queueTracksInWindow_, 
		&queueCurrGlobalHypotheses_);
	Hypothesis_PruningTrackWithGTP(
		nCurrentFrameIdx_, 
		stParam_.nMaxNumTracksInOptimization, 
		&queueTracksInWindow_, 
		&queuePtActiveTrees_);
	Hypothesis_RefreshHypotheses(queueCurrGlobalHypotheses_);

	/////////////////////////////////////////////////////////////////////////////
	// RESULT PACKING
	/////////////////////////////////////////////////////////////////////////////
	// measuring processing time
	timer_end = clock();
	processingTime = (double)(timer_end - timer_start) / CLOCKS_PER_SEC;
	dCurrentProcessingTime_ = processingTime;

	// packing current tracking result
	if (0 < queueCurrGlobalHypotheses_.size())
	{
		queueTracksInBestSolution_ = queueCurrGlobalHypotheses_.front().selectedTracks;
	}
	else
	{
		queueTracksInBestSolution_.clear();
	}

	// instance result
	CTrack3DResult currentResult = ResultWithTracks(&queueTracksInBestSolution_, frameIdx, processingTime);

	//// ID association
	//if (queueTrackingResult_.size() > 0)
	//{
	//	std::vector<float> arrMatchingCost;
	//	arrMatchingCost.resize(
	//		currentResult.object3DInfos.size() * queueTrackingResult_.back().object3DInfos.size(),
	//		std::numeric_limits<float>::infinity());

	//	int pos = 0;
	//	for (int i = 0; i < (int)currentResult.object3DInfos.size(); i++)
	//	{
	//		for (int j = 0; j < (int)queueTrackingResult_.back().object3DInfos.size(); j++)
	//		{
	//			arrMatchingCost[pos] = (currentResult.object3DInfos[i].recentPoints.back()
	//				- queueTrackingResult_.back().object3DInfos[j].recentPoints.back()).norm_L2();
	//			pos++;
	//		}
	//	}

	//	CHungarianMethod cHungarianMatcher;
	//	cHungarianMatcher.Initialize(arrMatchingCost, 
	//		(unsigned int)currentResult.object3DInfos.size(), 
	//		(unsigned int)queueTrackingResult_.back().object3DInfos.size());
	//	stMatchInfo *curMatchInfo = cHungarianMatcher.Match();
	//	int a = 0;
	//}

	queueTrackingResult_.push_back(currentResult);

	if (bVisualizeResult_) { VisualizeResult(frameIdx); }

	/////////////////////////////////////////////////////////////////////////////
	// WRAP-UP
	/////////////////////////////////////////////////////////////////////////////
	// memory clean-up
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		vecMatCurrentFrames_[camIdx].release(); 
	}

	// hypothesis backup
	queuePrevGlobalHypotheses_ = queueCurrGlobalHypotheses_;
	queueCurrGlobalHypotheses_.clear();

	// processing time
	queueProcessingTime_.push_back(dCurrentProcessingTime_);
	dCurrentSolvingTime_ = 0.0;

#ifdef HJ_PRINT_LOG_
	//// PRINT LOG
	//std::string strLog = "";
	//strLog += std::to_string(queueProcessingTime_.back()) + "\n";
	//hj::printLog(strLogFileName_.c_str(), strLog);
#endif

	return currentResult;
}


/************************************************************************
 Method Name: GetHuman3DBox
 Description:
	- find a corner points of 3D box for a given 3D point
 Input Arguments:
	- ptHeadCenter: 3D point of a center of the head
	- bodyWidth: a width of the box
	- camIdx: an index of the camera for display
 Return Values:
	- std::vector<hj::Point2D>: 3D tracking result from the processing
************************************************************************/
std::vector<hj::Point2D> CAssociator3D::GetHuman3DBox(
	hj::Point3D ptHeadCenter, 
	double bodyWidth, 
	int camIdx)
{
	std::vector<hj::Point2D> resultArray;
	double bodyWidthHalf = bodyWidth / 2.0;

	double offsetX[8] = { +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf };
	double offsetY[8] = { -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf };
	double z[8] = { ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, 0, 0, 0, 0 };

	unsigned int numPointsOnView = 0;
	for (unsigned int vertexIdx = 0; vertexIdx < 8; vertexIdx++)
	{
		hj::Point2D curPoint = this->WorldToImage(hj::Point3D(ptHeadCenter.x + offsetX[vertexIdx], ptHeadCenter.y + offsetY[vertexIdx], z[vertexIdx]), camIdx);
		resultArray.push_back(curPoint);
		if (curPoint.onView(stParam_.vecPCalibrationInfo[camIdx]->cCamModel.width(), stParam_.vecPCalibrationInfo[camIdx]->cCamModel.height())) { numPointsOnView++; }
	}

	// if there is no point can be seen from the view, clear a vector of point
	if (0 == numPointsOnView) { resultArray.clear(); }

	return resultArray;
}


/************************************************************************
 Method Name: WorldToImage
 Description:
	- convert a world point to a point on the image (3D -> 2D)
 Input Arguments:
	- point3D: 3D point
	- camIdx: an index of target camera
 Return Values:
	- hj::Point2D: image coordinate of the point (2D)
************************************************************************/
hj::Point2D CAssociator3D::WorldToImage(hj::Point3D point3D, int camIdx)
{
	return hj::WorldToImage(point3D, stParam_.vecPCalibrationInfo[camIdx]);
}


/************************************************************************
 Method Name: ImageToWorld
 Description:
	- convert an image point to 3D point with a specific height (2D -> 3D)
 Input Arguments:
	- point2D: a point on the image
	- z: height of reconstructed point
	- camIdx: an index of camera which a 2D point came from
 Return Values:
	- hj::Point2D: image coordinate of the point (2D)
************************************************************************/
hj::Point3D CAssociator3D::ImageToWorld(hj::Point2D point2D, double z, int camIdx)
{
	return hj::ImageToWorld(point2D, z, stParam_.vecPCalibrationInfo[camIdx]);;
}


/************************************************************************
 Method Name: CheckVisibility
 Description:m
	- return whether an input point can be seen at a specific camera
 Input Arguments:
	- testPoint: a 3D point
	- camIdx: camera index (for calibration information)
 Return Values:
	- bool: visible(true)/un-visible(false)
************************************************************************/
bool CAssociator3D::CheckVisibility(hj::Point3D testPoint, int camIdx, hj::Point2D *result2DPoint)
{
	hj::Point2D reprojectedPoint = this->WorldToImage(testPoint, camIdx);
	hj::Point2D reprojectedTopPoint = this->WorldToImage(
		hj::Point3D(testPoint.x, testPoint.y, stParam_.dDefaultHeight), camIdx);

	// pad in for detection probability, no pads for just finding the position of reprojected point
	double halfWidth = NULL == result2DPoint ? (reprojectedTopPoint - reprojectedPoint).norm_L2() / 6 : 0.0;

	if (NULL != result2DPoint) { *result2DPoint = reprojectedPoint; }
	if (reprojectedPoint.x < halfWidth 
		|| reprojectedPoint.x >= (double)stParam_.vecPCalibrationInfo[camIdx]->cCamModel.width() - halfWidth
		|| reprojectedPoint.y < halfWidth 
		|| reprojectedPoint.y >= (double)stParam_.vecPCalibrationInfo[camIdx]->cCamModel.height() - halfWidth)
	{
		return false;
	}
	return true;
}

/************************************************************************
 Method Name: CheckTrackletConnectivity
 Description:
	-
 Input Arguments:
	-
	-
 Return Values:
	-
************************************************************************/
bool CAssociator3D::CheckTrackletConnectivity(
	hj::Point3D endPoint, 
	hj::Point3D startPoint, 
	double sensitivity1, 
	double sensitivity2, 
	int timeGap)
{
	if (timeGap > 1) { return true; }
	double norm2 = (endPoint - startPoint).norm_L2();
	return norm2 <= std::max(
		stParam_.dCostTrackletLinkingMinDistance, 
		stParam_.dCalibrationError + stParam_.dDetectionError * (sensitivity1 + sensitivity2));
}

/************************************************************************
 Method Name: PointReconstruction
 Description:
	- Reconstruction 3D point with most recent points of 2D tracklets
 Input Arguments:
	- tracklet2Ds: 2D tracklets which are used to reconstruct 3D point
 Return Values:
	- CReconstruction: information structure of reconstructed point
************************************************************************/
CReconstruction CAssociator3D::PointReconstruction(CTrackletCombination &tracklet2Ds)
{
	CReconstruction resultReconstruction(nNumCameras_);
	resultReconstruction.bIsMeasurement     = false;
	resultReconstruction.point              = hj::Point3D(0.0, 0.0, 0.0);
	resultReconstruction.maxError           = 0.0;
	resultReconstruction.costReconstruction = DBL_MAX;
	resultReconstruction.costLink           = 0.0;
	resultReconstruction.velocity           = hj::Point3D(0.0, 0.0, 0.0);

	if (0 == tracklet2Ds.size()) { return resultReconstruction; }

	resultReconstruction.bIsMeasurement = true;
	resultReconstruction.tracklet2Ds    = tracklet2Ds;
	double fDistance = DBL_MAX;
	double maxError  = stParam_.dCalibrationError;
	//double maxError = 0;
	double probabilityReconstruction = 0.0;
	hj::Point2D curPoint(0.0, 0.0);

	//-------------------------------------------------
	// POINT RECONSTRUCTION
	//-------------------------------------------------
	switch (stParam_.nDetectionType)
	{
	case 0:
		// Full-body
		{
			std::vector<hj::Point2D_CamIdx> vecPointInfos;
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				CTracklet2D *curTracklet = resultReconstruction.tracklet2Ds.get(camIdx);
				if (NULL == curTracklet) { continue; }

				//// 0.8 shrink
				//hj::Rect curRect = curTracklet->rects.back();
				//curRect.y = curRect.y + 0.1 * curRect.h;
				//curRect.h *= 0.8;

				curPoint = curTracklet->rects.back().reconstructionPoint();
				vecPointInfos.push_back(hj::Point2D_CamIdx(curPoint, camIdx));
				maxError += stParam_.dDetectionError 
					* vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x);
				//maxError = std::max(maxError,  stParam_.dDetectionError * (double)vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x));

				// 3D position
				resultReconstruction.rawPoints.push_back(curTracklet->currentLocation3D);

				//// sensitivity
				//resultReconstruction.maxError += vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x);
			}
			//if (!stParam_.bConsiderSensitivity) { maxError = (double)stParam_.dMaxTrackletDistance / 2.0; }
			resultReconstruction.maxError = maxError;
			fDistance = this->NViewGroundingPointReconstruction(vecPointInfos, resultReconstruction.point);
		}
		break;
	default:
		// Head
		{
			std::vector<hj::Line3D> vecBackprojectionLines;
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				CTracklet2D *curTracklet = resultReconstruction.tracklet2Ds.get(camIdx);
				if (NULL == curTracklet) { continue; }

				curPoint = curTracklet->rects.back().reconstructionPoint();
				hj::Line3D curLine = curTracklet->backprojectionLines.back();

				vecBackprojectionLines.push_back(curLine);
				//maxError +=  stParam_.dDetectionError * vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x);
				maxError = std::max(maxError, stParam_.dDetectionError * (double)vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x));

				// 3D position
				resultReconstruction.rawPoints.push_back(curTracklet->currentLocation3D);

				//// sensitivity
				//resultReconstruction.maxError += vecMatProjectionSensitivity_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x);
			}
			if (!stParam_.bConsiderSensitivity) { maxError = stParam_.dMaxTrackletDistance / 2.0; }
			resultReconstruction.maxError = maxError;
			fDistance = this->NViewPointReconstruction(vecBackprojectionLines, resultReconstruction.point);
		}
		break;
	}
	resultReconstruction.smoothedPoint = resultReconstruction.point;

	if (2 > tracklet2Ds.size())
	{
		probabilityReconstruction = 0.5;
	}
	else
	{
		if (fDistance > maxError) { return resultReconstruction; }
		probabilityReconstruction = 1 == resultReconstruction.tracklet2Ds.size() ? 0.5 : 0.5 * hj::erfc(4.0 * fDistance / maxError - 2.0);
	}

	//-------------------------------------------------
	// DETECTION PROBABILITY
	//-------------------------------------------------
	double fDetectionProbabilityRatio = 1.0;
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		if (!this->CheckVisibility(resultReconstruction.point, camIdx)) { continue; }
		if (NULL == resultReconstruction.tracklet2Ds.get(camIdx))
		{
			// false negative
			fDetectionProbabilityRatio *= stParam_.dFNRate / (1 - stParam_.dFNRate);
			continue;
		}
		// positive
		fDetectionProbabilityRatio *= (1 - stParam_.dFPRate) / stParam_.dFPRate;
	}
	// reconstruction cost
	resultReconstruction.costReconstruction = log(1 - probabilityReconstruction) - log(probabilityReconstruction) - log(fDetectionProbabilityRatio);
	resultReconstruction.costSmoothedPoint = resultReconstruction.costReconstruction;

	return resultReconstruction;
}

/************************************************************************
Method Name: NViewPointReconstruction
Description:
- reconstruct 3D point with multiple back projection lines
Input Arguments:
- vecLines: vector of back projeciton lines
- outputPoint: (output) reconstructed point
Return Values:
- double: return distance between back projection lines and reconstructed point
************************************************************************/
double CAssociator3D::NViewPointReconstruction(std::vector<hj::Line3D> &vecLines, hj::Point3D &outputPoint)
{
	unsigned int numLines = (unsigned int)vecLines.size();
	if (numLines < 2)
	{
		outputPoint = vecLines[0].second;
		return stParam_.dMaxTrackletDistance / 2.0;
	}
	cv::Mat P = cv::Mat::zeros(3, 3, CV_64FC1);
	cv::Mat PP = cv::Mat::zeros(3, 3, CV_64FC1);
	cv::Mat A = cv::Mat::zeros(3, 3, CV_64FC1);
	cv::Mat b = cv::Mat::zeros(3, 1, CV_64FC1);
	cv::Mat resultPoint = cv::Mat::zeros(3, 1, CV_64FC1);
	cv::Mat I3x3 = cv::Mat::eye(3, 3, CV_64FC1);
	std::deque<cv::Mat> vecV;

	// reconstruction
	hj::Point3D vecDir;
	for (unsigned int lineIdx = 0; lineIdx < numLines; lineIdx++)
	{
		vecDir = vecLines[lineIdx].second - vecLines[lineIdx].first;
		vecDir /= vecDir.norm_L2();
		cv::Mat v = cv::Mat(vecDir.cv());
		vecV.push_back(v);
		cv::Mat s = cv::Mat(vecLines[lineIdx].first.cv());
		P = (v * v.t() - I3x3);
		PP = P.t() * P;
		A = A + PP;
		b = b + PP * s;
	}
	resultPoint = A.inv() * b;
	outputPoint.x = resultPoint.at<double>(0, 0);
	outputPoint.y = resultPoint.at<double>(1, 0);
	outputPoint.z = resultPoint.at<double>(2, 0);

	// error measurement
	double resultError = 0.0;
	double lambda = 0;
	cv::Mat diffVec = cv::Mat(3, 1, CV_64FC1);
	hj::Point3D diffPoint(0, 0, 0);
	for (unsigned int lineIdx = 0; lineIdx < numLines; lineIdx++)
	{
		cv::Mat s = cv::Mat(vecLines[lineIdx].first.cv());
		lambda = vecV[lineIdx].dot(resultPoint - s);
		diffVec = s + lambda * vecV[lineIdx] - resultPoint;
		diffPoint.x = diffVec.at<double>(0, 0);
		diffPoint.y = diffVec.at<double>(1, 0);
		diffPoint.z = diffVec.at<double>(2, 0);
		resultError += (1 / (double)numLines)*diffPoint.norm_L2();
	}

	return resultError;
}


/************************************************************************
Method Name: NViewGroundingPointReconstruction
Description:
- reconstruct 3D point with 2D grounding points (for PETS 9002 dataset)
Input Arguments:
- vecBottomCenterPoints: vector of bottom center point of each object
- outputPoint: (output) reconstructed point
Return Values:
- double: return distance between back projection lines and reconstructed point
************************************************************************/
double CAssociator3D::NViewGroundingPointReconstruction(std::vector<hj::Point2D_CamIdx> &vecPointInfos, hj::Point3D &outputPoint)
{
	double resultError = 0.0;
	double sumInvSensitivity = 0.0;
	unsigned int numPoints = (unsigned int)vecPointInfos.size();
	std::vector<hj::Point3D> vecReconstructedPoints;
	std::vector<double> vecInvSensitivity;
	vecReconstructedPoints.reserve(numPoints);

	outputPoint = hj::Point3D(0, 0, 0);
	for (unsigned int pointIdx = 0; pointIdx < numPoints; pointIdx++)
	{
		hj::Point2D curPoint = vecPointInfos[pointIdx].first;
		unsigned int curCamIdx = vecPointInfos[pointIdx].second;
		vecReconstructedPoints.push_back(this->ImageToWorld(curPoint, 0, curCamIdx));
		if (!stParam_.bConsiderSensitivity)
		{
			outputPoint += vecReconstructedPoints.back();
			continue;
		}

		int sensitivityX = std::min(std::max(0, (int)curPoint.x), vecMatProjectionSensitivity_[curCamIdx].cols);
		int sensitivityY = std::min(std::max(0, (int)curPoint.y), vecMatProjectionSensitivity_[curCamIdx].rows);
		double curSensitivity = vecMatProjectionSensitivity_[curCamIdx].at<float>(sensitivityY, sensitivityX);
		double invSensitivity = 1.0 / curSensitivity;

		vecInvSensitivity.push_back(invSensitivity);
		outputPoint += vecReconstructedPoints.back() * invSensitivity;
		sumInvSensitivity += invSensitivity;
	}
	if (stParam_.bConsiderSensitivity)
		outputPoint /= sumInvSensitivity;
	else
		outputPoint /= (double)numPoints;

	if (2 > numPoints)
	{
		if (stParam_.bConsiderSensitivity)
			return stParam_.dMaxSensitivityError;
		else
			return stParam_.dMaxTrackletDistance / 2.0;
	}
	for (unsigned int pointIdx = 0; pointIdx < numPoints; pointIdx++)
	{
		if (stParam_.bConsiderSensitivity)
			resultError += (1.0 / (double)numPoints) * vecInvSensitivity[pointIdx] * (outputPoint - vecReconstructedPoints[pointIdx]).norm_L2();
		else
			resultError += (1.0 / (double)numPoints) * (outputPoint - vecReconstructedPoints[pointIdx]).norm_L2();
	}

	return resultError;
}

/************************************************************************
Method Name: GetBackProjectionLine
Description:
- find two 3D points on the back-projection line correspondes to an input point
Input Arguments:
- point2D: image coordinate of the point
- camIdx: camera index (for calibration information)
Return Values:
- PSL_Line: a pair of 3D points have a height 2000 and 0
************************************************************************/
hj::Line3D CAssociator3D::GetBackProjectionLine(hj::Point2D point2D, int camIdx)
{
	hj::Point3D pointTop, pointBottom;
	pointTop = this->ImageToWorld(point2D, 2000, camIdx);
	pointBottom = this->ImageToWorld(point2D, 0, camIdx);
	return hj::Line3D(pointTop, pointBottom);
}

/************************************************************************
Method Name: GetDistanceFromBoundary
Description:
-
Input Arguments:
-
-
Return Values:
-
************************************************************************/
double CAssociator3D::GetDistanceFromBoundary(hj::Point3D point)
{
	hj::Point2D curPoint;
	double maxDistance = -100.0, curDistance = 0.0;
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		if (!CheckVisibility(point, camIdx, &curPoint)) { continue; }
		curDistance = vecMatDistanceFromBoundary_[camIdx].at<float>((int)curPoint.y, (int)curPoint.x);
		if (maxDistance < curDistance) { maxDistance = curDistance; }
	}
	return maxDistance;
}

/************************************************************************
Method Name: Tracklet2D_UpdateTracklets
Description:
- Update tracklets in each camera with information from 2D tracker
Input Arguments:
- curTrack2DResult: vector of 2D tracking results
- frameIdx: current frame index
Return Values:
- none
************************************************************************/
void CAssociator3D::Tracklet2D_UpdateTracklets(std::vector<CTrack2DResult> &curTrack2DResult, unsigned int frameIdx)
{
	bReceiveNewMeasurement_ = false;
	nNumTotalActive2DTracklet_ = 0;

	// check the validity of tracking result
	if (nNumCameras_ != curTrack2DResult.size()) { return; }

	/////////////////////////////////////////////////////////////////////
	// 2D TRACKLET UPDATE
	/////////////////////////////////////////////////////////////////////
	// find update infos (real updating will be done after this loop) and generate new tracklets
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		//this->m_matDetectionMap[camIdx] = cv::Mat::zeros(this->m_matDetectionMap[camIdx].rows, this->m_matDetectionMap[camIdx].cols, CV_64FC1);
		vecTracklet2DSet_[camIdx].newMeasurements.clear();
		if (frameIdx != curTrack2DResult[camIdx].frameIdx) { continue; }
		if (camIdx != curTrack2DResult[camIdx].camID) { continue; }

		unsigned int numObject = (unsigned int)curTrack2DResult[camIdx].object2DInfos.size();
		unsigned int numTracklet = (unsigned int)vecTracklet2DSet_[camIdx].activeTracklets.size();
		std::vector<CObject2DInfo*> tracklet2DUpdateInfos(numTracklet, NULL);

		//-------------------------------------------------
		// MATCHING AND GENERATING NEW 2D TRACKLET
		//-------------------------------------------------
		for (unsigned int objectIdx = 0; objectIdx < numObject; objectIdx++)
		{
			CObject2DInfo *curObject = &curTrack2DResult[camIdx].object2DInfos[objectIdx];

			// detection size crop
			curObject->box = curObject->box.cropWithSize(
				stParam_.vecPCalibrationInfo[camIdx]->cCamModel.width(), 
				stParam_.vecPCalibrationInfo[camIdx]->cCamModel.height());

			// set detection map
			//this->m_matDetectionMap[camIdx](curObject->box.cv()) = 1.0;

			// find appropriate 2D tracklet
			bool bNewTracklet2D = true;
			for (unsigned int tracklet2DIdx = 0; tracklet2DIdx < numTracklet; tracklet2DIdx++)
			{
				if (curObject->id == vecTracklet2DSet_[camIdx].activeTracklets[tracklet2DIdx]->id)
				{
					bNewTracklet2D = false;
					tracklet2DUpdateInfos[tracklet2DIdx] = curObject;
					break;
				}
			}
			if (!bNewTracklet2D) { continue; }

			// generate a new 2D tracklet
			CTracklet2D newTracklet(nNumCameras_);
			newTracklet.id = curObject->id;
			newTracklet.camIdx = camIdx;
			newTracklet.bActivated = true;
			newTracklet.rects.push_back(curObject->box);
			newTracklet.backprojectionLines.push_back(
				this->GetBackProjectionLine(curObject->box.bottomCenter(), camIdx));
			newTracklet.timeStart = frameIdx;
			newTracklet.timeEnd = frameIdx;
			newTracklet.duration = 1;

			// appearance
			cv::Rect cropRect = curObject->box.cv();
			if (cropRect.x < 0) { cropRect.x = 0; }
			if (cropRect.y < 0) { cropRect.y = 0; }
			if (cropRect.x + cropRect.width > vecMatCurrentFrames_[camIdx].cols)
				cropRect.width = vecMatCurrentFrames_[camIdx].cols - cropRect.x;
			if (cropRect.y + cropRect.height > vecMatCurrentFrames_[camIdx].rows)
				cropRect.height = vecMatCurrentFrames_[camIdx].rows - cropRect.y;

			cv::Mat patch = vecMatCurrentFrames_[camIdx](cropRect);
			newTracklet.RGBFeatureHead = GetRGBFeature(&patch, stParam_.nNumRGBHistogramBins);
			newTracklet.RGBFeatureTail = newTracklet.RGBFeatureHead.clone();

			// location in 3D
			newTracklet.currentLocation3D = this->ImageToWorld(curObject->box.bottomCenter(), 0.0, camIdx);

			// update tracklet info
			vecTracklet2DSet_[camIdx].tracklets.push_back(newTracklet);
			vecTracklet2DSet_[camIdx].activeTracklets.push_back(&vecTracklet2DSet_[camIdx].tracklets.back());
			vecTracklet2DSet_[camIdx].newMeasurements.push_back(&vecTracklet2DSet_[camIdx].tracklets.back());
			nNumTotalActive2DTracklet_++;
			bReceiveNewMeasurement_ = true;
		}

		//-------------------------------------------------
		// UPDATING EXISTING 2D TRACKLETS
		//-------------------------------------------------
		std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[camIdx].activeTracklets.begin();
		for (unsigned int objInfoIdx = 0; objInfoIdx < numTracklet; objInfoIdx++)
		{
			if (NULL == tracklet2DUpdateInfos[objInfoIdx])
			{
				if (!(*trackletIter)->bActivated)
				{
					trackletIter = vecTracklet2DSet_[camIdx].activeTracklets.erase(trackletIter);
					continue;
				}
				(*trackletIter)->bActivated = false;
				trackletIter++;
				continue;
			}

			// update tracklet
			CTracklet2D *curTracklet = *trackletIter;
			CObject2DInfo *curObject = tracklet2DUpdateInfos[objInfoIdx];
			curTracklet->bActivated = true;
			curTracklet->rects.push_back(curObject->box);
			curTracklet->backprojectionLines.push_back(
				this->GetBackProjectionLine(curObject->box.center(), camIdx));
			curTracklet->timeEnd = frameIdx;
			curTracklet->duration++;

			// appearance
			cv::Rect cropRect = curTracklet->rects.back().cv();
			if (cropRect.x < 0) { cropRect.x = 0; }
			if (cropRect.y < 0) { cropRect.y = 0; }
			if (cropRect.x + cropRect.width > vecMatCurrentFrames_[camIdx].cols)
				cropRect.width = vecMatCurrentFrames_[camIdx].cols - cropRect.x;
			if (cropRect.y + cropRect.height > vecMatCurrentFrames_[camIdx].rows)
				cropRect.height = vecMatCurrentFrames_[camIdx].rows - cropRect.y;

			cv::Mat patch = vecMatCurrentFrames_[camIdx](cropRect);
			curTracklet->RGBFeatureTail = GetRGBFeature(&patch, stParam_.nNumRGBHistogramBins);

			// location in 3D
			curTracklet->currentLocation3D = this->ImageToWorld(curObject->box.bottomCenter(), 0.0, camIdx);

			// association informations
			for (int subloopCamIdx = 0; subloopCamIdx < nNumCameras_; subloopCamIdx++)
			{
				curTracklet->bAssociableNewMeasurement[subloopCamIdx].clear();
			}

			// increase iterators
			trackletIter++;
		}
	}

	/////////////////////////////////////////////////////////////////////
	// ASSOCIATION CHECK
	/////////////////////////////////////////////////////////////////////
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		if (!bReceiveNewMeasurement_) { break; }

		for (std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[camIdx].activeTracklets.begin();
			trackletIter != vecTracklet2DSet_[camIdx].activeTracklets.end();
			trackletIter++)
		{
			hj::Line3D backProjectionLine1 = (*trackletIter)->backprojectionLines.back();
			for (int subloopCamIdx = 0; subloopCamIdx < nNumCameras_; subloopCamIdx++)
			{
				int numNewMeasurements = (int)vecTracklet2DSet_[subloopCamIdx].newMeasurements.size();
				(*trackletIter)->bAssociableNewMeasurement[subloopCamIdx].resize(numNewMeasurements, false);

				if (camIdx == subloopCamIdx) { continue; }
				for (int measurementIdx = 0; measurementIdx < numNewMeasurements; measurementIdx++)
				{
					CTracklet2D *curMeasurement = vecTracklet2DSet_[subloopCamIdx].newMeasurements[measurementIdx];
					hj::Line3D backProjectionLine2 = curMeasurement->backprojectionLines.back();

					hj::Point3D reconstructedPoint;
					std::vector<hj::Line3D> vecBackProjectionLines;
					vecBackProjectionLines.push_back(backProjectionLine1);
					vecBackProjectionLines.push_back(backProjectionLine2);

					double fDistance = this->NViewPointReconstruction(vecBackProjectionLines, reconstructedPoint);
					//double fDistance = this->StereoTrackletReconstruction(*trackletIter, curMeasurement, reconstructedPoints);

					if (stParam_.dMaxTrackletDistance >= fDistance)
					{
						(*trackletIter)->bAssociableNewMeasurement[subloopCamIdx][measurementIdx] = true;
					}
				}
			}
		}
	}
}

/************************************************************************
Method Name: GenerateTrackletCombinations
Description:
- Generate feasible 2D tracklet combinations
Input Arguments:
- vecBAssociationMap: feasible association map
- combination: current combination
- combinationQueue: queue for save combinations
- camIdx: current camera index
Return Values:
- none
************************************************************************/
void CAssociator3D::GenerateTrackletCombinations(
	std::vector<std::vector<bool>> &vecvecBAssociationMap,
	CTrackletCombination combination,
	std::deque<CTrackletCombination> &combinationQueue,
	int camIdx)
{
	if (camIdx >= nNumCameras_)
	{
		combinationQueue.push_back(combination);
		return;
	}

	if (NULL != combination.get(camIdx))
	{
		std::vector<std::vector<bool>> vecvecNewBAssociationMap(nNumCameras_);
		for (int subloopCamIdx = 0; subloopCamIdx < nNumCameras_; subloopCamIdx++)
		{
			vecvecNewBAssociationMap[subloopCamIdx] = vecvecBAssociationMap[subloopCamIdx];
			if (subloopCamIdx <= camIdx) { continue; }
			for (unsigned int mapIdx = 0; mapIdx < vecvecNewBAssociationMap[subloopCamIdx].size(); mapIdx++)
			{
				bool currentFlag = vecvecNewBAssociationMap[subloopCamIdx][mapIdx];
				bool measurementFlag = combination.get(camIdx)->bAssociableNewMeasurement[subloopCamIdx][mapIdx];
				vecvecNewBAssociationMap[subloopCamIdx][mapIdx] = currentFlag & measurementFlag;
			}
		}
		this->GenerateTrackletCombinations(vecvecNewBAssociationMap, combination, combinationQueue, camIdx + 1);
		return;
	}

	// null tracklet
	combination.set(camIdx, NULL);
	this->GenerateTrackletCombinations(vecvecBAssociationMap, combination, combinationQueue, camIdx + 1);

	for (int measurementIdx = 0; measurementIdx < vecTracklet2DSet_[camIdx].newMeasurements.size(); measurementIdx++)
	{
		if (!vecvecBAssociationMap[camIdx][measurementIdx]) { continue; }

		combination.set(camIdx, vecTracklet2DSet_[camIdx].newMeasurements[measurementIdx]);
		// AND operation of map
		std::vector<std::vector<bool>> vecvecNewBAssociationMap(nNumCameras_);
		for (int subloopCamIdx = 0; subloopCamIdx < nNumCameras_; subloopCamIdx++)
		{
			vecvecNewBAssociationMap[subloopCamIdx] = vecvecBAssociationMap[subloopCamIdx];
			if (subloopCamIdx <= camIdx) { continue; }
			for (int mapIdx = 0; mapIdx < vecvecNewBAssociationMap[subloopCamIdx].size(); mapIdx++)
			{
				bool currentFlag = vecvecNewBAssociationMap[subloopCamIdx][mapIdx];
				bool measurementFlag = vecTracklet2DSet_[camIdx].newMeasurements[measurementIdx]->bAssociableNewMeasurement[subloopCamIdx][mapIdx];
				vecvecNewBAssociationMap[subloopCamIdx][mapIdx] = currentFlag & measurementFlag;
			}
		}
		this->GenerateTrackletCombinations(vecvecNewBAssociationMap, combination, combinationQueue, camIdx + 1);
	}
}

/************************************************************************
Method Name: Track3D_Management
Description:
- updage established 3D tracks
Input Arguments:
- none
Return Values:
- none
************************************************************************/
void CAssociator3D::Track3D_Management(hj::TrackSet &outputSeedTracks, const unsigned int _frameIdx)
{
	outputSeedTracks.clear();

	// update 3D tracks
	this->Track3D_UpdateTracks(_frameIdx);

	// generate 3D tracks
	if (bReceiveNewMeasurement_)
	{
		this->Track3D_GenerateSeedTracks(outputSeedTracks, _frameIdx);
		this->Track3D_BranchTracks(&outputSeedTracks, _frameIdx);
	}

#if hj_PRINT_TRACKS == 1
	// track print
	hj::CreateDirectoryForWindows(std::string(TRACK_SAVE_PATH));
	char strTrackFileName[128];
	sprintf_s(strTrackFileName, "%s%04d.txt", TRACK_SAVE_PATH, _frameIdx);
	PrintTracks(queueTracksInWindow_, strTrackFileName, false);
#endif
}

/************************************************************************
Method Name: Track3D_UpdateTracks
Description:
- updage established 3D tracks
Input Arguments:
- none
Return Values:
- none
************************************************************************/
void CAssociator3D::Track3D_UpdateTracks(const unsigned int _frameIdx)
{
	std::deque<CTrack3D*> queueNewTracks;

	//---------------------------------------------------------
	// UPDATE ACTIVE TRACKS
	//---------------------------------------------------------	
	for (std::deque<CTrack3D*>::iterator trackIter = queueActiveTrack_.begin();
		trackIter != queueActiveTrack_.end();
		trackIter++)
	{
		if (!(*trackIter)->bValid) { continue; }
		CTrack3D *curTrack = *trackIter;

		// updating current tracklet information
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			if (NULL == curTrack->curTracklet2Ds.get(camIdx)) { continue; }
			if (!curTrack->curTracklet2Ds.get(camIdx)->bActivated)
			{
				if ((unsigned int)stParam_.nMinTrackletLength > curTrack->curTracklet2Ds.get(camIdx)->duration)
				{
					// invalidate track with 2D tracklet which has duration 1
					CTrackTree::SetValidityFlagInTrackBranch(curTrack, false);
					break;
				}
				curTrack->curTracklet2Ds.set(camIdx, NULL);
			}
			else
			{
				// time stamp
				curTrack->timeTrackletEnded[camIdx] = _frameIdx;
				curTrack->lastTrackletLocation3D[camIdx] = curTrack->curTracklet2Ds.get(camIdx)->currentLocation3D;
				hj::Point2D bottomCenter = curTrack->curTracklet2Ds.get(camIdx)->rects.back().bottomCenter();
				curTrack->lastTrackletSensitivity[camIdx] = vecMatProjectionSensitivity_[camIdx].at<float>((int)bottomCenter.y, (int)bottomCenter.x);
				curTrack->lastRGBFeature[camIdx].release();
				curTrack->lastRGBFeature[camIdx] = curTrack->curTracklet2Ds.get(camIdx)->RGBFeatureTail.clone();
			}
		}
		if (!curTrack->bValid) { continue; }

		// activation check and expiration
		if (0 == curTrack->curTracklet2Ds.size())
		{
			// de-activating
			curTrack->bActive = false;

			// cost update (with exit cost)
			std::vector<hj::Point3D> pointsIn3D;
			CReconstruction laCReconstruction = curTrack->reconstructions.back();
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				if (NULL == laCReconstruction.tracklet2Ds.get(camIdx)) { continue; }
				pointsIn3D.push_back(curTrack->lastTrackletLocation3D[camIdx]);
			}
			curTrack->costExit = std::min(
				stParam_.dCostExitMax, 
				-std::log(this->ComputeExitProbability(pointsIn3D, (int)curTrack->duration)));
			curTrack->costTotal = GetCost(curTrack);

			// move to time jump track list
			queuePausedTrack_.push_back(curTrack);
			continue;
		}

		// cost update (not essential, for just debuging)
		curTrack->costTotal = GetCost(curTrack);

		// point reconstruction and cost update
		CReconstruction newReconstruction = this->PointReconstruction(curTrack->curTracklet2Ds);
		newReconstruction.costSmoothedPoint = newReconstruction.costReconstruction;;
		newReconstruction.smoothedPoint = newReconstruction.point;

		//double curLinkProbability = ComputeLinkProbability(curTrack->reconstructions.back().point + curTrack->reconstructions.back().velocity, curReconstruction.point, 1);
		double curLinkProbability = ComputeLinkProbability(curTrack->reconstructions.back().point, newReconstruction.point);
		if (DBL_MAX == newReconstruction.costReconstruction || stParam_.dMinLinkingProbability > curLinkProbability)
		{
			// invalidate track
			curTrack->bValid = false;
			//this->SetValidityFlagInTrackBranch(curTrack, false);
			continue;
		}
		newReconstruction.costLink = -log(curLinkProbability);
		newReconstruction.velocity = newReconstruction.point - curTrack->reconstructions.back().smoothedPoint;

		curTrack->reconstructions.push_back(newReconstruction);
		curTrack->costReconstruction += newReconstruction.costSmoothedPoint;
		curTrack->costLink += newReconstruction.costLink;
		curTrack->duration++;
		curTrack->timeEnd++;

		// smoothing		
		int updateStartPos = curTrack->smoother.Insert(newReconstruction.point);
		bool trackValidity = true;
		double smoothedProbability = 0.0;
		for (int pos = updateStartPos; pos < curTrack->smoother.size(); pos++)
		{
			CReconstruction *curReconstruction = &curTrack->reconstructions[pos];
			if (MIN_SMOOTHING_LENGTH <= curTrack->duration)
			{
				curReconstruction->smoothedPoint = curTrack->smoother.GetResult(pos);
				//curTrack->smoothedTrajectory[pos] = curReconstruction->smoothedPoint;

				// update reconstruction cost
				curTrack->costReconstruction -= curReconstruction->costSmoothedPoint;
				smoothedProbability = ComputeReconstructionProbability(curReconstruction->smoothedPoint, &curReconstruction->rawPoints, &curReconstruction->tracklet2Ds, curReconstruction->maxError);
				if (0.0 == smoothedProbability)
				{
					trackValidity = false;
					break;
				}
				curReconstruction->costSmoothedPoint = -log(smoothedProbability);
				curTrack->costReconstruction += curReconstruction->costSmoothedPoint;

				if (0 == pos) { continue; }
				// update link cost
				curTrack->costLink -= curReconstruction->costLink;

				// TESTING
				smoothedProbability = ComputeLinkProbability(curTrack->reconstructions[pos - 1].smoothedPoint, curReconstruction->smoothedPoint);
				//smoothedProbability = ComputeLinkProbability(curTrack->reconstructions[pos-1].smoothedPoint + curTrack->reconstructions[pos-1].velocity, curReconstruction->smoothedPoint);

				if (0.0 == smoothedProbability)
				{
					trackValidity = false;
					break;
				}
				curReconstruction->costLink = -log(smoothedProbability);
				curTrack->costLink += curReconstruction->costLink;
			}
			else
			{
				curReconstruction->smoothedPoint = curReconstruction->point;
				if (0 == pos) { continue; }
			}

			// update velocity
			curReconstruction->velocity = curReconstruction->smoothedPoint - curTrack->reconstructions[pos - 1].smoothedPoint;
			if (stParam_.dMinMovingSpeed > curReconstruction->velocity.norm_L2()) 
			{ 
				curReconstruction->velocity = hj::Point3D(0.0, 0.0, 0.0); 
			}
		}
		curTrack->numOutpoint = 0;

		// increase iterator
		queueNewTracks.push_back(curTrack);
	}
	queueActiveTrack_ = queueNewTracks;
	queueNewTracks.clear();


	//---------------------------------------------------------
	// UPDATE DE-ACTIVATED TRACKS FOR TEMPORAL BRANCHING
	//---------------------------------------------------------
	std::deque<CTrack3D*> queueTerminatedTracks;
	for (std::deque<CTrack3D*>::iterator trackIter = queuePausedTrack_.begin();
		trackIter != queuePausedTrack_.end();
		trackIter++)
	{
		if (!(*trackIter)->bValid) { continue; }

		// handling expired track
		if ((*trackIter)->timeEnd + stParam_.nMaxTimeJump < _frameIdx)
		{
			// remove track instance for memory efficiency
			if (0.0 <= (*trackIter)->costTotal)
			{
				(*trackIter)->bValid = false;
			}
			else
			{
				// for logging
				queueTerminatedTracks.push_back(*trackIter);
			}
			continue;
		}

		// insert dummy reconstruction
		CReconstruction dummyReconstruction = this->PointReconstruction((*trackIter)->curTracklet2Ds);
		dummyReconstruction.bIsMeasurement = false;
		dummyReconstruction.point = (*trackIter)->reconstructions.back().smoothedPoint + (*trackIter)->reconstructions.back().velocity;
		dummyReconstruction.smoothedPoint = dummyReconstruction.point;
		dummyReconstruction.costLink = 0.0;
		dummyReconstruction.costReconstruction = 0.0;
		dummyReconstruction.costSmoothedPoint = 0.0;
		dummyReconstruction.maxError = 0.0;
		dummyReconstruction.velocity = (*trackIter)->reconstructions.back().velocity;
		(*trackIter)->reconstructions.push_back(dummyReconstruction);

		bool bPointOut = true;
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			if (this->CheckVisibility(dummyReconstruction.smoothedPoint, camIdx))
			{
				bPointOut = false;
				break;
			}
		}
		if (bPointOut) { (*trackIter)->numOutpoint++; }
		if (stParam_.nMaxOutpoint < (*trackIter)->numOutpoint)
		{
			// for logging
			queueTerminatedTracks.push_back(*trackIter);
			continue;
		}
		queueNewTracks.push_back(*trackIter);
	}
	queuePausedTrack_ = queueNewTracks;
	queueNewTracks.clear();
	queueTerminatedTracks.clear();

	//---------------------------------------------------------
	// MANAGE TRACKS IN PROCESSING WINDOW
	//---------------------------------------------------------
	for (std::deque<CTrack3D*>::iterator trackIter = queueTracksInWindow_.begin();
		trackIter != queueTracksInWindow_.end();
		trackIter++)
	{
		if (!(*trackIter)->bValid
			|| (*trackIter)->timeEnd + stParam_.nProcWindowSize <= _frameIdx)
		{
			continue;
		}

		queueNewTracks.push_back(*trackIter);
		//if ((*trackIter)->tree->maxGTProb < (*trackIter)->GTProb)
		//{
		//	(*trackIter)->tree->maxGTProb = (*trackIter)->GTProb;
		//}
	}
	queueTracksInWindow_ = queueNewTracks;
	queueNewTracks.clear();

	//---------------------------------------------------------
	// UPDATE TRACK TREES
	//---------------------------------------------------------
	// delete empty trees from the active tree list
	std::deque<CTrackTree*> queueNewActiveTrees;
	for (std::deque<CTrackTree*>::iterator treeIter = queuePtActiveTrees_.begin();
		treeIter != queuePtActiveTrees_.end();
		treeIter++)
	{
		// track validation
		double GTPSum = 0.0;
		std::deque<CTrack3D*> queueUpdated; // copy is faster than delete
		for (std::deque<CTrack3D*>::iterator trackIter = (*treeIter)->tracks.begin();
			trackIter != (*treeIter)->tracks.end();
			trackIter++)
		{
			GTPSum += (*trackIter)->GTProb;
			// reset fields for optimization
			(*trackIter)->BranchGTProb = 0.0;
			(*trackIter)->GTProb = 0.0;
			(*trackIter)->bCurrentBestSolution = false;

			//if (!(*trackIter)->bValid && (*trackIter)->timeGeneration + stParam_.nProcWindowSize > _frameIdx) { continue; }
			if (!(*trackIter)->bValid) { continue; }
			queueUpdated.push_back(*trackIter);
		}
		(*treeIter)->tracks = queueUpdated;
		//if (0 == (*treeIter)->tracks.size() || (0 == GTPSum && (*treeIter)->timeGeneration + stParam_.nNumFrameForConfirmation <= _frameIdx))
		if (0 == (*treeIter)->tracks.size())
		{
			(*treeIter)->bValid = false;
			continue;
		}
		queueNewActiveTrees.push_back(*treeIter);
	}
	queuePtActiveTrees_ = queueNewActiveTrees;
	queueNewActiveTrees.clear();

	// update unconfirmed trees
	std::deque<CTrackTree*> queueNewUnconfirmedTrees;
	for (std::deque<CTrackTree*>::iterator treeIter = queuePtUnconfirmedTrees_.begin();
		treeIter != queuePtUnconfirmedTrees_.end();
		treeIter++)
	{
		if (!(*treeIter)->bValid || (*treeIter)->bConfirmed) { continue; }
		if ((*treeIter)->timeGeneration + stParam_.nNumFramesForConfirmation <= _frameIdx)
		{
			(*treeIter)->bConfirmed = true;
			continue;
		}
		queueNewUnconfirmedTrees.push_back(*treeIter);
	}
	queuePtUnconfirmedTrees_ = queueNewUnconfirmedTrees;
	queueNewUnconfirmedTrees.clear();

	//---------------------------------------------------------
	// UPDATE HYPOTHESES
	//---------------------------------------------------------
	hj::HypothesisSet validHypothesesSet;
	for (hj::HypothesisSet::iterator hypothesisIter = queuePrevGlobalHypotheses_.begin();
		hypothesisIter != queuePrevGlobalHypotheses_.end();
		hypothesisIter++)
	{
		// validation
		for (size_t trackIdx = 0; trackIdx < (*hypothesisIter).selectedTracks.size(); trackIdx++)
		{
			if ((*hypothesisIter).selectedTracks[trackIdx]->bValid) { continue; }
			(*hypothesisIter).bValid = false;
			break;
		}
		if (!(*hypothesisIter).bValid) { continue; }

		// update related tracks
		hj::TrackSet newRelatedTrackSet;
		for (size_t trackIdx = 0; trackIdx < (*hypothesisIter).relatedTracks.size(); trackIdx++)
		{
			if (!(*hypothesisIter).relatedTracks[trackIdx]->bValid) { continue; }
			newRelatedTrackSet.push_back((*hypothesisIter).relatedTracks[trackIdx]);
		}
		(*hypothesisIter).relatedTracks = newRelatedTrackSet;
		validHypothesesSet.push_back(*hypothesisIter);
	}
	queuePrevGlobalHypotheses_ = validHypothesesSet;

	//---------------------------------------------------------
	// UPDATE TRACK LIST
	//---------------------------------------------------------
	// delete invalid tracks (in the list of track instances)
	for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
		trackIter != listTrack3D_.end();
		/* do in the loop */)
	{
		// delete invalid instance
		if (trackIter->bValid)
		{
			// clear new flag
			trackIter->bNewTrack = false;
			trackIter++;
			continue;
		}

		// delete the instance only when the tree is invalidated because of preserving the data structure
		if (trackIter->tree->bValid)
		{
			trackIter++;
			continue;
		}
		trackIter = listTrack3D_.erase(trackIter);
	}
}


/************************************************************************
Method Name: Track3D_GenerateSeedTracks
Description:
- Generate seeds of track tree with combinations of new measurements
Input Arguments:
- outputSeedTracks: (output) generated seed tracks
Return Values:
- none
************************************************************************/
void CAssociator3D::Track3D_GenerateSeedTracks(hj::TrackSet &outputSeedTracks, const unsigned int _frameIdx)
{
	// initialization
	outputSeedTracks.clear();

	//---------------------------------------------------------
	// FEASIBLE COMBINATIONS
	//---------------------------------------------------------
	std::deque<CTrackletCombination> queueSeeds;

	CTrackletCombination curCombination(nNumCameras_);
	std::vector<std::vector<bool>> vecvecNullAssociateMap(nNumCameras_);

	// null tracklet at the first camera
	vecvecNullAssociateMap[0].resize(vecTracklet2DSet_[0].newMeasurements.size(), false);
	for (int camIdx = 1; camIdx < nNumCameras_; camIdx++)
	{
		vecvecNullAssociateMap[camIdx].resize(vecTracklet2DSet_[camIdx].newMeasurements.size(), true);
	}
	this->GenerateTrackletCombinations(vecvecNullAssociateMap, curCombination, queueSeeds, 1);

	// real tracklets at the first camera
	for (std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[0].newMeasurements.begin();
		trackletIter != vecTracklet2DSet_[0].newMeasurements.end();
		trackletIter++)
	{
		curCombination.set(0, *trackletIter);
		this->GenerateTrackletCombinations((*trackletIter)->bAssociableNewMeasurement, curCombination, queueSeeds, 1);
	}

	// delete null track
	queueSeeds.erase(queueSeeds.begin());

	//---------------------------------------------------------
	// MAKE TRACKS AND TREES
	//---------------------------------------------------------
	for (unsigned int seedIdx = 0; seedIdx < queueSeeds.size(); seedIdx++)
	{
		// generate track
		CTrack3D newTrack(nNumCameras_);
		newTrack.Initialize(queueSeeds[seedIdx], nNewTrackID_, _frameIdx, NULL);

		// tracklet information
		newTrack.costRGB = 0.0;
		std::vector<hj::Point3D> pointsIn3D;
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			newTrack.timeTrackletEnded[camIdx] = 0;
			if (NULL == newTrack.curTracklet2Ds.get(camIdx)) { continue; }

			// ID
			newTrack.tracklet2DIDRecord[camIdx].push_back(newTrack.curTracklet2Ds.get(camIdx)->id);
			pointsIn3D.push_back(newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D);

			// others
			newTrack.timeTrackletEnded[camIdx]      = _frameIdx;
			newTrack.lastTrackletLocation3D[camIdx] = newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D;
			
			hj::Point2D bottomCenter = newTrack.curTracklet2Ds.get(camIdx)->rects.back().bottomCenter();
			newTrack.lastTrackletSensitivity[camIdx] = 
				vecMatProjectionSensitivity_[camIdx].at<float>((int)bottomCenter.y, (int)bottomCenter.x);
			newTrack.lastRGBFeature[camIdx] = 
				newTrack.curTracklet2Ds.get(camIdx)->RGBFeatureTail.clone();
		}

		// initiation cost
		newTrack.costEnter = bInitiationPenaltyFree_ ? 
			0 : std::min(stParam_.dCostEnterMax, -std::log(this->ComputeEnterProbability(pointsIn3D)));
		pointsIn3D.clear();

		// point reconstruction
		CReconstruction curReconstruction = this->PointReconstruction(newTrack.curTracklet2Ds);
		if (DBL_MAX == curReconstruction.costReconstruction) { continue; }
		newTrack.costReconstruction = curReconstruction.costReconstruction;
		newTrack.reconstructions.push_back(curReconstruction);

		// smoothing
		newTrack.smoother.SetQsets(&precomputedQsets);
		newTrack.smoother.Insert(curReconstruction.point);
		//newTrack.smoothedTrajectory.push_back(newTrack.smoother.GetResult(0));

		// generate a track instance
		listTrack3D_.push_back(newTrack);
		nNewTrackID_++;
		CTrack3D *curTrack = &listTrack3D_.back();

		// generate a new track tree
		CTrackTree newTree;
		newTree.Initialize(curTrack, nNewTreeID_++, _frameIdx, listTrackTree_);
		queuePtActiveTrees_.push_back(curTrack->tree);
		queuePtUnconfirmedTrees_.push_back(curTrack->tree);

		// insert to queues
		outputSeedTracks.push_back(curTrack);
	}
}


/************************************************************************
Method Name: Track3D_BranchTracks
Description:
- Generate branches of track tree with combinations of new measurements
Input Arguments:
- seedTracks: seed tracks
- outputBranchTracks: (output) track generated by branching step
Return Values:
- none
************************************************************************/
void CAssociator3D::Track3D_BranchTracks(hj::TrackSet *seedTracks, const unsigned int _frameIdx)
{
	//for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	//{
	//	// TODO: 브랜치 뜬 track을, tracklet 안에서 visual feature distance가 가장 작은 것부터 sorting, 점증적으로 cost 주기
	//}

	/////////////////////////////////////////////////////////////////////
	// SPATIAL BRANCHING
	/////////////////////////////////////////////////////////////////////
	hj::TrackSet queueBranchTracks;
	std::sort(queueActiveTrack_.begin(), queueActiveTrack_.end(), hjTrackGTPandLLDescend);
	for (std::deque<CTrack3D*>::iterator trackIter = queueActiveTrack_.begin();
		trackIter != queueActiveTrack_.end();
		trackIter++)
	{
		if (stParam_.bDoBranchCut && stParam_.nMaxNumTracksInOptimization <= queueBranchTracks.size()) { break; }

		//---------------------------------------------------------
		// FIND SPATIAL ASSOCIATION
		//---------------------------------------------------------	
		CTrackletCombination curCombination = (*trackIter)->curTracklet2Ds;

		// association map
		std::vector<std::vector<bool>> vecvecAssociationMap(nNumCameras_);
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			vecvecAssociationMap[camIdx].resize(vecTracklet2DSet_[camIdx].newMeasurements.size(), true);
		}
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			if (NULL == curCombination.get(camIdx)) { continue; }
			for (int subloopCamIdx = 0; subloopCamIdx < nNumCameras_; subloopCamIdx++)
			{
				for (int flagIdx = 0; flagIdx < vecTracklet2DSet_[subloopCamIdx].newMeasurements.size(); flagIdx++)
				{
					bool curTrackFlag = curCombination.get(camIdx)->bAssociableNewMeasurement[subloopCamIdx][flagIdx];
					bool curFlag = vecvecAssociationMap[subloopCamIdx][flagIdx];
					vecvecAssociationMap[subloopCamIdx][flagIdx] = curFlag && curTrackFlag;
				}
			}
		}

		// find feasible branches
		std::deque<CTrackletCombination> queueBranches;
		this->GenerateTrackletCombinations(vecvecAssociationMap, curCombination, queueBranches, 0);

		if (0 == queueBranches.size()) { continue; }
		if (1 == queueBranches.size() && queueBranches[0] == curCombination) { continue; }

		//---------------------------------------------------------
		// BRANCHING
		//---------------------------------------------------------
		for (std::deque<CTrackletCombination>::iterator branchIter = queueBranches.begin();
			branchIter != queueBranches.end();
			branchIter++)
		{
			if (*branchIter == curCombination) { continue; }
			CTrack3D *curTrack = *trackIter;

			// generate a new track with branching combination
			CTrack3D newTrack(nNumCameras_);
			newTrack.id = nNewTrackID_;
			newTrack.curTracklet2Ds = *branchIter;
			newTrack.bActive = true;
			newTrack.bValid = true;
			newTrack.tree = curTrack->tree;
			newTrack.parentTrack = curTrack;
			newTrack.childrenTrack.clear();
			newTrack.timeStart = curTrack->timeStart;
			newTrack.timeEnd = curTrack->timeEnd;
			newTrack.timeGeneration = _frameIdx;
			newTrack.duration = curTrack->duration;
			newTrack.bWasBestSolution = true;
			newTrack.bCurrentBestSolution = false;
			newTrack.bNewTrack = true;
			newTrack.GTProb = curTrack->GTProb;
			newTrack.costEnter = curTrack->costEnter;
			newTrack.costReconstruction = curTrack->costReconstruction;
			newTrack.costLink = curTrack->costLink;
			newTrack.costRGB = curTrack->costRGB;
			newTrack.costExit = 0.0;
			newTrack.costTotal = 0.0;

			// reconstruction
			CReconstruction newReconstruction = this->PointReconstruction(newTrack.curTracklet2Ds);
			if (DBL_MAX == newReconstruction.costReconstruction) { continue; }

			CReconstruction oldReconstruction = curTrack->reconstructions.back();
			CReconstruction preReconstruction = curTrack->reconstructions[(*trackIter)->reconstructions.size() - 2];
			//double curLinkProbability = ComputeLinkProbability(preReconstruction.point + preReconstruction.velocity, curReconstruction.point, 1);
			double curLinkProbability = ComputeLinkProbability(preReconstruction.point, newReconstruction.point);
			if (stParam_.dMinLinkingProbability > curLinkProbability) { continue; }

			newTrack.reconstructions = curTrack->reconstructions;
			newTrack.smoother = curTrack->smoother;

			newTrack.reconstructions.pop_back();
			newTrack.reconstructions.push_back(newReconstruction);
			newTrack.costReconstruction += newReconstruction.costReconstruction - oldReconstruction.costReconstruction;
			newTrack.costLink += -log(curLinkProbability) - oldReconstruction.costLink;

			// smoothing			
			newTrack.smoother.PopBack();
			int updateStartPos = newTrack.smoother.Insert(newReconstruction.point);
			bool trackValidity = true;
			double smoothedProbability;
			for (int pos = updateStartPos; pos < newTrack.smoother.size(); pos++)
			{
				CReconstruction *curReconstruction = &newTrack.reconstructions[pos];
				if (MIN_SMOOTHING_LENGTH <= curTrack->duration)
				{
					curReconstruction->smoothedPoint = newTrack.smoother.GetResult(pos);
					//newTrack.smoothedTrajectory[pos] = curReconstruction->smoothedPoint;

					// update reconstruction cost
					newTrack.costReconstruction -= curReconstruction->costSmoothedPoint;
					smoothedProbability = ComputeReconstructionProbability(curReconstruction->smoothedPoint, &curReconstruction->rawPoints, &curReconstruction->tracklet2Ds, curReconstruction->maxError);
					if (0.0 == smoothedProbability)
					{
						trackValidity = false;
						break;
					}
					curReconstruction->costSmoothedPoint = -log(smoothedProbability);
					newTrack.costReconstruction += curReconstruction->costSmoothedPoint;

					if (0 == pos) { continue; }
					// update link cost
					newTrack.costLink -= curReconstruction->costLink;

					// TESTING
					smoothedProbability = ComputeLinkProbability(newTrack.reconstructions[pos - 1].smoothedPoint, curReconstruction->smoothedPoint);
					//smoothedProbability = ComputeLinkProbability(newTrack.reconstructions[pos-1].smoothedPoint + newTrack.reconstructions[pos-1].velocity, curReconstruction->smoothedPoint);

					if (0.0 == smoothedProbability)
					{
						trackValidity = false;
						break;
					}
					curReconstruction->costLink = -log(smoothedProbability);
					newTrack.costLink += curReconstruction->costLink;
				}
				else
				{
					curReconstruction->smoothedPoint = curReconstruction->point;
					if (0 == pos) { continue; }
				}
				// update velocity
				curReconstruction->velocity = curReconstruction->smoothedPoint - newTrack.reconstructions[pos - 1].smoothedPoint;
				if (stParam_.dMinMovingSpeed > curReconstruction->velocity.norm_L2()) 
				{ 
					curReconstruction->velocity = hj::Point3D(0.0, 0.0, 0.0); 
				}
			}
			if (!trackValidity) { continue; }

			// copy 2D tracklet history and proecssig for clustering + appearance
			newTrack.costRGB = 0.0;
			std::deque<CTracklet2D*> queueNewlyInsertedTracklet2D;
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				// location in 3D
				newTrack.lastTrackletLocation3D[camIdx] = curTrack->lastTrackletLocation3D[camIdx];

				// appearance
				newTrack.lastRGBFeature[camIdx] = curTrack->lastRGBFeature[camIdx].clone();
				newTrack.timeTrackletEnded[camIdx] = curTrack->timeTrackletEnded[camIdx];

				// tracklet info
				newTrack.tracklet2DIDRecord[camIdx] = curTrack->tracklet2DIDRecord[camIdx];
				if (NULL == newTrack.curTracklet2Ds.get(camIdx)) { continue; }
				hj::Point2D bottomCenter = newTrack.curTracklet2Ds.get(camIdx)->rects.back().bottomCenter();
				double curTrackletSensitivity = vecMatProjectionSensitivity_[camIdx].at<float>((int)bottomCenter.y, (int)bottomCenter.x);
				if (0 == newTrack.tracklet2DIDRecord[camIdx].size() || newTrack.tracklet2DIDRecord[camIdx].back() != newTrack.curTracklet2Ds.get(camIdx)->id)
				{
					newTrack.tracklet2DIDRecord[camIdx].push_back(newTrack.curTracklet2Ds.get(camIdx)->id);
					queueNewlyInsertedTracklet2D.push_back(newTrack.curTracklet2Ds.get(camIdx));

					// when a newly inserted tracklet is not the most front tracklet
					if (newTrack.tracklet2DIDRecord[camIdx].size() > 1)
					{
						int timeGap = _frameIdx - newTrack.timeTrackletEnded[camIdx];

						// tracklet location
						trackValidity &= CheckTrackletConnectivity(
							newTrack.lastTrackletLocation3D[camIdx],
							newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D,
							curTrack->lastTrackletSensitivity[camIdx],
							curTrackletSensitivity,
							timeGap);

						// RGB cost
						newTrack.costRGB += ComputeRGBCost(
							&newTrack.lastRGBFeature[camIdx],
							&newTrack.curTracklet2Ds.get(camIdx)->RGBFeatureHead,
							timeGap);
					}
				}
				newTrack.lastRGBFeature[camIdx] = newTrack.curTracklet2Ds.get(camIdx)->RGBFeatureTail.clone();
				newTrack.lastTrackletLocation3D[camIdx] = newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D;
				newTrack.timeTrackletEnded[camIdx] = _frameIdx;
				newTrack.lastTrackletSensitivity[camIdx] = curTrackletSensitivity;
			}
			if (!trackValidity) { continue; }
			newTrack.costTotal = GetCost(&newTrack);

			// generate track instance
			listTrack3D_.push_back(newTrack);
			nNewTrackID_++;

			// insert to the track tree and related lists
			CTrack3D *branchTrack = &listTrack3D_.back();
			curTrack->tree->tracks.push_back(branchTrack);
			curTrack->childrenTrack.push_back(branchTrack);
			queueTracksInWindow_.push_back(branchTrack);

			// insert to global queue
			queueBranchTracks.push_back(branchTrack);
		}
	}

	// insert branches to the active track list
	size_t numBranches = queueActiveTrack_.size();
	queueActiveTrack_.insert(queueActiveTrack_.end(), queueBranchTracks.begin(), queueBranchTracks.end());

	/////////////////////////////////////////////////////////////////////
	// TEMPORAL BRANCHING
	/////////////////////////////////////////////////////////////////////
	int numTemporalBranch = 0;
	std::sort(queuePausedTrack_.begin(), queuePausedTrack_.end(), hjTrackGTPandLLDescend);
	for (std::deque<CTrack3D*>::iterator trackIter = queuePausedTrack_.begin();
		trackIter != queuePausedTrack_.end();
		trackIter++)
	{
		if (stParam_.bDoBranchCut && stParam_.nMaxNumTracksInOptimization <= numTemporalBranch) { break; }
		CTrack3D *curTrack = *trackIter;

		for (std::deque<CTrack3D*>::iterator seedTrackIter = seedTracks->begin();
			seedTrackIter != seedTracks->end();
			seedTrackIter++)
		{
			CTrack3D *seedTrack = *seedTrackIter;
			std::deque<CReconstruction> queueSeedReconstruction = (*seedTrackIter)->reconstructions;

			// link validation
			unsigned int timeGap = seedTrack->timeStart - curTrack->timeEnd;
			CReconstruction lastMeasurementReconstruction = curTrack->reconstructions[curTrack->duration - 1];
			double curLinkProbability = ComputeLinkProbability(
				lastMeasurementReconstruction.point, seedTrack->reconstructions.front().point, timeGap);
			if (stParam_.dMinLinkingProbability > curLinkProbability) { continue; }

			// generate a new track with branching combination
			CTrack3D newTrack(nNumCameras_);
			newTrack.id = nNewTrackID_;
			newTrack.curTracklet2Ds = seedTrack->curTracklet2Ds;
			newTrack.bActive = true;
			newTrack.bValid = true;
			newTrack.tree = curTrack->tree;
			newTrack.parentTrack = curTrack;
			newTrack.childrenTrack.clear();
			newTrack.timeStart = curTrack->timeStart;
			newTrack.timeEnd = seedTrack->timeEnd;
			newTrack.timeGeneration = _frameIdx;
			newTrack.duration = newTrack.timeEnd - newTrack.timeStart + 1;
			newTrack.bWasBestSolution = true;
			newTrack.bCurrentBestSolution = false;
			newTrack.bNewTrack = true;
			newTrack.GTProb = curTrack->GTProb;
			newTrack.costEnter = curTrack->costEnter;
			newTrack.costReconstruction = curTrack->costReconstruction;
			newTrack.costLink = curTrack->costLink;
			newTrack.costRGB = curTrack->costRGB;
			newTrack.costExit = 0.0;
			newTrack.costTotal = 0.0;

			// interpolation
			std::vector<CReconstruction> interpolatedReconstructions(timeGap, CReconstruction(nNumCameras_));
			std::vector<hj::Point3D> interpolatedPoints(timeGap);
			hj::Point3D delta = 
				(seedTrack->reconstructions.front().point - lastMeasurementReconstruction.point) / (double)timeGap;
			hj::Point3D lastPosition = lastMeasurementReconstruction.point;
			for (int pos = 0; pos < (int)timeGap - 1; pos++)
			{
				lastPosition += delta;
				interpolatedPoints[pos] = lastPosition;
				interpolatedReconstructions[pos].point = lastPosition;
				interpolatedReconstructions[pos].smoothedPoint = lastPosition;
				interpolatedReconstructions[pos].bIsMeasurement = false;
				interpolatedReconstructions[pos].maxError = 0.0;
				interpolatedReconstructions[pos].costReconstruction = 0.0;
				interpolatedReconstructions[pos].costSmoothedPoint = 0.0;
				interpolatedReconstructions[pos].costLink = 0.0;
			}
			interpolatedPoints.back() = seedTrack->reconstructions.front().point;
			interpolatedReconstructions.back() = seedTrack->reconstructions.front();

			// smoothing
			newTrack.reconstructions.insert(newTrack.reconstructions.begin(), curTrack->reconstructions.begin(), curTrack->reconstructions.begin() + curTrack->duration);
			newTrack.reconstructions.insert(newTrack.reconstructions.end(), interpolatedReconstructions.begin(), interpolatedReconstructions.end());
			newTrack.smoother = curTrack->smoother;
			int updateStartPos = newTrack.smoother.Insert(interpolatedPoints);
			bool trackValidity = true;
			double smoothedProbability = 0.0;
			for (int pos = updateStartPos; pos < newTrack.smoother.size(); pos++)
			{
				CReconstruction *curReconstruction = &newTrack.reconstructions[pos];
				if (MIN_SMOOTHING_LENGTH <= newTrack.duration)
				{
					curReconstruction->smoothedPoint = newTrack.smoother.GetResult(pos);
					//newTrack.smoothedTrajectory[pos] = curReconstruction->smoothedPoint;

					// update reconstruction cost
					newTrack.costReconstruction -= curReconstruction->costSmoothedPoint;
					smoothedProbability = ComputeReconstructionProbability(curReconstruction->smoothedPoint, &curReconstruction->rawPoints, &curReconstruction->tracklet2Ds, curReconstruction->maxError);
					if (0.0 == smoothedProbability)
					{
						trackValidity = false;
						break;
					}
					curReconstruction->costSmoothedPoint = -log(smoothedProbability);
					newTrack.costReconstruction += curReconstruction->costSmoothedPoint;

					if (0 == pos) { continue; }
					// update link cost
					newTrack.costLink -= curReconstruction->costLink;

					// TESTING
					smoothedProbability = ComputeLinkProbability(newTrack.reconstructions[pos - 1].smoothedPoint, curReconstruction->smoothedPoint);
					//smoothedProbability = ComputeLinkProbability(newTrack.reconstructions[pos-1].smoothedPoint + newTrack.reconstructions[pos-1].velocity, curReconstruction->smoothedPoint);

					if (0.0 == smoothedProbability)
					{
						trackValidity = false;
						break;
					}
					curReconstruction->costLink = -log(smoothedProbability);
					newTrack.costLink += curReconstruction->costLink;
				}
				else
				{
					curReconstruction->smoothedPoint = curReconstruction->point;
					if (0 == pos) { continue; }
				}

				// update velocity
				curReconstruction->velocity = curReconstruction->smoothedPoint - newTrack.reconstructions[pos - 1].smoothedPoint;
				if (stParam_.dMinMovingSpeed > curReconstruction->velocity.norm_L2()) 
				{ 
					curReconstruction->velocity = hj::Point3D(0.0, 0.0, 0.0); 
				}
			}
			if (!trackValidity) { continue; }

			// tracklet history
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				// location in 3D
				newTrack.lastTrackletLocation3D[camIdx] = curTrack->lastTrackletLocation3D[camIdx];

				// appearance
				newTrack.lastRGBFeature[camIdx] = curTrack->lastRGBFeature[camIdx].clone();
				newTrack.timeTrackletEnded[camIdx] = curTrack->timeTrackletEnded[camIdx];

				// tracklet info
				newTrack.tracklet2DIDRecord[camIdx] = curTrack->tracklet2DIDRecord[camIdx];
				if (NULL == seedTrack->curTracklet2Ds.get(camIdx)) { continue; }
				hj::Point2D bottomCenter = newTrack.curTracklet2Ds.get(camIdx)->rects.back().bottomCenter();
				double curTrackletSensitivity = vecMatProjectionSensitivity_[camIdx].at<float>((int)bottomCenter.y, (int)bottomCenter.x);
				if (0 == newTrack.tracklet2DIDRecord[camIdx].size() || newTrack.tracklet2DIDRecord[camIdx].back() != seedTrack->tracklet2DIDRecord[camIdx].back())
				{
					newTrack.tracklet2DIDRecord[camIdx].push_back(seedTrack->tracklet2DIDRecord[camIdx].back());

					// when a newly inserted tracklet is not the most front tracklet
					if (newTrack.tracklet2DIDRecord[camIdx].size() > 1)
					{
						int timeGap = _frameIdx - newTrack.timeTrackletEnded[camIdx];

						// tracklet location
						trackValidity &= CheckTrackletConnectivity(
							newTrack.lastTrackletLocation3D[camIdx],
							newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D,
							curTrack->lastTrackletSensitivity[camIdx],
							curTrackletSensitivity,
							timeGap);

						// RGB cost
						newTrack.costRGB += ComputeRGBCost(
							&newTrack.lastRGBFeature[camIdx],
							&newTrack.curTracklet2Ds.get(camIdx)->RGBFeatureHead,
							timeGap);
					}
				}
				newTrack.lastRGBFeature[camIdx] = newTrack.curTracklet2Ds.get(camIdx)->RGBFeatureTail.clone();
				newTrack.lastTrackletLocation3D[camIdx] = newTrack.curTracklet2Ds.get(camIdx)->currentLocation3D;
				newTrack.timeTrackletEnded[camIdx] = _frameIdx;
				newTrack.lastTrackletSensitivity[camIdx] = curTrackletSensitivity;
			}
			if (!trackValidity) { continue; }

			newTrack.costTotal = GetCost(&newTrack);

			// insert to the list of track instances
			listTrack3D_.push_back(newTrack);
			nNewTrackID_++;

			// insert to the track tree and related lists
			CTrack3D *branchTrack = &listTrack3D_.back();
			curTrack->tree->tracks.push_back(branchTrack);
			curTrack->childrenTrack.push_back(branchTrack);
			queueActiveTrack_.push_back(branchTrack);
			queueTracksInWindow_.push_back(branchTrack);
			numTemporalBranch++;
		}
	}

	// add seed tracks to active track list
	queueActiveTrack_.insert(queueActiveTrack_.end(), seedTracks->begin(), seedTracks->end());
	queueTracksInWindow_.insert(queueTracksInWindow_.end(), seedTracks->begin(), seedTracks->end());
}

/************************************************************************
Method Name: Track3D_GetWholeCandidateTracks
Description:
-
Input Arguments:
-
Return Values:
- hj_TrackSet: selected tracks
************************************************************************/
hj::TrackSet CAssociator3D::Track3D_GetWholeCandidateTracks(void)
{
	return queueTracksInWindow_;
}

/************************************************************************
Method Name: ComputeEnterProbability
Description:
- calculate entering probability
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeEnterProbability(std::vector<hj::Point3D> &points)
{
	double distanceFromBoundary = -100.0;
	for (int pointIdx = 0; pointIdx < points.size(); pointIdx++)
	{
		double curDistance = GetDistanceFromBoundary(points[pointIdx]);
		if (distanceFromBoundary < curDistance) { distanceFromBoundary = curDistance; }
	}
	if (distanceFromBoundary < 0) { return 1.0; }
	return distanceFromBoundary <= stParam_.dBoundaryDistance ? 
		1.0 : 
		stParam_.dEnterProbMax * 
		exp(-(double)(stParam_.dEnterProbDecayCoef * std::max(0.0, distanceFromBoundary - stParam_.dBoundaryDistance)));
}

/************************************************************************
Method Name: ComputeExitProbability
Description:
- calculate exit probability
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeExitProbability(std::vector<hj::Point3D> &points, int tracklength)
{
	double distanceFromBoundary = -100.0;
	for (int pointIdx = 0; pointIdx < points.size(); pointIdx++)
	{
		double curDistance = GetDistanceFromBoundary(points[pointIdx]);
		if (distanceFromBoundary < curDistance) { distanceFromBoundary = curDistance; }
	}
	if (distanceFromBoundary < 0) { return 1.0; }
	if (distanceFromBoundary < stParam_.dBoundaryDistance) { return stParam_.dExitProbMax; }

	double probability = stParam_.dExitProbMax;
	probability *= exp(-(double)(stParam_.dExitProbDistDecayCoef * 
		std::max(0.0, distanceFromBoundary - (double)stParam_.dBoundaryDistance)));	// decaying with distance
	probability *= exp(-(double)stParam_.dExitProbLengthDecayCoef * 
		std::max(0.0, (double)tracklength - (double)stParam_.nNumFramesForConfirmation));
	return probability;
}

/************************************************************************
Method Name: ComputeLinkProbability
Description:
- calculate linking probability
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeLinkProbability(hj::Point3D &prePoint, hj::Point3D &curPoint, unsigned int timeGap)
{
	double distance = (prePoint - curPoint).norm_L2();
	double maxDistance = stParam_.dMaxMovingSpeed * timeGap;
	return 0.5 * hj::erfc(4.0 * distance / maxDistance - 2.0);
}

/************************************************************************
Method Name: ComputeRGBCost
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeTrackletLinkCost(hj::Point3D preLocation, hj::Point3D curLocation, int timeGap)
{
	if (timeGap > 1) { return 0.0; }
	double norm2 = (preLocation - curLocation).norm_L2();
	return norm2 > stParam_.dCostTrackletLinkingMinDistance ? 
		stParam_.dCostTrackletLinkingCoef * (norm2 - stParam_.dCostTrackletLinkingMinDistance) : 0.0;	
}

/************************************************************************
Method Name: ComputeReconstructionProbability
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeReconstructionProbability(hj::Point3D point, std::vector<hj::Point3D> *rawPoints, CTrackletCombination *trackletCombination, double maxError)
{
	//-------------------------------------------------
	// RECONSTRUCTION PROBABILITY
	//-------------------------------------------------
	double probabilityReconstruction = 0.5;
	if (1 < rawPoints->size())
	{
		double fDistance = 0.0;
		for (int pointIdx = 0; pointIdx < rawPoints->size(); pointIdx++)
		{
			fDistance += (point - (*rawPoints)[pointIdx]).norm_L2();
		}
		fDistance /= (double)rawPoints->size();

		if (0.0 == maxError) { maxError = stParam_.bConsiderSensitivity ? stParam_.dMaxSensitivityError : (double)stParam_.dMaxTrackletDistance / 2.0; }
		if (fDistance > maxError) { return 0.0; }
		probabilityReconstruction = 0.5 * hj::erfc(4.0 * fDistance / maxError - 2.0);
	}

	//-------------------------------------------------
	// DETECTION PROBABILITY
	//-------------------------------------------------
	double fDetectionProbabilityRatio = 1.0;
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		if (!this->CheckVisibility(point, camIdx)) { continue; }
		if (NULL == trackletCombination->get(camIdx))
		{
			// false negative
			fDetectionProbabilityRatio *= stParam_.dFNRate / (1 - stParam_.dFNRate);
			continue;
		}
		// positive
		fDetectionProbabilityRatio *= (1 - stParam_.dFPRate) / stParam_.dFPRate;
	}
	return fDetectionProbabilityRatio * probabilityReconstruction / (1.0 - probabilityReconstruction);
}

/************************************************************************
Method Name: ComputeRGBCost
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::ComputeRGBCost(const cv::Mat *feature1, const cv::Mat *feature2, unsigned int timeGap)
{
	cv::Mat vecDiff = *feature1 - *feature2;
	vecDiff = vecDiff.t() * vecDiff;
	double norm2 = vecDiff.at<double>(0, 0);
	return norm2 > stParam_.dCostRGBMinDistance ? 
		stParam_.dCostRGBCoef * 
		std::exp(-stParam_.dCostRGBDecayCoef * (double)(timeGap - 1.0)) * (norm2 - stParam_.dCostRGBMinDistance) : 0.0;	
}

/************************************************************************
Method Name: CheckIncompatibility
Description:
- check incompatibility between two tracks
Input Arguments:
-
Return Values:
- whether they are incompatible or not
************************************************************************/
bool CAssociator3D::CheckIncompatibility(CTrack3D *track1, CTrack3D *track2, double maxTrackletDistance, double minProximity)
{
	assert(track1->curTracklet2Ds.numCameras() == track2->curTracklet2Ds.numCameras());

	bool bIncompatible = false;
	int  numCameras    = track1->curTracklet2Ds.numCameras();

	// check coupling
	if (track1->tree->id == track2->tree->id)
	{
		bIncompatible = true;
	}
	else
	{
		for (int camIdx = 0; camIdx < numCameras; camIdx++)
		{
			if (0 == track1->tracklet2DIDRecord[camIdx].size() || 0 == track2->tracklet2DIDRecord[camIdx].size()) { continue; }
			if (track1->tracklet2DIDRecord[camIdx][0] > track2->tracklet2DIDRecord[camIdx][0])
			{
				if (track1->tracklet2DIDRecord[camIdx][0] > track2->tracklet2DIDRecord[camIdx].back()) { continue; }
				for (std::deque<unsigned int>::iterator tracklet1Iter = track2->tracklet2DIDRecord[camIdx].begin();
					tracklet1Iter != track2->tracklet2DIDRecord[camIdx].end();
					tracklet1Iter++)
				{
					for (std::deque<unsigned int>::iterator tracklet2Iter = track1->tracklet2DIDRecord[camIdx].begin();
						tracklet2Iter != track1->tracklet2DIDRecord[camIdx].end();
						tracklet2Iter++)
					{
						if (*tracklet1Iter != *tracklet2Iter) { continue; }
						bIncompatible = true;
						break;
					}
					if (bIncompatible) { break; }
				}
			}
			else if (track1->tracklet2DIDRecord[camIdx][0] < track2->tracklet2DIDRecord[camIdx][0])
			{
				if (track1->tracklet2DIDRecord[camIdx].back() < track2->tracklet2DIDRecord[camIdx][0]) { continue; }
				for (std::deque<unsigned int>::iterator tracklet1Iter = track1->tracklet2DIDRecord[camIdx].begin();
					tracklet1Iter != track1->tracklet2DIDRecord[camIdx].end();
					tracklet1Iter++)
				{
					for (std::deque<unsigned int>::iterator tracklet2Iter = track2->tracklet2DIDRecord[camIdx].begin();
						tracklet2Iter != track2->tracklet2DIDRecord[camIdx].end();
						tracklet2Iter++)
					{
						if (*tracklet1Iter != *tracklet2Iter) { continue; }
						bIncompatible = true;
						break;
					}
					if (bIncompatible) { break; }
				}
			}
			else
			{
				bIncompatible = true;
				break;
			}
			if (bIncompatible) { break; }
		}
	}

	// check proximity and crossing
	if (!bIncompatible)
	{
		// check overlapping
		if (track1->timeEnd < track2->timeStart || track2->timeEnd < track1->timeStart) { return bIncompatible; }

		unsigned int timeStart = std::max(track1->timeStart, track2->timeStart);
		unsigned int timeEnd = std::min(track1->timeEnd, track2->timeEnd);
		unsigned int track1ReconIdx = timeStart - track1->timeStart;
		unsigned int track2ReconIdx = timeStart - track2->timeStart;
		unsigned int overlapLength = timeEnd - timeStart + 1;
		hj::Point3D *reconLocation1;
		hj::Point3D *reconLocation2;
		double distanceBetweenTracks = 0.0;
		for (unsigned int reconIdx = 0; reconIdx < overlapLength; reconIdx++)
		{
			reconLocation1 = &(track1->reconstructions[track1ReconIdx].point);
			reconLocation2 = &(track2->reconstructions[track2ReconIdx].point);
			track1ReconIdx++;
			track2ReconIdx++;
			// proximity
			distanceBetweenTracks = (*reconLocation1 - *reconLocation2).norm_L2();
			if (distanceBetweenTracks > maxTrackletDistance) { continue; }
			if (distanceBetweenTracks < minProximity) { return true; }
			if (reconIdx < overlapLength - 1)
			{
				// crossing
				if (hj::IsLineSegmentIntersect(
						hj::Line3D(*reconLocation1, track1->reconstructions[track1ReconIdx].point), 
						hj::Line3D(*reconLocation2, track2->reconstructions[track2ReconIdx].point)))
				{
					return true;
				}
			}
		}
	}
	return bIncompatible;
}

/************************************************************************
Method Name: CheckIncompatibility
Description:
- check incompatibility between two tracks
Input Arguments:
-
Return Values:
- whether they are incompatible or not
************************************************************************/
bool CAssociator3D::CheckIncompatibility(CTrackletCombination &combi1, CTrackletCombination &combi2)
{
	assert(combi1.numCameras() == combi2.numCameras());
	bool bIncompatible = false;
	for (int camIdx = 0; camIdx < combi1.numCameras(); camIdx++)
	{
		if (NULL == combi1.get(camIdx) || NULL == combi2.get(camIdx))
		{
			continue;
		}
		if (combi1.get(camIdx) == combi2.get(camIdx))
		{
			bIncompatible = true;
			break;
		}
	}

	return bIncompatible;
}

/************************************************************************
Method Name: GetRGBFeature
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
cv::Mat CAssociator3D::GetRGBFeature(const cv::Mat *patch, int numBins)
{
	std::vector<cv::Mat> channelImages(3);
	cv::split(*patch, channelImages);
	double numPixels = patch->rows * patch->cols;
	cv::Mat featureB = hj::histogram(channelImages[0], numBins);
	cv::Mat featureG = hj::histogram(channelImages[1], numBins);
	cv::Mat featureR = hj::histogram(channelImages[2], numBins);

	cv::vconcat(featureR, featureG, featureR);
	cv::vconcat(featureR, featureB, featureR);
	featureR = featureR / numPixels;

	return featureR;
}

/************************************************************************
Method Name: GetCost
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
double CAssociator3D::GetCost(CTrack3D *track)
{
	track->costReconstruction = 0.0;
	track->costLink = 0.0;
	for (int reconIdx = 0; reconIdx < track->reconstructions.size(); reconIdx++)
	{
		track->costReconstruction += track->reconstructions[reconIdx].costSmoothedPoint;
		track->costLink += track->reconstructions[reconIdx].costLink;
	}
	double resultCost = track->costEnter + track->costReconstruction + track->costLink + track->costRGB + track->costExit;
	//double resultCost = track->costEnter + track->costReconstruction + track->costLink + track->costExit;
	return resultCost;
}

/************************************************************************
Method Name: Hypothesis_UpdateHypotheses
Description:
-
Input Arguments:
- none
Return Values:
- none
************************************************************************/
void CAssociator3D::Hypothesis_UpdateHypotheses(hj::HypothesisSet &inoutUpdatedHypotheses, hj::TrackSet *newSeedTracks)
{
	// DEBUG
	nCountTrackInOptimization_ = 0;
	nCountUCTrackInOptimization_ = 0;

	vectorSolvingInfo_.clear();
	vectorSolvingInfo_.reserve(inoutUpdatedHypotheses.size());
	//---------------------------------------------------------
	// ADD NEW TRACKS TO RELATED TRACK QUEUE
	//---------------------------------------------------------
	hj::HypothesisSet newHypothesesSet;
	for (hj::HypothesisSet::iterator hypothesisIter = inoutUpdatedHypotheses.begin();
		hypothesisIter != inoutUpdatedHypotheses.end() && newHypothesesSet.size() < stParam_.nKBestSize;
		hypothesisIter++)
	{
		if (stParam_.nKBestSize <= newHypothesesSet.size()) { break; }
		//if (!(*hypothesisIter).bValid){ continue; }
		stGlobalHypothesis newGlobalHypothesis = (*hypothesisIter);
		stHypothesisSolvingInfo newInfo;

		// TODO...: cost <->solving relationship
		newInfo.rank = (int)newHypothesesSet.size();

		// add new branch tracks
		hj::TrackSet newRelatedTracks = (*hypothesisIter).relatedTracks;
		for (hj::TrackSet::iterator trackIter = (*hypothesisIter).relatedTracks.begin();
			trackIter != (*hypothesisIter).relatedTracks.end();
			trackIter++)
		{
			for (hj::TrackSet::iterator childTrackIter = (*trackIter)->childrenTrack.begin();
				childTrackIter != (*trackIter)->childrenTrack.end();
				childTrackIter++)
			{
				if ((*childTrackIter)->bNewTrack) { newRelatedTracks.push_back(*childTrackIter); }
			}
		}
		std::sort(newRelatedTracks.begin(), newRelatedTracks.end(), hjTrackGTPandLLDescend);
		if (newRelatedTracks.size() <= stParam_.nMaxNumTracksInOptimization)
		{
			newGlobalHypothesis.relatedTracks = newRelatedTracks;
		}
		else
		{
			newGlobalHypothesis.relatedTracks.insert(newGlobalHypothesis.relatedTracks.end(), newRelatedTracks.begin(), newRelatedTracks.begin() + stParam_.nMaxNumTracksInOptimization);
		}

		// add new seed tracks
		newGlobalHypothesis.relatedTracks.insert(newGlobalHypothesis.relatedTracks.end(), newSeedTracks->begin(), newSeedTracks->end());

		// save hypothesis
		newHypothesesSet.push_back(newGlobalHypothesis);

		// DEBUG
		nCountTrackInOptimization_ += (int)(*hypothesisIter).relatedTracks.size();
		for (hj::TrackSet::iterator trackIter = (*hypothesisIter).relatedTracks.begin();
			trackIter != (*hypothesisIter).relatedTracks.end();
			trackIter++)
		{
			if (!(*trackIter)->tree->bConfirmed) { nCountUCTrackInOptimization_++; }
		}
	}
	inoutUpdatedHypotheses = newHypothesesSet;
}

/************************************************************************
Method Name: Hypothesis_Formation
Description:
-
Input Arguments:
- none
Return Values:
- none
************************************************************************/
void CAssociator3D::Hypothesis_Formation(hj::HypothesisSet &outBranchHypotheses, hj::HypothesisSet *existingHypotheses)
{
	//---------------------------------------------------------
	// PROPAGATE GLOBAL HYPOTHESES
	//---------------------------------------------------------
	outBranchHypotheses.clear();
	if (0 == existingHypotheses->size())
	{
		this->Hypothesis_BranchHypotheses(outBranchHypotheses, &stParam_, &this->Track3D_GetWholeCandidateTracks());
	}
	else
	{
		// solve K times
		std::vector<HANDLE> vecHSolvingThreads(existingHypotheses->size());
		std::vector<stSolvingThreadParams> vecParams(existingHypotheses->size());
		std::vector<HypothesisSet> vecOutBranchHypotheses(existingHypotheses->size());

		std::sort(existingHypotheses->begin(), existingHypotheses->end(), hjSolutionLogLikelihoodDescendComparator);
		for (size_t hypothesisIdx = 0; hypothesisIdx < std::min(stParam_.nKBestSize, (int)existingHypotheses->size()); hypothesisIdx++)
		{
			//this->Hypothesis_BranchHypotheses(
			//	outBranchHypotheses,
			//	&stParam_,
			//	&(*existingHypotheses)[hypothesisIdx].relatedTracks,
			//	&(*existingHypotheses)[hypothesisIdx].selectedTracks);
			vecParams[hypothesisIdx].threadID              = (int)hypothesisIdx;
			vecParams[hypothesisIdx].outBranchHypotheses   = &vecOutBranchHypotheses[hypothesisIdx];
			vecParams[hypothesisIdx].params                = &stParam_;
			vecParams[hypothesisIdx].tracks                = &(*existingHypotheses)[hypothesisIdx].relatedTracks;
			vecParams[hypothesisIdx].initialSolutionTracks = &(*existingHypotheses)[hypothesisIdx].selectedTracks;
			
			vecHSolvingThreads[hypothesisIdx] = (HANDLE)_beginthreadex(
				0, 0, &AssociationSolvingWork, &vecParams[hypothesisIdx], 0, 0);
		}
		
		WaitForMultipleObjects((DWORD)vecHSolvingThreads.size(), vecHSolvingThreads.data(), true, INFINITE);
		for (size_t threadIdx = 0; threadIdx < vecHSolvingThreads.size(); threadIdx++)
		{
			// close thread handles
			CloseHandle(vecHSolvingThreads[threadIdx]);

			// merge results
			outBranchHypotheses.insert(outBranchHypotheses.end(),
				vecOutBranchHypotheses[threadIdx].begin(),
				vecOutBranchHypotheses[threadIdx].end());
		}
		vecHSolvingThreads.clear();
		vecParams.clear();
	}
	if (0 == outBranchHypotheses.size()) { return; }

	//---------------------------------------------------------
	// CALCULATE PROBABILITIES
	//---------------------------------------------------------
	double hypothesesTotalLoglikelihood = 0.0;
	for (size_t solutionIdx = 0; solutionIdx < outBranchHypotheses.size(); solutionIdx++)
	{
		hypothesesTotalLoglikelihood += outBranchHypotheses[solutionIdx].logLikelihood;
	}

	// probability of global hypothesis and global track probability
	for (size_t hypothesisIdx = 0; hypothesisIdx < outBranchHypotheses.size(); hypothesisIdx++)
	{
		outBranchHypotheses[hypothesisIdx].probability = outBranchHypotheses[hypothesisIdx].logLikelihood / hypothesesTotalLoglikelihood;
		for (size_t trackIdx = 0; trackIdx < outBranchHypotheses[hypothesisIdx].selectedTracks.size(); trackIdx++)
		{
			outBranchHypotheses[hypothesisIdx].selectedTracks[trackIdx]->GTProb += outBranchHypotheses[hypothesisIdx].probability;
		}
	}

	// sorting and confirmaiton of track tree
	std::sort(outBranchHypotheses.begin(), outBranchHypotheses.end(), hjSolutionLogLikelihoodDescendComparator);
}

/************************************************************************
Method Name: Hypothesis_BranchHypotheses
Description:
- construct graph with inserted tracks ans solve that graph
Input Arguments:
- tracks: related tracks
- initialselectedTracks: tracks in previous solution (or initial solution)
Return Values:
- none
************************************************************************/
#define hj_GRAPH_SOLUTION_DUPLICATION_RESOLUTION (1.0E-5)
void CAssociator3D::Hypothesis_BranchHypotheses(
	hj::HypothesisSet &outBranchHypotheses,
	hj::stParamAssociate3D *params,
	hj::TrackSet *tracks,
	hj::TrackSet *initialselectedTracks)
{
	if (0 == tracks->size()) { return; }

	//---------------------------------------------------------
	// GRAPH CONSTRUCTION
	//---------------------------------------------------------	
	hj::Graph *curGraph = new hj::Graph();
	hj::VertexSet vertexInGraph;
	hj::VertexSet vertexInInitialSolution;
	std::deque<hj::VertexSet> initialSolutionSet;
	std::deque<CTrack3D*> queueTracksInOptimization;

	for (size_t trackIdx = 0; trackIdx < tracks->size(); trackIdx++)
	{
		CTrack3D *curTrack = (*tracks)[trackIdx];
		if (!curTrack->bValid) { continue; }

		// cost	
		curTrack->costTotal = hj::CAssociator3D::GetCost(curTrack);

		if (0.0 < curTrack->costTotal) { continue; }

		// for miss penalty
		//double unpenaltyCost = 0.0;
		//for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		//{
		//	if (NULL == curTrack->curTracklet2Ds.get(camIdx)) { continue; }
		//	unpenaltyCost += log(stParam_.dFPRate/(1 - stParam_.dFPRate));
		//}

		// add vertex
		//curTrack->loglikelihood = -(curTrack->costTotal + unpenaltyCost);
		curTrack->loglikelihood = -curTrack->costTotal;
		hj::GraphVertex *curVertex = curGraph->AddVertex(curTrack->loglikelihood);
		vertexInGraph.push_back(curVertex);
		queueTracksInOptimization.push_back(curTrack);

		// find incompatibility constraints
		for (unsigned int compTrackIdx = 0; compTrackIdx < queueTracksInOptimization.size() - 1; compTrackIdx++)
		{
			if (hj::CAssociator3D::CheckIncompatibility(
				curTrack, 
				queueTracksInOptimization[compTrackIdx],
				params->dMaxMovingSpeed * 2.0,
				params->dMinTargetProximity))
			{
				continue; 
			}
			curGraph->AddEdge(vertexInGraph[compTrackIdx], vertexInGraph.back()); // for MCP
		}

		// initial solution related
		if (NULL == initialselectedTracks) { continue; }
		hj::TrackSet::iterator findIter = 
			std::find(initialselectedTracks->begin(), initialselectedTracks->end(), curTrack);
		if (initialselectedTracks->end() != findIter)
		{
			vertexInInitialSolution.push_back(curVertex);
		}
	}

	// validate initial solution
	if (0 < vertexInInitialSolution.size())
	{
		initialSolutionSet.push_back(vertexInInitialSolution);
	}

	if (0 == curGraph->Size())
	{
		delete curGraph;
		return;
	}

	//---------------------------------------------------------
	// GRAPH SOLVING
	//---------------------------------------------------------
	hj::CGraphSolver graphSolver;
	graphSolver.Initialize(curGraph, HJ_GRAPH_SOLVER_BLS);
	graphSolver.SetInitialSolutions(initialSolutionSet);
	hj::stGraphSolvingResult *curResult = graphSolver.Solve();

	//---------------------------------------------------------
	// SOLUTION CONVERSION AND DUPLICATION CHECK
	//---------------------------------------------------------
	for (size_t solutionIdx = 0; solutionIdx < curResult->vecSolutions.size(); solutionIdx++)
	{
		stGlobalHypothesis newHypothesis;
		newHypothesis.selectedTracks.clear();
		newHypothesis.relatedTracks = *tracks;
		newHypothesis.bValid = true;
		newHypothesis.logLikelihood = curResult->vecSolutions[solutionIdx].second;
		for (hj::VertexSet::iterator vertexIter = curResult->vecSolutions[solutionIdx].first.begin();
			vertexIter != curResult->vecSolutions[solutionIdx].first.end();
			vertexIter++)
		{
			newHypothesis.selectedTracks.push_back(queueTracksInOptimization[(*vertexIter)->id_]);
		}
		std::sort(newHypothesis.selectedTracks.begin(), newHypothesis.selectedTracks.end(), hjTrackIDAscendComparator);

		// duplication check
		bool bDuplicated = false;
		for (size_t compSolutionIdx = 0; compSolutionIdx < outBranchHypotheses.size(); compSolutionIdx++)
		{
			if (std::abs(outBranchHypotheses[compSolutionIdx].logLikelihood - newHypothesis.logLikelihood) > hj_GRAPH_SOLUTION_DUPLICATION_RESOLUTION)
			{
				continue;
			}

			if (outBranchHypotheses[compSolutionIdx].selectedTracks == newHypothesis.selectedTracks)
			{
				bDuplicated = true;
				break;
			}
		}
		if (bDuplicated) { continue; }

		outBranchHypotheses.push_back(newHypothesis);
	}

	//---------------------------------------------------------
	// WRAP-UP
	//---------------------------------------------------------
	delete curGraph;
}

/************************************************************************
Method Name: Hypothesis_PruningNScanBack
Description:
-
Input Arguments:
-
Return Values:
- none
************************************************************************/
void CAssociator3D::Hypothesis_PruningNScanBack(
	unsigned int nCurrentFrameIdx,
	unsigned int N,
	hj::TrackSet *tracksInWindow,
	std::deque<stGlobalHypothesis> *ptQueueHypothesis)
{
	if (NULL == ptQueueHypothesis) { return; }
	if (0 >= (*ptQueueHypothesis).size()) { return; }

	// DEBUG
	//PrintHypotheses(*ptQueueHypothesis, "data/hypotheses.txt", nCurrentFrameIdx);

	// track tree branch pruning
	int nTimeBranchPruning = (int)nCurrentFrameIdx - (int)N;
	hj::TrackSet *tracksInBestSolution = &(*ptQueueHypothesis).front().selectedTracks;
	CTrack3D *brachSeedTrack = NULL;
	for (int trackIdx = 0; trackIdx < tracksInBestSolution->size(); trackIdx++)
	{
		(*tracksInBestSolution)[trackIdx]->bCurrentBestSolution = true;

		// left unconfirmed tracks
		if ((*tracksInBestSolution)[trackIdx]->tree->timeGeneration + stParam_.nNumFramesForConfirmation > nCurrentFrameIdx) { continue; }

		// pruning
		brachSeedTrack = CTrackTree::FindOldestTrackInBranch((*tracksInBestSolution)[trackIdx], nTimeBranchPruning);
		if (NULL == brachSeedTrack->parentTrack) { continue; }
		for (int childIdx = 0; childIdx < brachSeedTrack->parentTrack->childrenTrack.size(); childIdx++)
		{
			if (brachSeedTrack->parentTrack->childrenTrack[childIdx] == brachSeedTrack) { continue; }
			CTrackTree::SetValidityFlagInTrackBranch(brachSeedTrack->parentTrack->childrenTrack[childIdx], false);
		}
	}
}

/************************************************************************
Method Name: Hypothesis_PruningKBest
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
void CAssociator3D::Hypothesis_PruningTrackWithGTP(unsigned int nCurrentFrameIdx, unsigned int nNumMaximumTrack, hj::TrackSet *tracksInWindow, std::deque<CTrackTree*> *queueActiveTrackTree)
{
	int numTrackLeft = 0;
	int numUCTrackLeft = 0;

	// GTP pruning with confirmed track tree
	for (int treeIdx = 0; treeIdx < queuePtActiveTrees_.size(); treeIdx++)
	{
		hj::TrackSet sortedTrackQueue = queuePtActiveTrees_[treeIdx]->tracks;
		std::sort(sortedTrackQueue.begin(), sortedTrackQueue.end(), hjTrackGTPandLLDescend);
		for (int trackIdx = stParam_.nNumFramesForConfirmation; trackIdx < sortedTrackQueue.size(); trackIdx++)
		{
			sortedTrackQueue[trackIdx]->bValid = false;
		}
	}

	// simple GTP pruning
	std::sort(tracksInWindow->begin(), tracksInWindow->end(), hjTrackGTPandLLDescend);
	for (hj::TrackSet::iterator trackIter = tracksInWindow->begin();
		trackIter != tracksInWindow->end();
		trackIter++)
	{
		if (!(*trackIter)->bValid) { continue; }
		if (!(*trackIter)->tree->bConfirmed)
		{
			numUCTrackLeft++;
			continue;
		}
		if (numTrackLeft < stParam_.nMaxNumTracksInOptimization && (*trackIter)->GTProb > 0.0)
		{
			numTrackLeft++;
			continue;
		}
		(*trackIter)->bValid = false;
	}

	// pruning unconfirmed tracks
	for (int treeIdx = 0; treeIdx < queuePtUnconfirmedTrees_.size(); treeIdx++)
	{
		hj::TrackSet sortedTrackQueue = queuePtUnconfirmedTrees_[treeIdx]->tracks;
		std::sort(sortedTrackQueue.begin(), sortedTrackQueue.end(), hjTrackGTPandLLDescend);
		for (int trackIdx = stParam_.nMaxNumTracksInUnconfirmedTree; trackIdx < sortedTrackQueue.size(); trackIdx++)
		{
			sortedTrackQueue[trackIdx]->bValid = false;
		}
	}
}

/************************************************************************
Method Name: Hypothesis_RefreshHypotheses
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
void CAssociator3D::Hypothesis_RefreshHypotheses(hj::HypothesisSet &inoutUpdatedHypotheses)
{
	// gathering tracks in unconfirmed track trees
	hj::TrackSet unconfirmedTracks;
	for (int treeIdx = 0; treeIdx < queuePtUnconfirmedTrees_.size(); treeIdx++)
	{
		for (hj::TrackSet::iterator trackIter = queuePtUnconfirmedTrees_[treeIdx]->tracks.begin();
			trackIter != queuePtUnconfirmedTrees_[treeIdx]->tracks.end();
			trackIter++)
		{
			if (!(*trackIter)->bValid) { continue; }
			unconfirmedTracks.push_back(*trackIter);
		}
	}

	// update hyptheses
	std::deque<stGlobalHypothesis> existingHypothesis = inoutUpdatedHypotheses;
	inoutUpdatedHypotheses.clear();
	stGlobalHypothesis *curHypotheis = NULL;
	bool bCurHypothesisValid = true;
	for (int hypothesisIdx = 0; hypothesisIdx < existingHypothesis.size(); hypothesisIdx++)
	{
		curHypotheis = &existingHypothesis[hypothesisIdx];

		// check validity
		bCurHypothesisValid = true;
		for (int trackIdx = 0; trackIdx < curHypotheis->selectedTracks.size(); trackIdx++)
		{
			if (curHypotheis->selectedTracks[trackIdx]->bValid) { continue; }
			bCurHypothesisValid = false;
			break;
		}
		if (!bCurHypothesisValid) { continue; }

		// update related track list
		hj::TrackSet newRelatedTracks = curHypotheis->selectedTracks;
		newRelatedTracks.insert(newRelatedTracks.end(), unconfirmedTracks.begin(), unconfirmedTracks.end());
		curHypotheis->relatedTracks = newRelatedTracks;

		// save valid hypothesis
		inoutUpdatedHypotheses.push_back(*curHypotheis);
	}
}

/************************************************************************
Method Name: ResultWithCandidateTracks
Description:
- generate track result with candidate tracks
Input Arguments:
-
Return Values:
- CTrack3DResult: a struct of result
************************************************************************/
CTrack3DResult CAssociator3D::ResultWithTracks(hj::TrackSet *trackSet, unsigned int nFrameIdx, double fProcessingTime)
{
	CTrack3DResult result3D;
	result3D.frameIdx = nFrameIdx;
	result3D.processingTime = fProcessingTime;
	if (0 == trackSet->size()) { return result3D; }

	for (hj::TrackSet::iterator trackIter = trackSet->begin();
		trackIter != trackSet->end();
		trackIter++)
	{
		if ((*trackIter)->timeEnd < nFrameIdx) { continue; }

		CTrack3D *curTrack = *trackIter;
		CObject3DInfo newObject(nNumCameras_);

		//---------------------------------------------------------
		// ID
		//---------------------------------------------------------
		newObject.id = curTrack->id;
		if (!stParam_.bShowTreeID)
		{
			// ID for tracking result
			bool bIDNotFound = true;
			for (size_t pairIdx = 0; pairIdx < queuePairTreeIDToTargetID_.size(); pairIdx++)
			{
				if (queuePairTreeIDToTargetID_[pairIdx].first == curTrack->tree->id)
				{
					newObject.id = queuePairTreeIDToTargetID_[pairIdx].second;
					bIDNotFound = false;
					break;
				}
			}
			if (bIDNotFound)
			{
				newObject.id = nNewTargetID_++;
				queuePairTreeIDToTargetID_.push_back(std::make_pair(curTrack->tree->id, newObject.id));
			}
		}

		//---------------------------------------------------------
		// TRAJECTORY AND BOX
		//---------------------------------------------------------
		int numPoint = 0;
		//int deferredLength = std::max((int)curTrack->timeEnd - (int)nFrameIdx, 0); 
		int deferredLength = curTrack->timeEnd - nFrameIdx;
		if (deferredLength > curTrack->reconstructions.size()) { continue; }
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++) { newObject.bVisibleInViews[camIdx] = false; }
		for (std::deque<CReconstruction>::reverse_iterator pointIter = curTrack->reconstructions.rbegin() + deferredLength;
			pointIter != curTrack->reconstructions.rend();
			pointIter++)
		{
			numPoint++;
			newObject.recentPoints.push_back((*pointIter).smoothedPoint);
			if (stParam_.nResultTrajectoryLength < numPoint) { break; }

			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				// trajectory
				hj::Point2D reprejectedPoint(0.0, 0.0);
				if (CheckVisibility((*pointIter).smoothedPoint, camIdx, &reprejectedPoint))
				{
					newObject.recentPoint2Ds[camIdx].push_back(reprejectedPoint);
				}

				// box
				if (1 < numPoint) { continue; }
				if (NULL == curTrack->curTracklet2Ds.get(camIdx))
				{
					hj::Rect curRect(0.0, 0.0, 0.0, 0.0);
					hj::Point3D topCenterInWorld((*pointIter).smoothedPoint);
					topCenterInWorld.z = 1700;
					hj::Point2D bottomCenterInImage = this->WorldToImage((*pointIter).smoothedPoint, camIdx);
					hj::Point2D topCenterInImage = this->WorldToImage(topCenterInWorld, camIdx);

					double height = (topCenterInImage - bottomCenterInImage).norm_L2();
					curRect.x = bottomCenterInImage.x - height * 2.5 / 17.0;
					curRect.y = topCenterInImage.y;
					curRect.w = height * 5.0 / 17.0;
					curRect.h = bottomCenterInImage.y - topCenterInImage.y;

					newObject.rectInViews[camIdx] = curRect;
					newObject.bVisibleInViews[camIdx] = false;
				}
				else
				{
					newObject.rectInViews[camIdx] = curTrack->curTracklet2Ds.get(camIdx)->rects.back();
					newObject.bVisibleInViews[camIdx] = true;
				}

				// 3D box
				//newObject.point3DBox[camIdx] = this->GetHuman3DBox((*pointIter).point, 500, camIdx);
			}
		}

		//---------------------------------------------------------
		// 2D DETECTIONS
		//---------------------------------------------------------
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			CTracklet2D *curTracklet = curTrack->reconstructions.back().tracklet2Ds.get(camIdx);
			hj::Point3D detectionLocation(0.0, 0.0, 0.0);
			if (NULL != curTracklet) { detectionLocation = this->ImageToWorld(curTracklet->rects.back().bottomCenter(), 0.0, camIdx); }
			newObject.curDetectionPosition.push_back(detectionLocation);
		}

		result3D.object3DInfos.push_back(newObject);
	}

	return result3D;
}

/************************************************************************
Method Name: PrintTracks
Description:
- print out information of tracks
Input Arguments:
- queueTracks: seletect tracks for print out
- strFilePathAndName: output file path and name
- bAppend: option for appending
Return Values:
- void
************************************************************************/
void CAssociator3D::PrintTracks(hj::TrackSet &queueTracks, char *strFilePathAndName, bool bAppend)
{
	if (0 == queueTracks.size()) { return; }
	try
	{
		FILE *fp;
		if (bAppend)
		{
			fopen_s(&fp, strFilePathAndName, "a");
		}
		else
		{
			fopen_s(&fp, strFilePathAndName, "w");
			fprintf_s(fp, "numCamera:%d\n", (int)nNumCameras_);
			fprintf_s(fp, "numTracks:%d\n", (int)queueTracks.size());
		}

		for (std::deque<CTrack3D*>::iterator trackIter = queueTracks.begin();
			trackIter != queueTracks.end();
			trackIter++)
		{
			CTrack3D *curTrack = *trackIter;
			fprintf_s(fp, "{\n\tid:%d\n\ttreeID:%d\n", (int)curTrack->id, (int)curTrack->tree->id);
			fprintf_s(fp, "\tnumReconstructions:%d\n\ttimeStart:%d\n\ttimeEnd:%d\n\ttimeGeneration:%d\n", (int)curTrack->reconstructions.size(), (int)curTrack->timeStart, (int)curTrack->timeEnd, (int)curTrack->timeGeneration);
			fprintf_s(fp, "\ttrackleIDs:\n\t{\n");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				if (0 == curTrack->tracklet2DIDRecord[camIdx].size())
				{
					fprintf(fp, "\t\tnumTracklet:0,[]\n");
					continue;
				}
				fprintf_s(fp, "\t\tnumTracklet:%d,[", (int)curTrack->tracklet2DIDRecord[camIdx].size());
				for (unsigned int trackletIDIdx = 0; trackletIDIdx < curTrack->tracklet2DIDRecord[camIdx].size(); trackletIDIdx++)
				{
					fprintf_s(fp, "%d", (int)curTrack->tracklet2DIDRecord[camIdx][trackletIDIdx]);
					if (curTrack->tracklet2DIDRecord[camIdx].size() - 1 != trackletIDIdx)
					{
						fprintf_s(fp, ",");
					}
					else
					{
						fprintf_s(fp, "]\n");
					}
				}
			}
			fprintf(fp, "\t}\n");
			fprintf_s(fp, "\ttotalCost:%e\n", curTrack->costTotal);
			fprintf_s(fp, "\treconstructionCost:%e\n", curTrack->costReconstruction);
			fprintf_s(fp, "\tlinkCost:%e\n", curTrack->costLink);
			fprintf_s(fp, "\tinitCost:%e\n", curTrack->costEnter);
			fprintf_s(fp, "\ttermCost:%e\n", curTrack->costExit);
			fprintf_s(fp, "\tRGBCost:%e\n", curTrack->costRGB);
			fprintf_s(fp, "\treconstructions:\n\t{\n");
			for (std::deque<CReconstruction>::iterator pointIter = curTrack->reconstructions.begin();
				pointIter != curTrack->reconstructions.end();
				pointIter++)
			{
				fprintf_s(fp, "\t\t%d:(%f,%f,%f),%e,%e,%e\n", (int)(*pointIter).bIsMeasurement, (*pointIter).point.x, (*pointIter).point.y, (*pointIter).point.z, (*pointIter).costReconstruction, (*pointIter).costSmoothedPoint, (*pointIter).costLink);
			}
			fprintf_s(fp, "\t}\n}\n");
		}

		fclose(fp);
	}
	catch (int nError)
	{
		printf("[ERROR](PrintTracks) cannot open file! error code %d\n", nError);
		return;
	}
}

/************************************************************************
Method Name: PrintHypotheses
Description:
- print out information of tracks
Input Arguments:
- queueTracks: seletect tracks for print out
- strFilePathAndName: output file path and name
- bAppend: option for appending
Return Values:
- void
************************************************************************/
void CAssociator3D::PrintHypotheses(hj::HypothesisSet &queueHypotheses, char *strFilePathAndName, unsigned int frameIdx)
{
	if (0 == queueHypotheses.size()) { return; }
	try
	{
		FILE *fp;
		std::string trackIDList;
		fopen_s(&fp, strFilePathAndName, "w");
		fprintf_s(fp, "frameIndex:%d\n", (int)frameIdx);
		fprintf_s(fp, "numHypotheses:%d\n", (int)queueHypotheses.size());
		int rank = 0;
		for (hj::HypothesisSet::iterator hypothesisIter = queueHypotheses.begin();
			hypothesisIter != queueHypotheses.end();
			hypothesisIter++, rank++)
		{
			fprintf_s(fp, "{\n\trank:%d\n", rank);
			fprintf_s(fp, "\tselectedTracks:%d,", (int)(*hypothesisIter).selectedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).selectedTracks) + "\n";
			fprintf_s(fp, trackIDList.c_str());
			fprintf_s(fp, "\n\t{\n");
			for (int trackIdx = 0; trackIdx < (*hypothesisIter).selectedTracks.size(); trackIdx++)
			{
				CTrack3D *curTrack = (*hypothesisIter).selectedTracks[trackIdx];
				fprintf_s(fp, "\t\t{id:%d,Total:%e,Recon:%e,Link:%e,Enr:%e,Ex:%e,RGB:%e}\n",
					(int)curTrack->id,
					curTrack->costTotal,
					curTrack->costReconstruction,
					curTrack->costLink,
					curTrack->costEnter,
					curTrack->costExit,
					curTrack->costRGB);
			}
			fprintf_s(fp, "\t}\n");

			fprintf_s(fp, "\trelatedTracks:%d,", (int)(*hypothesisIter).relatedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).relatedTracks) + "\n";
			fprintf_s(fp, trackIDList.c_str());

			fprintf_s(fp, "\tlogLikelihood:%e\n", (*hypothesisIter).logLikelihood);
			fprintf_s(fp, "\tprobability:%e\n", (*hypothesisIter).probability);
			fprintf_s(fp, "\tbValid:%d\n\t}\n", (int)(*hypothesisIter).bValid);
		}
		fclose(fp);
	}
	catch (int nError)
	{
		printf("[ERROR](PrintHypotheses) cannot open file! error code %d\n", nError);
		return;
	}
}

/************************************************************************
Method Name: PrintCurrentTrackTrees
Description:
- print out tree structure of current tracks
Input Arguments:
- strFilePathAndName: output file path and name
Return Values:
- void
************************************************************************/
void CAssociator3D::PrintCurrentTrackTrees(const char *strFilePath)
{
	try
	{
		FILE *fp;
		fopen_s(&fp, strFilePath, "w");
		int numTree = 0;

		// structure
		std::deque<stTrackInTreeInfo> queueNodes;
		for (std::list<CTrackTree>::iterator treeIter = listTrackTree_.begin();
			treeIter != listTrackTree_.end();
			treeIter++)
		{
			if (0 == (*treeIter).tracks.size()) { continue; }
			numTree++;

			CTrack3D* curTrack = (*treeIter).tracks.front();
			if (!curTrack->bValid) { continue; }
			stTrackInTreeInfo newInfo;
			newInfo.id = curTrack->id;
			newInfo.parentNode = 0;
			newInfo.timeGenerated = curTrack->timeGeneration;
			newInfo.GTP = (float)curTrack->GTProb;
			queueNodes.push_back(newInfo);
			CTrackTree::MakeTreeNodesWithChildren(curTrack->childrenTrack, (int)queueNodes.size(), queueNodes);
		}

		// write trees
		fprintf(fp, "numTrees:%d\n", numTree);
		fprintf(fp, "nodeLength:%d,{", (int)queueNodes.size());
		for (std::deque<stTrackInTreeInfo>::iterator idxIter = queueNodes.begin();
			idxIter != queueNodes.end();
			idxIter++)
		{
			fprintf(fp, "(%d,%d,%d,%f)", (*idxIter).id, (*idxIter).parentNode, (*idxIter).timeGenerated, (*idxIter).GTP);
			if (idxIter < queueNodes.end() - 1) { fprintf(fp, ","); }
		}
		fprintf(fp, "}");

		fclose(fp);
	}
	catch (int nError)
	{
		printf("[ERROR](PrintCurrentTrackTrees) cannot open file! error code %d\n", nError);
		return;
	}
}


/************************************************************************
Method Name: PrintResult
Description:
-
Input Arguments:
- none
Return Values:
- none
************************************************************************/
void CAssociator3D::PrintResult(const char *strFilepath, std::deque<CTrack3DResult> *queueResults)
{
	FILE *fp;
	try
	{
		fopen_s(&fp, strFilepath, "w");

		for (std::deque<CTrack3DResult>::iterator resultIter = queueResults->begin();
			resultIter != queueResults->end();
			resultIter++)
		{
			// first, count the visible objects
			int numObject = 0;
			for (std::vector<CObject3DInfo>::iterator objectIter = (*resultIter).object3DInfos.begin();
				objectIter != (*resultIter).object3DInfos.end();
				objectIter++)
			{
				if (0 == (*objectIter).recentPoints.size()) { continue; }
				numObject++;
			}

			fprintf_s(fp, "{\n\tframeIndex:%d\n\tnumObjects:%d\n", (int)(*resultIter).frameIdx, numObject);
			for (std::vector<CObject3DInfo>::iterator objectIter = (*resultIter).object3DInfos.begin();
				objectIter != (*resultIter).object3DInfos.end();
				objectIter++)
			{
				if (0 == (*objectIter).recentPoints.size()) { continue; }
				hj::Point3D curPoint = (*objectIter).recentPoints.front();
				fprintf_s(fp, "\t{id:%d,position:(%f,%f,%f)}\n", (int)(*objectIter).id, (float)curPoint.x, (float)curPoint.y, (float)curPoint.z);
			}
			fprintf_s(fp, "}\n");
		}

		fclose(fp);
	}
	catch (int nError)
	{
		printf("[ERROR](PrintResult) cannot open file! error code %d\n", nError);
		return;
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
void CAssociator3D::VisualizeResult(const unsigned int _frameIdx)
{
	if (!bVisualizeResult_) { return; }

	double rescale = 1.0;
	if (vecMatCurrentFrames_[0].cols > 960)
	{
		rescale = 960.0 / vecMatCurrentFrames_[0].cols;
	}

	/* boxes */
	for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
	{
		std::vector<hj::CObject3DInfo> vecObjects = queueTrackingResult_.back().object3DInfos;

		// resize
		if (1.0 > rescale)
		{
			cv::Size targetSize((int)((double)vecMatCurrentFrames_[camIdx].cols * rescale),
				(int)((double)vecMatCurrentFrames_[camIdx].rows * rescale));
			cv::resize(vecMatCurrentFrames_[camIdx], vecMatCurrentFrames_[camIdx], targetSize);
		}

		for (size_t trajectoryIdx = 0; trajectoryIdx < vecObjects.size(); trajectoryIdx++) 
		{
			hj::CObject3DInfo *curTrajectory = &vecObjects[trajectoryIdx];

			// drwa box with ID
			if (curTrajectory->bVisibleInViews[camIdx])
				hj::DrawBoxWithID(
					vecMatCurrentFrames_[camIdx], 
					curTrajectory->rectInViews[camIdx] * rescale,
					curTrajectory->id, 
					1, 
					1, 
					&vecColors_);
			else
				hj::DrawBoxWithID(
					vecMatCurrentFrames_[camIdx], 
					curTrajectory->rectInViews[camIdx] * rescale,
					curTrajectory->id, 
					0, 
					1, 
					&vecColors_);
			
			// draw recent locations
			std::vector<hj::Point2D> vecPoints = curTrajectory->recentPoint2Ds[camIdx];
			for (int pIdx = 0; pIdx < vecPoints.size(); pIdx++)
			{
				vecPoints[pIdx] *= rescale;
			}
			hj::DrawLine(
				vecMatCurrentFrames_[camIdx], 
				vecPoints,
				curTrajectory->id, 
				1, 
				&vecColors_);
		}

		// DEBUG
		//hj::Point2D ptOrigin2D = WorldToImage(hj::Point3D(.0, .0, .0), camIdx);
		//hj::Point2D ptX_500    = WorldToImage(hj::Point3D( 500.0,    .0, .0), camIdx); ptX_500 *= rescale;
		//hj::Point2D ptY_500    = WorldToImage(hj::Point3D(    .0, 500.0, .0), camIdx); ptY_500 *= rescale;
		//hj::Point2D ptX_N500   = WorldToImage(hj::Point3D(-500.0,    .0, .0), camIdx); ptX_N500 *= rescale;
		//hj::Point2D ptY_N500   = WorldToImage(hj::Point3D(    .0,-500.0, .0), camIdx); ptY_N500 *= rescale;
		// DEBUG: PETS with the first pedestrian
		//hj::Point2D ptOrigin2D = WorldToImage(hj::Point3D(-3937, -7478, .0), camIdx);
		//hj::Point2D ptX_500 = WorldToImage(hj::Point3D(-3937+2000.0, -7478, .0), camIdx); ptX_500 *= rescale;
		//hj::Point2D ptY_500 = WorldToImage(hj::Point3D(-3937, -7478+2000.0, .0), camIdx); ptY_500 *= rescale;
		//hj::Point2D ptX_N500 = WorldToImage(hj::Point3D(-3937-2000.0, -7478, .0), camIdx); ptX_N500 *= rescale;
		//hj::Point2D ptY_N500 = WorldToImage(hj::Point3D(-3937, -7478-2000.0, .0), camIdx); ptY_N500 *= rescale;
		//cv::line(
		//	vecMatCurrentFrames_[camIdx],
		//	ptX_500.cv(), ptX_N500.cv(),			
		//	cv::Scalar(255.0, .0, .0));
		//cv::line(
		//	vecMatCurrentFrames_[camIdx],
		//	ptY_500.cv(), ptY_N500.cv(),
		//	cv::Scalar(0.0, .0, 255.0));

		// DEBUG
/*		hj::Point2D GT1 = WorldToImage(hj::Point3D(-5946.9,   -6908.1, .0), camIdx); GT1 *= rescale;
		hj::Point2D GT2 = WorldToImage(hj::Point3D(-9367.1,   -6318.5, .0), camIdx); GT2 *= rescale;
		hj::Point2D GT3 = WorldToImage(hj::Point3D(-10100.0, -10400.0, .0), camIdx); GT3 *= rescale;
		hj::Point2D RES1 = WorldToImage(hj::Point3D(-4275.7,  -7568.4, .0), camIdx); RES1 *= rescale;
		hj::Point2D RES2 = WorldToImage(hj::Point3D(-9128.4, -12900.0, .0), camIdx); RES2 *= rescale;
		hj::Point2D RES3 = WorldToImage(hj::Point3D(-11500.0, -5644.3, .0), camIdx); RES3 *= rescale;
		cv::circle(vecMatCurrentFrames_[camIdx], GT1.cv(), 5, cv::Scalar(255, 0, 0), 2);
		cv::circle(vecMatCurrentFrames_[camIdx], GT2.cv(), 5, cv::Scalar(255, 0, 0), 2);
		cv::circle(vecMatCurrentFrames_[camIdx], GT3.cv(), 5, cv::Scalar(255, 0, 0), 2);
		cv::circle(vecMatCurrentFrames_[camIdx], RES1.cv(), 5, cv::Scalar(0, 0, 255), 2);
		cv::circle(vecMatCurrentFrames_[camIdx], RES2.cv(), 5, cv::Scalar(0, 0, 255), 2);
		cv::circle(vecMatCurrentFrames_[camIdx], RES3.cv(), 5, cv::Scalar(0, 0, 255), 2);	*/	
		
	}	
	matTrackingResult_ = hj::MakeMatTile(&vecMatCurrentFrames_, 2, 2);

	// writing frame info
	char strFrameInfo[100];
	sprintf_s(strFrameInfo, "Frame: %04d", _frameIdx);
	cv::rectangle(matTrackingResult_, cv::Rect(5, 2, 145, 22), cv::Scalar(0, 0, 0), CV_FILLED);
	cv::putText(matTrackingResult_, strFrameInfo, cv::Point(6, 20), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(255, 255, 255));
	// writing processing time
	sprintf_s(strFrameInfo, "processing time: %f", dCurrentProcessingTime_);
	cv::rectangle(matTrackingResult_, cv::Rect(5, 30, 175, 18), cv::Scalar(0, 0, 0), CV_FILLED);
	cv::putText(matTrackingResult_, strFrameInfo, cv::Point(6, 42), cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255, 255, 255));

	cv::namedWindow(strVisWindowName_);
	//cv::moveWindow(strVisWindowName_, 10, (int)(vecMatCurrentFrames_[0].rows) + 10);
	cv::moveWindow(strVisWindowName_, 10, 10);
	cv::imshow(strVisWindowName_, matTrackingResult_);
	cv::waitKey(1);
	matTrackingResult_.release();
}


/************************************************************************
 Method Name: SaveSnapshot
 Description:
	-
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CAssociator3D::SaveSnapshot(const char *strFilepath)
{
	FILE *fpTracklet2D, *fpTrack3D, *fpHypothesis, *fpResult, *fpInfo;
	char strFilename[128] = "";
	std::string trackIDList;
	try
	{
		sprintf_s(strFilename, "%ssnapshot_3D_tracklet.txt", strFilepath);
		fopen_s(&fpTracklet2D, strFilename, "w");
		fprintf_s(fpTracklet2D, "numCamera:%d\n", nNumCameras_);
		fprintf_s(fpTracklet2D, "frameIndex:%d\n\n", (int)nCurrentFrameIdx_);

		//---------------------------------------------------------
		// 2D TRACKLET RELATED
		//---------------------------------------------------------
		fprintf_s(fpTracklet2D, "vecTracklet2DSet:\n{\n");
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			fprintf_s(fpTracklet2D, "\t{\n");
			fprintf_s(fpTracklet2D, "\t\ttracklets:%d,\n\t\t{\n", (int)vecTracklet2DSet_[camIdx].tracklets.size());
			for (std::list<CTracklet2D>::iterator trackletIter = vecTracklet2DSet_[camIdx].tracklets.begin();
				trackletIter != vecTracklet2DSet_[camIdx].tracklets.end();
				trackletIter++)
			{
				CTracklet2D *curTracklet = &(*trackletIter);
				fprintf_s(fpTracklet2D, "\t\t\t{\n");
				fprintf_s(fpTracklet2D, "\t\t\t\tid:%d\n", (int)curTracklet->id);
				fprintf_s(fpTracklet2D, "\t\t\t\tcamIdx:%d\n", (int)curTracklet->camIdx);
				fprintf_s(fpTracklet2D, "\t\t\t\tbActivated:%d\n", (int)curTracklet->bActivated);
				fprintf_s(fpTracklet2D, "\t\t\t\ttimeStart:%d\n", (int)curTracklet->timeStart);
				fprintf_s(fpTracklet2D, "\t\t\t\ttimeEnd:%d\n", (int)curTracklet->timeEnd);
				fprintf_s(fpTracklet2D, "\t\t\t\tduration:%d\n", (int)curTracklet->duration);

				// rects
				fprintf_s(fpTracklet2D, "\t\t\t\trects:%d,\n\t\t\t\t{\n", (int)curTracklet->rects.size());
				for (int rectIdx = 0; rectIdx < curTracklet->rects.size(); rectIdx++)
				{
					hj::Rect curRect = curTracklet->rects[rectIdx];
					fprintf_s(fpTracklet2D, "\t\t\t\t\t(%f,%f,%f,%f)\n", (float)curRect.x, (float)curRect.y, (float)curRect.w, (float)curRect.h);
				}
				fprintf_s(fpTracklet2D, "\t\t\t\t}\n");

				// backprojectionLines
				fprintf_s(fpTracklet2D, "\t\t\t\tbackprojectionLines:%d,\n\t\t\t\t{\n", (int)curTracklet->backprojectionLines.size());
				for (int lineIdx = 0; lineIdx < curTracklet->backprojectionLines.size(); lineIdx++)
				{
					hj::Point3D pt1 = curTracklet->backprojectionLines[lineIdx].first;
					hj::Point3D pt2 = curTracklet->backprojectionLines[lineIdx].second;
					fprintf_s(fpTracklet2D, "\t\t\t\t\t(%f,%f,%f,%f,%f,%f)\n", (float)pt1.x, (float)pt1.y, (float)pt1.z, (float)pt2.x, (float)pt2.y, (float)pt2.z);
				}
				fprintf_s(fpTracklet2D, "\t\t\t\t}\n");

				// appearance
				fprintf_s(fpTracklet2D, "\t\t\t\tRGBFeatureHead:%d,(", curTracklet->RGBFeatureHead.rows);
				for (int idx = 0; idx < curTracklet->RGBFeatureHead.rows; idx++)
				{
					fprintf_s(fpTracklet2D, "%f", curTracklet->RGBFeatureHead.at<double>(idx, 0));
					if (idx < curTracklet->RGBFeatureHead.rows - 1) { fprintf_s(fpTracklet2D, ","); }
				}
				fprintf_s(fpTracklet2D, ")\n");
				fprintf_s(fpTracklet2D, "\t\t\t\tRGBFeatureTail:%d,(", curTracklet->RGBFeatureTail.rows);
				for (int idx = 0; idx < curTracklet->RGBFeatureTail.rows; idx++)
				{
					fprintf_s(fpTracklet2D, "%f", curTracklet->RGBFeatureTail.at<double>(idx, 0));
					if (idx < curTracklet->RGBFeatureTail.rows - 1) { fprintf_s(fpTracklet2D, ","); }
				}
				fprintf_s(fpTracklet2D, ")\n");

				// location in 3D
				fprintf_s(fpTracklet2D, "\t\t\t\tcurrentLocation3D:(%f,%f,%f)\n", curTracklet->currentLocation3D.x, curTracklet->currentLocation3D.y, curTracklet->currentLocation3D.z);

				// bAssociableNewMeasurement
				for (int compCamIdx = 0; compCamIdx < nNumCameras_; compCamIdx++)
				{
					fprintf_s(fpTracklet2D, "\t\t\t\tbAssociableNewMeasurement[%d]:%d,{", compCamIdx, (int)curTracklet->bAssociableNewMeasurement[compCamIdx].size());
					for (int flagIdx = 0; flagIdx < curTracklet->bAssociableNewMeasurement[compCamIdx].size(); flagIdx++)
					{
						fprintf_s(fpTracklet2D, "%d", (int)curTracklet->bAssociableNewMeasurement[compCamIdx][flagIdx]);
						if (flagIdx < curTracklet->bAssociableNewMeasurement[compCamIdx].size() - 1) { fprintf_s(fpTracklet2D, ","); }
					}
					fprintf_s(fpTracklet2D, "}\n");
				}
				fprintf_s(fpTracklet2D, "\t\t\t}\n");
			}

			fprintf_s(fpTracklet2D, "\t\tactiveTracklets:%d,{", (int)vecTracklet2DSet_[camIdx].activeTracklets.size());
			for (std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[camIdx].activeTracklets.begin();
				trackletIter != vecTracklet2DSet_[camIdx].activeTracklets.end();
				trackletIter++)
			{
				fprintf_s(fpTracklet2D, "%d", (int)(*trackletIter)->id);
				if (trackletIter != vecTracklet2DSet_[camIdx].activeTracklets.end() - 1) { fprintf_s(fpTracklet2D, ","); }
			}
			fprintf_s(fpTracklet2D, "}\n\t\tnewMeasurements:%d,{", (int)vecTracklet2DSet_[camIdx].newMeasurements.size());
			for (std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[camIdx].newMeasurements.begin();
				trackletIter != vecTracklet2DSet_[camIdx].newMeasurements.end();
				trackletIter++)
			{
				fprintf_s(fpTracklet2D, "%d", (int)(*trackletIter)->id);
				if (trackletIter != vecTracklet2DSet_[camIdx].newMeasurements.end() - 1) { fprintf_s(fpTracklet2D, ","); }
			}
			fprintf_s(fpTracklet2D, "}\n\t}\n");
		}
		fprintf_s(fpTracklet2D, "}\n");
		fprintf_s(fpTracklet2D, "nNumTotalActive2DTracklet:%d\n\n", (int)nNumTotalActive2DTracklet_);
		fclose(fpTracklet2D);

		//---------------------------------------------------------
		// 3D TRACK RELATED
		//---------------------------------------------------------
		sprintf_s(strFilename, "%ssnapshot_3D_track.txt", strFilepath);
		fopen_s(&fpTrack3D, strFilename, "w");
		fprintf_s(fpTrack3D, "numCamera:%d\n", nNumCameras_);
		fprintf_s(fpTrack3D, "frameIndex:%d\n\n", (int)nCurrentFrameIdx_);
		fprintf_s(fpTrack3D, "bReceiveNewMeasurement:%d\n", (int)bReceiveNewMeasurement_);
		fprintf_s(fpTrack3D, "bInitiationPenaltyFree:%d\n", (int)bInitiationPenaltyFree_);
		fprintf_s(fpTrack3D, "nNewTrackID:%d\n", (int)nNewTrackID_);
		fprintf_s(fpTrack3D, "nNewTreeID:%d\n", (int)nNewTreeID_);

		// track
		fprintf_s(fpTrack3D, "listTrack3D:%d,\n{\n", (int)listTrack3D_.size());
		for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
			trackIter != listTrack3D_.end();
			trackIter++)
		{
			CTrack3D *curTrack = &(*trackIter);
			fprintf_s(fpTrack3D, "\t{\n\t\tid:%d\n", (int)curTrack->id);

			// curTracklet2Ds
			fprintf_s(fpTrack3D, "\t\tcurTracklet2Ds:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				if (NULL == curTrack->curTracklet2Ds.get(camIdx))
				{
					fprintf_s(fpTrack3D, "-1");
				}
				else
				{
					fprintf_s(fpTrack3D, "%d", curTrack->curTracklet2Ds.get(camIdx)->id);
				}

				if (camIdx < nNumCameras_ - 1)
				{
					fprintf_s(fpTrack3D, ",");
				}
			}
			fprintf_s(fpTrack3D, "}\n");

			// tracklet2DIDs
			fprintf_s(fpTrack3D, "\t\ttrackleIDs:\n\t\t{\n");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				if (0 == curTrack->tracklet2DIDRecord[camIdx].size())
				{
					fprintf_s(fpTrack3D, "\t\t\tnumTracklet:0,{}\n");
					continue;
				}
				fprintf_s(fpTrack3D, "\t\t\tnumTracklet:%d,{", (int)curTrack->tracklet2DIDRecord[camIdx].size());
				for (unsigned int trackletIDIdx = 0; trackletIDIdx < curTrack->tracklet2DIDRecord[camIdx].size(); trackletIDIdx++)
				{
					fprintf_s(fpTrack3D, "%d", (int)curTrack->tracklet2DIDRecord[camIdx][trackletIDIdx]);
					if (curTrack->tracklet2DIDRecord[camIdx].size() - 1 != trackletIDIdx) { fprintf_s(fpTrack3D, ","); }
					else { fprintf_s(fpTrack3D, "}\n"); }
				}
			}
			fprintf_s(fpTrack3D, "\t\t}\n");

			fprintf_s(fpTrack3D, "\t\tbActive:%d\n", (int)curTrack->bActive);
			fprintf_s(fpTrack3D, "\t\tbValid:%d\n", (int)curTrack->bValid);
			fprintf_s(fpTrack3D, "\t\ttreeID:%d\n", (int)curTrack->tree->id);
			if (NULL == curTrack->parentTrack) { fprintf_s(fpTrack3D, "\t\tparentTrackID:-1\n"); }
			else { fprintf_s(fpTrack3D, "\t\tparentTrackID:%d\n", (int)curTrack->parentTrack->id); }

			// children tracks
			fprintf_s(fpTrack3D, "\t\tchildrenTrack:%d,", (int)curTrack->childrenTrack.size());
			trackIDList = hj::MakeTrackIDList(&curTrack->childrenTrack) + "\n";
			fprintf_s(fpTrack3D, trackIDList.c_str());

			// temporal information
			fprintf_s(fpTrack3D, "\t\ttimeStart:%d\n", (int)curTrack->timeStart);
			fprintf_s(fpTrack3D, "\t\ttimeEnd:%d\n", (int)curTrack->timeEnd);
			fprintf_s(fpTrack3D, "\t\ttimeGeneration:%d\n", (int)curTrack->timeGeneration);
			fprintf_s(fpTrack3D, "\t\tduration:%d\n", (int)curTrack->duration);

			// reconstrcution related
			fprintf(fpTrack3D, "\t\treconstructions:%d,\n\t\t{\n", (int)curTrack->reconstructions.size());
			for (std::deque<CReconstruction>::iterator pointIter = curTrack->reconstructions.begin();
				pointIter != curTrack->reconstructions.end();
				pointIter++)
			{
				fprintf_s(fpTrack3D, "\t\t\t%d,point:(%f,%f,%f),smoothedPoint:(%f,%f,%f),velocity:(%f,%f,%f),maxError:%e,costReconstruction:%e,costSmoothedPoint:%e,costLink:%e,", (int)(*pointIter).bIsMeasurement,
					(*pointIter).point.x, (*pointIter).point.y, (*pointIter).point.z,
					(*pointIter).smoothedPoint.x, (*pointIter).smoothedPoint.y, (*pointIter).smoothedPoint.z,
					(*pointIter).velocity.x, (*pointIter).velocity.y, (*pointIter).velocity.z,
					(*pointIter).maxError, (*pointIter).costReconstruction, (*pointIter).costSmoothedPoint, (*pointIter).costLink);
				fprintf_s(fpTrack3D, "tracklet2Ds:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					if (NULL == (*pointIter).tracklet2Ds.get(camIdx)) { fprintf_s(fpTrack3D, "-1"); }
					else { fprintf_s(fpTrack3D, "%d", (*pointIter).tracklet2Ds.get(camIdx)->id); }
					if (camIdx < nNumCameras_ - 1) { fprintf_s(fpTrack3D, ","); }
				}
				fprintf_s(fpTrack3D, "},");
				fprintf_s(fpTrack3D, "rawPoints:%d,{", (int)(*pointIter).rawPoints.size());
				for (int pointIdx = 0; pointIdx < (*pointIter).rawPoints.size(); pointIdx++)
				{
					fprintf_s(fpTrack3D, "(%f,%f,%f)", (*pointIter).rawPoints[pointIdx].x, (*pointIter).rawPoints[pointIdx].y, (*pointIter).rawPoints[pointIdx].z);
					if (pointIdx < (*pointIter).rawPoints.size() - 1) { fprintf_s(fpTrack3D, ","); }
				}
				fprintf_s(fpTrack3D, "}\n");
			}
			fprintf_s(fpTrack3D, "\t\t}\n");

			// point smoother
			fprintf_s(fpTrack3D, "\t\tsmoother:{span:%d,degree:%d}\n", SGS_DEFAULT_SPAN, SGS_DEFAULT_DEGREE);

			// cost
			fprintf_s(fpTrack3D, "\t\tcostTotal:%e\n", curTrack->costTotal);
			fprintf_s(fpTrack3D, "\t\tcostReconstruction:%e\n", curTrack->costReconstruction);
			fprintf_s(fpTrack3D, "\t\tcostLink:%e\n", curTrack->costLink);
			fprintf_s(fpTrack3D, "\t\tcostEnter:%e\n", curTrack->costEnter);
			fprintf_s(fpTrack3D, "\t\tcostExit:%e\n", curTrack->costExit);
			fprintf_s(fpTrack3D, "\t\tcostRGB:%e\n", curTrack->costRGB);

			// loglikelihood
			fprintf_s(fpTrack3D, "\t\tloglikelihood:%e\n", curTrack->loglikelihood);

			// GTP
			fprintf_s(fpTrack3D, "\t\tGTProb:%e\n", curTrack->GTProb);
			fprintf_s(fpTrack3D, "\t\tBranchGTProb:%e\n", curTrack->BranchGTProb);
			fprintf_s(fpTrack3D, "\t\tbWasBestSolution:%d\n", (int)curTrack->bWasBestSolution);
			fprintf_s(fpTrack3D, "\t\tbCurrentBestSolution:%d\n", (int)curTrack->bCurrentBestSolution);

			// HO-MHT
			fprintf_s(fpTrack3D, "\t\tbNewTrack:%d\n", (int)curTrack->bNewTrack);

			// tracklet related
			fprintf_s(fpTrack3D, "\t\ttimeTrackletEnded:(");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fprintf_s(fpTrack3D, "%d", (int)curTrack->timeTrackletEnded[camIdx]);
				if (camIdx < nNumCameras_ - 1) { fprintf_s(fpTrack3D, ","); }
			}
			fprintf_s(fpTrack3D, ")\n");
			fprintf_s(fpTrack3D, "\t\tlastTrackletLocation3D:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fprintf_s(fpTrack3D, "(%f,%f,%f)",
					curTrack->lastTrackletLocation3D[camIdx].x,
					curTrack->lastTrackletLocation3D[camIdx].y,
					curTrack->lastTrackletLocation3D[camIdx].z);
			}
			fprintf_s(fpTrack3D, "}\n");
			fprintf_s(fpTrack3D, "\t\tlastTrackletSensitivity:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fprintf_s(fpTrack3D, "%f", curTrack->lastTrackletSensitivity[camIdx]);
				if (camIdx < nNumCameras_ - 1) { fprintf_s(fpTrack3D, ","); }
			}
			fprintf_s(fpTrack3D, "}\n");
			fprintf_s(fpTrack3D, "\t\tlastRGBFeature:\n\t\t{\n");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fprintf_s(fpTrack3D, "\t\t\tlastRGBFeature[%d]:%d,(", camIdx, curTrack->lastRGBFeature[camIdx].rows);
				for (int idx = 0; idx < curTrack->lastRGBFeature[camIdx].rows; idx++)
				{
					fprintf_s(fpTrack3D, "%f", curTrack->lastRGBFeature[camIdx].at<double>(idx, 0));
					if (idx < curTrack->lastRGBFeature[camIdx].rows - 1) { fprintf_s(fpTrack3D, ","); }
				}
				fprintf_s(fpTrack3D, ")\n");
			}
			fprintf_s(fpTrack3D, "\t\t}\n");

			// termination
			fprintf_s(fpTrack3D, "\t\tbnumOutpoint:%d\n", (int)curTrack->numOutpoint);

			fprintf_s(fpTrack3D, "\t}\n");
		}
		fprintf_s(fpTrack3D, "}\n");

		// track tree
		fprintf_s(fpTrack3D, "listTrackTree:%d,\n{\n", (int)listTrackTree_.size());
		for (std::list<CTrackTree>::iterator treeIter = listTrackTree_.begin();
			treeIter != listTrackTree_.end();
			treeIter++)
		{
			CTrackTree *curTree = &(*treeIter);
			fprintf_s(fpTrack3D, "\t{\n\t\tid:%d\n", (int)curTree->id);
			fprintf_s(fpTrack3D, "\t\ttimeGeneration:%d\n", (int)curTree->timeGeneration);
			fprintf_s(fpTrack3D, "\t\tbValid:%d\n", (int)curTree->bValid);
			fprintf_s(fpTrack3D, "\t\tbConfirmed:%d\n", (int)curTree->bConfirmed);
			fprintf_s(fpTrack3D, "\t\ttracks:%d,", (int)curTree->tracks.size());
			trackIDList = hj::MakeTrackIDList(&curTree->tracks) + "\n";
			fprintf_s(fpTrack3D, trackIDList.c_str());
			fprintf_s(fpTrack3D, "\t}\n");
		}
		fprintf_s(fpTrack3D, "}\n");

		// new seed tracks
		fprintf_s(fpTrack3D, "queueNewSeedTracks:%d,", (int)queueNewSeedTracks_.size());
		trackIDList = hj::MakeTrackIDList(&queueNewSeedTracks_) + "\n";
		fprintf_s(fpTrack3D, trackIDList.c_str());

		// active tracks
		fprintf_s(fpTrack3D, "queueActiveTrack:%d,", (int)queueActiveTrack_.size());
		trackIDList = hj::MakeTrackIDList(&queueActiveTrack_) + "\n";
		fprintf_s(fpTrack3D, trackIDList.c_str());

		// puased tracks
		fprintf_s(fpTrack3D, "queuePausedTrack:%d,", (int)queuePausedTrack_.size());
		trackIDList = hj::MakeTrackIDList(&queuePausedTrack_) + "\n";
		fprintf_s(fpTrack3D, trackIDList.c_str());

		// tracks in window
		fprintf_s(fpTrack3D, "queueTracksInWindow:%d,", (int)queueTracksInWindow_.size());
		trackIDList = hj::MakeTrackIDList(&queueTracksInWindow_) + "\n";
		fprintf_s(fpTrack3D, trackIDList.c_str());

		// tracks in the best solution
		fprintf_s(fpTrack3D, "queueTracksInBestSolution:%d,", (int)queueTracksInBestSolution_.size());
		trackIDList = hj::MakeTrackIDList(&queueTracksInBestSolution_) + "\n";
		fprintf_s(fpTrack3D, trackIDList.c_str());

		// active trees
		fprintf_s(fpTrack3D, "queuePtActiveTrees:%d,{", (int)queuePtActiveTrees_.size());
		for (int treeIdx = 0; treeIdx < queuePtActiveTrees_.size(); treeIdx++)
		{
			fprintf_s(fpTrack3D, "%d", queuePtActiveTrees_[treeIdx]->id);
			if (treeIdx < queuePtActiveTrees_.size() - 1) { fprintf_s(fpTrack3D, ","); }
		}
		fprintf_s(fpTrack3D, "}\n");

		// unconfirmed trees
		fprintf_s(fpTrack3D, "queuePtUnconfirmedTrees:%d,{", (int)queuePtUnconfirmedTrees_.size());
		for (int treeIdx = 0; treeIdx < queuePtUnconfirmedTrees_.size(); treeIdx++)
		{
			fprintf_s(fpTrack3D, "%d", queuePtUnconfirmedTrees_[treeIdx]->id);
			if (treeIdx < queuePtUnconfirmedTrees_.size() - 1) { fprintf_s(fpTrack3D, ","); }
		}
		fprintf_s(fpTrack3D, "}\n");
		fclose(fpTrack3D);

		//---------------------------------------------------------
		// RESULTS
		//---------------------------------------------------------
		sprintf_s(strFilename, "%ssnapshot_3D_result.txt", strFilepath);
		fopen_s(&fpResult, strFilename, "w");
		fprintf_s(fpResult, "numCamera:%d\n", nNumCameras_);
		fprintf_s(fpResult, "frameIndex:%d\n\n", (int)nCurrentFrameIdx_);

		// instance result
		fprintf_s(fpResult, "queueTrackingResult:%d,\n{\n", (int)queueTrackingResult_.size());
		for (int resultIdx = 0; resultIdx < queueTrackingResult_.size(); resultIdx++)
		{
			CTrack3DResult *curResult = &queueTrackingResult_[resultIdx];
			fprintf_s(fpResult, "\t{\n\t\tframeIdx:%d\n", curResult->frameIdx);
			fprintf_s(fpResult, "\t\tprocessingTime:%f\n", curResult->processingTime);

			// object info
			fprintf_s(fpResult, "\t\tobjectInfo:%d,\n\t\t{\n", (int)curResult->object3DInfos.size());
			for (int objIdx = 0; objIdx < curResult->object3DInfos.size(); objIdx++)
			{
				CObject3DInfo *curObject = &curResult->object3DInfos[objIdx];
				fprintf_s(fpResult, "\t\t\t{\n\t\t\t\tid:%d\n", curObject->id);

				// points
				fprintf_s(fpResult, "\t\t\t\trecentPoints:%d,{", (int)curObject->recentPoints.size());
				for (int pointIdx = 0; pointIdx < curObject->recentPoints.size(); pointIdx++)
				{
					hj::Point3D curPoint = curObject->recentPoints[pointIdx];
					fprintf_s(fpResult, "(%f,%f,%f)", (float)curPoint.x, (float)curPoint.y, (float)curPoint.z);
					if (pointIdx < curObject->recentPoints.size() - 1) { fprintf_s(fpResult, ","); }
				}
				fprintf_s(fpResult, "}\n");

				// 2D points
				fprintf_s(fpResult, "\t\t\t\trecentPoint2Ds:\n\t\t\t\t{\n");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fprintf_s(fpResult, "\t\t\t\t\tcam%d:%d,{", camIdx, (int)curObject->recentPoint2Ds[camIdx].size());
					for (int pointIdx = 0; pointIdx < curObject->recentPoint2Ds[camIdx].size(); pointIdx++)
					{
						hj::Point2D curPoint = curObject->recentPoint2Ds[camIdx][pointIdx];
						fprintf_s(fpResult, "(%f,%f)", (float)curPoint.x, (float)curPoint.y);
						if (pointIdx < curObject->recentPoint2Ds[camIdx].size() - 1) { fprintf_s(fpResult, ","); }
					}
					fprintf_s(fpResult, "}\n");
				}
				fprintf_s(fpResult, "\t\t\t\t}\n");

				// 3D box points in each view
				fprintf_s(fpResult, "\t\t\t\tpoint3DBox:\n\t\t\t\t{\n");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fprintf_s(fpResult, "\t\t\t\t\tcam%d:%d,{", camIdx, (int)curObject->point3DBox[camIdx].size());
					for (int pointIdx = 0; pointIdx < curObject->point3DBox[camIdx].size(); pointIdx++)
					{
						hj::Point2D curPoint = curObject->point3DBox[camIdx][pointIdx];
						fprintf_s(fpResult, "(%f,%f)", (float)curPoint.x, (float)curPoint.y);
						if (pointIdx < curObject->recentPoint2Ds[camIdx].size() - 1) { fprintf_s(fpResult, ","); }
					}
					fprintf_s(fpResult, "}\n");
				}
				fprintf_s(fpResult, "\t\t\t\t}\n");

				// rects
				fprintf_s(fpResult, "\t\t\t\trectInViews:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					hj::Rect curRect = curObject->rectInViews[camIdx];
					fprintf_s(fpResult, "(%f,%f,%f,%f)", (float)curRect.x, (float)curRect.y, (float)curRect.w, (float)curRect.h);
					if (camIdx < nNumCameras_ - 1) { fprintf_s(fpResult, ","); }
				}
				fprintf_s(fpResult, "}\n");

				// visibility
				fprintf_s(fpResult, "\t\t\t\tbVisibleInViews:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fprintf_s(fpResult, "%d", (int)curObject->bVisibleInViews[camIdx]);
					if (camIdx < nNumCameras_ - 1) { fprintf_s(fpResult, ","); }
				}
				fprintf_s(fpResult, "}\n\t\t\t}\n");
			}
			fprintf_s(fpResult, "\t\t}\n\t}\n");
		}
		fprintf_s(fpResult, "}\n");
		fclose(fpResult);

		//---------------------------------------------------------
		// HYPOTHESES
		//---------------------------------------------------------
		sprintf_s(strFilename, "%ssnapshot_3D_hypotheses.txt", strFilepath);
		fopen_s(&fpHypothesis, strFilename, "w");
		fprintf_s(fpHypothesis, "numCamera:%d\n", nNumCameras_);
		fprintf_s(fpHypothesis, "frameIndex:%d\n\n", (int)nCurrentFrameIdx_);

		// queuePrevGlobalHypotheses
		fprintf_s(fpHypothesis, "queuePrevGlobalHypotheses:%d,\n{\n", (int)queuePrevGlobalHypotheses_.size());
		for (hj::HypothesisSet::iterator hypothesisIter = queuePrevGlobalHypotheses_.begin();
			hypothesisIter != queuePrevGlobalHypotheses_.end();
			hypothesisIter++)
		{
			fprintf_s(fpHypothesis, "\t{\n\t\tselectedTracks:%d,", (int)(*hypothesisIter).selectedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).selectedTracks) + "\n";
			fprintf_s(fpHypothesis, trackIDList.c_str());

			fprintf_s(fpHypothesis, "\t\trelatedTracks:%d,", (int)(*hypothesisIter).relatedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).relatedTracks) + "\n";
			fprintf_s(fpHypothesis, trackIDList.c_str());

			fprintf_s(fpHypothesis, "\t\tlogLikelihood:%e\n", (*hypothesisIter).logLikelihood);
			fprintf_s(fpHypothesis, "\t\tprobability:%e\n", (*hypothesisIter).probability);
			fprintf_s(fpHypothesis, "\t\tbValid:%d\n\t}\n", (int)(*hypothesisIter).bValid);
		}
		fprintf_s(fpHypothesis, "}\n");

		// queueCurrGlobalHypotheses
		fprintf_s(fpHypothesis, "queueCurrGlobalHypotheses:%d,\n{\n", (int)queueCurrGlobalHypotheses_.size());
		for (hj::HypothesisSet::iterator hypothesisIter = queueCurrGlobalHypotheses_.begin();
			hypothesisIter != queueCurrGlobalHypotheses_.end();
			hypothesisIter++)
		{
			fprintf_s(fpHypothesis, "\t{\n\t\tselectedTracks:%d,", (int)(*hypothesisIter).selectedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).selectedTracks) + "\n";
			fprintf_s(fpHypothesis, trackIDList.c_str());

			fprintf_s(fpHypothesis, "\t\trelatedTracks:%d,", (int)(*hypothesisIter).relatedTracks.size());
			trackIDList = hj::MakeTrackIDList(&(*hypothesisIter).relatedTracks) + "\n";
			fprintf_s(fpHypothesis, trackIDList.c_str());

			fprintf_s(fpHypothesis, "\t\tlogLikelihood:%e\n", (*hypothesisIter).logLikelihood);
			fprintf_s(fpHypothesis, "\t\tprobability:%e\n", (*hypothesisIter).probability);
			fprintf_s(fpHypothesis, "\t\tbValid:%d\n\t}\n", (int)(*hypothesisIter).bValid);
		}
		fprintf_s(fpHypothesis, "}\n");
		fclose(fpHypothesis);

		//---------------------------------------------------------
		// VISUALIZATION
		//---------------------------------------------------------
		sprintf_s(strFilename, "%ssnapshot_3D_info.txt", strFilepath);
		fopen_s(&fpInfo, strFilename, "w");
		fprintf_s(fpInfo, "numCamera:%d\n", nNumCameras_);
		fprintf_s(fpInfo, "frameIndex:%d\n\n", (int)nCurrentFrameIdx_);
		fprintf_s(fpInfo, "nNewVisualizationID:%d\n", (int)nNewTargetID_);
		fprintf_s(fpInfo, "queuePairTreeIDToVisualizationID:%d,{", (int)queuePairTreeIDToTargetID_.size());
		for (int pairIdx = 0; pairIdx < queuePairTreeIDToTargetID_.size(); pairIdx++)
		{
			fprintf_s(fpInfo, "(%d,%d)", queuePairTreeIDToTargetID_[pairIdx].first, queuePairTreeIDToTargetID_[pairIdx].second);
		}
		fprintf_s(fpInfo, "}\n");

		fprintf_s(fpInfo, "()()\n");
		fprintf_s(fpInfo, "('')\n");

		fclose(fpInfo);
	}
	catch (int nError)
	{
		printf("[ERROR](SaveSnapShot) cannot open file! error code %d\n", nError);
		return;
	}
}

/************************************************************************
Method Name: LoadSnapshot
Description:
-
Input Arguments:
- none
Return Values:
-
************************************************************************/
bool CAssociator3D::LoadSnapshot(const char *strFilepath)
{
	FILE *fpTracklet, *fpTrack, *fpHypothesis, *fpResult, *fpInfo;
	char strFilename[128] = "";
	int readingInt = 0;
	float readingFloat = 0.0;

	std::deque<std::pair<CTrack3D*, unsigned int>> treeIDPair;
	std::deque<std::pair<CTrack3D*, std::deque<unsigned int>>> childrenTrackIDPair;

	try
	{
		// file open
		sprintf_s(strFilename, "%ssnapshot_3D_tracklet.txt", strFilepath);
		fopen_s(&fpTracklet, strFilename, "r");
		if (NULL == fpTracklet) { return false; }
		fscanf_s(fpTracklet, "numCamera:%d\n", &readingInt);
		assert(nNumCameras_ == readingInt);
		fscanf_s(fpTracklet, "frameIndex:%d\n\n", &readingInt);
		nCurrentFrameIdx_ = (unsigned int)readingInt;

		//---------------------------------------------------------
		// 2D TRACKLET RELATED
		//---------------------------------------------------------
		fscanf_s(fpTracklet, "vecTracklet2DSet:\n{\n");
		for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
		{
			int numTrackletSet = 0;
			fscanf_s(fpTracklet, "\t{\n\t\ttracklets:%d,\n\t\t{\n", &numTrackletSet);
			vecTracklet2DSet_[camIdx].tracklets.clear();
			for (int trackletIdx = 0; trackletIdx < numTrackletSet; trackletIdx++)
			{
				CTracklet2D newTracklet(nNumCameras_);

				fscanf_s(fpTracklet, "\t\t\t{\n");
				fscanf_s(fpTracklet, "\t\t\t\tid:%d\n", &readingInt);
				newTracklet.id = (unsigned int)readingInt;
				fscanf_s(fpTracklet, "\t\t\t\tcamIdx:%d\n", &readingInt);
				newTracklet.camIdx = (unsigned int)readingInt;
				fscanf_s(fpTracklet, "\t\t\t\tbActivated:%d\n", &readingInt);
				newTracklet.bActivated = 0 < readingInt ? true : false;
				fscanf_s(fpTracklet, "\t\t\t\ttimeStart:%d\n", &readingInt);
				newTracklet.timeStart = (unsigned int)readingInt;
				fscanf_s(fpTracklet, "\t\t\t\ttimeEnd:%d\n", &readingInt);
				newTracklet.timeEnd = (unsigned int)readingInt;
				fscanf_s(fpTracklet, "\t\t\t\tduration:%d\n", &readingInt);
				newTracklet.duration = (unsigned int)readingInt;

				// rects
				int numRect = 0;
				fscanf_s(fpTracklet, "\t\t\t\trects:%d,\n\t\t\t\t{\n", &numRect);
				for (int rectIdx = 0; rectIdx < numRect; rectIdx++)
				{
					float x, y, w, h;
					fscanf_s(fpTracklet, "\t\t\t\t\t(%f,%f,%f,%f)\n", &x, &y, &w, &h);
					newTracklet.rects.push_back(hj::Rect((double)x, (double)y, (double)w, (double)h));
				}
				fscanf_s(fpTracklet, "\t\t\t\t}\n");

				// backprojectionLines
				int numLine = 0;
				fscanf_s(fpTracklet, "\t\t\t\tbackprojectionLines:%d,\n\t\t\t\t{\n", &numLine);
				for (int lineIdx = 0; lineIdx < numLine; lineIdx++)
				{
					float x1, y1, z1, x2, y2, z2;
					fscanf_s(fpTracklet, "\t\t\t\t\t(%f,%f,%f,%f,%f,%f)\n", &x1, &y1, &z1, &x2, &y2, &z2);
					newTracklet.backprojectionLines.push_back(std::make_pair(hj::Point3D((double)x1, (double)y1, (double)z1), hj::Point3D((double)x2, (double)y2, (double)z2)));
				}
				fscanf_s(fpTracklet, "\t\t\t\t}\n");

				// appearance
				fscanf_s(fpTracklet, "\t\t\t\tRGBFeatureHead:%d,(", &readingInt);
				newTracklet.RGBFeatureHead = cv::Mat::zeros(readingInt, 1, CV_64FC1);
				for (int idx = 0; idx < readingInt; idx++)
				{
					fscanf_s(fpTracklet, "%f", &readingFloat);
					newTracklet.RGBFeatureHead.at<double>(idx, 0) = (double)readingFloat;
					if (idx < readingInt - 1) { fscanf_s(fpTracklet, ","); }
				}
				fscanf_s(fpTracklet, ")\n");
				fscanf_s(fpTracklet, "\t\t\t\tRGBFeatureTail:%d,(", &readingInt);
				newTracklet.RGBFeatureTail = cv::Mat::zeros(readingInt, 1, CV_64FC1);
				for (int idx = 0; idx < readingInt; idx++)
				{
					fscanf_s(fpTracklet, "%f", &readingFloat);
					newTracklet.RGBFeatureTail.at<double>(idx, 0) = (double)readingFloat;
					if (idx < readingInt - 1) { fscanf_s(fpTracklet, ","); }
				}
				fscanf_s(fpTracklet, ")\n");

				// location in 3D
				fscanf_s(fpTracklet, "\t\t\t\tcurrentLocation3D:(%lf,%lf,%lf)\n", 
					&newTracklet.currentLocation3D.x,
					&newTracklet.currentLocation3D.y, 
					&newTracklet.currentLocation3D.z);

				// bAssociableNewMeasurement
				for (int compCamIdx = 0; compCamIdx < nNumCameras_; compCamIdx++)
				{
					int numFlags = 0;
					fscanf_s(fpTracklet, "\t\t\t\tbAssociableNewMeasurement[%d]:%d,{", &readingInt, &numFlags);
					for (int flagIdx = 0; flagIdx < numFlags; flagIdx++)
					{
						fscanf_s(fpTracklet, "%d", &readingInt);
						bool curFlag = 0 < readingInt ? true : false;
						newTracklet.bAssociableNewMeasurement[compCamIdx].push_back(curFlag);
						if (flagIdx < numFlags - 1) { fscanf_s(fpTracklet, ","); }
					}
					fscanf_s(fpTracklet, "}\n");
				}
				fscanf_s(fpTracklet, "\t\t\t}\n");
				vecTracklet2DSet_[camIdx].tracklets.push_back(newTracklet);
			}

			int numActiveTracklet = 0;
			vecTracklet2DSet_[camIdx].activeTracklets.clear();
			fscanf_s(fpTracklet, "\t\tactiveTracklets:%d,{", &numActiveTracklet);
			for (int trackletIdx = 0; trackletIdx < numActiveTracklet; trackletIdx++)
			{
				fscanf_s(fpTracklet, "%d", &readingInt);
				// search active tracklet
				for (std::list<CTracklet2D>::iterator trackletIter = vecTracklet2DSet_[camIdx].tracklets.begin();
					trackletIter != vecTracklet2DSet_[camIdx].tracklets.end();
					trackletIter++)
				{
					if ((*trackletIter).id != (unsigned int)readingInt) { continue; }
					vecTracklet2DSet_[camIdx].activeTracklets.push_back(&(*trackletIter));
					break;
				}
				if (trackletIdx < numActiveTracklet - 1) { fscanf_s(fpTracklet, ","); }
			}

			int numNewTracklet = 0;
			fscanf_s(fpTracklet, "}\n\t\tnewMeasurements:%d,{", &numNewTracklet);
			for (int trackletIdx = 0; trackletIdx < numNewTracklet; trackletIdx++)
			{
				fscanf_s(fpTracklet, "%d", &readingInt);
				// search active tracklet
				for (std::list<CTracklet2D>::iterator trackletIter = vecTracklet2DSet_[camIdx].tracklets.begin();
					trackletIter != vecTracklet2DSet_[camIdx].tracklets.end();
					trackletIter++)
				{
					if ((*trackletIter).id != (unsigned int)readingInt) { continue; }
					vecTracklet2DSet_[camIdx].newMeasurements.push_back(&(*trackletIter));
					break;
				}
				if (trackletIdx < numNewTracklet - 1) { fscanf_s(fpTracklet, ","); }
			}
			fscanf_s(fpTracklet, "}\n\t}\n");
		}
		fscanf_s(fpTracklet, "}\n");
		fscanf_s(fpTracklet, "nNumTotalActive2DTracklet:%d\n\n", &readingInt);
		nNumTotalActive2DTracklet_ = (unsigned int)readingInt;

		fclose(fpTracklet);

		//---------------------------------------------------------
		// 3D TRACK RELATED
		//---------------------------------------------------------
		// file open
		sprintf_s(strFilename, "%ssnapshot_3D_track.txt", strFilepath);
		fopen_s(&fpTrack, strFilename, "r");
		if (NULL == fpTrack) { return false; }
		fscanf_s(fpTrack, "numCamera:%d\n", &readingInt);
		assert(nNumCameras_ == readingInt);
		fscanf_s(fpTrack, "frameIndex:%d\n\n", &readingInt);
		nCurrentFrameIdx_ = (unsigned int)readingInt;

		fscanf_s(fpTrack, "bReceiveNewMeasurement:%d\n", &readingInt);
		bReceiveNewMeasurement_ = 0 < readingInt ? true : false;
		fscanf_s(fpTrack, "bInitiationPenaltyFree:%d\n", &readingInt);
		bInitiationPenaltyFree_ = 0 < readingInt ? true : false;
		fscanf_s(fpTrack, "nNewTrackID:%d\n", &readingInt);
		nNewTrackID_ = (unsigned int)readingInt;
		fscanf_s(fpTrack, "nNewTreeID:%d\n", &readingInt);
		nNewTreeID_ = (unsigned int)readingInt;

		// track
		int numTrack = 0;
		listTrack3D_.clear();
		fscanf_s(fpTrack, "listTrack3D:%d,\n{\n", &numTrack);
		for (int trackIdx = 0; trackIdx < numTrack; trackIdx++)
		{
			CTrack3D newTrack(nNumCameras_);
			int parentTrackID = 0, treeID = 0;
			std::deque<unsigned int> childrenTrackID;
			fscanf_s(fpTrack, "\t{\n\t\tid:%d\n", &readingInt);
			newTrack.id = (unsigned int)readingInt;

			// curTracklet2Ds
			fscanf_s(fpTrack, "\t\tcurTracklet2Ds:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fscanf_s(fpTrack, "%d", &readingInt);
				if (-1 == readingInt)
				{
					newTrack.curTracklet2Ds.set(camIdx, NULL);
				}
				else
				{
					for (std::deque<CTracklet2D*>::iterator trackletIter = vecTracklet2DSet_[camIdx].activeTracklets.begin();
						trackletIter != vecTracklet2DSet_[camIdx].activeTracklets.end();
						trackletIter++)
					{
						if ((*trackletIter)->id != (unsigned int)readingInt) { continue; }
						newTrack.curTracklet2Ds.set(camIdx, *trackletIter);
						break;
					}
				}
				if (camIdx < nNumCameras_ - 1) { fscanf_s(fpTrack, ","); }
			}
			fscanf_s(fpTrack, "}\n");

			// tracklet2DIDs
			fscanf_s(fpTrack, "\t\ttrackleIDs:\n\t\t{\n");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				int numTracklet = 0;
				fscanf_s(fpTrack, "\t\t\tnumTracklet:%d,{", &numTracklet);

				for (int trackletIdx = 0; trackletIdx < numTracklet; trackletIdx++)
				{
					fscanf_s(fpTrack, "%d", &readingInt);
					newTrack.tracklet2DIDRecord[camIdx].push_back((unsigned int)readingInt);
					if (trackletIdx < numTracklet - 1) { fscanf_s(fpTrack, ","); }
				}

				fscanf_s(fpTrack, "}\n");
			}
			fscanf_s(fpTrack, "\t\t}\n");

			fscanf_s(fpTrack, "\t\tbActive:%d\n", &readingInt);
			newTrack.bActive = 0 < readingInt ? true : false;
			fscanf_s(fpTrack, "\t\tbValid:%d\n", &readingInt);
			newTrack.bValid = 0 < readingInt ? true : false;
			fscanf_s(fpTrack, "\t\ttreeID:%d\n", &treeID);
			fscanf_s(fpTrack, "\t\tparentTrackID:%d\n", &parentTrackID);
			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)parentTrackID != (*trackIter).id) { continue; }
				newTrack.parentTrack = &(*trackIter);
				break;
			}

			// children tracks
			int numChildrenTrack = 0;
			fscanf_s(fpTrack, "\t\tchildrenTrack:%d,{", &numChildrenTrack);
			for (int childTrackIdx = 0; childTrackIdx < numChildrenTrack; childTrackIdx++)
			{
				fscanf_s(fpTrack, "%d", &readingInt);
				childrenTrackID.push_back((unsigned int)readingInt);
				if (childTrackIdx < numChildrenTrack - 1) { fscanf_s(fpTrack, ","); }
			}
			fscanf_s(fpTrack, "}\n");

			// temporal information
			fscanf_s(fpTrack, "\t\ttimeStart:%d\n", &readingInt);
			newTrack.timeStart = (unsigned int)readingInt;
			fscanf_s(fpTrack, "\t\ttimeEnd:%d\n", &readingInt);
			newTrack.timeEnd = (unsigned int)readingInt;
			fscanf_s(fpTrack, "\t\ttimeGeneration:%d\n", &readingInt);
			newTrack.timeGeneration = (unsigned int)readingInt;
			fscanf_s(fpTrack, "\t\tduration:%d\n", &readingInt);
			newTrack.duration = (unsigned int)readingInt;

			// reconstrcution related
			int numReconstruction = 0;
			fscanf_s(fpTrack, "\t\treconstructions:%d,\n\t\t{\n", &numReconstruction);
			std::deque<hj::Point3D> points(newTrack.duration);
			std::deque<hj::Point3D> smoothedPoints(newTrack.duration);
			for (int pointIdx = 0; pointIdx < numReconstruction; pointIdx++)
			{
				CReconstruction newReconstruction(nNumCameras_);
				float x, y, z, sx, sy, sz, vx, vy, vz, maxError, costReconstruction, costSmoothedPoint, costLink;
				fscanf_s(fpTrack, "\t\t\t%d,point:(%f,%f,%f),smoothedPoint:(%f,%f,%f),velocity:(%f,%f,%f),maxError:%e,costReconstruction:%e,costSmoothedPoint:%e,costLink:%e,",
					&readingInt, &x, &y, &z, &sx, &sy, &sz, &vx, &vy, &vz, &maxError, &costReconstruction, &costSmoothedPoint, &costLink);
				newReconstruction.bIsMeasurement = 0 < readingInt ? true : false;
				newReconstruction.point = hj::Point3D((double)x, (double)y, (double)z);
				newReconstruction.smoothedPoint = hj::Point3D((double)sx, (double)sy, (double)sz);
				newReconstruction.velocity = hj::Point3D((double)vx, (double)vy, (double)vz);
				newReconstruction.maxError = (double)maxError;
				newReconstruction.costLink = (double)costLink;
				newReconstruction.costReconstruction = (double)costReconstruction;
				newReconstruction.costSmoothedPoint = (double)costSmoothedPoint;

				if (pointIdx < (int)newTrack.duration)
				{
					points[pointIdx].x = x;
					points[pointIdx].y = y;
					points[pointIdx].z = z;
					smoothedPoints[pointIdx].x = sx;
					smoothedPoints[pointIdx].y = sy;
					smoothedPoints[pointIdx].z = sz;
				}

				fscanf_s(fpTrack, "tracklet2Ds:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fscanf_s(fpTrack, "%d", &readingInt);

					if (0 > readingInt)
					{
						newReconstruction.tracklet2Ds.set(camIdx, NULL);
					}
					else
					{
						for (std::list<CTracklet2D>::iterator trackletIter = vecTracklet2DSet_[camIdx].tracklets.begin();
							trackletIter != vecTracklet2DSet_[camIdx].tracklets.end();
							trackletIter++)
						{
							if ((*trackletIter).id != (unsigned int)readingInt) { continue; }
							newReconstruction.tracklet2Ds.set(camIdx, &(*trackletIter));
							break;
						}
					}
					if (camIdx < nNumCameras_ - 1) { fscanf_s(fpTrack, ","); }
				}
				fscanf_s(fpTrack, "},");
				fscanf_s(fpTrack, "rawPoints:%d,{", &readingInt);
				for (int pointIdx = 0; pointIdx < readingInt; pointIdx++)
				{
					fscanf_s(fpTrack, "(%f,%f,%f)", &x, &y, &z);
					newReconstruction.rawPoints.push_back(hj::Point3D((double)x, (double)y, (double)z));
					if (pointIdx < readingInt - 1) { fscanf_s(fpTrack, ","); }
				}
				fscanf_s(fpTrack, "}\n");

				newTrack.reconstructions.push_back(newReconstruction);
			}
			fscanf_s(fpTrack, "\t\t}\n");

			// point smoother
			int span, degree;
			fscanf_s(fpTrack, "\t\tsmoother:{span:%d,degree:%d}\n", &span, &degree);
			newTrack.smoother.SetSmoother(points, smoothedPoints, span, degree);
			newTrack.smoother.SetQsets(&precomputedQsets);

			// cost
			fscanf_s(fpTrack, "\t\tcostTotal:%e\n", &readingFloat);
			newTrack.costTotal = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tcostReconstruction:%e\n", &readingFloat);
			newTrack.costReconstruction = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tcostLink:%e\n", &readingFloat);
			newTrack.costLink = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tcostEnter:%e\n", &readingFloat);
			newTrack.costEnter = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tcostExit:%e\n", &readingFloat);
			newTrack.costExit = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tcostRGB:%e\n", &readingFloat);
			newTrack.costRGB = (double)readingFloat;

			// loglikelihood
			fscanf_s(fpTrack, "\t\tloglikelihood:%e\n", &readingFloat);
			newTrack.loglikelihood = (double)readingFloat;

			// GTP
			fscanf_s(fpTrack, "\t\tGTProb:%e\n", &readingFloat);
			newTrack.GTProb = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tBranchGTProb:%e\n", &readingFloat);
			newTrack.BranchGTProb = (double)readingFloat;
			fscanf_s(fpTrack, "\t\tbWasBestSolution:%d\n", &readingInt);
			newTrack.bWasBestSolution = 0 < readingInt ? true : false;
			fscanf_s(fpTrack, "\t\tbCurrentBestSolution:%d\n", &readingInt);
			newTrack.bCurrentBestSolution = 0 < readingInt ? true : false;

			// HO-MHT
			fscanf_s(fpTrack, "\t\tbNewTrack:%d\n", &readingInt);
			newTrack.bNewTrack = 0 < readingInt ? true : false;

			// tracklet related
			fscanf_s(fpTrack, "\t\ttimeTrackletEnded:(");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fscanf_s(fpTrack, "%d", &readingInt);
				newTrack.timeTrackletEnded[camIdx] = (unsigned int)readingInt;
				if (camIdx < nNumCameras_ - 1) { fscanf_s(fpTrack, ","); }
			}
			fscanf_s(fpTrack, ")\n");
			fscanf_s(fpTrack, "\t\tlastTrackletLocation3D:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				float tx = 0, ty = 0, tz = 0;
				fscanf_s(fpTrack, "(%f,%f,%f)", &tx, &ty, &tz);
				newTrack.lastTrackletLocation3D[camIdx] = hj::Point3D((double)tx, (double)ty, (double)tz);
			}
			fscanf_s(fpTrack, "}\n");
			fscanf_s(fpTrack, "\t\tlastTrackletSensitivity:{");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				fscanf_s(fpTrack, "%f", &readingFloat);
				newTrack.lastTrackletSensitivity[camIdx] = (double)readingFloat;
				if (camIdx < nNumCameras_ - 1) { fscanf_s(fpTrack, ","); }
			}
			fscanf_s(fpTrack, "}\n");
			fscanf_s(fpTrack, "\t\tlastRGBFeature:\n\t\t{\n");
			for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
			{
				int dimFeature = 0;
				fscanf_s(fpTrack, "\t\t\tlastRGBFeature[%d]:%d,(", &readingInt, &dimFeature);
				newTrack.lastRGBFeature[camIdx] = cv::Mat::zeros(dimFeature, 1, CV_64FC1);
				for (int idx = 0; idx < dimFeature; idx++)
				{
					fscanf_s(fpTrack, "%f", &readingFloat);
					newTrack.lastRGBFeature[camIdx].at<double>(idx, 0) = (double)readingFloat;
					if (idx < dimFeature - 1) { fscanf_s(fpTrack, ","); }
				}
				fscanf_s(fpTrack, ")\n");
			}
			fscanf_s(fpTrack, "\t\t}\n");

			// termination
			fscanf_s(fpTrack, "\t\tbnumOutpoint:%d\n", &readingInt);
			newTrack.numOutpoint = (unsigned int)readingInt;

			// track info end
			fscanf_s(fpTrack, "\t}\n");

			// generate instance
			listTrack3D_.push_back(newTrack);
			treeIDPair.push_back(std::make_pair(&listTrack3D_.back(), (unsigned int)treeID));
			childrenTrackIDPair.push_back(std::make_pair(&listTrack3D_.back(), childrenTrackID));
		}
		fscanf_s(fpTrack, "}\n");

		// children track
		for (int queueIdx = 0; queueIdx < childrenTrackIDPair.size(); queueIdx++)
		{
			CTrack3D *curTrack = childrenTrackIDPair[queueIdx].first;
			for (int idIdx = 0; idIdx < childrenTrackIDPair[queueIdx].second.size(); idIdx++)
			{
				unsigned int curChildID = childrenTrackIDPair[queueIdx].second[idIdx];
				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if (curChildID != (*trackIter).id) { continue; }
					curTrack->childrenTrack.push_back(&(*trackIter));
					break;
				}
			}
		}

		// track tree
		int numTree = 0;
		listTrackTree_.clear();
		fscanf_s(fpTrack, "listTrackTree:%d,\n{\n", &numTree);
		for (int treeIdx = 0; treeIdx < numTree; treeIdx++)
		{
			CTrackTree newTree;

			fscanf_s(fpTrack, "\t{\n\t\tid:%d\n", &readingInt);
			newTree.id = (unsigned int)readingInt;
			fscanf_s(fpTrack, "\t\ttimeGeneration:%d\n", &readingInt);
			newTree.timeGeneration = (unsigned int)readingInt;
			fscanf_s(fpTrack, "\t\tbValid:%d\n", &readingInt);
			newTree.bValid = 0 < readingInt ? true : false;
			fscanf_s(fpTrack, "\t\tbConfirmed:%d\n", &readingInt);
			newTree.bConfirmed = 0 < readingInt ? true : false;
			int numTracks = 0;
			fscanf_s(fpTrack, "\t\ttracks:%d,{", &numTracks);
			for (int trackIdx = 0; trackIdx < numTracks; trackIdx++)
			{
				fscanf_s(fpTrack, "%d", &readingInt);
				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if ((unsigned int)readingInt != (*trackIter).id) { continue; }
					newTree.tracks.push_back(&(*trackIter));
					break;
				}
				if (trackIdx < numTracks - 1) { fscanf_s(fpTrack, ","); }
			}
			fscanf_s(fpTrack, "}\n");
			fscanf_s(fpTrack, "\t}\n");
			listTrackTree_.push_back(newTree);
		}
		fscanf_s(fpTrack, "}\n");

		// set track's tree pointers
		for (int pairIdx = 0; pairIdx < treeIDPair.size(); pairIdx++)
		{
			CTrack3D *curTrack = treeIDPair[pairIdx].first;
			for (std::list<CTrackTree>::iterator treeIter = listTrackTree_.begin();
				treeIter != listTrackTree_.end();
				treeIter++)
			{
				if (treeIDPair[pairIdx].second != (*treeIter).id) { continue; }
				curTrack->tree = &(*treeIter);
				break;
			}
		}

		// new seed tracks
		int numSeedTracks = 0;
		queueNewSeedTracks_.clear();
		fscanf_s(fpTrack, "queueNewSeedTracks:%d,{", &numSeedTracks);
		for (int trackIdx = 0; trackIdx < numSeedTracks; trackIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (trackIdx < numSeedTracks - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)readingInt != (*trackIter).id) { continue; }
				queueNewSeedTracks_.push_back(&(*trackIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// active tracks
		int numActiveTracks = 0;
		queueActiveTrack_.clear();
		fscanf_s(fpTrack, "queueActiveTrack:%d,{", &numActiveTracks);
		for (int trackIdx = 0; trackIdx < numActiveTracks; trackIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (trackIdx < numActiveTracks - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)readingInt != (*trackIter).id) { continue; }
				queueActiveTrack_.push_back(&(*trackIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// puased tracks
		int numPuasedTracks = 0;
		queuePausedTrack_.clear();
		fscanf_s(fpTrack, "queuePausedTrack:%d,{", &numPuasedTracks);
		for (int trackIdx = 0; trackIdx < numPuasedTracks; trackIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (trackIdx < numPuasedTracks - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)readingInt != (*trackIter).id) { continue; }
				queuePausedTrack_.push_back(&(*trackIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// tracks in window
		int numInWindowTracks = 0;
		queueTracksInWindow_.clear();
		fscanf_s(fpTrack, "queueTracksInWindow:%d,{", &numInWindowTracks);
		for (int trackIdx = 0; trackIdx < numInWindowTracks; trackIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (trackIdx < numInWindowTracks - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)readingInt != (*trackIter).id) { continue; }
				queueTracksInWindow_.push_back(&(*trackIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// tracks in the best solution
		int numBSTracks = 0;
		queueTracksInBestSolution_.clear();
		fscanf_s(fpTrack, "queueTracksInBestSolution:%d,{", &numBSTracks);
		for (int trackIdx = 0; trackIdx < numBSTracks; trackIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (trackIdx < numBSTracks - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
				trackIter != listTrack3D_.end();
				trackIter++)
			{
				if ((unsigned int)readingInt != (*trackIter).id) { continue; }
				queueTracksInBestSolution_.push_back(&(*trackIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// active trees
		int numActiveTrees = 0;
		queuePtActiveTrees_.clear();
		fscanf_s(fpTrack, "queuePtActiveTrees:%d,{", &numActiveTrees);
		for (int treeIdx = 0; treeIdx < numActiveTrees; treeIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (treeIdx < numActiveTrees - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrackTree>::iterator treeIter = listTrackTree_.begin();
				treeIter != listTrackTree_.end();
				treeIter++)
			{
				if ((unsigned int)readingInt != (*treeIter).id) { continue; }
				queuePtActiveTrees_.push_back(&(*treeIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");

		// unconfirmed trees
		int numUCTrees = 0;
		queuePtUnconfirmedTrees_.clear();
		fscanf_s(fpTrack, "queuePtUnconfirmedTrees:%d,{", &numUCTrees);
		for (int treeIdx = 0; treeIdx < numUCTrees; treeIdx++)
		{
			fscanf_s(fpTrack, "%d", &readingInt);
			if (treeIdx < numUCTrees - 1) { fscanf_s(fpTrack, ","); }

			for (std::list<CTrackTree>::iterator treeIter = listTrackTree_.begin();
				treeIter != listTrackTree_.end();
				treeIter++)
			{
				if ((unsigned int)readingInt != (*treeIter).id) { continue; }
				queuePtUnconfirmedTrees_.push_back(&(*treeIter));
				break;
			}
		}
		fscanf_s(fpTrack, "}\n");
		fclose(fpTrack);

		//---------------------------------------------------------
		// RESULTS
		//---------------------------------------------------------
		// file open
		sprintf_s(strFilename, "%ssnapshot_3D_result.txt", strFilepath);
		fopen_s(&fpResult, strFilename, "r");
		if (NULL == fpResult) { return false; }
		fscanf_s(fpResult, "numCamera:%d\n", &readingInt);
		assert(nNumCameras_ == readingInt);
		fscanf_s(fpResult, "frameIndex:%d\n\n", &readingInt);
		nCurrentFrameIdx_ = (unsigned int)readingInt;

		// instance result
		int numTrackingResults = 0;
		queueTrackingResult_.clear();
		fscanf_s(fpResult, "queueTrackingResult:%d,\n{\n", &numTrackingResults);
		for (int resultIdx = 0; resultIdx < numTrackingResults; resultIdx++)
		{
			CTrack3DResult curResult;
			fscanf_s(fpResult, "\t{\n\t\tframeIdx:%d\n", &readingInt);
			curResult.frameIdx = (unsigned int)readingInt;
			fscanf_s(fpResult, "\t\tprocessingTime:%f\n", &readingFloat);
			curResult.processingTime = (double)readingFloat;

			// object info
			int numObjects = 0;
			fscanf_s(fpResult, "\t\tobjectInfo:%d,\n\t\t{\n", &numObjects);
			for (int objIdx = 0; objIdx < numObjects; objIdx++)
			{
				CObject3DInfo curObject(nNumCameras_);
				fscanf_s(fpResult, "\t\t\t{\n\t\t\t\tid:%d\n", &readingInt);
				curObject.id = (unsigned int)readingInt;

				// points
				int numPoints = 0;
				fscanf_s(fpResult, "\t\t\t\trecentPoints:%d,{", &numPoints);
				for (int pointIdx = 0; pointIdx < numPoints; pointIdx++)
				{
					float x, y, z;
					fscanf_s(fpResult, "(%f,%f,%f)", &x, &y, &z);
					if (pointIdx < numPoints - 1) { fscanf_s(fpResult, ","); }
					curObject.recentPoints.push_back(hj::Point3D((double)x, (double)y, (double)z));
				}
				fscanf_s(fpResult, "}\n");

				// 2D points
				fscanf_s(fpResult, "\t\t\t\trecentPoint2Ds:\n\t\t\t\t{\n");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fscanf_s(fpResult, "\t\t\t\t\tcam%d:%d,{", &readingInt, &numPoints);
					for (int pointIdx = 0; pointIdx < numPoints; pointIdx++)
					{
						float x, y;
						fscanf_s(fpResult, "(%f,%f)", &x, &y);
						if (pointIdx < numPoints - 1) { fscanf_s(fpResult, ","); }
						curObject.recentPoint2Ds[camIdx].push_back(hj::Point2D((double)x, (double)y));
					}
					fscanf_s(fpResult, "}\n");
				}
				fscanf_s(fpResult, "\t\t\t\t}\n");

				// 3D box points in each view
				fscanf_s(fpResult, "\t\t\t\tpoint3DBox:\n\t\t\t\t{\n");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fscanf_s(fpResult, "\t\t\t\t\tcam%d:%d,{", &readingInt, &numPoints);
					for (int pointIdx = 0; pointIdx < numPoints; pointIdx++)
					{
						float x, y;
						fscanf_s(fpResult, "(%f,%f)", &x, &y);
						if (pointIdx < numPoints - 1) { fscanf_s(fpResult, ","); }
						curObject.point3DBox[camIdx].push_back(hj::Point2D((double)x, (double)y));
					}
					fscanf_s(fpResult, "}\n");
				}
				fscanf_s(fpResult, "\t\t\t\t}\n");

				// rects
				fscanf_s(fpResult, "\t\t\t\trectInViews:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					float x, y, w, h;
					fscanf_s(fpResult, "(%f,%f,%f,%f)", &x, &y, &w, &h);
					if (camIdx < nNumCameras_ - 1) { fscanf_s(fpResult, ","); }
					curObject.rectInViews[camIdx] = hj::Rect((double)x, (double)y, (double)w, (double)h);
				}
				fscanf_s(fpResult, "}\n");

				// visibility
				fscanf_s(fpResult, "\t\t\t\tbVisibleInViews:{");
				for (int camIdx = 0; camIdx < nNumCameras_; camIdx++)
				{
					fscanf_s(fpResult, "%d", &readingInt);
					if (camIdx < nNumCameras_ - 1) { fscanf_s(fpResult, ","); }
					curObject.bVisibleInViews[camIdx] = 0 < readingInt ? true : false;
				}
				fscanf_s(fpResult, "}\n\t\t\t}\n");

				curResult.object3DInfos.push_back(curObject);
			}
			fscanf_s(fpResult, "\t\t}\n\t}\n");

			queueTrackingResult_.push_back(curResult);
		}
		fscanf_s(fpResult, "}\n");
		fclose(fpResult);

		//---------------------------------------------------------
		// HYPOTHESES
		//---------------------------------------------------------
		// file open
		sprintf_s(strFilename, "%ssnapshot_3D_hypotheses.txt", strFilepath);
		fopen_s(&fpHypothesis, strFilename, "r");
		if (NULL == fpHypothesis) { return false; }
		fscanf_s(fpHypothesis, "numCamera:%d\n", &readingInt);
		assert(nNumCameras_ == readingInt);
		fscanf_s(fpHypothesis, "frameIndex:%d\n\n", &readingInt);
		nCurrentFrameIdx_ = (unsigned int)readingInt;

		// queuePrevGlobalHypotheses
		int numPrevGH = 0;
		queuePrevGlobalHypotheses_.clear();
		fscanf_s(fpHypothesis, "queuePrevGlobalHypotheses:%d,\n{\n", &numPrevGH);
		for (int hypothesisIdx = 0; hypothesisIdx < numPrevGH; hypothesisIdx++)
		{
			stGlobalHypothesis newHypothesis;

			int numSelectedTracks = 0;
			fscanf_s(fpHypothesis, "\t{\n\t\tselectedTracks:%d,{", &numSelectedTracks);
			for (int trackIdx = 0; trackIdx < numSelectedTracks; trackIdx++)
			{
				fscanf_s(fpHypothesis, "%d", &readingInt);
				if (trackIdx < numSelectedTracks - 1) { fscanf_s(fpHypothesis, ","); }

				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if ((unsigned int)readingInt != (*trackIter).id) { continue; }
					newHypothesis.selectedTracks.push_back(&(*trackIter));
					break;
				}
			}
			fscanf_s(fpHypothesis, "}\n");

			int numRelatedTracks = 0;
			fscanf_s(fpHypothesis, "\t\trelatedTracks:%d,{", &numRelatedTracks);
			for (int trackIdx = 0; trackIdx < numRelatedTracks; trackIdx++)
			{
				fscanf_s(fpHypothesis, "%d", &readingInt);
				if (trackIdx < numRelatedTracks - 1) { fscanf_s(fpHypothesis, ","); }

				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if ((unsigned int)readingInt != (*trackIter).id) { continue; }
					newHypothesis.relatedTracks.push_back(&(*trackIter));
					break;
				}
			}
			fscanf_s(fpHypothesis, "}\n");
			fscanf_s(fpHypothesis, "\t\tlogLikelihood:%e\n", &readingFloat);
			newHypothesis.logLikelihood = (double)readingFloat;
			fscanf_s(fpHypothesis, "\t\tprobability:%e\n", &readingFloat);
			newHypothesis.probability = (double)readingFloat;
			fscanf_s(fpHypothesis, "\t\tbValid:%d\n\t}\n", &readingInt);
			newHypothesis.bValid = 0 < readingInt ? true : false;

			queuePrevGlobalHypotheses_.push_back(newHypothesis);

		}
		fscanf_s(fpHypothesis, "}\n");

		// queueCurrGlobalHypotheses
		int numCurrGH = 0;
		queueCurrGlobalHypotheses_.clear();
		fscanf_s(fpHypothesis, "queueCurrGlobalHypotheses:%d,\n{\n", &numCurrGH);
		for (int hypothesisIdx = 0; hypothesisIdx < numCurrGH; hypothesisIdx++)
		{
			stGlobalHypothesis newHypothesis;

			int numSelectedTracks = 0;
			fscanf_s(fpHypothesis, "\t{\n\t\tselectedTracks:%d,{", &numSelectedTracks);
			for (int trackIdx = 0; trackIdx < numSelectedTracks; trackIdx++)
			{
				fscanf_s(fpHypothesis, "%d", &readingInt);
				if (trackIdx < numSelectedTracks - 1) { fscanf_s(fpHypothesis, ","); }

				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if ((unsigned int)readingInt != (*trackIter).id) { continue; }
					newHypothesis.selectedTracks.push_back(&(*trackIter));
					break;
				}
			}
			fscanf_s(fpHypothesis, "}\n");

			int numRelatedTracks = 0;
			fscanf_s(fpHypothesis, "\t\trelatedTracks:%d,{", &numRelatedTracks);
			for (int trackIdx = 0; trackIdx < numRelatedTracks; trackIdx++)
			{
				fscanf_s(fpHypothesis, "%d", &readingInt);
				if (trackIdx < numRelatedTracks - 1) { fscanf_s(fpHypothesis, ","); }

				for (std::list<CTrack3D>::iterator trackIter = listTrack3D_.begin();
					trackIter != listTrack3D_.end();
					trackIter++)
				{
					if ((unsigned int)readingInt != (*trackIter).id) { continue; }
					newHypothesis.relatedTracks.push_back(&(*trackIter));
					break;
				}
			}
			fscanf_s(fpHypothesis, "}\n");
			fscanf_s(fpHypothesis, "\t\tlogLikelihood:%e\n", &readingFloat);
			newHypothesis.logLikelihood = (double)readingFloat;
			fscanf_s(fpHypothesis, "\t\tprobability:%e\n", &readingFloat);
			newHypothesis.probability = (double)readingFloat;
			fscanf_s(fpHypothesis, "\t\tbValid:%d\n\t}\n", &readingInt);
			newHypothesis.bValid = 0 < readingInt ? true : false;

			queueCurrGlobalHypotheses_.push_back(newHypothesis);
		}
		fscanf_s(fpHypothesis, "}\n");
		fclose(fpHypothesis);

		//---------------------------------------------------------
		// VISUALIZATION
		//---------------------------------------------------------
		// file open
		sprintf_s(strFilename, "%ssnapshot_3D_info.txt", strFilepath);
		fopen_s(&fpInfo, strFilename, "r");
		if (NULL == fpInfo) { return false; }
		fscanf_s(fpInfo, "numCamera:%d\n", &readingInt);
		assert(nNumCameras_ == readingInt);
		fscanf_s(fpInfo, "frameIndex:%d\n\n", &readingInt);
		nCurrentFrameIdx_ = (unsigned int)readingInt;

		fscanf_s(fpInfo, "nNewVisualizationID:%d\n", &readingInt);
		nNewTargetID_ = (unsigned int)readingInt;

		int numPair = 0;
		queuePairTreeIDToTargetID_.clear();
		fscanf_s(fpInfo, "queuePairTreeIDToVisualizationID:%d,{", &numPair);
		for (int pairIdx = 0; pairIdx < numPair; pairIdx++)
		{
			int id1, id2;
			fscanf_s(fpInfo, "(%d,%d)", &id1, &id2);
			queuePairTreeIDToTargetID_.push_back(std::make_pair(id1, id2));
		}
		fscanf_s(fpInfo, "}\n");

		fclose(fpInfo);
	}
	catch (int nError)
	{
		printf("[ERROR](LoadSnapShot) cannot open file! error code %d\n", nError);
		return false;
	}

	return true;
}

/************************************************************************
Method Name: IndexCombination
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
std::deque<std::vector<unsigned int>> CAssociator3D::IndexCombination(std::deque<std::deque<unsigned int>> &inputIndexDoubleArray, size_t curLevel, std::deque<std::vector<unsigned int>> curCombination)
{
	if (0 == curLevel)
	{
		std::deque<std::vector<unsigned int>> newCombinations;
		for (unsigned int indexIdx = 0; indexIdx < (unsigned int)inputIndexDoubleArray[0].size(); indexIdx++)
		{
			std::vector<unsigned int> curIndex(inputIndexDoubleArray.size(), 0);
			curIndex[0] = indexIdx;
			newCombinations.push_back(curIndex);
		}
		return IndexCombination(inputIndexDoubleArray, 1, newCombinations);
	}
	else if (inputIndexDoubleArray.size() <= curLevel)
	{
		return curCombination;
	}

	std::deque<std::vector<unsigned int>> newCombinations;
	for (unsigned int indexIdx = 0; indexIdx < (unsigned int)inputIndexDoubleArray[curLevel].size(); indexIdx++)
	{
		for (size_t combiIdx = 0; combiIdx < curCombination.size(); combiIdx++)
		{
			std::vector<unsigned int> newIndex = curCombination[combiIdx];
			newIndex[curLevel] = indexIdx;
			newCombinations.push_back(newIndex);
		}
	}
	return IndexCombination(inputIndexDoubleArray, curLevel + 1, newCombinations);
}


}

//()()
//('')HAANJU.YOO

