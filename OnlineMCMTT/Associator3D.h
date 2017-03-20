/**************************************************************************
* Title        : CAssociator3D
* Author       : Haanju Yoo
* Initial Date : 2014.01.22 (ver. 0.9)
* Version Num. : 1.0 (since 2016.11.04)
* Description  :
*	Generate 3D trajectories by spatial-temporal association between 
	tracklets.
**************************************************************************/


#pragma once

#include "GraphSolver.h"
#include "types.hpp"
#include "haanju_3D.hpp"


namespace hj
{

/////////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////////

//struct stGraphSolution
//{
//	TrackSet tracks;
//	double logLikelihood; // = weight sum
//	double probability;
//	//double weightSum;
//	bool bValid;
//};

struct stGlobalHypothesis
{
	TrackSet selectedTracks;
	TrackSet relatedTracks;
	double logLikelihood; // = weight sum
	double probability;
	bool bValid;
};
typedef std::deque<stGlobalHypothesis> HypothesisSet;
typedef std::pair<unsigned int, unsigned int> PAIR_UINT;

// for GRAPH SOLVING ANALYSIS
struct stHypothesisSolvingInfo
{
	int rank;
	double costInitialSolution;
	double costAdded;
	double reviousSolution;
	stGraphSolvingResult solvingResult;
};

/////////////////////////////////////////////////////////////////////////
// CLASS DECLARATION
/////////////////////////////////////////////////////////////////////////
class CAssociator3D
{
	//////////////////////////////////////////////////////////////////////////
	// METHODS
	//////////////////////////////////////////////////////////////////////////
public:
	CAssociator3D(void);
	~CAssociator3D(void);
	void Initialize(stParamAssociate3D &_stParams);
	void Finalize(void);
	CTrack3DResult Run(std::vector<CTrack2DResult> &curTrack2DResult, std::vector<cv::Mat> curFrames, int frameIdx);

	// Visualization related
	std::vector<hj::Point2D> GetHuman3DBox(hj::Point3D ptHeadCenter, double bodyWidth, int camIdx);

private:
	//----------------------------------------------------------------
	// 3D geometry related
	//----------------------------------------------------------------
	hj::Point2D WorldToImage(hj::Point3D point3D, int camIdx);
	hj::Point3D ImageToWorld(hj::Point2D point2D, double z, int camIdx);

	bool CheckVisibility(hj::Point3D testPoint, int camIdx, hj::Point2D *result2DPoint = NULL);
	bool CheckTrackletConnectivity(hj::Point3D endPoint, hj::Point3D startPoint, double sensitivity1, double sensitivity2, int timeGap);

	CReconstruction PointReconstruction(CTrackletCombination &tracklet2Ds);
	double NViewPointReconstruction(std::vector<hj::Line3D> &vecLines, hj::Point3D &outputPoint);
	double NViewGroundingPointReconstruction(std::vector<hj::Point2D_CamIdx> &vecPointInfos, hj::Point3D &outputPoint);
	hj::Line3D GetBackProjectionLine(hj::Point2D point2D, int camIdx);
	double GetDistanceFromBoundary(hj::Point3D point);

	//----------------------------------------------------------------
	// 2D tracklet related
	//----------------------------------------------------------------
	void Tracklet2D_UpdateTracklets(std::vector<CTrack2DResult> &curTrack2DResult, unsigned int frameIdx);
	void GenerateTrackletCombinations(std::vector<std::vector<bool>> &vecvecBAssociationMap, CTrackletCombination combination, std::deque<CTrackletCombination> &combinationQueue, int camIdx);

	//----------------------------------------------------------------
	// 3D track related
	//----------------------------------------------------------------	
	void Track3D_Management(TrackSet &outputSeedTracks, const unsigned int _frameIdx);
	void Track3D_UpdateTracks(const unsigned int _frameIdx);
	void Track3D_GenerateSeedTracks(TrackSet &outputSeedTracks, const unsigned int _frameIdx);
	void Track3D_BranchTracks(TrackSet *seedTracks, const unsigned int _frameIdx);
	TrackSet Track3D_GetWholeCandidateTracks(void);

	//void Track3D_SolveHOMHT(void);
	//void Track3D_Pruning_GTP(void);
	//void Track3D_Pruning_KBest(void);
	//void Track3D_RepairDataStructure(void);

	// cost calculation
	double ComputeEnterProbability(std::vector<hj::Point3D> &points);
	double ComputeExitProbability(std::vector<hj::Point3D> &points, int trackLength);
	double ComputeLinkProbability(hj::Point3D &prePoint, hj::Point3D &curPoint, unsigned int timeGap = 1);
	double ComputeRGBCost(const cv::Mat *feature1, const cv::Mat *feature2, unsigned int timeGap);
	double ComputeTrackletLinkCost(hj::Point3D preLocation, hj::Point3D curLocation, int timeGap);
	double ComputeReconstructionProbability(hj::Point3D point, std::vector<hj::Point3D> *rawPoints, CTrackletCombination *trackletCombination, double maxError = 0.0);

	// miscellaneous
public:
	static bool CheckIncompatibility(CTrack3D *track1, CTrack3D *track2, double maxTrackletDistance, double minProximity);
	static bool CheckIncompatibility(CTrackletCombination &combi1, CTrackletCombination &combi2);
	static cv::Mat GetRGBFeature(const cv::Mat *patch, int numBins);
	static double GetCost(CTrack3D *track);
private:

	//----------------------------------------------------------------
	// Hypothesis related
	//----------------------------------------------------------------	
	void Hypothesis_UpdateHypotheses(hj::HypothesisSet &inoutUpdatedHypotheses, TrackSet *newSeedTracks);
	void Hypothesis_Formation(hj::HypothesisSet &outBranchHypotheses, hj::HypothesisSet *existingHypotheses);
public:
	static void Hypothesis_BranchHypotheses(hj::HypothesisSet &outBranchHypotheses, stParamAssociate3D *params, TrackSet *tracks, TrackSet *initialSolutionTracks = NULL);
private:
	void Hypothesis_PruningNScanBack(unsigned int nCurrentFrameIdx, unsigned int N, TrackSet *tracksInWindow, std::deque<stGlobalHypothesis> *ptQueueHypothesis = NULL);
	void Hypothesis_PruningTrackWithGTP(unsigned int nCurrentFrameIdx, unsigned int nNumMaximumTrack, TrackSet *tracksInWindow, std::deque<CTrackTree*> *queueActiveTrackTree);
	void Hypothesis_RefreshHypotheses(hj::HypothesisSet &inoutUpdatedHypotheses);

	//----------------------------------------------------------------
	// interface
	//----------------------------------------------------------------
	CTrack3DResult ResultWithTracks(TrackSet *trackSet, unsigned int nFrameIdx, double fProcessingTime = 0);
	void PrintTracks(hj::TrackSet &queueTracks, char *strFilePathAndName, bool bAppend);
	void PrintHypotheses(hj::HypothesisSet &queueHypotheses, char *strFilePathAndName, unsigned int frameIdx);
	void PrintCurrentTrackTrees(const char *strFilePath);
	void PrintResult(const char *strFilepath, std::deque<CTrack3DResult> *queueResults);

	//----------------------------------------------------------------
	// ETC
	//----------------------------------------------------------------
	static std::deque<std::vector<unsigned int>> IndexCombination(std::deque<std::deque<unsigned int>> &inputIndexDoubleArray, size_t curLevel, std::deque<std::vector<unsigned int>> curCombination);
	void VisualizeResult(const unsigned int _frameIdx);
	void SaveSnapshot(const char *strFilepath);
	bool LoadSnapshot(const char *strFilepath);

	//////////////////////////////////////////////////////////////////////////
	// VARIABLES
	//////////////////////////////////////////////////////////////////////////
private:
	stParamAssociate3D stParam_;	
	bool bInit_;
	bool bSnapshotReaded_;
	int  nNumCameras_;
	
	std::vector<cv::Mat> vecMatProjectionSensitivity_;
	std::vector<cv::Mat> vecMatDistanceFromBoundary_;

	char strDatasetPath_[128];

	// frame related
	std::vector<cv::Mat> vecMatCurrentFrames_;
	unsigned int nCurrentFrameIdx_;
	unsigned int nNumFramesForProc_;
	unsigned int nCountForPenalty_;
	unsigned int nLastPrintedDeferredResultFrameIdx_;
	unsigned int nLastPrintedInstantResultFrameIdx_;
	double dCurrentProcessingTime_;
	double dCurrentSolvingTime_;

	// logging related
	std::string strTime_;
	std::string strLogFileName_;
	std::string strTrackLogFileName_;
	std::deque<double> queueProcessingTime_;

	//----------------------------------------------------------------
	// 2D tracklet
	//----------------------------------------------------------------
	std::vector<stTracklet2DSet> vecTracklet2DSet_;
	unsigned int nNumTotalActive2DTracklet_;

	//----------------------------------------------------------------
	// 3D track
	//----------------------------------------------------------------
	bool bReceiveNewMeasurement_;
	bool bInitiationPenaltyFree_;
	unsigned int nNewTrackID_;
	unsigned int nNewTreeID_;
	std::list<CTrackTree> listTrackTree_;
	std::list<CTrack3D> listTrack3D_;

	TrackSet queueNewSeedTracks_;
	TrackSet queueActiveTrack_;
	TrackSet queuePausedTrack_;
	TrackSet queueTracksInWindow_;
	TrackSet queueTracksInBestSolution_;

	std::deque<CTrackTree*> queuePtActiveTrees_;
	std::deque<CTrackTree*> queuePtUnconfirmedTrees_;

	// for result saving
	std::deque<CTrack3DResult> queueTrackingResult_;

	// optimization related
	//hj::CGraphSolver  cGraphSolver_;
	hj::HypothesisSet queuePrevGlobalHypotheses_;
	hj::HypothesisSet queueCurrGlobalHypotheses_;
	std::vector<stHypothesisSolvingInfo> vectorSolvingInfo_;

	// for tracking result
	unsigned int nNewTargetID_;
	std::deque<PAIR_UINT> queuePairTreeIDToTargetID_;

	// evaluation
	//std::vector<std::pair<CEvaluator, int>> *vecEvaluator_;

	// for debugging
	int nCountTrackInOptimization_;
	int nCountUCTrackInOptimization_;

	//----------------------------------------------------------------
	// Visualization
	//----------------------------------------------------------------
	bool        bVisualizeResult_;
	cv::Mat     matTrackingResult_;
	std::string strVisWindowName_;
	std::vector<cv::Scalar> vecColors_;
};

}

//()()
//('')HAANJU.YOO

