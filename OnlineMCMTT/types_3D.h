/******************************************************************************
* Title        : types_3D
* Author       : Haanju Yoo
* Initial Date : 2013.08.29 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  : Basic types for 3D association.
******************************************************************************
*               .__                           __.
*                \ `\~~---..---~~~~~~--.---~~| /
*                 `~-.   `                   .~         _____
*                     ~.                .--~~    .---~~~    /
*                      / .-.      .-.      |  <~~        __/
*                     |  |_|      |_|       \  \     .--'
*                    /-.      -       .-.    |  \_   \_
*                    \-'   -..-..-    `-'    |    \__  \_
*                     `.                     |     _/  _/
*                     ~-                .,-\   _/  _/
*                      /                 -~~~~\ /_  /_
*                     |               /   |    \  \_  \_
*                     |   /          /   /      | _/  _/
*                     |  |          |   /    .,-|/  _/
*                     )__/           \_/    -~~~| _/
*                       \                      /  \
*                        |           |        /_---`
*                        \    .______|      ./
*                        (   /        \    /
*                        `--'          /__/
*
******************************************************************************/

#ifndef __HAANJU_TYPES_3D_H__
#define __HAANJU_TYPES_3D_H__

#include "Point3DSmoother.h"

namespace hj
{

typedef std::pair<Point2D, unsigned int> Point2D_CamIdx;

class CTracklet2D
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructor
	CTracklet2D(int numCameras)
	{	
		id = 0;
		camIdx = 0;
		numCam = numCameras;
		bActivated = false;

		timeStart = 0;
		timeEnd = 0;
		duration = 0;

		bAssociableNewMeasurement.resize(numCam);
	}

	// desctructor
	~CTracklet2D()
	{
		if (!RGBFeatureHead.empty()) { RGBFeatureHead.release(); }
		if (!RGBFeatureTail.empty()) { RGBFeatureTail.release(); }
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	unsigned int camIdx;
	unsigned int numCam;
	bool bActivated;

	// spatial information
	std::deque<Rect> rects;
	std::deque<Line3D> backprojectionLines;

	// temporal information
	unsigned int timeStart;
	unsigned int timeEnd;
	unsigned int duration;

	// appearance 
	cv::Mat RGBFeatureHead;
	cv::Mat RGBFeatureTail;

	// location in 3D
	Point3D currentLocation3D;

	// matching related
	std::vector<std::vector<bool>> bAssociableNewMeasurement;
};


struct stTracklet2DSet
{
	std::list<CTracklet2D>   tracklets;
	std::deque<CTracklet2D*> activeTracklets; // for fast searching
	std::deque<CTracklet2D*> newMeasurements; // for generating seeds
};


class CTrackletCombination
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CTrackletCombination(int numCameras)
	{
		numCam = numCameras;
		tracklets.resize(numCam, NULL);
		numTracklets = 0;
	}

	// operator
	CTrackletCombination& operator=(const CTrackletCombination &a);
	bool operator==(const CTrackletCombination &a);

	// methods
	void set(unsigned int camIdx, CTracklet2D *tracklet);
	CTracklet2D* get(unsigned int camIdx);
	void print(void);
	bool checkCoupling(CTrackletCombination &a);
	unsigned int size(void);
	int numCameras(void) { return numCam; }

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
private:
	std::vector<CTracklet2D*> tracklets;
	unsigned int numTracklets;
	int numCam;
};


class CReconstruction
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CReconstruction(int numCameras) : tracklet2Ds(numCameras)
	{
		numCam = numCameras;
		bIsMeasurement = false;
		point = smoothedPoint = velocity = Point3D(0.0, 0.0, 0.0);
		maxError = 0.0;
		costReconstruction = costSmoothedPoint = costLink = DBL_MAX;
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	int numCam;
	bool bIsMeasurement;	
	CTrackletCombination tracklet2Ds;
	std::vector<Point3D> rawPoints;
	Point3D point;
	Point3D smoothedPoint;
	Point3D velocity;
	double maxError;
	double costReconstruction;
	double costSmoothedPoint;
	double costLink;
};

class CTrackTree;
class CTrack3D
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CTrack3D(int numCameras);
	~CTrack3D();

	void Initialize(
		CTrackletCombination &trackletCombination, 
		unsigned int id,
		unsigned int timeGeneration, 
		CTrack3D *parentTrack = NULL);
	void RemoveFromTree();
	//void SetKalmanFilter(PSN_Point3D &initialPoint);

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	CTrackletCombination curTracklet2Ds;
	std::vector<std::deque<unsigned int>> tracklet2DIDRecord;
	bool bActive; // for update and branching
	bool bValid;  // for deletion

	// for tree
	CTrackTree *tree;
	CTrack3D *parentTrack;
	std::deque<CTrack3D*> childrenTrack;

	// temporal information
	unsigned int timeStart;
	unsigned int timeEnd;
	unsigned int timeGeneration;
	unsigned int duration;

	// reconstruction related
	std::deque<CReconstruction> reconstructions;

	// smoothing related
	CPoint3DSmoother smoother;

	// cost
	double costTotal;
	double costReconstruction;
	double costLink;
	double costEnter;
	double costExit;
	double costRGB;

	// loglikelihood
	double loglikelihood;

	// global track probability
	double GTProb;
	double BranchGTProb;
	bool bWasBestSolution;
	bool bCurrentBestSolution;

	// HO-HMT
	bool bNewTrack;

	// tracklet related	
	std::vector<unsigned int> timeTrackletEnded;
	std::vector<Point3D> lastTrackletLocation3D;
	std::vector<double> lastTrackletSensitivity;
	std::vector<cv::Mat> lastRGBFeature;

	// termination
	int numOutpoint; // count tpoints
};

typedef std::deque<CTrack3D*> TrackSet;
struct stTrackInTreeInfo
{
	int id;
	int parentNode;
	int timeGenerated;
	float GTP;
};

class CTrackTree
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CTrackTree();
	~CTrackTree();

	void Initialize(CTrack3D *seedTrack, unsigned int id, unsigned int timeGeneration, std::list<CTrackTree> &treeList);
	void ResetGlobalTrackProbInTree();

	CTrack3D* FindPruningPoint(unsigned int timeWindowStart, CTrack3D *rootOfBranch = NULL);
	static bool CheckBranchContainsBestSolution(CTrack3D *rootOfBranch);
	static void SetValidityFlagInTrackBranch(CTrack3D* rootOfBranch, bool bValid);
	static void GetTracksInBranch(CTrack3D* rootOfBranch, std::deque<CTrack3D*> &queueOutput);
	static void MakeTreeNodesWithChildren(std::deque<CTrack3D*> queueChildrenTracks, const int parentNodeIdx, std::deque<stTrackInTreeInfo> &outQueueNodes);
	static double GTProbOfBrach(CTrack3D *rootOfBranch);
	static double MaxGTProbOfBrach(CTrack3D *rootOfBranch);
	static void InvalidateBranchWithMinGTProb(CTrack3D *rootOfBranch, double minGTProb);
	static CTrack3D* FindMaxGTProbBranch(CTrack3D* branchSeedTrack, size_t timeIndex);
	static CTrack3D* FindOldestTrackInBranch(CTrack3D *trackInBranch, int nMostPreviousFrameIdx);
	//static bool CheckConnectivityOfTrees(TrackTree *tree1, TrackTree *tree2);

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	unsigned int timeGeneration;
	bool bValid;
	bool bConfirmed;	// for confirmation
	std::deque<CTrack3D*> tracks; // seed at the first

	//unsigned int numMeasurements;
	//std::deque<stTracklet2DInfo> tracklet2Ds[NUM_CAM];

	// pruning related
	//double maxGTProb;
};


}


#endif

//()()
//('')HAANJU.YOO


