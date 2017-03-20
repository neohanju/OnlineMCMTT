#include "types_3D.h"

namespace hj
{

/////////////////////////////////////////////////////////////////////////
// CTrackletCombination MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////
CTrackletCombination& CTrackletCombination::operator=(const CTrackletCombination &a)
{
	this->tracklets = a.tracklets;	
	this->numTracklets = a.numTracklets;
	return *this;
}

bool CTrackletCombination::operator==(const CTrackletCombination &a)
{
	if (this->numTracklets != a.numTracklets) { return false; }
	if (this->tracklets.size() != a.tracklets.size()) { return false; }
	for (unsigned int camIdx = 0; camIdx < this->tracklets.size(); camIdx++)
	{
		if (this->tracklets[camIdx] != a.tracklets[camIdx]) { return false; }
	}
	return true;
}


/************************************************************************
 Method Name: set
 Description:
	- set the tracklet in the combination
 Input Arguments:
	- camIdx: target camera index
	- tracklet: target tracklet
 Return Values:
	- void
************************************************************************/
void CTrackletCombination::set(unsigned int camIdx, CTracklet2D *tracklet)
{
	if (this->tracklets[camIdx] == tracklet) { return; }
	if (NULL == this->tracklets[camIdx])
	{
		numTracklets++;
	}
	else if (NULL == tracklet)
	{
		numTracklets--;
	}
	this->tracklets[camIdx] = tracklet;
}

/************************************************************************
 Method Name: get
 Description:
	- get the tracklet in the combination
 Input Arguments:
	- camIdx: target camera index
 Return Values:
	- tracklet of camIdx in the combination
************************************************************************/
CTracklet2D* CTrackletCombination::get(unsigned int camIdx)
{
	return this->tracklets[camIdx];
}

/************************************************************************
 Method Name: print
 Description:
	- print current tracklet combination into the console
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CTrackletCombination::print()
{
	printf("[");
	for (unsigned int camIdx = 0; camIdx < this->tracklets.size(); camIdx++)
	{
		if (NULL != tracklets[camIdx])
		{
			printf("%d", tracklets[camIdx]->id);
		}
		else
		{
			printf("x");
		}

		if (camIdx < this->tracklets.size() - 1)
		{
			printf(",");
		}
		else
		{
			printf("]\n");
		}
	}
}

/************************************************************************
 Method Name: checkCoupling
 Description:
	- check wheather tracklets are compatible or not
 Input Arguments:
	- compCombination: the target of comparison
 Return Values:
	- true: incompatible/ false: compatible
************************************************************************/
bool CTrackletCombination::checkCoupling(CTrackletCombination &compCombination)
{
	for (unsigned int camIdx = 0; 
		camIdx < std::min(this->tracklets.size(), compCombination.tracklets.size()); 
		camIdx++)
	{
		if (this->tracklets[camIdx] == compCombination.tracklets[camIdx] && NULL != this->tracklets[camIdx])
		{
			return true;
		}
	}
	return false;
}

/************************************************************************
 Method Name: size
 Description:
	- return the number of tracklets in the combination
 Input Arguments:
	- none
 Return Values:
	- number of tracklets in the combination
************************************************************************/
unsigned int CTrackletCombination::size()
{
	return this->numTracklets;
}


/////////////////////////////////////////////////////////////////////////
// Track3D MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////
CTrack3D::CTrack3D(int numCameras)
	: id(0)
	, bActive(true)
	, bValid(true)
	, tree(NULL)
	, parentTrack(NULL)
	, timeStart(0)
	, timeEnd(0)
	, timeGeneration(0)
	, duration(1)
	, costTotal(0.0)
	, costReconstruction(0.0)
	, costLink(0.0)
	, costEnter(0.0)
	, costExit(0.0)
	, loglikelihood(0.0)
	, GTProb(0.0)
	, BranchGTProb(0.0)
	, bWasBestSolution(true)
	, bCurrentBestSolution(false)
	, bNewTrack(true)
	, curTracklet2Ds(numCameras)
{
	tracklet2DIDRecord.resize(numCameras);
	timeTrackletEnded.resize(numCameras, 0);
	lastTrackletLocation3D.resize(numCameras, Point3D(0.0, 0.0, 0.0));
	lastTrackletSensitivity.resize(numCameras, 0.0);
	lastRGBFeature.resize(numCameras);
}

CTrack3D::~CTrack3D()
{

}

/************************************************************************
 Method Name: Initialize
 Description:
	- initialize the current track
 Input Arguments:
	- trackletCombination: current tracklet combination of the current track
	- id: its ID
	- timeGeneration: birth time of the track (not equal to starting time
	  because it can be the child track, which is generated much later than
	  its starting time)
	- parentTrack: its parent track. It can be NULL for seed tracks
 Return Values:
	- none
************************************************************************/
void CTrack3D::Initialize(CTrackletCombination &trackletCombination, unsigned int id, unsigned int timeGeneration, CTrack3D *parentTrack)
{
	this->id = id;
	this->curTracklet2Ds = trackletCombination;
	for (int camIdx = 0; camIdx < this->curTracklet2Ds.numCameras(); camIdx++)
	{
		this->timeTrackletEnded[camIdx] = 0;
		this->lastTrackletLocation3D[camIdx] = Point3D(0.0, 0.0, 0.0);
		this->lastTrackletSensitivity[camIdx] = 0;
	}
	this->numOutpoint = 0;
	this->timeEnd = timeGeneration;
	this->timeGeneration = timeGeneration;
	if (NULL == parentTrack)
	{
		this->timeStart = timeGeneration;		
	}
	else
	{
		this->timeStart = parentTrack->timeStart;
		this->tree = parentTrack->tree;
		this->parentTrack = parentTrack;
		this->costEnter = parentTrack->costEnter;
		this->loglikelihood = parentTrack->loglikelihood;
	}
	this->duration = this->timeEnd - this->timeStart + 1;	
}

/************************************************************************
 Method Name: RemoveFromTree
 Description:
	- remove track from the tree. children track's parent track will be
	  modified.
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CTrack3D::RemoveFromTree()
{
	// find valid parent (for children's adoption)
	CTrack3D *newParentTrack = this->parentTrack;
	while (NULL != newParentTrack && !newParentTrack->bValid) { newParentTrack = newParentTrack->parentTrack; }

	// remove from children's parent pointer and adopt them to new parent
	for (std::deque<CTrack3D*>::iterator childTrackIter = this->childrenTrack.begin();
		childTrackIter != this->childrenTrack.end();
		childTrackIter++)
	{
		(*childTrackIter)->parentTrack = newParentTrack;
		if (NULL != (*childTrackIter)->parentTrack)
		{
			(*childTrackIter)->parentTrack->childrenTrack.push_back(*childTrackIter);
		}
	}

	// remove from parent's children list
	if (NULL != this->parentTrack)
	{
		for (std::deque<CTrack3D*>::iterator childTrackIter = this->parentTrack->childrenTrack.begin();
			childTrackIter != this->parentTrack->childrenTrack.end();
			childTrackIter++)
		{
			if ((*childTrackIter)->id != this->id) { continue; }
			this->parentTrack->childrenTrack.erase(childTrackIter);
			break;
		}
	}
}

//#define KALMAN_PROCESSNOISE_SIG (1.0E-5)
//#define KALMAN_MEASUREMENTNOISE_SIG (1.0E-5)
//#define KALMAN_POSTERROR_COV (0.1)
//#define KALMAN_CONFIDENCE_LEVEN (9)
//void Track3D::SetKalmanFilter(PSN_Point3D &initialPoint)
//{
//	this->KF.init(6, 3, 0);
//	this->KFMeasurement = cv::Mat(3, 1, CV_32FC1);
//
//	cv::setIdentity(this->KF.transitionMatrix); // [1,0,0,1,0,0; 0,1,0,0,1,0; 0,0,1,0,0,1; 0,0,0,1,0,0, ...]
//	this->KF.transitionMatrix.at<float>(0, 3) = 1.0f;
//	this->KF.transitionMatrix.at<float>(1, 4) = 1.0f;
//	this->KF.transitionMatrix.at<float>(2, 5) = 1.0f;
//
//	cv::setIdentity(this->KF.measurementMatrix);
//	cv::setIdentity(this->KF.processNoiseCov, cv::Scalar::all(KALMAN_PROCESSNOISE_SIG));
//	cv::setIdentity(this->KF.measurementNoiseCov, cv::Scalar::all(KALMAN_MEASUREMENTNOISE_SIG));
//	cv::setIdentity(this->KF.errorCovPost, cv::Scalar::all(KALMAN_POSTERROR_COV));
//
//	this->KF.statePost.at<float>(0, 0) = (float)initialPoint.x;
//	this->KF.statePost.at<float>(1, 0) = (float)initialPoint.y;
//	this->KF.statePost.at<float>(2, 0) = (float)initialPoint.z;
//	this->KF.statePost.at<float>(3, 0) = 0.0f;
//	this->KF.statePost.at<float>(4, 0) = 0.0f;
//	this->KF.statePost.at<float>(5, 0) = 0.0f;
//	cv::Mat curKFPrediction = this->KF.predict();
//}


/////////////////////////////////////////////////////////////////////////
// CTrackTree MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////
CTrackTree::CTrackTree()
	: id(0)
	, timeGeneration(0)
	, bValid(true)
	, bConfirmed(false)
	//	, numMeasurements(0)
	//	, maxGTProb(0.0)
{
}

CTrackTree::~CTrackTree()
{
}

/************************************************************************
 Method Name: Initialize
 Description:
	- initialize the current track tree with the seed track. also save the
	  instance of the current track tree into the list container.
 Input Arguments:
	- seedTrack: seed track of the curren track tree
	- id: its ID
	- timeGeneration: track tree's birth time
	- treeList: saving container for instance of the current track tree
 Return Values:
	- none
************************************************************************/
void CTrackTree::Initialize(CTrack3D *seedTrack, unsigned int id, unsigned int timeGeneration, std::list<CTrackTree> &treeList)
{
	this->id = id;
	this->timeGeneration = timeGeneration;
	this->tracks.push_back(seedTrack);
	treeList.push_back(*this);
	seedTrack->tree = &treeList.back();
}

/************************************************************************
 Method Name: ResetGlobalTrackProbInTree
 Description:
	- reset GTP in the entire tracks in the current track tree
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CTrackTree::ResetGlobalTrackProbInTree()
{
	for (std::deque<CTrack3D*>::iterator trackIter = this->tracks.begin();
		trackIter != this->tracks.end();
		trackIter++)
	{
		(*trackIter)->GTProb = 0.0;
	}
}

/************************************************************************
 Method Name: FindPruningPoint
 Description:
	- find the lowest track in the branch that is the pruing point. children
	  of it will be used in pruing
 Input Arguments:
	- timeWindowStart: where the pruning is starting
	- rootOfBranch: initial searching track. search is going to dive into the
	  children of this track
 Return Values:
	- track which has the children will be condiered in the pruning step
************************************************************************/
CTrack3D* CTrackTree::FindPruningPoint(unsigned int timeWindowStart, CTrack3D *rootOfBranch)
{
	if (NULL == rootOfBranch)
	{
		if (0 == this->tracks.size()) { return NULL; }
		rootOfBranch = this->tracks[0];
	}
	if (rootOfBranch->timeGeneration >= timeWindowStart) { return NULL; }

	// if more than one child placed inside of processing window, then the current node is the pruning point
	CTrack3D *curPruningPoint = NULL;
	for (TrackSet::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		curPruningPoint = FindPruningPoint(timeWindowStart, *trackIter);
		if (NULL == curPruningPoint) { break; }
	}

	if (NULL == curPruningPoint) { return rootOfBranch; }
	return curPruningPoint;
}

/************************************************************************
 Method Name: CheckBranchContainsBestSolution
 Description:
	- check wheather the current branch has a track in the best solution or not
 Input Arguments:
	- rootOfBranch: the root of current branch
 Return Values:
	- wheather the current branch has a track in the best solution or not
************************************************************************/
bool CTrackTree::CheckBranchContainsBestSolution(CTrack3D *rootOfBranch)
{
	if (rootOfBranch->bCurrentBestSolution) { return true; }

	for (TrackSet::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		if (CheckBranchContainsBestSolution(*trackIter)) { return true; }
	}
	return false;
}

/************************************************************************
 Method Name: SetValidityFlagInTrackBranch
 Description:
	- Recursively set validity flags to 'bValid'
 Input Arguments:
	- rootOfBranch: the root of current branch
	- bValid: the value of validity flag
 Return Values:
	- none
************************************************************************/
void CTrackTree::SetValidityFlagInTrackBranch(CTrack3D* rootOfBranch, bool bValid)
{
	rootOfBranch->bValid = bValid;
	for (std::deque<CTrack3D*>::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		SetValidityFlagInTrackBranch(*trackIter, bValid);
	}
}

/************************************************************************
 Method Name: GetTracksInBranch
 Description:
	- Recursively gather pointers of tracks in descendants
 Input Arguments:
	- rootOfBranch: the root of current branch
	- queueOutput: the pointer vector of tracks in the branch
 Return Values:
	- none
************************************************************************/
void CTrackTree::GetTracksInBranch(CTrack3D* rootOfBranch, std::deque<CTrack3D*> &queueOutput)
{
	queueOutput.push_back(rootOfBranch);
	for (std::deque<CTrack3D*>::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		GetTracksInBranch(*trackIter, queueOutput);
	}
}

/************************************************************************
 Method Name: MakeTreeNodesWithChildren
 Description:
	- Recursively generate track node information
 Input Arguments:
	- queueChildrenTracks: queue of children tracks
	- parentNodeIdx: parent node index
	- outQueueNodes: output node index queue
 Return Values:
	- none
************************************************************************/
void CTrackTree::MakeTreeNodesWithChildren(std::deque<CTrack3D*> queueChildrenTracks, const int parentNodeIdx, std::deque<stTrackInTreeInfo> &outQueueNodes)
{
	for (std::deque<CTrack3D*>::iterator trackIter = queueChildrenTracks.begin();
		trackIter != queueChildrenTracks.end();
		trackIter++)
	{
		stTrackInTreeInfo newInfo;
		newInfo.id = (*trackIter)->id;
		newInfo.parentNode = parentNodeIdx;
		newInfo.timeGenerated = (*trackIter)->timeGeneration;
		newInfo.GTP = (float)(*trackIter)->GTProb;
		//if ((*trackIter)->bValid) { outQueueNodes.push_back(newInfo); }
		outQueueNodes.push_back(newInfo);
		MakeTreeNodesWithChildren((*trackIter)->childrenTrack, (int)outQueueNodes.size(), outQueueNodes);
	}
}

/************************************************************************
 Method Name: GTProbOfBrach
 Description:
	- Recursively sum global track probabilities of track branch
 Input Arguments:
	- rootOfBranch: the root of current branch
 Return Values:
	- sum of global track probability
************************************************************************/
double CTrackTree::GTProbOfBrach(CTrack3D *rootOfBranch)
{
	double GTProb = rootOfBranch->GTProb;
	for (std::deque<CTrack3D*>::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		GTProb += GTProbOfBrach(*trackIter);
	}
	return GTProb;
}

/************************************************************************
 Method Name: MaxGTProbOfBrach
 Description:
	- find the maximum global track probability in the branch and set
	  that value into 'BranchGTProb' property
 Input Arguments:
	- rootOfBranch: the seed track of the branch
 Return Values:
	- maximum value of global track probability in the branch
************************************************************************/
double CTrackTree::MaxGTProbOfBrach(CTrack3D *rootOfBranch)
{
	double MaxGTProb = rootOfBranch->GTProb;
	for (std::deque<CTrack3D*>::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		double curGTProb = MaxGTProbOfBrach(*trackIter);
		if (MaxGTProb < curGTProb) { MaxGTProb = curGTProb; }
	}
	rootOfBranch->BranchGTProb = MaxGTProb;
	return MaxGTProb;
}

/************************************************************************
 Method Name: InvalidateBranchWithMinGTProb
 Description:
	- invalidate the tracks in the branch which have smaller GTP then minGTProb
 Input Arguments:
	- rootOfBranch: the seed track of the branch
	- minGTProb: GTP threshold
 Return Values:
	- none
************************************************************************/
void CTrackTree::InvalidateBranchWithMinGTProb(CTrack3D *rootOfBranch, double minGTProb)
{
	if (rootOfBranch->GTProb < minGTProb) { rootOfBranch->bValid = false; }
	for (std::deque<CTrack3D*>::iterator trackIter = rootOfBranch->childrenTrack.begin();
		trackIter != rootOfBranch->childrenTrack.end();
		trackIter++)
	{
		InvalidateBranchWithMinGTProb(*trackIter, minGTProb);
	}
}

/************************************************************************
 Method Name: FindMaxGTProbBranch
 Description:
	- find the track having the maximum GTB int the branch
 Input Arguments:
	- queueChildrenTracks: queue of children tracks
 Return Values:
	- none
************************************************************************/
CTrack3D* CTrackTree::FindMaxGTProbBranch(CTrack3D* branchSeedTrack, size_t timeIndex)
{
	if (branchSeedTrack->timeGeneration >= timeIndex) { return NULL; }

	CTrack3D* maxGTProbChild = NULL;
	for (std::deque<CTrack3D*>::iterator trackIter = branchSeedTrack->childrenTrack.begin();
		trackIter != branchSeedTrack->childrenTrack.end();
		trackIter++)
	{
		CTrack3D* curGTProbChild = FindMaxGTProbBranch((*trackIter), timeIndex);
		if (NULL == curGTProbChild) { continue; }
		if (NULL == maxGTProbChild || curGTProbChild > maxGTProbChild) { maxGTProbChild = curGTProbChild; }
	}
	return NULL == maxGTProbChild ? branchSeedTrack : maxGTProbChild;
}

/************************************************************************
 Method Name: FindOldestTrackInBranch
 Description:
	- find the oldest (but after some frame index 't') track in the branch
 Input Arguments:
	- trackInBranch: tracks in the current branch
	- nMostPreviousFrameIdx: 't' in the description
 Return Values:
	- the oldest track in the branch
************************************************************************/
CTrack3D* CTrackTree::FindOldestTrackInBranch(CTrack3D *trackInBranch, int nMostPreviousFrameIdx)
{
	CTrack3D *oldestTrack = trackInBranch;
	while (true)
	{
		if (NULL == oldestTrack->parentTrack) { break; }
		if (nMostPreviousFrameIdx >= (int)oldestTrack->parentTrack->timeGeneration) { break; }
		oldestTrack = oldestTrack->parentTrack;
	}
	return oldestTrack;
}

} // end of namespace 'haaanju'


//()()
//('')HAANJU.YOO


