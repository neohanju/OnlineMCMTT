#include "MainController.h"
#include "haanju_fileIO.hpp"
#include "haanju_string.hpp"
#include "haanju_misc.hpp"

static volatile bool gArrGDTThreadRun[MAX_NUM_SAME_THREAD];   // Grab + Detection + Tracking Thread
static volatile bool gAssociationThreadRun;                   // Association Thread
static volatile bool gGUIThreadRun;                 // Visualization Thread


//----------------------------------------------------------------
// GRAB/DETECTION/TRACKING THREAD
//----------------------------------------------------------------
struct stGDTThreadParams
{
	int                 nCamIdx;
	int                 nMaxNumGrabFail;
	bool                bDetect;
	std::string         strDetectionDir;
	hj::DETECTION_TYPE  nDetectionType;	
	CMainController*    pMainController;
	hj::CDataManager*   pDataManager;
	hj::CFrameGrabber*  pGrabber;
	CDetectorCrosstalk* pDetector;
	hj::CSCMTTracker*   pTracker;
};
static stGDTThreadParams gStGDTThreadParams[MAX_NUM_SAME_THREAD];

unsigned int __stdcall GDTWork(void *data)
{
	stGDTThreadParams *pParams = (stGDTThreadParams*)data;
	
	int     curCamIdx = pParams->nCamIdx;
	int     nFrameIndex;
	cv::Mat matFrame;
	std::vector<hj::CDetection> vecDetections;
	hj::CTrack2DResult cTrackResult;
	bool bStreamAlive = true;

	hj::printf_debug("Thread for view %d is started\n", curCamIdx);
	int cntFail = 0;

	while (gArrGDTThreadRun[curCamIdx] && pParams->pDataManager->GetRunFlag())
	{
		/* frame grabbing */
		switch (pParams->pGrabber->GrabFrame())
		{
		case HJ_GR_NORMAL:
			cntFail = 0;
			break;
		case HJ_GR_DATASET_ENDED:
			bStreamAlive = false;
			hj::printf_debug("  View %d: reach the end of the dataset\n", curCamIdx);
			break;
		default:
			hj::printf_debug("  View %d: Fail to grab a frame\n", curCamIdx);
			cntFail++;
			break;
		}		

		/* grabbing failure */
		if (cntFail > pParams->nMaxNumGrabFail || !bStreamAlive)
		{
			gArrGDTThreadRun[curCamIdx] = false;
			pParams->pMainController->TerminateProcess(curCamIdx);
			break;
		}

		matFrame    = pParams->pGrabber->GetFrame();
		nFrameIndex = pParams->pGrabber->GetFrameIndex();

		if (matFrame.empty()) { continue; }

		hj::printf_debug("  View %d: Frame is grabbed (frame no.%04d)\n", 
			curCamIdx, nFrameIndex);


		/* object detection */		
		if (pParams->bDetect)
		{
			vecDetections = pParams->pDetector->Detect(matFrame, nFrameIndex);
		}
		else
		{			
			std::string strFilePath 
				= pParams->strDetectionDir + "/" + hj::FormattedString("%04d.txt", nFrameIndex);
			vecDetections 
				= hj::ReadDetectionResultWithTxt(
					strFilePath, 
					pParams->nDetectionType);
		}
		hj::printf_debug("  View %d: Objects are detected\n", 
			curCamIdx);


		/* tracklet generation (2D tracking) */
		cTrackResult = pParams->pTracker->Track(vecDetections, matFrame, nFrameIndex);
		hj::printf_debug("  View %d: Tracklets are generated\n",
			curCamIdx);


		/* save results to data manager */
		while (pParams->pDataManager->IsFrameBufferFull(curCamIdx)
			|| pParams->pDataManager->IsTrack2DBufferFull(curCamIdx))
		{
			hj::printf_debug("  >> suspend GDT thread no.%d because of buffer full\n", curCamIdx);
			//SuspendThread(GetCurrentThread());
			::Sleep(3);
		}
		pParams->pDataManager->SetFrameImage(curCamIdx, matFrame, nFrameIndex);
		pParams->pDataManager->SetTrack2DResult(curCamIdx, cTrackResult);

		/* request association work */
		pParams->pMainController->WakeupAssociationThread(GDT_THREAD, nFrameIndex);
		hj::printf_debug("  >> suspend GDT thread no.%d\n", curCamIdx);
		SuspendThread(GetCurrentThread());		

		/* wrap-up */
		matFrame.release();
		vecDetections.clear();

		::Sleep(10);
	}

	pParams->pGrabber->Finalize();
	if (pParams->bDetect) { pParams->pDetector->Finalize(); }
	pParams->pTracker->Finalize();

	hj::printf_debug("Thread for view %d is terminated\n", curCamIdx);

	return 0;
}


//----------------------------------------------------------------
// ASSOCIATION THREAD
//----------------------------------------------------------------
struct stAssociationThreadParams
{	
	int nNumCams;	
	CMainController*   pMainController;
	hj::CDataManager*  pDataManager;
	hj::CAssociator3D* pAssociator;
	hj::CEvaluator*    pEvaluator;
};
static stAssociationThreadParams gStAssociationThreadParam;

unsigned int __stdcall AssociationWork(void *data)
{
	stAssociationThreadParams *pParams = (stAssociationThreadParams*)data;

	std::vector<cv::Mat> vecMatFrames(pParams->nNumCams);
	std::vector<hj::CTrack2DResult> vecTrack2DResults(pParams->nNumCams);
	unsigned int nFrameIdx = 0;
	unsigned int nFrameIdxRead = 0;
	bool bEvaluate = pParams->pDataManager->GetEvaluateFlag();

	hj::printf_debug("Thread for association is started\n");
	while (gAssociationThreadRun && pParams->pDataManager->GetRunFlag())
	{
		for (int camIdx = 0; camIdx < pParams->nNumCams; camIdx++)
		{
			while (!pParams->pDataManager->IsFrameBufferFull(camIdx)
				|| !pParams->pDataManager->IsTrack2DBufferFull(camIdx))
			{
				if (!pParams->pDataManager->GetRunFlag()) { break; }
				hj::printf_debug("  >> sleep association thread\n");
				::Sleep(5);
			}
			
			/* load inputs */
			pParams->pDataManager->GetFrameImage(camIdx, vecMatFrames[camIdx], &nFrameIdxRead);
			pParams->pDataManager->GetTrack2DResult(camIdx, &vecTrack2DResults[camIdx]);
			
			// check frame indices syncronization
			if (0 == camIdx)
			{
				nFrameIdx = nFrameIdxRead;
			}
			else if (nFrameIdx != nFrameIdxRead)
			{
				hj::printf_debug("  >> [ERROR] frame indices mismatch\n");
				gAssociationThreadRun = false;
				pParams->pMainController->TerminateProcess(-1); // assign temporary ID for 'A' thread
				break;
			}
		}

		// do association
		hj::CTrack3DResult curResult = 
			pParams->pAssociator->Run(vecTrack2DResults, vecMatFrames, nFrameIdx);

		// evaluation
		if (bEvaluate)
		{
			pParams->pEvaluator->SetResult(curResult);
		}

		/* request new frame/tracklets */
		pParams->pMainController->WakeupGDTThreads();

		bool bTerminated = true;
		for (int i = 0; i < MAX_NUM_SAME_THREAD; i++)
		{
			if (!gArrGDTThreadRun[i]) {	continue; }
			bTerminated = false;
			break;
		}
		if (bTerminated)
		{
			hj::printf_debug("  >> terminate association thread\n");
			break;
		}
		SuspendThread(GetCurrentThread());
		hj::printf_debug("  >> suspend association thread\n");
	}

	// evaluate
	if (bEvaluate)
	{		
		pParams->pEvaluator->Evaluate();
		pParams->pEvaluator->PrintResultToConsole();
		pParams->pEvaluator->PrintResultToFile();
		pParams->pEvaluator->Finalize();
	}

	pParams->pAssociator->Finalize();

	hj::printf_debug("Thread for association is terminated\n");	

	return 0;
}


/////////////////////////////////////////////////////////////////////////
// CMainController MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////
CMainController::CMainController()
	: bInit_(false)
	, bSystemRun_(false)	
	, bDetect_(false)
{
}


CMainController::~CMainController()
{
}


bool CMainController::Initialize(std::string _strParamXMLPath)
{
	if (bInit_) { return false; }

	/* read parameters */
	if (!cDataManager_.Initialize(_strParamXMLPath))
	{
		return false;
	}
	bDetect_   = cDataManager_.GetDetectFlag();	
	bEvaluate_ = cDataManager_.GetEvaluateFlag();
	
	/* view related */
	vecCameraIDs_ = cDataManager_.GetCameraIDs();
	numCameras_   = (int)vecCameraIDs_.size();
	vecInputFrames_.resize(numCameras_);

	/* instantiate each modules */	
	vecFrameGrabbers_.resize(numCameras_);	
	vecDetectors_.resize(numCameras_);
	vecMultiTracker2Ds_.resize(numCameras_);
	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{
		// frame grabbers
		vecFrameGrabbers_[camIdx].Initialize(
			camIdx, cDataManager_.GetFrameGrabberParams(camIdx));

		// detectors
		if (bDetect_)
		{
			vecDetectors_[camIdx].Initialize(
				cDataManager_.GetDetect2DParams(camIdx));
		}

		// 2D trackers
		vecMultiTracker2Ds_[camIdx].Initialize(
			camIdx, cDataManager_.GetTrack2DParams(camIdx));
	}

	// associator
	cAssociator3D_.Initialize(cDataManager_.GetAssociate3DParams());

	// evaluator
	if (bEvaluate_)
	{
		cEvaluator_.Initialize(cDataManager_.GetEvaluatorParams());
	}

	/* thread related */
	for (int flagIdx = 0; flagIdx < MAX_NUM_SAME_THREAD; flagIdx++)
	{
		gArrGDTThreadRun[flagIdx] = true;
	}
	gAssociationThreadRun   = true;
	gGUIThreadRun           = true;

	// grabber + detector + tracker thread
	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{
		hj::stParamFrameGrabber grabberParams = cDataManager_.GetFrameGrabberParams(camIdx);
		hj::stParamDetect2D detectorParams = cDataManager_.GetDetect2DParams(camIdx);

		gStGDTThreadParams[camIdx].nCamIdx         = camIdx;
		gStGDTThreadParams[camIdx].nMaxNumGrabFail = hj::GRABBING == grabberParams.nInputSource ? 5 : 1;
		gStGDTThreadParams[camIdx].bDetect         = bDetect_;		
		gStGDTThreadParams[camIdx].strDetectionDir = detectorParams.strDetectionDir;
		gStGDTThreadParams[camIdx].nDetectionType  = detectorParams.nDetectionType;
		gStGDTThreadParams[camIdx].pMainController = this;
		gStGDTThreadParams[camIdx].pDataManager    = &cDataManager_;
		gStGDTThreadParams[camIdx].pGrabber        = &vecFrameGrabbers_[camIdx];
		gStGDTThreadParams[camIdx].pDetector       = &vecDetectors_[camIdx];  // TODO: replace this with a real detector
		gStGDTThreadParams[camIdx].pTracker        = &vecMultiTracker2Ds_[camIdx];

		vecHGDTThreads_[camIdx] = (HANDLE)_beginthreadex(
			0, 0, &GDTWork, &gStGDTThreadParams[camIdx], CREATE_SUSPENDED, 0);
	}
	nNumGDTResults_ = 0;

	// associator thread
	gStAssociationThreadParam.nNumCams        = numCameras_;
	gStAssociationThreadParam.pMainController = this;
	gStAssociationThreadParam.pDataManager    = &cDataManager_;
	gStAssociationThreadParam.pAssociator     = &cAssociator3D_;
	gStAssociationThreadParam.pEvaluator      = &cEvaluator_;	
	hAssoicationThread_ = (HANDLE)_beginthreadex(
		0, 0, &AssociationWork, &gStAssociationThreadParam, CREATE_SUSPENDED, 0);
	
	bInit_      = true;
	bSystemRun_ = true;

	cDataManager_.SetRunFlag(bSystemRun_);

	// locks
	InitializeSRWLock(&lockGDT_);
	InitializeSRWLock(&lockFrameIdx_);

	return true;
}


void CMainController::Finalize()
{
	if (!bInit_) { return; }

	/* Grab/Detect/Track related */
	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{
		// terminate thread
		gArrGDTThreadRun[camIdx] = false;
		if (NULL != vecHGDTThreads_[camIdx])
		{
			WaitForSingleObject(vecHGDTThreads_[camIdx], INFINITE);
			CloseHandle(vecHGDTThreads_[camIdx]);
		}
	}

	// wait until the association thread is terminated
	gAssociationThreadRun = false;
	if (NULL != hAssoicationThread_)
	{
		WaitForSingleObject(hAssoicationThread_, INFINITE);
		CloseHandle(hAssoicationThread_);
	}
}

void CMainController::Reset()
{

}


void CMainController::Run()
{
	// resume Grab/Detect/Track threads
	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{
		ResumeThread(vecHGDTThreads_[camIdx]);
	}

	// wait until the association thread is terminated
	WaitForSingleObject(hAssoicationThread_, INFINITE);
}


bool CMainController::TerminateProcess(int _nThreadID)
{
	hj::printf_debug("  >> Overall process is terminated by thread %d\n", _nThreadID);
	cDataManager_.SetRunFlag(false);

	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{		
		gArrGDTThreadRun[camIdx] = false;
		hj::printf_debug("  >> resume GDT thread no.%d\n", camIdx);
		ResumeThread(vecHGDTThreads_[camIdx]);
	}
	
	gAssociationThreadRun = false;
	ResumeThread(hAssoicationThread_);

	gGUIThreadRun = false;

	return true;
}


bool CMainController::WakeupGDTThreads()
{
	for (int camIdx = 0; camIdx < numCameras_; camIdx++)
	{
		if (!gArrGDTThreadRun[camIdx]) { return false; }
		hj::printf_debug("  >> resume GDT thread no.%d\n", camIdx);
		ResumeThread(vecHGDTThreads_[camIdx]);
	}
	return true;
}


bool CMainController::WakeupAssociationThread(HJ_THREAD_TYPE _threadType, unsigned int _frameIdx)
{
	bool bGDTFull = false;	
	bool bTerminate = false;

	//------------------------------------------		
	AcquireSRWLockExclusive(&lockGDT_);		
	//------------------------------------------
	if (0 == nNumGDTResults_)
	{
		nFrameIdx_ = _frameIdx;
	}
	else if (nFrameIdx_ != _frameIdx)
	{
		hj::printf_debug("  >> [ERROR] Frame indices are not matched!\n");
		bTerminate = true;
	}
	//------------------------------------------	
	ReleaseSRWLockExclusive(&lockGDT_);
	//------------------------------------------

	if (bTerminate)
	{
		TerminateProcess(-2);  // assign temporary ID for main controller
		hj::printf_debug("  >> resume association thread to terminate\n");
		ResumeThread(hAssoicationThread_);
		return true;
	}

	switch (_threadType)
	{
	case GDT_THREAD:
		//------------------------------------------
		AcquireSRWLockExclusive(&lockGDT_);
		//------------------------------------------
		nNumGDTResults_++;
		if (nNumGDTResults_ >= numCameras_) { bGDTFull = true; }
		//------------------------------------------
		ReleaseSRWLockExclusive(&lockGDT_);
		//------------------------------------------
		break;
	default:
		break;
	}
	
	if (bGDTFull)
	{
		//------------------------------------------
		AcquireSRWLockExclusive(&lockGDT_);
		//------------------------------------------
		nNumGDTResults_ = 0;
		//------------------------------------------
		ReleaseSRWLockExclusive(&lockGDT_);
		//------------------------------------------

		hj::printf_debug("  >> resume association thread\n");
		ResumeThread(hAssoicationThread_);
		return true;
	}
	return false;
}


//()()
//('')HAANJU.YOO

