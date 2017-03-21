/******************************************************************************
* Title        : CMainController
* Author       : Haanju Yoo
* Initial Date : 2016.10.03
* Version Num. : 0.9
* Description  : managing overall tracking procedure (in a single thread)
******************************************************************************/

#pragma once

#include <atomic>
#include <process.h>
#include "DataManager.h"
#include "FrameGrabber.h"
#include "SCMTTracker.h"
#include "Associator3D.h"
#include "Evaluator.h"

typedef enum { GDT_THREAD = 0, A_THREAD, GUI_THREAD } HJ_THREAD_TYPE;

class CMainController
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CMainController();
	~CMainController();

	bool Initialize(std::string _strParamXMLPath);
	void Finalize();
	void Reset();
	void Run();

	bool TerminateProcess(int _nThreadID);

	bool WakeupGDTThreads();
	bool WakeupAssociationThread(HJ_THREAD_TYPE _threadType, unsigned int _frameIdx);	

private:

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:

private:
	hj::CDataManager                cDataManager_;	
	std::vector<hj::CFrameGrabber>  vecFrameGrabbers_;	
	std::vector<hj::CSCMTTracker>   vecMultiTracker2Ds_;
	hj::CAssociator3D               cAssociator3D_;	
	hj::CEvaluator                  cEvaluator_;

	volatile bool bInit_;
	volatile bool bSystemRun_;
	int numCameras_;
	std::vector<int>     vecCameraIDs_;
	std::vector<cv::Mat> vecInputFrames_;
	bool bDetect_;
	bool bEvaluate_;

	HANDLE vecHGDTThreads_[MAX_NUM_SAME_THREAD];
	HANDLE hAssoicationThread_;	
	std::atomic<int> nNumGDTResults_;
	std::atomic<unsigned int> nFrameIdx_;
	SRWLOCK lockGDT_, lockFrameIdx_;
};

//()()
//('')HAANJU.YOO

