/******************************************************************************
* Title        : CDataManager
* Author       : Haanju Yoo
* Initial Date : 2016.09.20
* Version Num. : 0.9
* Description  : managing global variables with semaphore
******************************************************************************/

#pragma once

#include <process.h>
#include <windows.h>
#include "types.hpp"

namespace hj
{

class CDataManager
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CDataManager();	
	~CDataManager();

	bool Initialize(std::string _strParameterFilePath);
	bool Finalize();
	bool UpdateParameters(std::string _strParameterFilePath);

	bool Ready() { return bDataReady_; }
	bool IsFrameBufferFull(const int _camIdx);
	bool IsDetectionBufferFull(const int _camIdx);
	bool IsTrack2DBufferFull(const int _camIdx);
	bool IsTrack3DBufferFull();
	bool IsGUIBufferFull();

	/* set result */
	bool SetFrameImage(const int _camIdx, const cv::Mat _matImage, unsigned int _frameIdx);
	bool SetDetectionResult(const int _camIdx, const hj::DetectionSet _vecDetections);
	bool SetTrack2DResult(const int _camIdx, const hj::CTrack2DResult _cTrack2DResult);
	bool SetTrack3DResult(const hj::CTrack3DResult _cTrack3DResult);
	bool SetGUIAssignResult(const std::vector<std::pair<int, int>> _vecAssignedIDtoRealID);

	/* set misc. */
	void SetRunFlag(bool _flag);

	/* get parameters */
	hj::stParamFrameGrabber  GetFrameGrabberParams(const int _camIdx);
	hj::stParamDetect2D      GetDetect2DParams(const int _camIdx);
	hj::stParamTrack2D       GetTrack2DParams(const int _camIdx);
	hj::stParamAssociate3D   GetAssociate3DParams();	
	hj::stParamEvaluator     GetEvaluatorParams();

	/* get results */
	bool GetFrameImage(const int _camIdx, cv::Mat &_receiver, unsigned int *_frameIdx = NULL);
	bool GetDetectionResult(const int _camIdx, hj::DetectionSet *_receiver);
	bool GetTrack2DResult(const int _camIdx, hj::CTrack2DResult *_receiver);
	bool GetTrack3DResult(hj::CTrack3DResult *_receiver);
	bool GetGUIAssignResult(std::vector<std::pair<int, int>> *_receiver);

	/* read functions */
	std::vector<int> GetCameraIDs();
	//hj::CCalibrationInfo* GetCalibInfoPt(int _camIdx);

	bool GetDetectFlag();	
	bool GetEvaluateFlag();
	bool GetRunFlag();

private:
	bool ReadParametersFromXML(std::string _strFilepath);
	bool ReadCalibrationInfos(std::string _strDir, std::vector<int> _vecCamIDs);	
	bool ReadProjectionSensitivity(cv::Mat &_matSensitivity, std::string _strFilepath);
	bool ReadDistanceFromBoundary(cv::Mat &_matDistance, std::string _strFilepath);

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
private:
	bool bDataReady_;

	hj::stViewInformation  stViewInfos_;
	std::vector<hj::stParamFrameGrabber>  vecFrameGrabberParams_;
	std::vector<hj::stParamDetect2D>      vecDetect2DParams_;
	std::vector<hj::stParamTrack2D>       vecTrack2DParams_;
	hj::stParamAssociate3D                associate3DParams_;	
	hj::stParamEvaluator                  evaluatorParams_;
	SRWLOCK lockParamAccess_;

	bool bRealtimeOperation_;
	bool bRecord_;
	bool bDetect_;	
	bool bEvaluate_;

	std::string strInputSource_;
	std::string strDatasetPath_;
	std::string strDetectorPath_;
	std::string strResultPath_;	
	double dFrameRate_;

	/* calibration related */
	std::string strCalibrationPath_;
	std::vector<hj::CCalibrationInfo>      vecCalibrationInfo_;	

	/* grabber related */	
	std::vector<cv::Mat>               vecMatInputFrameBuffer_;
	std::vector<SRWLOCK>               vecLockFrameImage_;
	std::vector<unsigned int>          vecNFrameIndex_;
	volatile bool                      arrFrameBufferFull_[MAX_NUM_SAME_THREAD];

	/* detector related */
	std::vector<hj::DetectionSet>      vecvecDetectionResult_;
	std::vector<SRWLOCK>               vecLockDetectionResult_;
	volatile bool                      arrDetectionBufferFull_[MAX_NUM_SAME_THREAD];

	/* 2D tracker related */
	std::vector<hj::CTrack2DResult>    vecCTrack2DResult_;
	std::vector<SRWLOCK>               vecLockTrack2DResult_;
	volatile bool                      arrTrack2DBufferFull_[MAX_NUM_SAME_THREAD];
	
	/* 3D associator related */
	hj::CTrack3DResult                 cTrack3DResult_;
	SRWLOCK                            lockTrack3DResult_;
	volatile bool                      bTrack3DBufferFull_;

	/* GUI */
	std::vector<std::pair<int, int>>   vecAssignedIDtoRealID_;
	SRWLOCK                            lockGUI_;
	volatile bool                      bGUIBufferFull_;

	/* overall operation */
	volatile bool bSystemRun_;
	SRWLOCK       lockRunFlag_;

};

}

//()()
//('')HAANJU.YOO


