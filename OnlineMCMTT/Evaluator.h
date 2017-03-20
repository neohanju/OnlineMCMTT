#pragma once

#include "types_3D.h"

namespace hj
{

typedef std::pair<unsigned int, hj::Point3D> pointInfo;
typedef std::deque<pointInfo> pointInfoSet;

struct stEvaluationResult
{
	double fMOTA;
	double fMOTP;
	double fMOTAL;
	double fRecall;
	double fPrecision;
	double fMissTargetPerGroundTruth;
	double fFalseAlarmPerGroundTruth;
	double fFalseAlarmPerFrame;
	int nMissed;
	int nFalsePositives;
	int nIDSwitch;
	int nMostTracked;
	int nPartilalyTracked;
	int nMostLost;
	int nFragments;
};

class CEvaluator
{
public:
	CEvaluator(void);
	~CEvaluator(void);

	void Initialize(stParamEvaluator _stParams);
	void Finalize(void);

	void SetResult(hj::CTrack3DResult &trackResult);
	void SetResult(hj::TrackSet &trackSet, unsigned int timeIdx);	
	void LoadResultFromText(std::string strFilepath);
	void Evaluate(void);
	//stEvaluationResult EvaluateWithCrop(double cropMargin);

	stEvaluationResult* GetEvaluationResult(void) { return &m_stResult; }
	void PrintResultToConsole();
	void PrintResultToFile(void);
	void PrintResultToFile(const char *strFilepathAndName);
	void PrintResultMatrix(const char *strFilepathAndName);
	std::string PrintResultToString();

private:
	bool bInit;
	int m_nNumObj;
	int m_nNumTime;
	int m_nSavedResult;

	cv::Mat matXgt;
	cv::Mat matYgt;
	cv::Mat matX;
	cv::Mat matY;

	hj::Rect m_rectCropZone;
	hj::Rect m_rectCropZoneMargin;
	hj::Rect m_rectInnerCropZone;

	std::deque<unsigned int> m_queueID;
	std::vector<pointInfoSet> m_queueSavedResult;

	stParamEvaluator m_stParams;
	stEvaluationResult m_stResult;
};

}

//()()
//('')HAANJU.YOO


