/******************************************************************************
* Title        : CGUIManager
* Author       : Haanju Yoo
* Initial Date : 2016.12.16
* Version Num. : 0.9
* Description  : managing graphical interface
******************************************************************************/

#pragma once
#include "types.hpp"
#include <opencv2\highgui\highgui.hpp>

#define NOMINMAX
#include <windows.h>

namespace hj
{

class CGUIManager
{
public:
	CGUIManager();
	~CGUIManager();

	bool Initialize(double _xmin, double _xmax, double _ymin, double _ymax);
	bool Finalize();

	void SetResult(hj::CTrackLidarResult _result);
	void DrawResult();
	void HandleEvent(int evt, int x, int y, int flags);
	hj::CTrackLidarResult GetUserDefinedResult(hj::CTrackLidarResult _result);
	std::vector<std::pair<int, int>> GetAssignedIDs() { return vecPairAssignedID2TrackID_; }

private:
	void DrawGrids(cv::Mat _targetImage);
	void DrawIDButtons(cv::Mat _targetImage);
	void HandleButtonInput(const char _pressedKey);

private:
	// scaling info
	bool   bInit_;
	double dXmin_, dXmax_, dYmin_, dYmax_;
	double dImageWidth_, dImageHeight_;
	double dInputRescale_;
	double dInputRestore_;

	// vis. info
	hj::Point2D ptOrigin_;
	hj::Point2D ptX1000_, ptXN1000_, ptY1000_, ptYN1000_;
	std::vector<std::pair<hj::Point2D, hj::Point2D>> vecXGrids_;
	std::vector<std::pair<hj::Point2D, hj::Point2D>> vecYGrids_;

	// draw
	bool bFirstDraw_;
	std::vector<cv::Scalar> vecColors_;

	// GUI related
	int nNumIDs_;
	int nSelectedIDIdx_;
	char pressedKey_;
	std::vector<bool> vecIDSelected_;
	std::vector<int>  vecID_;
	std::vector<hj::Rect> vecIDRects_;
	std::string strWindowName_;
	SRWLOCK lockEventHandle_;

	// tracking info
	hj::CTrackLidarResult recentResult_;
	std::vector<std::pair<int, int>> vecPairAssignedID2TrackID_; // Assigned ID, real ID
	std::vector<hj::CObjectLidarInfo> vecAssignedIDResultBackup_;
	SRWLOCK lockTrackingResult_;
};

}


//()()
//('')HAANJU.YOO

