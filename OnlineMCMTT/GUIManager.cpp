#include "GUIManager.h"
#include "haanju_visualize.hpp"
#include "haanju_misc.hpp"

// display
#define DISPLAY_IMG_WIDTH     (800)
#define DISPLAY_IMG_HEIGHT    (800)
#define DISPLAY_NUMBER_HEIGHT (30)

// matching
#define CLICK_DISTANCE (300)

const cv::Scalar gridColor(102, 62, 0);
const cv::Scalar bgColor(51, 31, 0);

void onMouse(int evt, int x, int y, int flags, void* param)
{
	hj::CGUIManager *p = (hj::CGUIManager *)param;
	p->HandleEvent(evt, x, y, flags);
}


namespace hj
{

CGUIManager::CGUIManager()
	: bInit_(false)
	, bFirstDraw_(true)
	, dInputRescale_(1.0)
	, dInputRestore_(1.0)
	, nNumIDs_(10)
	, nSelectedIDIdx_(-1)
	, strWindowName_("Tracking Pannel")
{
}


CGUIManager::~CGUIManager()
{
	Finalize();
}


bool CGUIManager::Initialize(double _xmin, double _xmax, double _ymin, double _ymax)
{
	if (bInit_) { return false; }

	vecColors_ = GenerateColors(300);

	// set range
	dXmin_ = _xmin;
	dXmax_ = _xmax;
	dYmin_ = _ymin;
	dYmax_ = _ymax;
	dImageWidth_  = dXmax_ - dXmin_ + 1;
	dImageHeight_ = dYmax_ - dYmin_ + 1;

	// expend to be a squaure
	if (dImageWidth_ < dImageHeight_)
	{
		double xCenter = dImageWidth_ * 0.5 + dXmin_;
		dXmin_ = xCenter - dImageHeight_ * 0.5;
		dXmax_ = xCenter + dImageHeight_ * 0.5;
		dImageWidth_ = dImageHeight_;
	}
	else
	{
		double yCenter = dImageHeight_ * 0.5 + dYmin_;
		dYmin_ = yCenter - dImageWidth_ * 0.5;
		dYmax_ = yCenter + dImageWidth_ * 0.5;
		dImageHeight_ = dImageWidth_;
	}	

	// calculate rescaling factors
	dInputRescale_ = (double)DISPLAY_IMG_WIDTH / dImageWidth_;
	dInputRestore_ = 1.0 / dInputRescale_;

	//cv::Mat matDisplay =
	//	cv::Mat::zeros((int)(DISPLAY_IMG_HEIGHT + DISPLAY_NUMBER_HEIGHT),
	//	(int)DISPLAY_IMG_WIDTH, CV_8UC3);

	ptOrigin_ = hj::Point2D(-dXmin_ * dInputRescale_, -dYmin_ * dInputRescale_);
	ptX1000_  = hj::Point2D((-dXmin_ + 1000) * dInputRescale_, -dYmin_ * dInputRescale_);
	ptXN1000_ = hj::Point2D((-dXmin_ - 1000) * dInputRescale_, -dYmin_ * dInputRescale_);
	ptY1000_  = hj::Point2D(-dXmin_ * dInputRescale_, (-1000 - dYmin_) * dInputRescale_);
	ptYN1000_ = hj::Point2D(-dXmin_ * dInputRescale_, (+1000 - dYmin_) * dInputRescale_);

	//// draw cooridnates
	//cv::line(matDisplay, ptXN1000_.cv(), ptX1000_.cv(), cv::Scalar(0, 255, 0));
	//cv::line(matDisplay, ptYN1000_.cv(), ptY1000_.cv(), cv::Scalar(255, 0, 0));

	// set grid
	int gridSize = 1000;
	std::vector<double> vecGridXPts;
	for (int x = (int)dXmin_; x < (int)dXmax_; x++)
	{
		if (0 == x % gridSize)
		{
			double xCoord = ((double)x - dXmin_) * dInputRescale_;
			std::pair<hj::Point2D, hj::Point2D> gridPair(
				hj::Point2D(xCoord, 0), hj::Point2D(xCoord, DISPLAY_IMG_HEIGHT-1));
			vecXGrids_.push_back(gridPair);
		}
	}
	std::vector<double> vecGridYPts;
	for (int y = (int)dYmin_; y < (int)dYmax_; y++)
	{
		if (0 == y % gridSize)
		{
			double yCoord = ((double)y - dYmin_) * dInputRescale_;
			std::pair<hj::Point2D, hj::Point2D> gridPair(
				hj::Point2D(0, yCoord), hj::Point2D(DISPLAY_IMG_WIDTH-1, yCoord));
			vecYGrids_.push_back(gridPair);
		}
	}
	//DrawGrids(matDisplay);


	// ID button related
	nNumIDs_ = 10;	
	vecIDSelected_.resize(nNumIDs_, false);
	vecIDRects_.resize(nNumIDs_, hj::Rect(0.0, 0.0, 0.0, 0.0));
	vecID_.resize(nNumIDs_, 1);

	double rectInterval = (double)DISPLAY_IMG_WIDTH / nNumIDs_;

	for (int i = 0; i < nNumIDs_; i++)
	{
		vecIDRects_[i].x = i * rectInterval;
		vecIDRects_[i].y = DISPLAY_IMG_HEIGHT + 1.0;
		vecIDRects_[i].w = rectInterval;
		vecIDRects_[i].h = DISPLAY_NUMBER_HEIGHT;

		vecID_[i] = i+1;
	}
	//DrawIDButtons(matDisplay);

	// ID mapping related
	vecPairAssignedID2TrackID_.resize(nNumIDs_);
	vecAssignedIDResultBackup_.resize(nNumIDs_);
	for (int i = 0; i < nNumIDs_; i++)
	{
		vecPairAssignedID2TrackID_[i] = std::make_pair(vecID_[i], -1);
		vecAssignedIDResultBackup_[i].id = vecID_[i];
	}

	//cv::imshow("Tracking Pannel", matDisplay);
	//cv::moveWindow("Tracking Pannel", 1100, 0);
	//cv::waitKey(1);

	bFirstDraw_ = true;

	InitializeSRWLock(&lockEventHandle_);
	InitializeSRWLock(&lockTrackingResult_);

	bInit_ = true;
	return true;
}


bool CGUIManager::Finalize()
{
	if (!bInit_) { return false; }

	cv::destroyWindow(strWindowName_.c_str());

	bInit_ = false;
	return true;
}


void CGUIManager::SetResult(hj::CTrackLidarResult _result)
{ 
	//------------------------------------------
	AcquireSRWLockExclusive(&lockTrackingResult_);
	//------------------------------------------
	recentResult_ = _result; 
	//------------------------------------------
	ReleaseSRWLockExclusive(&lockTrackingResult_);
	//------------------------------------------
}


void CGUIManager::DrawResult()
{	
	cv::Mat matDisplay =
		cv::Mat::zeros((int)(DISPLAY_IMG_HEIGHT + DISPLAY_NUMBER_HEIGHT), 
		(int)DISPLAY_IMG_WIDTH, CV_8UC3);
	matDisplay = bgColor;

	DrawGrids(matDisplay);
	DrawIDButtons(matDisplay);

	// draw cooridnates
	cv::line(matDisplay, ptXN1000_.cv(), ptX1000_.cv(), cv::Scalar(0, 255, 0));
	cv::line(matDisplay, ptYN1000_.cv(), ptY1000_.cv(), cv::Scalar(0, 255, 0));

	// rescale
	for (int objIdx = 0; objIdx < recentResult_.objectLidarInfos.size();
		objIdx++)
	{
		hj::Point2D scannedPoint;
		scannedPoint.x = recentResult_.objectLidarInfos[objIdx].location.x;
		scannedPoint.y = -recentResult_.objectLidarInfos[objIdx].location.y;

		//if (!cropZone.contain(scannedPoint)) { continue; }
		scannedPoint.x -= dXmin_;
		scannedPoint.y -= dYmin_;

		hj::Point2D rescaledPoint = scannedPoint * dInputRescale_;

		bool bSelectedOne = false;
		int  displayID    = recentResult_.objectLidarInfos[objIdx].id;
		for (int idx = 0; idx < vecPairAssignedID2TrackID_.size(); idx++)
		{
			if (vecPairAssignedID2TrackID_[idx].second == displayID)
			{
				displayID = vecPairAssignedID2TrackID_[idx].first;
				bSelectedOne = true;
				break;
			}
		}
		cv::Scalar curColor = bSelectedOne? vecColors_[displayID] : cv::Scalar(100.0, 100.0, 100.0);

		cv::circle(matDisplay, rescaledPoint.cv(), (int)(CLICK_DISTANCE * dInputRescale_), curColor);
		cv::putText(matDisplay, std::to_string(displayID),
			cv::Point((int)rescaledPoint.x - 1, (int)rescaledPoint.y - 2),
			cv::FONT_HERSHEY_SIMPLEX, 0.3, curColor);
		
	}
	cv::imshow(strWindowName_.c_str(), matDisplay);
	if (bFirstDraw_)
	{
		cv::moveWindow(strWindowName_.c_str(), 800, 0);
		
		// mouse callback
		cv::setMouseCallback(strWindowName_.c_str(), onMouse, this);
		bFirstDraw_ = false;
	}
	const char pressedKey = cv::waitKey(30);
	HandleButtonInput(pressedKey);
}

void CGUIManager::HandleEvent(int _event, int _x, int _y, int _flags)
{
	hj::Point2D clickedPoint((double)_x, (double)_y);

	//------------------------------------------
	AcquireSRWLockShared(&lockTrackingResult_);
	//------------------------------------------
	hj::CTrackLidarResult curResult = recentResult_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockTrackingResult_);
	//------------------------------------------

	//hj::printf_debug("%f, %f\n", clickedPoint.x, clickedPoint.y);
	

	////------------------------------------------
	//AcquireSRWLockExclusive(&lockEventHandle_);
	////------------------------------------------
	if (_event == CV_EVENT_LBUTTONDOWN)
	{
	}
	else if (_event == CV_EVENT_LBUTTONUP)
	{
		// check wheather ID clicked or not
		if (_y > DISPLAY_IMG_HEIGHT)
		{
			int nSelectedIDIdx = -1;
			for (int i = 0; i < nNumIDs_; i++)
			{
				if (vecIDRects_[i].contain(clickedPoint))
				{
					vecIDSelected_[i] = true;
					nSelectedIDIdx = i;
				}
				else
				{
					vecIDSelected_[i] = false;
				}
			}

			if (0 <= nSelectedIDIdx)
			{
				nSelectedIDIdx_ = nSelectedIDIdx;
			}
			else
			{
				nSelectedIDIdx_ = -1;
			}
		}
		else if (0 <= nSelectedIDIdx_)
		{
			double minDistance = CLICK_DISTANCE;
			int    matchedID = -1;

			hj::Point3D clickedLocation(0.0, 0.0, 0.0);
			clickedLocation.x = _x * dInputRestore_ + dXmin_;
			clickedLocation.y = -(_y * dInputRestore_ + dYmin_);

			for (int i = 0; i < curResult.objectLidarInfos.size(); i++)
			{
				double curDistance =
					(curResult.objectLidarInfos[i].location - clickedLocation).norm_L2();

				if (minDistance > curDistance)
				{
					minDistance = curDistance;
					matchedID = curResult.objectLidarInfos[i].id;
				}
			}

			if (minDistance < CLICK_DISTANCE)
			{
				vecPairAssignedID2TrackID_[nSelectedIDIdx_].second = matchedID;
			}
		}
		//else
		//{
		//	// TEMPORAL (manually select the coordinates)
		//	for (int i = 0; i < nNumIDs_; i++)
		//	{
		//		if (!vecIDSelected_[i])
		//		{
		//			continue;
		//		}
		//		vecAssignedIDResultBackup_[i].location.x =   _x * dInputRestore_ + dXmin_;
		//		vecAssignedIDResultBackup_[i].location.y = -(_y * dInputRestore_ + dYmin_);

		//		break;
		//	}
		//}
	}
	else if (_event == CV_EVENT_MOUSEMOVE)
	{
	}
	else if (_event == CV_EVENT_RBUTTONDOWN)
	{
	}
	else if (_event == CV_EVENT_RBUTTONUP)
	{
	}
	else if (_event == CV_EVENT_MBUTTONDOWN)
	{
	}
	else if (_event == CV_EVENT_MBUTTONUP)
	{
	}
	else if (_event == CV_EVENT_LBUTTONDBLCLK)
	{
	}
	else if (_event == CV_EVENT_RBUTTONDBLCLK)
	{
	}


	//if (_flags & CV_EVENT_FLAG_LBUTTON)
	//{
	//	cout << "\tCV_EVENT_FLAG_LBUTTON" << endl;
	//}
	//if (_flags & CV_EVENT_FLAG_RBUTTON)
	//{
	//	cout << "\tCV_EVENT_FLAG_RBUTTON" << endl;
	//}
	//if (_flags & CV_EVENT_FLAG_MBUTTON)
	//{
	//	cout << "\tCV_EVENT_FLAG_MBUTTON" << endl;
	//}
	//if (_flags & CV_EVENT_FLAG_CTRLKEY)
	//{
	//	cout << "\tCV_EVENT_FLAG_CTRLKEY" << endl;
	//}
	//if (_flags & CV_EVENT_FLAG_SHIFTKEY)
	//{
	//	cout << "\tCV_EVENT_FLAG_SHIFTKEY" << endl;
	//}
	//if (_flags & CV_EVENT_FLAG_ALTKEY)
	//{
	//	cout << "\tCV_EVENT_FLAG_ALTKEY" << endl;
	//}

	////------------------------------------------
	//ReleaseSRWLockShared(&lockEventHandle_);
	////------------------------------------------
}


hj::CTrackLidarResult CGUIManager::GetUserDefinedResult(hj::CTrackLidarResult _result)
{
	hj::CTrackLidarResult newResult;
	newResult.frameIdx  = _result.frameIdx;
	newResult.timeStamp = _result.timeStamp;

	for (int i = 0; i < vecPairAssignedID2TrackID_.size(); i++)
	{
		if (0 > vecPairAssignedID2TrackID_[i].second)
		{
			// not assigned
			continue;
		}

		for (int findIdx = 0; findIdx < _result.objectLidarInfos.size(); findIdx++)
		{
			if (_result.objectLidarInfos[findIdx].id !=
				vecPairAssignedID2TrackID_[i].second)
			{
				continue;
			}
			//newResult.objectLidarInfos.push_back(_result.objectLidarInfos[findIdx]);
			//newResult.objectLidarInfos.back().id = vecPairAssignedID2TrackID_[i].first;

			// update location
			vecAssignedIDResultBackup_[i].location = _result.objectLidarInfos[findIdx].location;
			vecAssignedIDResultBackup_[i].frameIdx = _result.frameIdx;
		}
		newResult.objectLidarInfos.push_back(vecAssignedIDResultBackup_[i]);
	}

	return newResult;
}


void CGUIManager::DrawGrids(cv::Mat _targetImage)
{
	for (int i = 0; i < vecXGrids_.size(); i++)
	{
		cv::line(_targetImage, 
			vecXGrids_[i].first.cv(), vecXGrids_[i].second.cv(), 
			gridColor);
	}
	for (int i = 0; i < vecYGrids_.size(); i++)
	{
		cv::line(_targetImage, 
			vecYGrids_[i].first.cv(), vecYGrids_[i].second.cv(), 
			gridColor);
	}
}


void CGUIManager::DrawIDButtons(cv::Mat _targetImage)
{
	for (int i = 0; i < nNumIDs_; i++)
	{
		if (vecIDSelected_[i])
		{
			cv::rectangle(_targetImage, vecIDRects_[i].cv(), gridColor, CV_FILLED);
			cv::putText(_targetImage, std::to_string(vecID_[i]),
				cv::Point((int)vecIDRects_[i].center().x - 10, (int)vecIDRects_[i].center().y + 8),
				0, 0.8, cv::Scalar(255, 255, 255));
		}
		else
		{
			cv::rectangle(_targetImage, vecIDRects_[i].cv(), gridColor);
			cv::putText(_targetImage, std::to_string(vecID_[i]),
				cv::Point((int)vecIDRects_[i].center().x - 10, (int)vecIDRects_[i].center().y + 8),
				0, 0.8, gridColor);
		}
		
	}
}


void CGUIManager::HandleButtonInput(const char _pressedKey)
{
	////------------------------------------------
	//AcquireSRWLockExclusive(&lockEventHandle_);
	////------------------------------------------

	if (-1 == _pressedKey)
	{
		// do nothing
	}
	else if (32 == _pressedKey)
	{
		// hot key for ID switch between 1 and 2
		hj::CObjectLidarInfo tempObject = vecAssignedIDResultBackup_[0];
		int tempID = vecPairAssignedID2TrackID_[0].second;

		vecPairAssignedID2TrackID_[0].second = vecPairAssignedID2TrackID_[1].second;
		vecAssignedIDResultBackup_[0] = vecAssignedIDResultBackup_[1];
		vecAssignedIDResultBackup_[0].id = vecPairAssignedID2TrackID_[0].first;

		vecPairAssignedID2TrackID_[1].second = tempID;
		vecAssignedIDResultBackup_[1] = tempObject;
		vecAssignedIDResultBackup_[1].id = vecPairAssignedID2TrackID_[1].first;
	}
	else
	{
		int nSelectedIDIdx = -1;
		for (int i = 0; i < vecID_.size(); i++)
		{
			if (48 + vecID_[i] == _pressedKey)
			{
				vecIDSelected_[i] = true;
				nSelectedIDIdx = i;
			}
			else
			{
				vecIDSelected_[i] = false;
			}
		}
		if (-1 < nSelectedIDIdx)
		{
			nSelectedIDIdx_ = nSelectedIDIdx;
		}
	}
	////------------------------------------------
	//ReleaseSRWLockShared(&lockEventHandle_);
	////------------------------------------------
}


}


//()()
//('')HAANJU.YOO
