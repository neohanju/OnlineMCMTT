#include "opencv2\highgui\highgui.hpp"
#include "haanju_string.hpp"
#include "haanju_fileIO.hpp"
#include "haanju_misc.hpp"
#include "FrameGrabber.h"
#include <iostream>

//using namespace mvIMPACT::acquire;
const unsigned int MV_TIMEOUT_MS = 500;

////-----------------------------------------------------------------------------
#ifdef HJ_MATRIX_VISION

//-----------------------------------------------------------------------------
// Start the acquisition manually if this was requested(this is to prepare the driver for data capture and tell the device to start streaming data)
inline void manuallyStartAcquisitionIfNeeded(mvIMPACT::acquire::Device* pDev, const mvIMPACT::acquire::FunctionInterface& fi)
//-----------------------------------------------------------------------------
{
	if (pDev->acquisitionStartStopBehaviour.read() == mvIMPACT::acquire::assbUser)
	{
		const mvIMPACT::acquire::TDMR_ERROR result = static_cast<mvIMPACT::acquire::TDMR_ERROR>(fi.acquisitionStart());
		if (result != mvIMPACT::acquire::DMR_NO_ERROR)
		{
			std::cout << "'FunctionInterface.acquisitionStart' returned with an unexpected result: " 
				      << result << "(" << mvIMPACT::acquire::ImpactAcquireException::getErrorCodeAsString(result) 
				      << ")" << std::endl;
		}
	}
}

//-----------------------------------------------------------------------------
// Stop the acquisition manually if this was requested
inline void manuallyStopAcquisitionIfNeeded(mvIMPACT::acquire::Device* pDev, const mvIMPACT::acquire::FunctionInterface& fi)
//-----------------------------------------------------------------------------
{
	if (pDev->acquisitionStartStopBehaviour.read() == mvIMPACT::acquire::assbUser)
	{
		const mvIMPACT::acquire::TDMR_ERROR result = static_cast<mvIMPACT::acquire::TDMR_ERROR>(fi.acquisitionStop());
		if (result != mvIMPACT::acquire::DMR_NO_ERROR)
		{
			std::cout << "'FunctionInterface.acquisitionStop' returned with an unexpected result: " 
				      << result << "(" << mvIMPACT::acquire::ImpactAcquireException::getErrorCodeAsString(result)
				      << ")" << std::endl;
		}
	}
}

#endif
////-----------------------------------------------------------------------------

namespace hj
{

CFrameGrabber::CFrameGrabber()
	: bInit_(false)
	, bRealtimeOperation_(false)
	, bTerminate_(false)
	, nCamID_(0)
	, nCurrentFrameIndex_(0)
#ifdef HJ_MATRIX_VISION
	, pMVDev_(NULL)
	, pMVFunctionInterface_(NULL)
	, pMVRequest_(NULL)
	, nMVOpenCVDataType_(CV_8UC1)
#endif
{
}


CFrameGrabber::~CFrameGrabber()
{
	Finalize();
}


bool CFrameGrabber::Initialize(int _nID, stParamFrameGrabber &_stParams)
{
	SetParameters(_stParams);

	if (bInit_)
	{
		Finalize();
	}
	nCamID_ = _nID;

	if (hj::GRABBING == stParams_.nInputSource)
	{
		nInputType_ = HJ_CAM_MV_GIGE;
#ifndef HJ_MATRIX_VISION
		printf("[ERROR] This system does not have the MATRIX VISION's grabbing SDK\n");
		return false;
#else
		// MATRIX VISION's Giga Ethernet Camara
		/* check device availabilty */
		int numDevice = cMVDeviceManager_.deviceCount();
		if (nCamID_ > numDevice)
		{
			printf("[ERROR] Not enough number of devices: %d\n", nCamID_);
			return false;
		}
		pMVDev_ = cMVDeviceManager_[nCamID_];
		if (!pMVDev_)
		{
			printf("[ERROR] Unalbe to find the device: %d\n", nCamID_);
			return false;
		}

		// get a pointer and open it
		try
		{
			pMVDev_->open();
			assert(NULL != pMVDev_);
			nDeviceID_ = nCamID_;
			pMVFunctionInterface_ = new mvIMPACT::acquire::FunctionInterface(pMVDev_);
		}
		catch (mvIMPACT::acquire::ImpactAcquireException &e)
		{
			// this e.g. might happen if the same device is already opened in another process...
			printf("[ERROR] An error occurred while opening the device(error code: %d)\n",
				e.getErrorCode());
			return false;
		}

		TDMR_ERROR result = DMR_NO_ERROR;
		while (DMR_NO_ERROR == (result = static_cast<TDMR_ERROR>(pMVFunctionInterface_->imageRequestSingle()))) 
		{};

		if (result != DEV_NO_FREE_REQUEST_AVAILABLE)
		{
			std::cout << "'FunctionInterface.imageRequestSingle' returned with an unexpected result: "
				<< result << "(" << mvIMPACT::acquire::ImpactAcquireException::getErrorCodeAsString(result)
				<< ")" << std::endl;
		}
		manuallyStartAcquisitionIfNeeded(pMVDev_, *pMVFunctionInterface_);
		pMVRequest_ = 0;
		pMVPreviousRequest_ = 0;
#endif
	}
	else
	{
		// for the approperiate operation, 'nCurrentFrameIndex' has to be subtraced by one
		nCurrentFrameIndex_ = stParams_.nStartFrameIndex - 1;

		/* set image foler path */
		strImageFoler_ = stParams_.strInputDir;

		if (hj::PILSNU == stParams_.nInputSource || hj::PETS == stParams_.nInputSource)
		{
			nInputType_ = HJ_READ_IMAGE;
		}
		else
		{
			nInputType_ = HJ_READ_IMAGE_TIME_FORMAT;
			filePos_    = -1;
			if (!hj::GetFileList(strImageFoler_, "frame_*.jpg", vecStrFileNames_))
			{
				hj::printf_debug("[ERROR] There is no matched file\n");
				return false;
			}
			std::sort(vecStrFileNames_.begin(), vecStrFileNames_.end());
			prevTime_ = hj::GetTimeFromFileName(vecStrFileNames_[0], "frame_", ".jpg");
		}		
	}

	bRealtimeOperation_ = stParams_.bDoRealtimeOperation;
	bInit_ = true;
	bTerminate_ = false;
	
	return true;
}


bool CFrameGrabber::Finalize()
{
	if (!bInit_)
	{
		return true;
	}

#ifdef HJ_MATRIX_VISION
	if (NULL != pMVDev_)
	{		
		manuallyStopAcquisitionIfNeeded(pMVDev_, *pMVFunctionInterface_);

		if (pMVPreviousRequest_)
		{
			pMVPreviousRequest_->unlock();
			pMVPreviousRequest_ = 0;
		}

		if (NULL != pMVFunctionInterface_)
		{
			pMVFunctionInterface_->imageRequestReset(0, 0);

			// extract and unlock all requests that are now returned as 'aborted'
			int requestNr = INVALID_ID;
			while ((requestNr = pMVFunctionInterface_->imageRequestWaitFor(0)) >= 0)
			{
				pMVRequest_ = pMVFunctionInterface_->getRequest(requestNr);
				std::cout << "Request " << requestNr << " did return with status "
					<< pMVRequest_->requestResult.readS() << std::endl;
				pMVRequest_->unlock();
			}
			pMVRequest_ = 0;

			delete pMVFunctionInterface_;
			pMVFunctionInterface_ = NULL;
		}

		pMVDev_->close();
		pMVDev_ = NULL;
	}
#endif

	if (!matFrameBuffer_.empty()) { matFrameBuffer_.release(); }

	return true;
}

bool CFrameGrabber::Reset()
{
	if (bInit_)
	{
		Finalize();
	}
	Initialize(nCamID_, stParams_);

	return true;
}


bool CFrameGrabber::SetParameters(stParamFrameGrabber &_stParams)
{
	stParams_ = _stParams;
	return true;
}


HJ_GRAB_RESULT CFrameGrabber::GrabFrame()
{
	assert(bInit_);
	
	bool bGrabbed = false;
	HJ_GRAB_RESULT resultCode = HJ_GR_NORMAL;

	if (!matFrame_.empty())
	{
		matFramePrev_ = matFrame_.clone();
		matFrame_.release(); 
	}

	nCurrentFrameIndex_++;
	if (HJ_READ_IMAGE == nInputType_)
	{
		bGrabbed = GrabImage(nCurrentFrameIndex_);
		if (!bGrabbed) { resultCode = HJ_GR_DATASET_ENDED; }
	}
	else if (HJ_READ_IMAGE_TIME_FORMAT == nInputType_)
	{
		long timeGap = 0;
		time_t timeCurrGrab = clock();
		if (nCurrentFrameIndex_ > stParams_.nStartFrameIndex)
		{
			timeGap = (long)(timeCurrGrab - timePrevGrab_);			
		}
		timePrevGrab_ = timeCurrGrab;

		//bGrabbed = GrabImage(timeGap);
		bGrabbed = GrabImage((long)10);
		if (!bGrabbed) { resultCode = HJ_GR_UNKNOWN_ERROR; }
	}
	else if (HJ_READ_VIDEO == nInputType_)
	{
		bGrabbed = GrabVideoFrame();
		if (!bGrabbed) { resultCode = HJ_GR_UNKNOWN_ERROR; }
	}
	else if (HJ_CAM_MV_GIGE == nInputType_)
	{		
		bGrabbed = GrabMVGECamera();
		if (!bGrabbed) { resultCode = HJ_GR_HARDWARE_FAIL; }
	}

	if (matFrame_.empty() || !bGrabbed)
	{
		if (!matFramePrev_.empty())
		{
			matFrame_ = matFramePrev_.clone();
		}
	}

	return resultCode;
}


cv::Mat CFrameGrabber::GetFrame()
{
	assert(bInit_);
	return matFrame_;
}


bool CFrameGrabber::GrabImage(int _nFrameIndex)
{
	if (stParams_.bExistFrameRange && _nFrameIndex > stParams_.nEndFrameIndex)
	{
		return false;
	}

	char strImagePath[300];
	sprintf_s(strImagePath, sizeof(strImagePath),
		"%s/frame_%04d.jpg", strImageFoler_.c_str(), _nFrameIndex);
	matFrame_ = cv::imread(strImagePath, cv::IMREAD_COLOR);
	if (matFrame_.empty())
	{
		return false;
	}
	return true;
}


bool CFrameGrabber::GrabImage(long _nTimeGap)
{
	// convert _nTimeGap to formatted time gap
	long timeTemp = _nTimeGap;
	long miliseconds = 0, seconds = 0, minutes = 0, hours = 0;
	long fileTime = 0;

	if (0 <= filePos_ && 0 < _nTimeGap && bRealtimeOperation_)
	{
		timeTemp = _nTimeGap;
		// ms
		miliseconds = timeTemp % 1000;
		timeTemp = timeTemp / 1000;
		// sec
		seconds = timeTemp % 60;
		timeTemp = timeTemp / 60;
		// min
		minutes = timeTemp % 60;
		timeTemp = timeTemp / 60;
		// hours
		hours = timeTemp % 24;

		long formattedTimeGap = miliseconds;
		formattedTimeGap += seconds * 1000;
		formattedTimeGap += minutes * 100000;
		formattedTimeGap += hours * 10000000;

		// find nearest (right side) file
		long targetTime = prevTime_ + formattedTimeGap;
		bool bFound = false;
		for (filePos_++; filePos_ < (int)vecStrFileNames_.size(); filePos_++)
		{
			fileTime = hj::GetTimeFromFileName(vecStrFileNames_[filePos_], "frame_", ".jpg");
			if (fileTime > targetTime)
			{
				bFound = true;
				break;
			}
		}
		if (!bFound) { return false; }
	}
	else
	{
		filePos_++;
		if (filePos_ >= (int)vecStrFileNames_.size())
		{
			return false;
		}
		fileTime = hj::GetTimeFromFileName(vecStrFileNames_[filePos_], "frame_", ".jpg");
	}
	
	char strImagePath[300];
	sprintf_s(strImagePath, sizeof(strImagePath),
		"%s/%s", strImageFoler_.c_str(), vecStrFileNames_[filePos_].c_str());
	matFrame_ = cv::imread(strImagePath, cv::IMREAD_COLOR);
	if (matFrame_.empty())
	{
		return false;
	}
	prevTime_ = fileTime;
	return true;
}


bool CFrameGrabber::GrabVideoFrame()
{
	return true;
}


bool CFrameGrabber::GrabMVGECamera()
{	
#ifdef HJ_MATRIX_VISION

	const unsigned int timeout_ms = 500;
	int requestNr = INVALID_ID;

	// for real-time grabbing
	do
	{
		// get a queued result from the default capture queue without waiting(just check if there are pending results)
		requestNr = pMVFunctionInterface_->imageRequestWaitFor(0);
		if (pMVFunctionInterface_->isRequestNrValid(requestNr))
		{
			pMVRequest_ = pMVFunctionInterface_->getRequest(requestNr);
			if (pMVRequest_)
			{
				// discard an outdated buffer
				pMVRequest_->unlock();
				
				// request next
				pMVFunctionInterface_->imageRequestSingle();
			}
			pMVRequest_ = 0;
		}

	} while (pMVFunctionInterface_->isRequestNrValid(requestNr));

	if (!pMVRequest_)
	{
		// there was no image in the queue we have to wait for results from the default capture queue
		requestNr   = pMVFunctionInterface_->imageRequestWaitFor(timeout_ms);
		pMVRequest_ = pMVFunctionInterface_->isRequestNrValid(requestNr) 
			? pMVFunctionInterface_->getRequest(requestNr) : 0;
	}

	if (pMVRequest_)
	{
		if (pMVRequest_->isOK())
		{
			SetMVGECamImage(pMVRequest_);
		}
		else
		{
			std::cout << "Error: " << pMVRequest_->requestResult.readS() << std::endl;
		}

		if (pMVPreviousRequest_)
		{
			pMVPreviousRequest_->unlock();
		}
		pMVPreviousRequest_ = pMVRequest_;
		pMVRequest_         = 0;

		// send a new image request into the capture queue
		pMVFunctionInterface_->imageRequestSingle();
	}
	else
	{
		// If the error code is -2119(DEV_WAIT_FOR_REQUEST_FAILED), the documentation will provide
		// additional information under TDMR_ERROR in the interface reference
		std::cout << "imageRequestWaitFor failed ("
				    << requestNr << ", " << ImpactAcquireException::getErrorCodeAsString(requestNr) << ")"
				    << ", timeout value too small?" << std::endl;
	}	
#endif

	return true;
}

#ifdef HJ_MATRIX_VISION
bool CFrameGrabber::SetMVGECamImage(mvIMPACT::acquire::Request* _pRequest)
{
	bool bFormatFound = true;
	if (matFrameBuffer_.empty())
	{
		// allocate buffer		
		switch (_pRequest->imagePixelFormat.read())
		{
		case ibpfMono8:
			nMVOpenCVDataType_ = CV_8UC1;
			break;
		case ibpfMono10:
		case ibpfMono12:
		case ibpfMono14:
		case ibpfMono16:
			nMVOpenCVDataType_ = CV_16UC1;
			break;
		case ibpfMono32:
			nMVOpenCVDataType_ = CV_32SC1;
			break;
		case ibpfBGR888Packed:
		case ibpfRGB888Packed:
			nMVOpenCVDataType_ = CV_8UC3;
			break;
		case ibpfRGBx888Packed:
			nMVOpenCVDataType_ = CV_8UC4;
			break;
		case ibpfRGB101010Packed:
		case ibpfRGB121212Packed:
		case ibpfRGB141414Packed:
		case ibpfRGB161616Packed:
			nMVOpenCVDataType_ = CV_16UC3;
			break;
		case ibpfMono12Packed_V1:
		case ibpfMono12Packed_V2:
		case ibpfBGR101010Packed_V2:
		case ibpfRGB888Planar:
		case ibpfRGBx888Planar:
		case ibpfYUV422Packed:
		case ibpfYUV422_10Packed:
		case ibpfYUV422_UYVYPacked:
		case ibpfYUV422_UYVY_10Packed:
		case ibpfYUV422Planar:
		case ibpfYUV444Packed:
		case ibpfYUV444_10Packed:
		case ibpfYUV444_UYVPacked:
		case ibpfYUV444_UYV_10Packed:
		case ibpfYUV444Planar:
		case ibpfYUV411_UYYVYY_Packed:
		default:
			printf("[ERROR] Don't know how to render this pixel format (%s)\
									in OpenCV! Select another one e.g. by writing to \
									mvIMPACT::acquire::ImageDestination::pixelFormat!\n",
				_pRequest->imagePixelFormat.readS().c_str());
			bFormatFound = false;
			break;
		}

		if (!bFormatFound) { return false; }
	}

	matFrameBuffer_ = cv::Mat(
		cv::Size(_pRequest->imageWidth.read(), _pRequest->imageHeight.read()),
		nMVOpenCVDataType_,
		_pRequest->imageData.read(),
		_pRequest->imageLinePitch.read());

	cv::cvtColor(matFrameBuffer_, matFrame_, CV_RGBA2RGB);

	return true;
}
#endif


}



//()()
//('')HAANJU.YOO


