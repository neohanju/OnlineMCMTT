#include "DataManager.h"
#include "pugixml.hpp"
#include "haanju_string.hpp"
#include "haanju_misc.hpp"

/* INPUT SOURCES */
#define PARAM_TAG_GRABBING "GRABBING"
#define PARAM_TAG_PILSNU   "PILSNU"
#define PARAM_TAG_DKU      "DKU"

/* NODE NAMES */
#define PARAM_TAG_PARAMETERS     "Parameters"
#define PARAM_TAG_OPERATION      "Operation"
#define PARAM_TAG_SETTING        "Setting"
#define PARAM_TAG_CAMERA         "Camera"
#define PARAM_TAG_DATASET        "Dataset"
#define PARAM_TAG_CALIBRATION    "Calibration"
#define PARAM_TAG_DETECTION      "Detection"
#define PARAM_TAG_TRACKING_2D    "Tracking2D"
#define PARAM_TAG_ASSOCIATION_3D "Association3D"
#define PARAM_TAG_RESULT         "Result"
#define PARAM_TAG_EVALUATION     "Evaluation"

/* ATTRIBUTES NAMES */
// Operation
#define PARAM_TAG_INPUT_SOURCE              "input_source"
#define PARAM_TAG_REALTIME_OPERATION_TAG    "do_realtime_op"
#define PARAM_TAG_RECORD_FLAG               "do_recording"
#define PARAM_TAG_DETECT_FLAG               "do_detection"
#define PARAM_TAG_EVALUATION_FLAG           "do_evaluation"
#define PARAM_TAG_DETECTOR_PATH             "detector_path"
// Setting
#define PARAM_TAG_SETTING_NAME "name"
// Camera
#define PARAM_TAG_CAMERA_IDS "camera_ids"
#define PARAM_TAG_FRAME_RATE "frame_rate"
// Dataset
#define PARAM_TAG_PATH        "path"
#define PARAM_TAG_START_FRAME "start_frame"
#define PARAM_TAG_END_FRAME   "end_frame"
// Detection
#define PARAM_TAG_DETECTION_TYPE            "type"
#define PARAM_TAG_NMS_RATIO_1               "nsm_overlap_ratio_1"
#define PARAM_TAG_NMS_RATIO_2               "nsm_overlap_ratio_2"
#define PARAM_TAG_THRESHOLD                 "threshold"
// Tracking2D
#define PARAM_TAG_IMAGE_RESCALE             "image_rescale"
#define PARAM_TAG_DET_MIN_HEIGHT            "detection_min_height"
#define PARAM_TAG_DET_MAX_HEIGHT            "detection_max_height"
#define PARAM_TAG_MAX_TRACKLET_LENGTH       "max_tracklet_length"
#define PARAM_TAG_MIN_NUM_FEATURES          "min_num_features"
#define PARAM_TAG_MAX_NUM_FEATURES          "max_num_features"
#define PARAM_TAG_BACK_TRACKING_LENGTH      "back_tracking_length"
#define PARAM_TAG_FEATURE_TRACK_WIN_RATIO   "feature_track_window_size_ratio"
#define PARAM_TAG_MAX_BOX_DISTANCE          "max_box_distance"
#define PARAM_TAG_MIN_BOX_OVERLAP_RATIO     "min_box_overlap_ratio"
#define PARAM_TAG_MAX_BOX_CENTER_DIFF_RATIO "max_box_center_diff_ratio"
#define PARAM_TAG_MIN_OPT_FLOW_MAJ_RATIO    "min_optical_flow_majority_ratio"
#define PARAM_TAG_MAX_DETECTION_DISTANCE    "max_detection_distance"
#define PARAM_TAG_MAX_HEIGHT_DIFFERENCE     "max_height_difference"
#define PARAM_TAG_VISUALIZE                 "visualize"
// Association3D
#define PARAM_TAG_PROC_WINDOW_SIZE                   "proc_window_size"
#define PARAM_TAG_K_BEST_SIZE                        "K_best_size"
#define PARAM_TAG_MAX_NUM_TRACKS_IN_OPTIMIZATION     "max_num_tracks_in_optimization"
#define PARAM_TAG_MAX_NUM_TRACKS_IN_UNCONFIRMED_TREE "maX_num_tracks_in_unconfirmed_tree"
#define PARAM_TAG_MAX_NUM_TRACKS_IN_CONFIRMED_TREE   "max_num_tracks_in_confirmed_tree"
#define PARAM_TAG_NUM_FRAMES_FOR_CONFIRMATION        "num_frames_for_confirmation"
#define PARAM_TAG_DO_BRANCH_CUT                      "do_branch_cut"
#define PARAM_TAG_MIN_TRACKLET_LENGTH                "min_tracklet_length"
#define PARAM_TAG_MAX_TRACKLET_DISTANCE              "max_tracklet_distance"
#define PARAM_TAG_CONSIDER_SENSITIVITY               "consider_sensitivity"
#define PARAM_TAG_MAX_SENSITIVITY_ERROR              "max_sensitivity_error"
#define PARAM_TAG_MIN_TARGET_PROXIMITY               "min_target_proximity"
#define PARAM_TAG_DEFAULT_HEIGHT                     "default_height"
#define PARAM_TAG_CLUSTERING_DISTANCE                "clustering_distance"
#define PARAM_TAG_MIN_LINKING_PROBABILITY            "min_linking_probability"
#define PARAM_TAG_MAX_TIME_JUMP                      "max_time_jump"
#define PARAM_TAG_COST_RGB_MIN_DISTANCE              "cost_rgb_min_distance"
#define PARAM_TAG_COST_RGB_COEF                      "cost_rgb_coef"
#define PARAM_TAG_COST_RGB_DECAY_COEF                "cost_rgb_decay_coef"
#define PARAM_TAG_COST_TRACKLET_LINKING_MIN_DISTANCE "cost_tracklet_linking_min_distance"
#define PARAM_TAG_COST_TRACKLET_LINKING_COEF         "cost_tracklet_linking_coef"
#define PARAM_TAG_MIN_CONSTRUCT_PROBABILITY          "min_construct_probability"
#define PARAM_TAG_FP_RATE                            "fp_rate"
#define PARAM_TAG_FN_RATE                            "fn_rate"
#define PARAM_TAG_ENTER_PENALTY_FREE_LENGTH          "enter_penalty_free_length"
#define PARAM_TAG_BOUNDARY_DISTANCE                  "boundary_distance"
#define PARAM_TAG_ENTER_PROB_MAX                     "enter_prob_max"
#define PARAM_TAG_ENTER_PROB_DECAY_COEF              "enter_prob_decay_coef"
#define PARAM_TAG_EXIT_PROB_MAX                      "exit_prob_max"
#define PARAM_TAG_EXIT_PROB_DIST_DECAY_COEF          "exit_prob_dist_decay_coef"
#define PARAM_TAG_EXIT_PROB_LENGTH_DECAY_COEF        "exit_prob_length_decay_coef"
#define PARAM_TAG_COST_ENTER_MAX                     "cost_enter_max"
#define PARAM_TAG_COST_EXIT_MAX                      "cost_exit_max"
#define PARAM_TAG_MAX_OUTPOINT                       "max_outpoint"
#define PARAM_TAG_DETECTION_ERROR                    "detection_error"
#define PARAM_TAG_CALIBRATION_ERROR                  "calibration_error"
#define PARAM_TAG_KALMAN_PROC_NOISE_SIGMA            "kalman_proc_noise_sigma"
#define PARAM_TAG_KALMAN_MEASURE_NOISE_SIGMA         "kalman_measure_noise_sigma"
#define PARAM_TAG_KALMAN_POST_ERROR_COV              "kalmna_post_error_cov"
#define PARAM_TAG_KALMAN_CONFIDENCE_LEVEL            "kalman_confidence_level"
#define PARAM_TAG_VELOCITY_LEARNING_RATE             "velocity_learning_rate"
//#define PARAM_TAG_FRAME_RATE "frame_rate"
#define PARAM_TAG_MAX_MOVING_SPEED         "max_moving_speed"
#define PARAM_TAG_MIN_MOVING_SPEED         "min_moving_speed"
#define PARAM_TAG_IMAGE_PATCH_WIDTH        "image_patch_width"
#define PARAM_TAG_NUM_RGB_HISTOGRAM_BINS   "num_rgb_histogram_bins"
#define PARAM_TAG_SOLVER_TIME_LIMIT        "solver_time_limit"
#define PARAM_TAG_SOLVER_MAX_ITER          "solver_max_iter"
#define PARAM_TAG_RESULT_TRAJECTORY_LENGTH "result_trajectory_length"
#define PARAM_TAG_SHOW_TREE_ID             "show_tree_id"
// Evaluation
#define PARAM_TAG_CROP_XMIN "crop_xmin"
#define PARAM_TAG_CROP_XMAX "crop_xmax"
#define PARAM_TAG_CROP_YMIN "crop_ymin"
#define PARAM_TAG_CROP_YMAX "crop_ymax"


namespace hj
{

CDataManager::CDataManager()
	: bDataReady_(false)
	, bSystemRun_(false)
{
}


CDataManager::~CDataManager()
{
}


bool CDataManager::Initialize(std::string _strParameterFilePath)
{
	/* default values for each parameters */
	bRealtimeOperation_ = false;
	bDataReady_ = false;
	bRecord_ = false;
	bDetect_ = true;	
	strCalibrationPath_ = "";
	strDatasetPath_     = "";
	strDetectorPath_    = "";
	strResultPath_      = "";

	InitializeSRWLock(&lockParamAccess_);	
	InitializeSRWLock(&lockTrack3DResult_);
	InitializeSRWLock(&lockGUI_);

	std::fill(std::begin(arrFrameBufferFull_), std::end(arrFrameBufferFull_), false);	
	std::fill(std::begin(arrDetectionBufferFull_), std::end(arrDetectionBufferFull_), false);
	std::fill(std::begin(arrTrack2DBufferFull_), std::end(arrTrack2DBufferFull_), false);
	bTrack3DBufferFull_ = false;
	bGUIBufferFull_ = false;

	/* set parameters and calibration infos */
	UpdateParameters(_strParameterFilePath);

	bSystemRun_ = false;
	InitializeSRWLock(&lockRunFlag_);

	return bDataReady_;
}


bool CDataManager::Finalize()
{
	// make sure to release slim redears-writers lock
	AcquireSRWLockExclusive(&lockParamAccess_);
	ReleaseSRWLockExclusive(&lockParamAccess_);

	vecFrameGrabberParams_.clear();
	vecDetect2DParams_.clear();
	vecTrack2DParams_.clear();
	
	vecCalibrationInfo_.clear();
	for (int camIdx = 0; camIdx < stViewInfos_.nNumCameras; camIdx++)
	{
		//------------------------------------------
		AcquireSRWLockExclusive(&vecLockFrameImage_[camIdx]);
		//------------------------------------------
		if (!vecMatInputFrameBuffer_[camIdx].empty())
		{
			vecMatInputFrameBuffer_[camIdx].release();
		}
		//------------------------------------------
		ReleaseSRWLockExclusive(&vecLockFrameImage_[camIdx]);
		//------------------------------------------
	}

	bDataReady_ = false;
	return !bDataReady_;
}


bool CDataManager::UpdateParameters(std::string _strParameterFilePath)
{
	bDataReady_  = ReadParametersFromXML(_strParameterFilePath);
	bDataReady_ &= ReadCalibrationInfos(strCalibrationPath_, stViewInfos_.vecCamIDs);
	return bDataReady_;
}


bool CDataManager::IsFrameBufferFull(const int _camIdx)
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------
	bool returnValue = arrFrameBufferFull_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------	

	return returnValue;
}


bool CDataManager::IsDetectionBufferFull(const int _camIdx)
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------
	bool returnValue = arrDetectionBufferFull_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------

	return returnValue;
}


bool CDataManager::IsTrack2DBufferFull(const int _camIdx)
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------
	bool returnValue = arrTrack2DBufferFull_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------
	
	return returnValue;
}


bool CDataManager::IsTrack3DBufferFull()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockTrack3DResult_);
	//------------------------------------------
	bool returnValue = bTrack3DBufferFull_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockTrack3DResult_);
	//------------------------------------------

	return returnValue;
}


bool CDataManager::IsGUIBufferFull()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockGUI_);
	//------------------------------------------
	bool returnValue = bGUIBufferFull_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockGUI_);
	//------------------------------------------

	return returnValue;
}


bool CDataManager::SetFrameImage(const int _camIdx, const cv::Mat _matImage, unsigned int _frameIdx)
{
	assert(_camIdx < (int)vecFrameGrabberParams_.size() && _camIdx < (int)vecNFrameIndex_.size());
	if (arrFrameBufferFull_[_camIdx]) { return false; }

	bool result = true;	
	//------------------------------------------
	AcquireSRWLockExclusive(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------
	try
	{
		if (_matImage.size() != vecMatInputFrameBuffer_[_camIdx].size()
			|| _matImage.type() != vecMatInputFrameBuffer_[_camIdx].type())
		{
			if (!vecMatInputFrameBuffer_[_camIdx].empty()) { vecMatInputFrameBuffer_[_camIdx].release(); }
			vecMatInputFrameBuffer_[_camIdx].create(_matImage.size(), _matImage.type());
		}
		_matImage.copyTo(vecMatInputFrameBuffer_[_camIdx]);
		vecNFrameIndex_[_camIdx] = _frameIdx;
		arrFrameBufferFull_[_camIdx] = true;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] SetFrameImage: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------	

	return result;
}


bool CDataManager::SetDetectionResult(
	const int _camIdx, 
	const hj::DetectionSet _vecDetections)
{
	assert(_camIdx < (int)vecvecDetectionResult_.size());
	if (arrDetectionBufferFull_[_camIdx]) { return false; }

	bool result = true;
	//------------------------------------------
	AcquireSRWLockExclusive(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------
	try
	{
		vecvecDetectionResult_[_camIdx]  = _vecDetections;
		arrDetectionBufferFull_[_camIdx] = true;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] SetDetectionResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------

	return result;
}


bool CDataManager::SetTrack2DResult(const int _camIdx, const hj::CTrack2DResult _cTrack2DResult)
{
	assert(_camIdx < (int)vecCTrack2DResult_.size());
	if (arrTrack2DBufferFull_[_camIdx]) { return false; }

	bool result = true;	
	//------------------------------------------
	AcquireSRWLockExclusive(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------
	try
	{
		vecCTrack2DResult_[_camIdx]    = _cTrack2DResult;
		arrTrack2DBufferFull_[_camIdx] = true;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] SetTrack2DResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------	

	return result;
}


bool CDataManager::SetTrack3DResult(const hj::CTrack3DResult _cTrack3DResult)
{
	if (bTrack3DBufferFull_) { return false; }

	bool result = true;
	//------------------------------------------
	AcquireSRWLockExclusive(&lockTrack3DResult_);
	//------------------------------------------
	try
	{
		cTrack3DResult_     = _cTrack3DResult;
		bTrack3DBufferFull_ = true;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] SetTrack3DResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&lockTrack3DResult_);
	//------------------------------------------

	return result;
}


bool CDataManager::SetGUIAssignResult(const std::vector<std::pair<int, int>> _vecAssignedIDtoRealID)
{
	if (bTrack3DBufferFull_) { return false; }

	bool result = true;
	//------------------------------------------
	AcquireSRWLockExclusive(&lockGUI_);
	//------------------------------------------
	try
	{
		vecAssignedIDtoRealID_ = _vecAssignedIDtoRealID;
		bGUIBufferFull_ = true;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] SetGUIAssignResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&lockGUI_);
	//------------------------------------------

	return result;
}

void CDataManager::SetRunFlag(bool _flag)
{
	//------------------------------------------
	AcquireSRWLockExclusive(&lockRunFlag_);
	//------------------------------------------
	bSystemRun_ = _flag;
	//------------------------------------------
	ReleaseSRWLockExclusive(&lockRunFlag_);
	//------------------------------------------
}


hj::stParamFrameGrabber CDataManager::GetFrameGrabberParams(const int _camIdx)
{
	assert(bDataReady_ && _camIdx < (int)vecFrameGrabberParams_.size());
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	hj::stParamFrameGrabber returnValue = vecFrameGrabberParams_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------

	return returnValue;
}


hj::stParamDetect2D CDataManager::GetDetect2DParams(const int _camIdx)
{
	assert(bDataReady_ && _camIdx < (int)vecDetect2DParams_.size());
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	hj::stParamDetect2D returnValue = vecDetect2DParams_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------

	return returnValue;
}


hj::stParamTrack2D CDataManager::GetTrack2DParams(const int _camIdx)
{
	assert(bDataReady_ && _camIdx < (int)vecTrack2DParams_.size());
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	hj::stParamTrack2D returnValue = vecTrack2DParams_[_camIdx];
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------

	return returnValue;
}


hj::stParamAssociate3D CDataManager::GetAssociate3DParams()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	hj::stParamAssociate3D returnValue = associate3DParams_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------

	return returnValue;
}


hj::stParamEvaluator CDataManager::GetEvaluatorParams()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	hj::stParamEvaluator returnValue = evaluatorParams_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------

	return returnValue;
}


bool CDataManager::GetFrameImage(const int _camIdx, cv::Mat &_receiver, unsigned int *_frameIdx)
{
	assert(_camIdx < (int)vecMatInputFrameBuffer_.size());

	bool result = true;
	//------------------------------------------
	AcquireSRWLockShared(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------
	try
	{
		if (_receiver.size() != vecMatInputFrameBuffer_[_camIdx].size()
			|| _receiver.type() != vecMatInputFrameBuffer_[_camIdx].type())
		{
			if (!_receiver.empty()) { _receiver.release(); }
			_receiver.create(
				vecMatInputFrameBuffer_[_camIdx].size(), 
				vecMatInputFrameBuffer_[_camIdx].type());
		}
		vecMatInputFrameBuffer_[_camIdx].copyTo(_receiver);
		arrFrameBufferFull_[_camIdx] = false;
		*_frameIdx = vecNFrameIndex_[_camIdx];
	}
	catch (int e)
	{
		hj::printf_debug("[Error] GetFrameImage: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockFrameImage_[_camIdx]);
	//------------------------------------------
	return result;
}


bool CDataManager::GetDetectionResult(const int _camIdx, hj::DetectionSet *_receiver)
{
	assert(_camIdx < (int)vecvecDetectionResult_.size());

	bool result = true;
	//------------------------------------------
	AcquireSRWLockShared(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------
	try
	{
		*_receiver = vecvecDetectionResult_[_camIdx];
		arrDetectionBufferFull_[_camIdx] = false;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] GetDetectionResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockDetectionResult_[_camIdx]);
	//------------------------------------------
	return result;
}


bool CDataManager::GetTrack2DResult(const int _camIdx, hj::CTrack2DResult *_receiver)
{
	assert(_camIdx < (int)vecCTrack2DResult_.size());
	bool result = true;
	//------------------------------------------
	AcquireSRWLockShared(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------
	try
	{
		*_receiver = vecCTrack2DResult_[_camIdx];
		arrTrack2DBufferFull_[_camIdx] = false;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] GetTrack2DResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockShared(&vecLockTrack2DResult_[_camIdx]);
	//------------------------------------------
	return result;
}


bool CDataManager::GetTrack3DResult(hj::CTrack3DResult *_receiver)
{
	bool result = true;
	//------------------------------------------
	AcquireSRWLockShared(&lockTrack3DResult_);
	//------------------------------------------
	try
	{
		*_receiver = cTrack3DResult_;
		bTrack3DBufferFull_ = false;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] GetTrack3DResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockShared(&lockTrack3DResult_);
	//------------------------------------------
	return result;
}


bool CDataManager::GetGUIAssignResult(std::vector<std::pair<int, int>> *_receiver)
{
	bool result = true;
	//------------------------------------------
	AcquireSRWLockShared(&lockGUI_);
	//------------------------------------------
	try
	{
		*_receiver = vecAssignedIDtoRealID_;
		bGUIBufferFull_ = false;
	}
	catch (int e)
	{
		hj::printf_debug("[Error] GetGUIAssignResult: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockShared(&lockGUI_);
	//------------------------------------------
	return result;
}


std::vector<int> CDataManager::GetCameraIDs()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	std::vector<int> returnValue = stViewInfos_.vecCamIDs;
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	return returnValue;
}


bool CDataManager::GetDetectFlag()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	bool returnValue = bDetect_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	return returnValue;
}


bool CDataManager::GetEvaluateFlag()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	bool returnValue = bEvaluate_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockParamAccess_);
	//------------------------------------------
	return returnValue;
}


bool CDataManager::GetRunFlag()
{
	assert(bDataReady_);
	//------------------------------------------
	AcquireSRWLockShared(&lockRunFlag_);
	//------------------------------------------
	bool returnValue = bSystemRun_;
	//------------------------------------------
	ReleaseSRWLockShared(&lockRunFlag_);
	//------------------------------------------
	return returnValue;
}


bool CDataManager::ReadParametersFromXML(std::string _strFilepath)
{	
	hj::stViewInformation    stReadViewInfos;
	hj::stParamFrameGrabber  stReadParamFrameGrabber;
	hj::stParamDetect2D      stReadParamDetect2D;
	hj::stParamTrack2D       stReadParamTrack2D;
	hj::stParamAssociate3D   stReadParamAssociate3D;
	hj::stParamEvaluator     stReadParamEvaluate;

	//----------------------------------------------------------------
	// XML FILE OPEN
	//----------------------------------------------------------------
	pugi::xml_document     doc;
	pugi::xml_parse_result xmlResult = doc.load_file(_strFilepath.c_str());
	if (!xmlResult)
	{
		hj::printf_debug("[Error] Cannot open XML file from '%s'\n", _strFilepath.c_str());
		return false;
	}

	pugi::xml_node rootNode = doc.child(PARAM_TAG_PARAMETERS);
	if (!rootNode)
	{ 
		hj::printf_debug("[Error] Cannot find root node from '%s'\n", _strFilepath.c_str());
		return false; 
	}

	//----------------------------------------------------------------
	// OPERATION DATA
	//----------------------------------------------------------------
	pugi::xml_node operationNode = rootNode.child(PARAM_TAG_OPERATION);
	if (!operationNode) { return false; }
	
	strInputSource_ = "";
	for (pugi::xml_attribute_iterator attrIter = operationNode.attributes_begin();
		attrIter != operationNode.attributes_end();
		attrIter++)
	{
		if (0 == std::string(PARAM_TAG_INPUT_SOURCE).compare(attrIter->name()))
		{
			strInputSource_ = attrIter->as_string();
			if (0 == std::string(PARAM_TAG_GRABBING).compare(strInputSource_))
			{
				// default input is camera
				stReadParamFrameGrabber.nInputSource = hj::GRABBING;
				stReadParamFrameGrabber.bExistFrameRange = false;
			}
			else if (0 == std::string(PARAM_TAG_PILSNU).compare(strInputSource_))
			{
				stReadParamFrameGrabber.nInputSource = hj::PILSNU;
				stReadParamFrameGrabber.bExistFrameRange = true;
			}
			else // if (std::string::npos != strInputSource_.find(std::string(PARAM_TAG_DKU)))
			{
				stReadParamFrameGrabber.nInputSource = hj::DKU;
				stReadParamFrameGrabber.bExistFrameRange = true;
			}
		}
		else if (!std::string(PARAM_TAG_REALTIME_OPERATION_TAG).compare(attrIter->name()))
		{
			bRealtimeOperation_ = attrIter->as_bool();
		}
		else if (!std::string(PARAM_TAG_RECORD_FLAG).compare(attrIter->name()))
		{
			bRecord_ = attrIter->as_bool();
		}
		else if (!std::string(PARAM_TAG_DETECT_FLAG).compare(attrIter->name()))
		{
			bDetect_ = attrIter->as_bool();
		}
		else if (!std::string(PARAM_TAG_EVALUATION_FLAG).compare(attrIter->name()))
		{
			bEvaluate_ = attrIter->as_bool();
		}
		else if (!std::string(PARAM_TAG_DETECTOR_PATH).compare(attrIter->name()))
		{
			strDetectorPath_ = attrIter->as_string();
		}		
	}

	//----------------------------------------------------------------
	// INPUT SETTING
	//----------------------------------------------------------------
	for (pugi::xml_node settingNode = rootNode.child(PARAM_TAG_SETTING);
		settingNode;
		settingNode = settingNode.next_sibling(PARAM_TAG_SETTING))
	{
		if (0 != strInputSource_.compare(settingNode.attribute(PARAM_TAG_SETTING_NAME).as_string()))
		{
			continue;
		}

		for (pugi::xml_node_iterator nodeIter = settingNode.children().begin();
			nodeIter != settingNode.children().end();
			nodeIter++)
		{			
			/* Dataset */
			if (0 == std::string(PARAM_TAG_DATASET).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_PATH).compare(attrIter->name()))
					{
						strDatasetPath_ = attrIter->as_string();
					}
					else if (0 == std::string(PARAM_TAG_START_FRAME).compare(attrIter->name()))
					{
						stReadParamFrameGrabber.nStartFrameIndex = attrIter->as_int();						
					}
					else if (0 == std::string(PARAM_TAG_END_FRAME).compare(attrIter->name()))
					{
						stReadParamFrameGrabber.nEndFrameIndex = attrIter->as_int();						
					}
				}
			}
			/* Camera */
			else if (0 == std::string(PARAM_TAG_CAMERA).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_CAMERA_IDS).compare(attrIter->name()))
					{
						// camera number and ids
						stReadViewInfos.nNumCameras =
							hj::StringDelimite(attrIter->as_string(), ',', stReadViewInfos.vecCamIDs);
					}
					else if (0 == std::string(PARAM_TAG_FRAME_RATE).compare(attrIter->name()))
					{
						dFrameRate_ = attrIter->as_double();
					}
				}
			}			
			/* Calibration */
			else if (0 == std::string(PARAM_TAG_CALIBRATION).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_PATH).compare(attrIter->name()))
					{
						strCalibrationPath_ = attrIter->as_string();
						if (0 != std::string("").compare(strDatasetPath_))
						{
							strCalibrationPath_ = strDatasetPath_ + "/" + strCalibrationPath_;
						}
					}
				}
			}
			/* Detection */
			else if (0 == std::string(PARAM_TAG_DETECTION).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_DETECTION_TYPE).compare(attrIter->name()))
					{
						std::string strDetectionType = attrIter->as_string();
						if (0 == strDetectionType.compare("head"))
						{
							stReadParamDetect2D.nDetectionType = hj::HEAD;
						}
						else if (0 == strDetectionType.compare("fullbody"))
						{
							stReadParamDetect2D.nDetectionType = hj::FULLBODY;
						}
					}
					else if (0 == std::string(PARAM_TAG_PATH).compare(attrIter->name()))
					{
						stReadParamDetect2D.strDetectionDir = strDatasetPath_ + "/" + attrIter->as_string();
					}
					else if (0 == std::string(PARAM_TAG_IMAGE_RESCALE).compare(attrIter->name()))
					{
						stReadParamDetect2D.dImageRescale = attrIter->as_double();
						stReadParamDetect2D.dImageRescaleRecover = 1.0 / stReadParamDetect2D.dImageRescale;
					}
					else if (0 == std::string(PARAM_TAG_THRESHOLD).compare(attrIter->name()))
					{
						stReadParamDetect2D.nDetectionThreshold = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_NMS_RATIO_1).compare(attrIter->name()))
					{
						stReadParamDetect2D.dNMSOverlapRatio_1 = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_NMS_RATIO_2).compare(attrIter->name()))
					{
						stReadParamDetect2D.dNMSOverlapRatio_2 = attrIter->as_double();
					}
				}
			}
			/* Tracking 2D */
			else if (0 == std::string(PARAM_TAG_TRACKING_2D).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_VISUALIZE).compare(attrIter->name()))
					{
						stReadParamTrack2D.bVisualize = attrIter->as_bool();
					}
					else if (0 == std::string(PARAM_TAG_IMAGE_RESCALE).compare(attrIter->name()))
					{
						stReadParamTrack2D.dImageRescale = attrIter->as_double();
						stReadParamTrack2D.dImageRescaleRecover = 1.0 / stReadParamTrack2D.dImageRescale;
					}
					else if (0 == std::string(PARAM_TAG_DET_MIN_HEIGHT).compare(attrIter->name()))
					{
						stReadParamTrack2D.dDetectionMinHeight = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_DET_MAX_HEIGHT).compare(attrIter->name()))
					{
						stReadParamTrack2D.dDetectionMaxHeight = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_TRACKLET_LENGTH).compare(attrIter->name()))
					{
						stReadParamTrack2D.nMaxTrackletLength = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MIN_NUM_FEATURES).compare(attrIter->name()))
					{
						stReadParamTrack2D.nMinNumFeatures = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MAX_NUM_FEATURES).compare(attrIter->name()))
					{
						stReadParamTrack2D.nMaxNumFeatures = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_BACK_TRACKING_LENGTH).compare(attrIter->name()))
					{
						stReadParamTrack2D.nBackTrackingLength = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_FEATURE_TRACK_WIN_RATIO).compare(attrIter->name()))
					{
						stReadParamTrack2D.dFeatureTrackWindowSizeRatio = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_BOX_DISTANCE).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMaxBoxDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MIN_BOX_OVERLAP_RATIO).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMinBoxOverlapRatio = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_BOX_CENTER_DIFF_RATIO).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMaxBoxCenterDiffRatio = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MIN_OPT_FLOW_MAJ_RATIO).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMinOpticalFlowMajorityRatio = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_DETECTION_DISTANCE).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMaxDetectionDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_HEIGHT_DIFFERENCE).compare(attrIter->name()))
					{
						stReadParamTrack2D.dMaxHeightDifference = attrIter->as_double();
					}
				}
			}
			/* Association 3D */
			else if (0 == std::string(PARAM_TAG_ASSOCIATION_3D).compare(nodeIter->name()))
			{
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_VISUALIZE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.bVisualize = attrIter->as_bool();
					}
					else if (0 == std::string(PARAM_TAG_PROC_WINDOW_SIZE).compare(attrIter->name()))
					{						
						stReadParamAssociate3D.nProcWindowSize = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_K_BEST_SIZE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nKBestSize = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MAX_NUM_TRACKS_IN_OPTIMIZATION).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMaxNumTracksInOptimization = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MAX_NUM_TRACKS_IN_UNCONFIRMED_TREE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMaxNumTracksInUnconfirmedTree = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MAX_NUM_TRACKS_IN_CONFIRMED_TREE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMaxNumTracksInConfirmedTree = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_NUM_FRAMES_FOR_CONFIRMATION).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nNumFramesForConfirmation = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_DO_BRANCH_CUT).compare(attrIter->name()))
					{
						stReadParamAssociate3D.bDoBranchCut = attrIter->as_bool();
					}
					else if (0 == std::string(PARAM_TAG_MIN_TRACKLET_LENGTH).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMinTrackletLength = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_MAX_TRACKLET_DISTANCE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMaxTrackletDistance = attrIter->as_double();					
					}
					else if (0 == std::string(PARAM_TAG_CONSIDER_SENSITIVITY).compare(attrIter->name()))
					{						
						stReadParamAssociate3D.bConsiderSensitivity = attrIter->as_bool();						
					}
					else if (0 == std::string(PARAM_TAG_MAX_SENSITIVITY_ERROR).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMaxSensitivityError = attrIter->as_double();						
					}
					else if (0 == std::string(PARAM_TAG_MIN_TARGET_PROXIMITY).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMinTargetProximity = attrIter->as_double();						
					}
					else if (0 == std::string(PARAM_TAG_DEFAULT_HEIGHT).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dDefaultHeight = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_CLUSTERING_DISTANCE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dClusteringDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MIN_LINKING_PROBABILITY).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMinLinkingProbability = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_TIME_JUMP).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMaxTimeJump = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_COST_RGB_MIN_DISTANCE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostRGBMinDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_RGB_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostRGBCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_RGB_DECAY_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostRGBDecayCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_TRACKLET_LINKING_MIN_DISTANCE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostTrackletLinkingMinDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_TRACKLET_LINKING_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostTrackletLinkingCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MIN_CONSTRUCT_PROBABILITY).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMinConstructProbability = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_FP_RATE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dFPRate = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_FN_RATE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dFNRate = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_ENTER_PENALTY_FREE_LENGTH).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nEnterPenaltyFreeLength = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_BOUNDARY_DISTANCE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dBoundaryDistance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_ENTER_PROB_MAX).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dEnterProbMax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_ENTER_PROB_DECAY_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dEnterProbDecayCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_EXIT_PROB_MAX).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dExitProbMax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_EXIT_PROB_DIST_DECAY_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dExitProbDistDecayCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_EXIT_PROB_LENGTH_DECAY_COEF).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dExitProbLengthDecayCoef = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_ENTER_MAX).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostEnterMax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_COST_EXIT_MAX).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCostExitMax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_OUTPOINT).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nMaxOutpoint = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_DETECTION_ERROR).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dDetectionError = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_CALIBRATION_ERROR).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dCalibrationError = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_KALMAN_PROC_NOISE_SIGMA).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dKalmanProcessNoiseSigma = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_KALMAN_MEASURE_NOISE_SIGMA).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dKalmanMeasurementNoiseSigma = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_KALMAN_POST_ERROR_COV).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dKalmanPostErrorCovariance = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_KALMAN_CONFIDENCE_LEVEL).compare(attrIter->name()))
					{		
						stReadParamAssociate3D.nKalmanConfidenceLevel = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_VELOCITY_LEARNING_RATE).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dVelocityLearningRate = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MAX_MOVING_SPEED).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMaxMovingSpeed = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_MIN_MOVING_SPEED).compare(attrIter->name()))
					{
						stReadParamAssociate3D.dMinMovingSpeed = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_IMAGE_PATCH_WIDTH).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nImagePatchWidth = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_NUM_RGB_HISTOGRAM_BINS).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nNumRGBHistogramBins = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_SOLVER_TIME_LIMIT).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nSolverTimeLimit = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_SOLVER_MAX_ITER).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nSolverMaxIter = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_RESULT_TRAJECTORY_LENGTH).compare(attrIter->name()))
					{
						stReadParamAssociate3D.nResultTrajectoryLength = attrIter->as_int();
					}
					else if (0 == std::string(PARAM_TAG_SHOW_TREE_ID).compare(attrIter->name()))
					{						
						stReadParamAssociate3D.bShowTreeID = attrIter->as_bool();					
					}					
				}
			}/* Evaluation */
			else if (0 == std::string(PARAM_TAG_EVALUATION).compare(nodeIter->name()))
			{
				double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
				for (pugi::xml_attribute_iterator attrIter = nodeIter->attributes_begin();
					attrIter != nodeIter->attributes_end();
					attrIter++)
				{
					if (0 == std::string(PARAM_TAG_CROP_XMIN).compare(attrIter->name()))
					{
						xmin = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_CROP_XMAX).compare(attrIter->name()))
					{
						xmax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_CROP_YMIN).compare(attrIter->name()))
					{
						ymin = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_CROP_YMAX).compare(attrIter->name()))
					{
						ymax = attrIter->as_double();
					}
					else if (0 == std::string(PARAM_TAG_PATH).compare(attrIter->name()))
					{
						stReadParamEvaluate.strGroundTruthPath = attrIter->as_string();
					}
				}
				stReadParamEvaluate.cropZone =
					hj::Rect(xmin, ymin, xmax - xmin + 1, ymax - ymin + 1);
			}
		}
	}


	//----------------------------------------------------------------
	// SET VALUES
	//----------------------------------------------------------------	
	//------------------------------------------
	AcquireSRWLockExclusive(&lockParamAccess_);
	//------------------------------------------

	/* operation variables */
	stViewInfos_  = stReadViewInfos;	

	/* parameters */	
	vecCalibrationInfo_.resize(stViewInfos_.nNumCameras);
	vecFrameGrabberParams_.resize(stViewInfos_.nNumCameras);
	vecNFrameIndex_.resize(stViewInfos_.nNumCameras, 0);
	vecLockFrameImage_.resize(stViewInfos_.nNumCameras);
	vecDetect2DParams_.resize(stViewInfos_.nNumCameras);
	vecLockDetectionResult_.resize(stViewInfos_.nNumCameras);
	vecTrack2DParams_.resize(stViewInfos_.nNumCameras);
	vecLockTrack2DResult_.resize(stViewInfos_.nNumCameras);
	
	/* view */
	for (int camIdx = 0; camIdx < stViewInfos_.nNumCameras; camIdx++)
	{
		/* frame grabber */
		vecFrameGrabberParams_[camIdx] = stReadParamFrameGrabber;
		vecFrameGrabberParams_[camIdx].bDoRealtimeOperation = bRealtimeOperation_;
		vecFrameGrabberParams_[camIdx].nCamIndex   = camIdx;
		vecFrameGrabberParams_[camIdx].strInputDir = strDatasetPath_ +
			"/" + hj::FormattedString("View_%03d", stViewInfos_.vecCamIDs[camIdx]);
		InitializeSRWLock(&vecLockFrameImage_[camIdx]);		

		/* detector */
		vecDetect2DParams_[camIdx] = stReadParamDetect2D;
		vecDetect2DParams_[camIdx].nCamIndex       = camIdx;
		vecDetect2DParams_[camIdx].stViewInfo      = stViewInfos_;
		vecDetect2DParams_[camIdx].strDetectorDir  = strDetectorPath_;
		vecDetect2DParams_[camIdx].strDetectionDir += "/" +
			hj::FormattedString("View_%03d", stViewInfos_.vecCamIDs[camIdx]);
		vecDetect2DParams_[camIdx].pCalibrationInfo = &vecCalibrationInfo_[camIdx];
		InitializeSRWLock(&vecLockDetectionResult_[camIdx]);

		/* tracker */
		vecTrack2DParams_[camIdx] = stReadParamTrack2D;
		vecTrack2DParams_[camIdx].nCamIndex = camIdx;
		vecTrack2DParams_[camIdx].pCalibrationInfo = &vecCalibrationInfo_[camIdx];
		InitializeSRWLock(&vecLockTrack2DResult_[camIdx]);
	}	

	/* associator */
	associate3DParams_ = stReadParamAssociate3D;
	associate3DParams_.stViewInfo     = stViewInfos_;	
	associate3DParams_.dFrameRate     = dFrameRate_;
	associate3DParams_.nDetectionType = stReadParamDetect2D.nDetectionType;
	associate3DParams_.vecPCalibrationInfo.resize(stViewInfos_.nNumCameras);
	for (int camIdx = 0; camIdx < stViewInfos_.nNumCameras; camIdx++)
	{
		associate3DParams_.vecPCalibrationInfo[camIdx] = &vecCalibrationInfo_[camIdx];
	}

	/* evaluator */
	evaluatorParams_.strGroundTruthPath = strDatasetPath_ + "/groundTruth/ground_truth_XY.txt";	

	/* buffers */
	vecMatInputFrameBuffer_.resize(stViewInfos_.nNumCameras);
	vecvecDetectionResult_.resize(stViewInfos_.nNumCameras);
	vecCTrack2DResult_.resize(stViewInfos_.nNumCameras);
	for (int camIdx = 0; camIdx < stViewInfos_.nNumCameras; camIdx++)
	{
		InitializeSRWLock(&vecLockFrameImage_[camIdx]);
		InitializeSRWLock(&vecLockDetectionResult_[camIdx]);
		InitializeSRWLock(&vecLockTrack2DResult_[camIdx]);
	}	

	//------------------------------------------
	ReleaseSRWLockExclusive(&lockParamAccess_);
	//------------------------------------------

	return true;
}


bool CDataManager::ReadCalibrationInfos(std::string _strDir, std::vector<int> _vecCamIDs)
{
	assert(vecCalibrationInfo_.size() == _vecCamIDs.size());

	//------------------------------------------
	AcquireSRWLockExclusive(&lockParamAccess_);
	//------------------------------------------
	bool result = true;
	try
	{		
		for (int camIdx = 0; camIdx < (int)_vecCamIDs.size(); camIdx++)
		{
			vecCalibrationInfo_[camIdx].nCamIdx = camIdx;

			//----------------------------------------------------
			// READ CALIBRATION INFORMATION
			//----------------------------------------------------
			hj::printf_debug("Reading calibration information of camera %d\n", _vecCamIDs[camIdx]);

			std::string strFilePath = _strDir + "/" +
				hj::FormattedString("View_%03d.xml", _vecCamIDs[camIdx]);
			if (vecCalibrationInfo_[camIdx].cCamModel.fromXml(strFilePath))
			{
				hj::printf_debug("done ...\n");
				hj::printf_debug(" Camera position : (%f, %f, %f)\n",
					vecCalibrationInfo_[camIdx].cCamModel.cposx(),
					vecCalibrationInfo_[camIdx].cCamModel.cposy(),
					vecCalibrationInfo_[camIdx].cCamModel.cposz());
			}
			else
			{
				hj::printf_debug("fail!!\n");
				return false;
			}

			//----------------------------------------------------
			// READ PROJECTION SENSITIVITY MATRIX
			//----------------------------------------------------
			hj::printf_debug(" Read projection sensitivity : ");
			strFilePath = _strDir + "/" +
				hj::FormattedString("ProjectionSensitivity_View%03d.txt", _vecCamIDs[camIdx]);
			if (!vecCalibrationInfo_[camIdx].matProjectionSensitivity.empty())
			{
				vecCalibrationInfo_[camIdx].matProjectionSensitivity.release();
			}
			if (this->ReadProjectionSensitivity(
				vecCalibrationInfo_[camIdx].matProjectionSensitivity, strFilePath))
			{
				hj::printf_debug("success\n");
			}
			else
			{
				hj::printf_debug("fail\n");
			}

			//----------------------------------------------------
			// READ DISTANCE FROM BOUNDARY MATRIX
			//----------------------------------------------------	
			hj::printf_debug(" Read distance from boundary : ");
			strFilePath = _strDir + "/" +
				hj::FormattedString("DistanceFromBoundary_View%03d.txt", _vecCamIDs[camIdx]);
			if (!vecCalibrationInfo_[camIdx].matDistanceFromBoundary.empty())
			{
				vecCalibrationInfo_[camIdx].matDistanceFromBoundary.release();
			}
			if (this->ReadDistanceFromBoundary(
				vecCalibrationInfo_[camIdx].matDistanceFromBoundary, strFilePath))
			{
				hj::printf_debug("success\n");
			}
			else
			{
				hj::printf_debug("fail\n");
			}
		}
	}
	catch (int e)
	{
		hj::printf_debug("[Error] Cannot read calibration information with error code: %d\n", e);
		result = false;
	}
	//------------------------------------------
	ReleaseSRWLockExclusive(&lockParamAccess_);
	//------------------------------------------

	return result;
}


/************************************************************************
 Method Name: ReadProjectionSensitivity
 Description:
	- read sensitivity of projection matrix
 Input Arguments:
	- matSensitivity: map of sensitivity of a projection matrix
	- strFilepath   : path of target file
 Return Values:
	- whether file read successfully
************************************************************************/
bool CDataManager::ReadProjectionSensitivity(cv::Mat &_matSensitivity, std::string _strFilepath)
{
	FILE *fp;
	int numRows = 0, numCols = 0;
	float curSensitivity = 0.0;

	try
	{
		fopen_s(&fp, _strFilepath.c_str(), "r");
		fscanf_s(fp, "row:%d,col:%d\n", &numRows, &numCols);
		_matSensitivity = cv::Mat::zeros(numRows, numCols, CV_32FC1);

		for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
		{
			for (int colIdx = 0; colIdx < numCols; colIdx++)
			{
				fscanf_s(fp, "%f,", &curSensitivity);
				_matSensitivity.at<float>(rowIdx, colIdx) = curSensitivity;
			}
			fscanf_s(fp, "\n");
		}

		fclose(fp);
	}
	catch (int e)
	{
		hj::printf_debug("[ERROR] Cannot open projection sensitivity matrix from '%s'! error code %d\n", 
			_strFilepath.c_str(), e);
		return false;
	}
	return true;
}


/************************************************************************
 Method Name: ReadDistanceFromBoundary
 Description:
	- read matrix of distance from boundary
 Input Arguments:
	- matDistance: map of distance from boundary
	- strFilepath: path of target file
 Return Values:
	- whether file read successfully
************************************************************************/
bool CDataManager::ReadDistanceFromBoundary(cv::Mat &_matDistance, std::string _strFilepath)
{
	FILE *fp;
	int numRows = 0, numCols = 0;
	float curDistance = 0.0;

	try
	{
		fopen_s(&fp, _strFilepath.c_str(), "r");
		fscanf_s(fp, "row:%d,col:%d\n", &numRows, &numCols);
		_matDistance = cv::Mat::zeros(numRows, numCols, CV_32FC1);

		for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
		{
			for (int colIdx = 0; colIdx < numCols; colIdx++)
			{
				fscanf_s(fp, "%f,", &curDistance);
				_matDistance.at<float>(rowIdx, colIdx) = curDistance;
			}
			fscanf_s(fp, "\n");
		}

		fclose(fp);
		
	}
	catch (int e)
	{
		hj::printf_debug("[ERROR] Cannot open distance from boundary matrix from '%s'! error code %d\n", 
			_strFilepath.c_str(), e);
		return false;
	}
	return true;
}


}

//()()
//('')HAANJU.YOO
