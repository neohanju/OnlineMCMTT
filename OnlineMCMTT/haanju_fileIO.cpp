#include <Windows.h>

#include "haanju_fileIO.hpp"


/************************************************************************
 Method Name: IsDirectoryExists
 Description:
	- Check wheter the directory exists or not
 Input Arguments:
	- dirName: directory path
 Return Values:
	- true: exists / false: non-exists
************************************************************************/
bool hj::CreateDirectoryForWindows(const std::string &dirName)
{
	std::wstring wideStrDirName = L"";
	wideStrDirName.assign(dirName.begin(), dirName.end());
	if (CreateDirectory(wideStrDirName.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()) { return true; }
	return false;
}


/************************************************************************
 Method Name: GetFileList
 Description:
	- Get file name list with '_fileFormat' at '_dirPath'
 Input Arguments:
	- _dirPath   : directory path
	- _fileFormat: filename format used in 'dir' command
	- _outputVecFileNameList: name list of found files
 Return Values:
	- true: exists / false: non-exists
************************************************************************/
bool hj::GetFileList(const std::string _dirPath, const std::string _fileFormat, std::vector<std::string> &_outputVecFileNameList)
{
	HANDLE dir;
	WIN32_FIND_DATA fileData;
	std::string strDirFormat_ = _dirPath + "/" + _fileFormat;
#ifdef UNICODE
	std::wstring strDirFormat = L"";
	strDirFormat.assign(strDirFormat_.begin(), strDirFormat_.end());
#else
	std::string strDirFormat(strDirFormat_);
#endif

	// check existance
	if (INVALID_HANDLE_VALUE == (dir = FindFirstFile(strDirFormat.c_str(), &fileData)))
	{
		return false;
	}

	// read files
	_outputVecFileNameList.clear();
	do
	{
#ifdef UNICODE
		std::wstring strFileName_w(fileData.cFileName);
		std::string  strFileName = "";
		strFileName.assign(strFileName_w.begin(), strFileName_w.end());
#else
		std::string strFileName(fileData.cFileName);
#endif
		_outputVecFileNameList.push_back(strFileName);
	} while (0 != FindNextFile(dir, &fileData) || ERROR_NO_MORE_FILES != GetLastError());

	return true;
}


/************************************************************************
 Method Name: printLog
 Description:
	- print out log file
 Input Arguments:
	- filename: file path
	- strLog  : log string
 Return Values:
	- none
************************************************************************/
void hj::printLog(const char *filename, std::string strLog)
{
	try
	{
		FILE *fp;
		fopen_s(&fp, filename, "a");
		fprintf(fp, strLog.c_str());
		fclose(fp);
	}
	catch (DWORD dwError)
	{
		printf("[ERROR] cannot open logging file! error code %d\n", dwError);
		return;
	}
}


//************************************************************************
// Method Name: GrabFrame
// Description: 
//	- frame grabbing
// Input Arguments:
//	- strDatasetPath: path of the dataset
//	- camID: camera ID
//	- frameIdx: frame index
// Return Values:
//	- grabbed frame in cv::Mat container
//************************************************************************/
//cv::Mat hj::GrabFrame(std::string strDatasetPath, unsigned int camID, unsigned int frameIdx)
//{
//	char inputFilePath[300];
//	if (PSN_INPUT_TYPE)	
//	{
//		sprintf_s(inputFilePath, sizeof(inputFilePath), "%s\\View_%03d\\frame_%04d.jpg", strDatasetPath.c_str(), camID, frameIdx);						
//	} 
//	else 
//	{
//		sprintf_s(inputFilePath, sizeof(inputFilePath), "%s\\%d_%d.jpg", strDatasetPath.c_str(), camID, frameIdx);						
//	}	
//	return cv::imread(inputFilePath, cv::IMREAD_COLOR);
//}


/************************************************************************
 Method Name: ReadDetectionResultWithTxt
 Description:
	-
 Input Arguments:
	-
 Return Values:
	- Track3D*:
************************************************************************/
std::vector<hj::CDetection> hj::ReadDetectionResultWithTxt(std::string _strFilePath, DETECTION_TYPE _detectionType)
{
	std::vector<hj::CDetection> vec_result;
	int num_detection = 0;
	float x, y, w, h, temp;

	FILE *fid;
	try {		
		fopen_s(&fid, _strFilePath.c_str(), "r");
		if (NULL == fid) { return vec_result; }

		switch (_detectionType)
		{
		case hj::PARTS:	// Parts
				//// read # of detections
				//fscanf_s(fid, "numBoxes:%d\n", &num_detection);
				//vec_result.reserve(num_detection);

				//// read box infos
				//for (int detect_idx = 0; detect_idx < num_detection; detect_idx++)
				//{
				//	fscanf_s(fid, "{\n\tROOT:{%f,%f,%f,%f}\n", &x, &y, &w, &h);
				//	CDetection cur_detection;
				//	cur_detection.box = Rect((double)x, (double)y, (double)w, (double)h);

				//	// read part info
				//	cur_detection.vecPartBoxes.reserve(8);
				//	for (unsigned int partIdx = 0; partIdx < NUM_DETECTION_PART; partIdx++)
				//	{
				//		char strPartName[20];
				//		sprintf_s(strPartName, "\t%s:", DETCTION_PART_NAME[partIdx].c_str());
				//		fscanf_s(fid, strPartName);
				//		fscanf_s(fid, "{%f,%f,%f,%f}\n", &x, &y, &w, &h);
				//		PSN_Rect partBox((double)x, (double)y, (double)w, (double)h);
				//		cur_detection.vecPartBoxes.push_back(partBox);
				//	}
				//	fscanf_s(fid, "}\n");
				//	vec_result.push_back(cur_detection);
				//}			
			break;
		default: // Full-body or head
			// read # of detections
			fscanf_s(fid, "%d\n", &num_detection);
			vec_result.reserve(num_detection);

			// read box infos
			for (int detect_idx = 0; detect_idx < num_detection; detect_idx++)
			{
				fscanf_s(fid, "%f %f %f %f %f %f\n", &temp, &temp, &w, &h, &x, &y);
				CDetection cur_detection;
				cur_detection.box = Rect((double)x, (double)y, (double)w, (double)h);
				//curDetection.vecPartBoxes.reserve(8);
				vec_result.push_back(cur_detection);
			}
			break;
		}
		fclose(fid);
	}
	catch (DWORD dwError) {
		printf("[ERROR] file open error with detection result reading: %d\n", dwError);
	}
	return vec_result;
}


/************************************************************************
 Method Name: Read2DTrackResultWithTxt
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
std::vector<hj::CTrack2DResult> hj::Read2DTrackResultWithTxt(std::string strDatasetPath, unsigned int frameIdx, std::vector<unsigned int> vecCamIDs)
{
	std::vector<hj::CTrack2DResult> resultSet;
	unsigned int inCamIdx = 0, inFrameIdx = 0;
	unsigned int numModel = 0;
	char textFilePath[128];

	for (unsigned int camIdx = 0; camIdx < vecCamIDs.size(); camIdx++)
	{
		hj::CTrack2DResult curResult;
		curResult.camID = camIdx;
		curResult.frameIdx = frameIdx;

		FILE *fid;

		sprintf_s(textFilePath, sizeof(textFilePath), "%s\\trackletInput\\trackResult\\track2D_result_cam%d_frame%04d.txt", strDatasetPath.c_str(), vecCamIDs[camIdx], frameIdx);
		fopen_s(&fid, textFilePath, "r");
		fscanf_s(fid, "camIdx:%d\nframeIdx:%d\nnumModel:%d\n", &inCamIdx, &inFrameIdx, &numModel);
		for (unsigned int modelIdx = 0; modelIdx < numModel; modelIdx++)
		{
			hj::CObject2DInfo curObject;
			int x, y, w, h;

			fscanf_s(fid, "{id:%d,box:{%d,%d,%d,%d}}\n", &curObject.id, &x, &y, &w, &h);
			curObject.id--;
			curObject.box.x = (double)x;
			curObject.box.y = (double)y;
			curObject.box.w = (double)w;
			curObject.box.h = (double)h;

			curResult.object2DInfos.push_back(curObject);
		}

		resultSet.push_back(curResult);
		fclose(fid);
	}

	return resultSet;
}


/************************************************************************
 Method Name: Read2DTrackResultWithTxt
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
hj::CTrack2DResult hj::Read2DTrackResultWithTxt(std::string strDataPath, unsigned int camID, unsigned int frameIdx)
{
	char strFilename[128];
	sprintf_s(strFilename, "track2D_result_cam%d_frame%04d.txt", (int)camID, (int)frameIdx);
	std::string strFilePathAndName = strDataPath + std::string(strFilename);

	hj::CTrack2DResult result;
	result.camID = camID;
	result.frameIdx = frameIdx;

	int readInt1 = 0, readInt2 = 0;
	float readFloat = 0.0f;
	float x = 0.0f, y = 0.0f, w = 0.0f, h = 0.0f;

	FILE *fp;
	try
	{
		fopen_s(&fp, strFilePathAndName.c_str(), "r");

		// frame infos
		fscanf_s(fp, "camIdx:%d\nframeIdx:%d\n", &readInt1, &readInt2);

		// object infos
		int numObj = 0;
		fscanf_s(fp, "numObjectInfos:%d{\n", &numObj);
		result.object2DInfos.reserve(numObj);

		for (size_t objIdx = 0; objIdx < numObj; objIdx++)
		{
			hj::CObject2DInfo curObject;

			fscanf_s(fp, "\t{\n");
			////////////////
			fscanf_s(fp, "\t\tid:%d\n", &readInt1);						curObject.id = (unsigned int)readInt1;
			fscanf_s(fp, "\t\tbox:(%f,%f,%f,%f)\n", &x, &y, &w, &h);	curObject.box = hj::Rect((double)x, (double)y, (double)w, (double)h);
			fscanf_s(fp, "\t\thead:(%f,%f,%f,%f)\n", &x, &y, &w, &h);	curObject.head = hj::Rect((double)x, (double)y, (double)w, (double)h);
			fscanf_s(fp, "\t\tscore:%f\n", &readFloat);					curObject.score = (double)readFloat;

			fscanf_s(fp, "\t\tfeaturePointsPrev:%d,{", &readInt1);		curObject.prevFeatures.reserve(readInt1);
			for (size_t fIdx = 0; fIdx < readInt1; fIdx++)
			{
				fscanf_s(fp, "(%f,%f)", &x, &y);
				if (readInt1 > fIdx + 1) { fscanf_s(fp, ","); }
				curObject.prevFeatures.push_back(cv::Point2f(x, y));
			}
			fscanf_s(fp, "}\n");
			fscanf_s(fp, "\t\tfeaturePointsCurr:%d,{", &readInt1);		curObject.currFeatures.reserve(readInt1);
			for (size_t fIdx = 0; fIdx < readInt1; fIdx++)
			{
				fscanf_s(fp, "(%f,%f)", &x, &y);
				if (readInt1 > fIdx + 1) { fscanf_s(fp, ","); }
				curObject.currFeatures.push_back(cv::Point2f(x, y));
			}
			fscanf_s(fp, "}\n");
			////////////////
			fscanf_s(fp, "\t}\n");

			result.object2DInfos.push_back(curObject);
		}
		fscanf_s(fp, "}\n");

		// detection rects	
		fscanf_s(fp, "detectionRects:%d,{", &readInt1);
		result.vecDetectionRects.reserve(readInt1);
		for (size_t rectIdx = 0; rectIdx < readInt1; rectIdx++)
		{
			fscanf_s(fp, "(%f,%f,%f,%f)", &x, &y, &w, &h);
			result.vecDetectionRects.push_back(hj::Rect((double)x, (double)y, (double)w, (double)h));
			if (readInt1 > rectIdx + 1) { fscanf_s(fp, ","); }
		}
		fscanf_s(fp, "}\n");

		// tracker rects
		fscanf_s(fp, "trackerRects:%d,{", &readInt1);
		result.vecTrackerRects.reserve(readInt1);
		for (size_t rectIdx = 0; rectIdx < readInt1; rectIdx++)
		{
			fscanf_s(fp, "(%f,%f,%f,%f)", &x, &y, &w, &h);
			result.vecTrackerRects.push_back(hj::Rect((double)x, (double)y, (double)w, (double)h));
			if (readInt1 > rectIdx + 1) { fscanf_s(fp, ","); }
		}
		fscanf_s(fp, "}\n");

		fclose(fp);
	}
	catch (DWORD dwError)
	{
		printf("[ERROR](FilePrintResult) cannot open file! error code %d\n", dwError);
		return result;
	}

	return result;
}


//()()
//('')HAANJU.YOO


