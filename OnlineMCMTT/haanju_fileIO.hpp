/******************************************************************************
* Title        : haanju_fileIO
* Author       : Haanju Yoo
* Initial Date : 2013.08.29 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  : contains file input/output operation functions.
******************************************************************************
*               .__                           __.
*                \ `\~~---..---~~~~~~--.---~~| /
*                 `~-.   `                   .~         _____
*                     ~.                .--~~    .---~~~    /
*                      / .-.      .-.      |  <~~        __/
*                     |  |_|      |_|       \  \     .--'
*                    /-.      -       .-.    |  \_   \_
*                    \-'   -..-..-    `-'    |    \__  \_
*                     `.                     |     _/  _/
*                     ~-                .,-\   _/  _/
*                      /                 -~~~~\ /_  /_
*                     |               /   |    \  \_  \_
*                     |   /          /   /      | _/  _/
*                     |  |          |   /    .,-|/  _/
*                     )__/           \_/    -~~~| _/
*                       \                      /  \
*                        |           |        /_---`
*                        \    .______|      ./
*                        (   /        \    /
*                        `--'          /__/
*
******************************************************************************/

#ifndef __HAANJU_FILEIO_HPP__
#define __HAANJU_FILEIO_HPP__

#include "types_3D.h"

namespace hj
{

bool CreateDirectoryForWindows(const std::string &dirName);
bool GetFileList(const std::string _dirPath, const std::string _fileFormat, std::vector<std::string> &_outputVecFileNameList);
void printLog(const char *filename, std::string strLog);
std::string MakeTrackIDList(TrackSet *tracks);
std::vector<CDetection> ReadDetectionResultWithTxt(std::string _strFilePath, DETECTION_TYPE _detectionType = FULLBODY);
std::vector<CTrack2DResult> Read2DTrackResultWithTxt(std::string strDatasetPath, unsigned int frameIdx, std::vector<unsigned int> vecCamIDs);
CTrack2DResult Read2DTrackResultWithTxt(std::string strDataPath, unsigned int camID, unsigned int frameIdx);

}


#endif

//()()
//('')HAANJU.YOO


