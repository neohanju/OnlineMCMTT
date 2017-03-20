/******************************************************************************
* Title        : haanju_matOperations
* Author       : Haanju Yoo
* Initial Date : 2013.08.29 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  : contains matrix operations.
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

#ifndef __HAANJU_MATOPERATIONS_HPP__
#define __HAANJU_MATOPERATIONS_HPP__

#include <opencv2\core.hpp>

namespace hj
{

void appendRow(cv::Mat &dstMat, cv::Mat row);
void appendCol(cv::Mat &dstMat, cv::Mat col);

template<typename _Tp> _Tp MatTotalSum(cv::Mat &inputMat)
{
	_Tp resultSum = 0;
	for (int rowIdx = 0; rowIdx < inputMat.rows; rowIdx++)
	{
		for (int colIdx = 0; colIdx < inputMat.cols; colIdx++)
		{
			resultSum += inputMat.at<_Tp>(rowIdx, colIdx);
		}
	}
	return resultSum;
}

template<typename _Tp> bool MatLowerThan(cv::Mat &inputMat, _Tp compValue)
{
	//_Tp maxValue = std::max_element(inputMat.begin(), inputMat.end());
	//for(int rowIdx = 0; rowIdx < inputMat.rows; rowIdx++)
	//{
	//	for(int colIdx = 0; colIdx < inputMat.cols; colIdx++)
	//	{
	//		if(inputMat.at<_Tp>(rowIdx, colIdx) >= compValue)
	//		{
	//			return false;
	//		}
	//	}
	//}
	return (std::max_element(inputMat.begin(), inputMat.end()) < compValue) ? true : false;
}

template<typename _Tp> bool MatContainLowerThan(cv::Mat &inputMat, _Tp compValue)
{
	for (int rowIdx = 0; rowIdx < inputMat.rows; rowIdx++)
	{
		for (int colIdx = 0; colIdx < inputMat.cols; colIdx++)
		{
			if (inputMat.at<_Tp>(rowIdx, colIdx) < compValue)
			{
				return true;
			}
		}
	}
	return false;
}

template<typename _Tp> std::vector<_Tp> mat2vec_C1(cv::Mat &inputMat)
{
	std::vector<_Tp> vecResult;

	for (int rowIdx = 0; rowIdx < inputMat.rows; rowIdx++)
	{
		for (int colIdx = 0; colIdx < inputMat.cols; colIdx++)
		{
			vecResult.push_back(inputMat.at<_Tp>(rowIdx, colIdx));
		}
	}

	return vecResult;
}

}


#endif

//()()
//('')HAANJU.YOO


