/******************************************************************************
* Title        : haanju_math
* Author       : Haanju Yoo
* Initial Date : 2013.08.29 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  : contains basic math operations.
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

#ifndef __HAANJU_MATH_HPP__
#define __HAANJU_MATH_HPP__

#include <opencv2\core.hpp>
#include "types.hpp"

namespace hj
{

std::vector<std::vector<unsigned int>> nchoosek(int n, int k);
double  erf(double x);
double  erfc(double x);
cv::Mat histogram(cv::Mat singleChannelImage, int numBin);
double  NSSD(cv::Mat matPatch1, cv::Mat matPatch2, cv::InputArray mask = cv::noArray());
bool    IsLineSegmentIntersect(Line3D &line1, Line3D &line2);
double  Triangulation(Line3D &line1, Line3D &line2, Point3D &midPoint3D);
double  EstimateBoxHeight(const hj::Rect box, hj::CCalibrationInfo &calibInfo, double z = 0.0, hj::Point3D *location3D = NULL);
hj::Point2D PolarToEuclidean(double _distance, double _degree);
double  MinDistanceBetweenPointSets(std::vector<hj::Point2D> &_vecPoints1, std::vector<hj::Point2D> &_vecPoints2);


/************************************************************************
 Method Name: Median
 Description:
	- find the median value of the vector
 Input Arguments:
	- vector: target vector
 Return Values:
	- the median value
************************************************************************/
template<typename T>T hj::Median(std::vector<T> vector)
{
	assert(0 < vector.size());
	std::sort(vector.begin(), vector.end());
	size_t medianIndex = vector.size() >> 1;
	if (0 == vector.size() % 2)
	{
		return (vector[medianIndex - 1] + vector[medianIndex]) / 2.0;
	}
	return vector[medianIndex];
}

}

#endif

//()()
//('')HAANJU.YOO


