#include <assert.h>
#include "Point3DSmoother.h"

namespace hj
{

CPoint3DSmoother::CPoint3DSmoother(void)
{
}

CPoint3DSmoother::~CPoint3DSmoother(void)
{
}

int CPoint3DSmoother::Reset(std::vector<Point3D> &points)
{
	// insertion	
	std::vector<double> pointsX(points.size(), 0.0), pointsY(points.size(), 0.0), pointsZ(points.size(), 0.0);
	for (int pointIdx = 0; pointIdx < points.size(); pointIdx++)
	{
		pointsX[pointIdx] = points[pointIdx].x;
		pointsY[pointIdx] = points[pointIdx].y;
		pointsZ[pointIdx] = points[pointIdx].z;
	}
	smootherX_.Reset(pointsX);
	smootherY_.Reset(pointsY);
	smootherZ_.Reset(pointsZ);

	// smoothing
	this->Update(0, (int)points.size());
	return 0;
}

/************************************************************************
 Method Name: Insert
 Description:
	- insert a single point into the smoothing sequence
 Input Arguments:
	- point: point being inserted
 Return Values:
	- the number of points that are modified by smoothing (from the back
	  of the sequence)
************************************************************************/
int CPoint3DSmoother::Insert(Point3D &point)
{
	// insertion
	int refreshPos = smootherX_.Insert(point.x);
	smootherY_.Insert(point.y);
	smootherZ_.Insert(point.z);

	// smoothing
	this->Update(refreshPos, 1);
	return refreshPos;
}

/************************************************************************
 Method Name: Insert
 Description:
	- insert muliple points into the smoothing sequence
 Input Arguments:
	- points: points being inserted
 Return Values:
	- the number of points that are modified by smoothing (from the back
	  of the sequence
************************************************************************/
int CPoint3DSmoother::Insert(std::vector<Point3D> &points)
{
	// insertion
	std::vector<double> pointsX(points.size(), 0.0), pointsY(points.size(), 0.0), pointsZ(points.size(), 0.0);
	for (int pointIdx = 0; pointIdx < points.size(); pointIdx++)
	{
		pointsX[pointIdx] = points[pointIdx].x;
		pointsY[pointIdx] = points[pointIdx].y;
		pointsZ[pointIdx] = points[pointIdx].z;
	}
	int refreshPos = smootherX_.Insert(pointsX);
	smootherY_.Insert(pointsY);
	smootherZ_.Insert(pointsZ);

	// smoothing
	this->Update(refreshPos, (int)points.size());
	return refreshPos;
}

/************************************************************************
 Method Name: ReplaceBack
 Description:
	- replace the last point of the sequence
 Input Arguments:
	- point: point for replacement
 Return Values:
	- the number of points that are modified by smoothing (from the back
	  of the sequence
************************************************************************/
int CPoint3DSmoother::ReplaceBack(Point3D &point)
{
	// replacement
	int refreshPos = smootherX_.ReplaceBack(point.x);
	smootherY_.ReplaceBack(point.y);
	smootherZ_.ReplaceBack(point.z);
	smoothedPoints_.pop_back();

	// smoothing
	this->Update(refreshPos, 1);
	return refreshPos;
}

/************************************************************************
 Method Name: SetQsets
 Description:
	- Set precomputed Q sets for smoothing
 Input Arguments:
	- Qsets: precomputed Q sets
 Return Values:
	- none
************************************************************************/
void CPoint3DSmoother::SetQsets(std::vector<Qset> *Qsets)
{
	smootherX_.SetPrecomputedQsets(Qsets);
	smootherY_.SetPrecomputedQsets(Qsets);
	smootherZ_.SetPrecomputedQsets(Qsets);
}

/************************************************************************
 Method Name: PopBack
 Description:
	- pop the last point from smoothed sequence
 Input Arguments:
	- none
 Return Values:
	- none
************************************************************************/
void CPoint3DSmoother::PopBack(void)
{
	assert(0 < smoothedPoints_.size());
	smoothedPoints_.pop_back();
	smootherX_.PopBack();
	smootherY_.PopBack();
	smootherZ_.PopBack();
}

/************************************************************************
 Method Name: SetSmoother
 Description:
	- set the smoother by copying other smoother
 Input Arguments:
	- data: original sequence
	- smoothedPoints: smoothed sequence
	- span: window size of smoother
	- degree: degree for smoother
 Return Values:
	- none
************************************************************************/
void CPoint3DSmoother::SetSmoother(std::deque<Point3D> &data, std::deque<Point3D> &smoothedPoints, int span, int degree)
{
	smoothedPoints_ = smoothedPoints;
	size_ = smoothedPoints_.size();
	back_ = &smoothedPoints_.back();
	std::deque<double> pointX, pointY, pointZ, smoothX, smoothY, smoothZ;
	for (int pos = 0; pos < smoothedPoints_.size(); pos++)
	{
		pointX.push_back(data[pos].x);
		pointY.push_back(data[pos].y);
		pointZ.push_back(data[pos].z);
		smoothX.push_back(smoothedPoints_[pos].x);
		smoothY.push_back(smoothedPoints_[pos].y);
		smoothZ.push_back(smoothedPoints_[pos].z);
	}
	smootherX_.SetSmoother(pointX, smoothX, span, degree);
	smootherY_.SetSmoother(pointY, smoothY, span, degree);
	smootherZ_.SetSmoother(pointZ, smoothZ, span, degree);
}

/************************************************************************
 Method Name: GetResult
 Description:
	- get the specific point of smoothed sequence
 Input Arguments:
	- pos: the position of target point
 Return Values:
	- target point
************************************************************************/
Point3D CPoint3DSmoother::GetResult(int pos)
{
	assert(pos < smoothedPoints_.size());
	return smoothedPoints_[pos];
}

/************************************************************************
 Method Name: GetResults
 Description:
	- get points in the specific range of smoothed sequence
 Input Arguments:
	- startPos: starting position of the range
	- endPos: ending position of the range
 Return Values:
	- points sequence in the range
************************************************************************/
std::vector<Point3D> CPoint3DSmoother::GetResults(int startPos, int endPos)
{
	if (endPos < 0) { endPos = (int)smoothedPoints_.size(); }
	assert(endPos <= smoothedPoints_.size() && startPos <= endPos);
	std::vector<Point3D> results;
	results.insert(results.end(), smoothedPoints_.begin() + startPos, smoothedPoints_.begin() + endPos);
	return results;
}

/************************************************************************
 Method Name: GetSmoother
 Description:
	- backup the current smoother (by input arguments)
 Input Arguments:
	- data: backup of the original sequence
	- smoothedPoints: backup of the smoothed sequence
	- span: backup of current smoother's window size
	- degree: backup of the degree of current smoother
 Return Values:
	- none
************************************************************************/
void CPoint3DSmoother::GetSmoother(std::deque<Point3D> &data, std::deque<Point3D> &smoothedPoints, int &span, int &degree)
{
	smoothedPoints = smoothedPoints_;
	std::deque<double> pointX, pointY, pointZ, smoothX, smoothY, smoothZ;
	smootherX_.GetSmoother(pointX, smoothX, span, degree);
	smootherY_.GetSmoother(pointY, smoothY, span, degree);
	smootherZ_.GetSmoother(pointZ, smoothZ, span, degree);
	data.resize(pointX.size());
	for (int pos = 0; pos < smoothedPoints_.size(); pos++)
	{
		data[pos].x = pointX[pos];
		data[pos].y = pointY[pos];
		data[pos].z = pointZ[pos];
	}
}

/************************************************************************
 Method Name: Update
 Description:
	- do smoothing in the continuous sub interval of the sequence
 Input Arguments:
	- refresPos: the starting position of the interval
	- numPoints: the size of the interval
 Return Values:
	- none
************************************************************************/
void CPoint3DSmoother::Update(int refreshPos, int numPoints)
{
	int endPos = (int)smoothedPoints_.size() + numPoints;
	smoothedPoints_.erase(smoothedPoints_.begin() + refreshPos, smoothedPoints_.end());
	for (int pos = refreshPos; pos < endPos; pos++)
	{
		smoothedPoints_.push_back(Point3D(smootherX_.GetResult(pos), smootherY_.GetResult(pos), smootherZ_.GetResult(pos)));
	}
	size_ = smoothedPoints_.size();
	back_ = &smoothedPoints_.back();
}

}

//()()
//('')HAANJU.YOO
