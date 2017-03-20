/******************************************************************************
* Title        : Point3DSmoother
* Author       : Haanju Yoo
* Initial Date : 2015.05.04 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  :
*	Smoothing 3D point via 'SGSmoother' class instance at each coordinate.
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

#pragma once

#include "types.hpp"
#include "SGSmoother.h"

namespace hj
{

class CPoint3DSmoother
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CPoint3DSmoother();
	~CPoint3DSmoother();

	// setter
	int  Reset(std::vector<Point3D> &points);
	int  Insert(Point3D &point);
	int  Insert(std::vector<Point3D> &points);
	int  ReplaceBack(Point3D &point);
	void SetQsets(std::vector<Qset> *Qsets);
	void PopBack(void);
	void SetSmoother(std::deque<Point3D> &data, std::deque<Point3D> &smoothedPoints, int span, int degree);

	// getter
	Point3D GetResult(int pos);
	std::vector<Point3D> GetResults(int startPos, int endPos = -1);
	size_t size(void) const { return size_; }
	Point3D back(void) const { return *back_; }
	void GetSmoother(std::deque<Point3D> &data, std::deque<Point3D> &smoothedPoints, int &span, int &degree);

private:
	void Update(int refreshPos, int numPoints);

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
private:
	size_t size_;
	Point3D *back_;
	CSGSmoother smootherX_, smootherY_, smootherZ_;
	std::deque<Point3D> smoothedPoints_;
};

}

//()()
//('')HAANJU.YOO

