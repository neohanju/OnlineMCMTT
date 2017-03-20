#include "haanju_3D.hpp"

namespace hj
{

hj::Point2D WorldToImage(hj::Point3D point3D, hj::CCalibrationInfo *pCalibInfo)
{
	hj::Point2D resultPoint2D(0, 0);
	pCalibInfo->cCamModel.worldToImage(
		point3D.x, point3D.y, point3D.z, resultPoint2D.x, resultPoint2D.y);
	resultPoint2D.x--;
	resultPoint2D.y--;
	return resultPoint2D;
}

hj::Point3D ImageToWorld(hj::Point2D point2D, double z, hj::CCalibrationInfo *pCalibInfo)
{
	hj::Point3D resultPoint3D(0, 0, 0);
	resultPoint3D.z = z;
	pCalibInfo->cCamModel.imageToWorld(
		point2D.x + 1.0, point2D.y + 1.0, z, resultPoint3D.x, resultPoint3D.y);
	//	stParam_.pCalibrationInfo->cCamModel.imageToWorld(point2D.x, point2D.y, z, resultPoint3D.x, resultPoint3D.y);
	return resultPoint3D;
}

std::vector<hj::Point2D> GetHuman3DBox(
	hj::Point3D ptHeadCenter,
	double bodyWidth,
	hj::CCalibrationInfo *pCalibInfo)
{
	std::vector<hj::Point2D> resultArray;
	double bodyWidthHalf = bodyWidth / 2.0;

	double offsetX[8] = { +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf };
	double offsetY[8] = { -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf, -bodyWidthHalf, +bodyWidthHalf, +bodyWidthHalf, -bodyWidthHalf };
	double z[8] = { ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, ptHeadCenter.z + bodyWidthHalf, 0, 0, 0, 0 };

	unsigned int numPointsOnView = 0;
	for (unsigned int vertexIdx = 0; vertexIdx < 8; vertexIdx++)
	{
		hj::Point2D curPoint = WorldToImage(hj::Point3D(ptHeadCenter.x + offsetX[vertexIdx], ptHeadCenter.y + offsetY[vertexIdx], z[vertexIdx]), pCalibInfo);
		resultArray.push_back(curPoint);
		if (curPoint.onView(pCalibInfo->cCamModel.width(), pCalibInfo->cCamModel.height())) { numPointsOnView++; }
	}

	// if there is no point can be seen from the view, clear a vector of point
	if (0 == numPointsOnView) { resultArray.clear(); }

	return resultArray;
}

}