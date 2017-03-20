/******************************************************************************
* Title        : haanju_3D
* Author       : Haanju Yoo
* Initial Date : 2016.12.20 (ver. 0.9)
* Version Num. : 0.9
* Description  : contains operations related with 3D.
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

#include "types_3D.h"

namespace hj
{
	hj::Point2D WorldToImage(hj::Point3D point3D, hj::CCalibrationInfo *pCalibInfo);
	hj::Point3D ImageToWorld(hj::Point2D point2D, double z, hj::CCalibrationInfo *pCalibInfo);
	std::vector<hj::Point2D> GetHuman3DBox(
		hj::Point3D ptHeadCenter,
		double bodyWidth,
		hj::CCalibrationInfo *pCalibInfo);
}

//()()
//('')HAANJU.YOO
