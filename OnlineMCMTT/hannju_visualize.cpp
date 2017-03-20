#include <vector>
#include "haanju_visualize.hpp"


/************************************************************************
 Method Name: GenerateColors
 Description:
	- generate distinguishable color set
 Input Arguments:
	- numColor: the number of colors needed
 Return Values:
	- vector of RGB color coordinates
************************************************************************/
std::vector<cv::Scalar> hj::GenerateColors(unsigned int numColor)
{
	// refer: http://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
	//# use golden ratio
	//golden_ratio_conjugate = 0.618033988749895
	//h = rand # use random start value
	//gen_html {
	//  h += golden_ratio_conjugate
	//  h %= 1
	//  hsv_to_rgb(h, 0.5, 0.95)
	//}

	double golden_ratio_conjugate = 0.618033988749895;
	//double hVal = (double)std::rand()/(INT_MAX);
	double hVal = 0.0;
	std::vector<cv::Scalar> resultColors;
	resultColors.reserve(numColor);
	for (unsigned int colorIdx = 0; colorIdx < numColor; colorIdx++)
	{
		hVal += golden_ratio_conjugate;
		hVal = std::fmod(hVal, 1.0);
		resultColors.push_back(hj::hsv2rgb(hVal, 0.5, 0.95));
	}
	return resultColors;
}


/************************************************************************
 Method Name: hsv2rgb
 Description:
	- HSV -> RGB
 Input Arguments:
	- h: hue
	- s: saturation
	- v: value
 Return Values:
	- RGB color coordinates
************************************************************************/
cv::Scalar hj::hsv2rgb(double h, double s, double v)
{
	// refer: http://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
	//# HSV values in [0..1[
	//# returns [r, g, b] values from 0 to 255
	//def hsv_to_rgb(h, s, v)
	//  h_i = (h*6).to_i
	//  f = h*6 - h_i
	//  p = v * (1 - s)
	//  q = v * (1 - f*s)
	//  t = v * (1 - (1 - f) * s)
	//  r, g, b = v, t, p if h_i==0
	//  r, g, b = q, v, p if h_i==1
	//  r, g, b = p, v, t if h_i==2
	//  r, g, b = p, q, v if h_i==3
	//  r, g, b = t, p, v if h_i==4
	//  r, g, b = v, p, q if h_i==5
	//  [(r*256).to_i, (g*256).to_i, (b*256).to_i]
	//end

	int h_i = (int)(h * 6);
	double f = h * 6 - (double)h_i;
	double p = v * (1 - s);
	double q = v * (1 - f * s);
	double t = v * (1 - (1 - f) * s);
	double r, g, b;
	switch (h_i)
	{
	case 0: r = v; g = t; b = p; break;
	case 1: r = q; g = v; b = p; break;
	case 2: r = p; g = v; b = t; break;
	case 3: r = p; g = q; b = v; break;
	case 4: r = t; g = p; b = v; break;
	case 5: r = v; g = p; b = q; break;
	default:
		break;
	}

	return cv::Scalar((int)(r * 255), (int)(g * 255), (int)(b * 255));
}


/************************************************************************
 Method Name: getColorByID
 Description:
	- get distinguishable color by its ID. ID can be larger than the size
	  of distinguishable colors
 Input Arguments:
	- vecColors: distinguishable color space
	- nID: its ID
 Return Values:
	- RGB color coordinates
************************************************************************/
cv::Scalar hj::getColorByID(unsigned int nID, std::vector<cv::Scalar> *vecColors)
{
	if (NULL == vecColors) { return cv::Scalar(255, 255, 255); }
	unsigned int colorIdx = nID % vecColors->size();
	return (*vecColors)[colorIdx];
}


/************************************************************************
 Method Name: DrawBoxWithID
 Description:
	- draw box and
 Input Arguments:
	- imageFrame: frame for drawing
	- curRect: box for drawing
	- nID: box's ID
	- lineStyle: 0: thin / 1: thick
	- fontSize: 0: small / 1:big
	- vecColors: a distinguishable color space. when it is set by NULL, draw
	  items with white color
 Return Values:
	- none
************************************************************************/
void hj::DrawBoxWithID(cv::Mat &imageFrame, hj::Rect box, unsigned int nID, int lineStyle, int fontSize, std::vector<cv::Scalar> *vecColors)
{
	// get label length
	unsigned int labelLength = nID > 0 ? 0 : 1;
	unsigned int tempLabel = nID;
	while (tempLabel > 0)
	{
		tempLabel /= 10;
		labelLength++;
	}

	cv::Scalar curColor = hj::getColorByID(nID, vecColors);
	if (0 == fontSize)
	{
		cv::rectangle(imageFrame, box.cv(), curColor, 1);
		cv::rectangle(imageFrame, cv::Rect((int)box.x, (int)box.y - 10, 7 * labelLength, 14), curColor, CV_FILLED);
		cv::putText(imageFrame, std::to_string(nID), cv::Point((int)box.x, (int)box.y - 1), cv::FONT_HERSHEY_SIMPLEX, 0.3, cv::Scalar(0, 0, 0));
	}
	else
	{
		cv::rectangle(imageFrame, box.cv(), curColor, 1 + lineStyle);
		cv::putText(imageFrame, std::to_string(nID), cv::Point((int)box.x, (int)box.y + 40), cv::FONT_HERSHEY_SIMPLEX, 0.5, curColor);
	}
}


/************************************************************************
 Method Name: Draw3DBoxWithID
 Description:
	- draw 3D box for the target on the frame with its ID
 Input Arguments:
	- imageFrame: the frame for drawing
	- pointArray: pre-projected corners of 3D box
	- nID: box's ID
	- vecColors: a distinguishable color space. when it is set by NULL, draw
	  items with white color
 Return Values:
	- none
************************************************************************/
void hj::Draw3DBoxWithID(cv::Mat &imageFrame, std::vector<hj::Point2D> &pointArray, unsigned int nID, std::vector<cv::Scalar> *vecColors)
{
	// point 1 to 4 : roof box (clock-wise)
	// point 5 to 8 : bottom box (clock-wise, point 1 must be at a vertically upside of point 5)

	if (pointArray.size() < 8) { return; }
	cv::Scalar curColor = hj::getColorByID(nID, vecColors);

	// draw roof and bottom box
	std::vector<cv::Point> points;
	unsigned int drawOrder[10] = { 0, 1, 2, 3, 0, 4, 5, 6, 7, 4 };
	for (unsigned int pointIdx = 0; pointIdx < 10; pointIdx++)
	{
		points.push_back(pointArray[drawOrder[pointIdx]].cv());
	}

	const cv::Point *pts = (const cv::Point*)cv::Mat(points).data;
	int npts = cv::Mat(points).rows;
	cv::polylines(imageFrame, &pts, &npts, 1, false, curColor);

	// draw rest part of rectangle
	for (unsigned int pointIdx = 1; pointIdx < 4; pointIdx++)
	{
		cv::line(imageFrame, pointArray[pointIdx].cv(), pointArray[pointIdx + 4].cv(), curColor);
	}

	// get label length
	unsigned int labelLength = nID > 0 ? 0 : 1;
	unsigned int tempLabel = nID;
	while (tempLabel > 0)
	{
		tempLabel /= 10;
		labelLength++;
	}

	// draw label
	cv::rectangle(imageFrame, cv::Rect((int)pointArray[0].x, (int)pointArray[0].y - 10, 7 * labelLength, 14), curColor, CV_FILLED);
	cv::putText(imageFrame, std::to_string(nID), cv::Point((int)pointArray[0].x, (int)pointArray[0].y - 1), cv::FONT_HERSHEY_SIMPLEX, 0.3, cv::Scalar(255, 255, 255));

}


/************************************************************************
 Method Name: DrawTriangleWithID
 Description:
	- For topview, draw triangle indicates the location of target
 Input Arguments:
	- imageFrame: the frame for drawing
	- point: location of the center of triangle
	- nID: its ID
	- vecColors: a distinguishable color space. when it is set by NULL, draw
	  items with white color
 Return Values:
	- none
************************************************************************/
void hj::DrawTriangleWithID(cv::Mat &imageFrame, hj::Point2D &point, unsigned int nID, std::vector<cv::Scalar> *vecColors)
{
	double size = 20;
	cv::Scalar curColor = hj::getColorByID(nID, vecColors);
	// draw five points
	double cosd30xSize = size * 0.866025403784439;
	double sind30xSize = size * 0.5;
	std::vector<cv::Point> points;
	hj::Point2D topPoint(point.x, point.y - size);
	hj::Point2D leftPoint(point.x - cosd30xSize, point.y + sind30xSize);
	hj::Point2D rightPoint(point.x + cosd30xSize, point.y + sind30xSize);
	points.push_back(topPoint.cv());
	points.push_back(leftPoint.cv());
	points.push_back(rightPoint.cv());
	points.push_back(topPoint.cv());
	points.push_back(leftPoint.cv()); // for handling thick line

	const cv::Point *pts = (const cv::Point*)cv::Mat(points).data;
	int npts = cv::Mat(points).rows;
	cv::polylines(imageFrame, &pts, &npts, 1, false, curColor);
	cv::putText(imageFrame, std::to_string(nID), cv::Point((int)point.x, (int)point.y - 2), cv::FONT_HERSHEY_SIMPLEX, 0.3, cv::Scalar(255, 255, 255));
}


/************************************************************************
 Method Name: DrawLine
 Description:
	- draw polylines on the frame
 Input Arguments:
	- imageFrame: the frame for drawing
	- pointArray: vetices of line
	- nID: its ID (for color)
	- lineThickness: 0: thin / 1: thick
	- vecColors: a distinguishable color space. when it is set by NULL, draw
	  items with white color
 Return Values:
	- none
************************************************************************/
void hj::DrawLine(cv::Mat &imageFrame, std::vector<hj::Point2D> &pointArray, unsigned int nID, int lineThickness, std::vector<cv::Scalar> *vecColors)
{
	cv::Scalar curColor = hj::getColorByID(nID, vecColors);
	for (int pointIdx = 0; pointIdx < (int)pointArray.size() - 1; pointIdx++)
	{
		cv::line(imageFrame, pointArray[pointIdx].cv(), pointArray[pointIdx + 1].cv(), curColor, 1 + lineThickness);
	}
}


/************************************************************************
 Method Name: MakeMatTile
 Description:
	- Make one result image for multi-camera input for visualization
 Input Arguments:
	- imageArray:
	- numRows: number of rows in tile image
	- numCols: number of colums in tile image
 Return Values:
	- cv::Mat: result image
************************************************************************/
cv::Mat hj::MakeMatTile(std::vector<cv::Mat> *imageArray, unsigned int numRows, unsigned int numCols)
{
	// column first arrangement
	unsigned int numImage = (unsigned int)imageArray->size();
	unsigned int acturalNumCols = numImage < numCols ? numImage : numCols;
	unsigned int acturalNumRows = (0 == numImage % numCols) ? numImage / numCols : numImage / numCols + 1;

	// find maximum size image
	cv::Size maxSize = (*imageArray)[0].size();
	for (unsigned int imageIdx = 0; imageIdx < numImage; imageIdx++)
	{
		if ((*imageArray)[imageIdx].rows > maxSize.height) { maxSize.height = (*imageArray)[imageIdx].rows; }
		if ((*imageArray)[imageIdx].cols > maxSize.width) { maxSize.width = (*imageArray)[imageIdx].cols; }
	}

	// make augmenting matrix
	std::vector<cv::Mat> augMats;
	cv::Mat baseMat;
	baseMat = 1 == (*imageArray)[0].channels() ? cv::Mat::zeros(maxSize, CV_8UC1) : cv::Mat::zeros(maxSize, CV_8UC3);
	unsigned int numAugMats = acturalNumRows * acturalNumCols;
	for (unsigned int imageIdx = 0; imageIdx < numAugMats; imageIdx++)
	{
		cv::Mat augMat = baseMat.clone();
		if (imageIdx < numImage)
		{
			(*imageArray)[imageIdx].copyTo(augMat(cv::Rect(0, 0, (*imageArray)[imageIdx].cols, (*imageArray)[imageIdx].rows)));
		}
		augMats.push_back(augMat);
	}

	// matrix concatenation
	cv::Mat hConcatMat;
	cv::Mat resultMat;
	for (unsigned int rowIdx = 0; rowIdx < acturalNumRows; rowIdx++)
	{
		unsigned int startIdx = rowIdx * acturalNumCols;
		for (unsigned int colIdx = 0; colIdx < acturalNumCols; colIdx++)
		{
			if (0 == colIdx)
			{
				hConcatMat = augMats[startIdx].clone();
				continue;
			}
			cv::hconcat(hConcatMat, augMats[startIdx + colIdx], hConcatMat);
		}
		if (0 == rowIdx)
		{
			resultMat = hConcatMat.clone();
			continue;
		}
		cv::vconcat(resultMat, hConcatMat, resultMat);
	}

	return resultMat;
}


/************************************************************************
 Method Name: PointRescale
 Description:
	- 
 Input Arguments:
	- 
 Return Values:
	- 
************************************************************************/
hj::Point2D hj::PointRescale(hj::Point2D _point, double _rescale, double _xmin, double _ymin)
{
	hj::Point2D resultPoint = _point;

	resultPoint.x -= _xmin;
	resultPoint.y -= _ymin;
	resultPoint *= _rescale;

	return resultPoint;
}

//()()
//('')HAANJU.YOO


