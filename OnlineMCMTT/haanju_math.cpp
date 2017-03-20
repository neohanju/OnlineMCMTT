#include <vector>
#include <algorithm>
#include <opencv2\core.hpp>
#include "types.hpp"
#include "haanju_math.hpp"


/////////////////////////////////////////////////////////////////////////
// PREDEFINE
/////////////////////////////////////////////////////////////////////////

// reference for erf function: http://math.stackexchange.com/questions/263216/error-function-erf-with-better-precision
static const double haanju_tiny = 1e-300,
haanju_half = 5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
haanju_one = 1.00000000000000000000e+00,  /* 0x3FF00000, 0x00000000 */
haanju_two = 2.00000000000000000000e+00,  /* 0x40000000, 0x00000000 */
										  /* haanju_c = (float)0.84506291151 */
haanju_erx = 8.45062911510467529297e-01,  /* 0x3FEB0AC1, 0x60000000 */

// Coefficients for approximation to erf on [0,0.84375]
haanju_efx = 1.28379167095512586316e-01,  /* 0x3FC06EBA, 0x8214DB69 */
haanju_efx8 = 1.02703333676410069053e+00, /* 0x3FF06EBA, 0x8214DB69 */
haanju_pp0 = 1.28379167095512558561e-01,  /* 0x3FC06EBA, 0x8214DB68 */
haanju_pp1 = -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
haanju_pp2 = -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
haanju_pp3 = -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
haanju_pp4 = -2.37630166566501626084e-05, /* 0xBEF8EAD6, 0x120016AC */
haanju_qq1 = 3.97917223959155352819e-01,  /* 0x3FD97779, 0xCDDADC09 */
haanju_qq2 = 6.50222499887672944485e-02,  /* 0x3FB0A54C, 0x5536CEBA */
haanju_qq3 = 5.08130628187576562776e-03,  /* 0x3F74D022, 0xC4D36B0F */
haanju_qq4 = 1.32494738004321644526e-04,  /* 0x3F215DC9, 0x221C1A10 */
haanju_qq5 = -3.96022827877536812320e-06, /* 0xBED09C43, 0x42A26120 */

// Coefficients for approximation to erf in [0.84375,1.25]
haanju_pa0 = -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
haanju_pa1 = 4.14856118683748331666e-01,  /* 0x3FDA8D00, 0xAD92B34D */
haanju_pa2 = -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
haanju_pa3 = 3.18346619901161753674e-01,  /* 0x3FD45FCA, 0x805120E4 */
haanju_pa4 = -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
haanju_pa5 = 3.54783043256182359371e-02,  /* 0x3FA22A36, 0x599795EB */
haanju_pa6 = -2.16637559486879084300e-03, /* 0xBF61BF38, 0x0A96073F */
haanju_qa1 = 1.06420880400844228286e-01,  /* 0x3FBB3E66, 0x18EEE323 */
haanju_qa2 = 5.40397917702171048937e-01,  /* 0x3FE14AF0, 0x92EB6F33 */
haanju_qa3 = 7.18286544141962662868e-02,  /* 0x3FB2635C, 0xD99FE9A7 */
haanju_qa4 = 1.26171219808761642112e-01,  /* 0x3FC02660, 0xE763351F */
haanju_qa5 = 1.36370839120290507362e-02,  /* 0x3F8BEDC2, 0x6B51DD1C */
haanju_qa6 = 1.19844998467991074170e-02,  /* 0x3F888B54, 0x5735151D */

// Coefficients for approximation to erfc in [1.25,1/0.35]
haanju_ra0 = -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
haanju_ra1 = -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
haanju_ra2 = -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
haanju_ra3 = -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
haanju_ra4 = -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
haanju_ra5 = -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
haanju_ra6 = -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
haanju_ra7 = -9.81432934416914548592e+00, /* 0xC023A0EF, 0xC69AC25C */
haanju_sa1 = 1.96512716674392571292e+01,  /* 0x4033A6B9, 0xBD707687 */
haanju_sa2 = 1.37657754143519042600e+02,  /* 0x4061350C, 0x526AE721 */
haanju_sa3 = 4.34565877475229228821e+02,  /* 0x407B290D, 0xD58A1A71 */
haanju_sa4 = 6.45387271733267880336e+02,  /* 0x40842B19, 0x21EC2868 */
haanju_sa5 = 4.29008140027567833386e+02,  /* 0x407AD021, 0x57700314 */
haanju_sa6 = 1.08635005541779435134e+02,  /* 0x405B28A3, 0xEE48AE2C */
haanju_sa7 = 6.57024977031928170135e+00,  /* 0x401A47EF, 0x8E484A93 */
haanju_sa8 = -6.04244152148580987438e-02, /* 0xBFAEEFF2, 0xEE749A62 */

// Coefficients for approximation to erfc in [1/.35,28]
haanju_rb0 = -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
haanju_rb1 = -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
haanju_rb2 = -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
haanju_rb3 = -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
haanju_rb4 = -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
haanju_rb5 = -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
haanju_rb6 = -4.83519191608651397019e+02, /* 0xC07E384E, 0x9BDC383F */
haanju_sb1 = 3.03380607434824582924e+01,  /* 0x403E568B, 0x261D5190 */
haanju_sb2 = 3.25792512996573918826e+02,  /* 0x40745CAE, 0x221B9F0A */
haanju_sb3 = 1.53672958608443695994e+03,  /* 0x409802EB, 0x189D5118 */
haanju_sb4 = 3.19985821950859553908e+03,  /* 0x40A8FFB7, 0x688C246A */
haanju_sb5 = 2.55305040643316442583e+03,  /* 0x40A3F219, 0xCEDF3BE6 */
haanju_sb6 = 4.74528541206955367215e+02,  /* 0x407DA874, 0xE79FE763 */
haanju_sb7 = -2.24409524465858183362e+01; /* 0xC03670E2, 0x42712D62 */


/////////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////////

/************************************************************************
 Method Name: nchoosek
 Description:
	- generate combinations consist of k elements from integer set {0, ..., n-1}
 Input Arguments:
	- n: choose from {0, ..., n-1}
	- k: the number of elements for output combinations
 Return Values:
	- result combinations
************************************************************************/
std::vector<std::vector<unsigned int>> hj::nchoosek(int n, int k)
{
	std::vector<std::vector<unsigned int>> outputCombinations;
	if (n < k || n <= 0) { return outputCombinations; }

	std::vector<bool> v(n);
	std::vector<unsigned int> curCombination;

	std::fill(v.begin() + k, v.end(), true);
	do
	{
		curCombination.clear();
		curCombination.reserve(k);
		for (int idx = 0; idx < n; ++idx)
		{
			if (!v[idx]) { curCombination.push_back(unsigned int(idx)); }
		}
		outputCombinations.push_back(curCombination);
	} while (std::next_permutation(v.begin(), v.end()));
	return outputCombinations;
}

/************************************************************************
 Method Name: erf
 Description:
	- error function
 Input Arguments:
	- x: any real value
 Return Values:
	- error
************************************************************************/
double hj::erf(double x)
{
	bool bFastErf = false;

	if (bFastErf)
	{
		// FAST VERSION
		// constants
		double a1 = 0.254829592;
		double a2 = -0.284496736;
		double a3 = 1.421413741;
		double a4 = -1.453152027;
		double a5 = 1.061405429;
		double p = 0.3275911;

		// Save the sign of x
		int sign = 1;
		if (x < 0) { sign = -1; }
		x = fabs(x);

		// A&S formula 7.1.26
		double t = 1.0 / (1.0 + p*x);
		double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

		return sign*y;
	}

	// ACCURATE VERSION
	int n0, hx, ix, i;
	double R, S, P, Q, s, y, z, r;
	n0 = ((*(int*)&haanju_one) >> 29) ^ 1;
	hx = *(n0 + (int*)&x);
	ix = hx & 0x7fffffff;
	if (ix >= 0x7ff00000)
	{
		/* erf(nan)=nan */
		i = ((unsigned)hx >> 31) << 1;
		return (double)(1 - i) + haanju_one / x; /* erf(+-inf)=+-1 */
	}

	if (ix < 0x3feb0000)
	{
		/* |x|<0.84375 */
		if (ix < 0x3e300000)
		{
			/* |x|<2**-28 */
			if (ix < 0x00800000) { return 0.125 * (8.0 * x + haanju_efx8 * x); } /*avoid underflow */
			return x + haanju_efx * x;
		}
		z = x*x;
		r = haanju_pp0 + z * (haanju_pp1 + z * (haanju_pp2 + z * (haanju_pp3 + z * haanju_pp4)));
		s = haanju_one + z * (haanju_qq1 + z * (haanju_qq2 + z * (haanju_qq3 + z * (haanju_qq4 + z * haanju_qq5))));
		y = r / s;
		return x + x*y;
	}

	if (ix < 0x3ff40000)
	{
		/* 0.84375 <= |x| < 1.25 */
		s = fabs(x) - haanju_one;
		P = haanju_pa0 + s*(haanju_pa1 + s*(haanju_pa2 + s*(haanju_pa3 + s*(haanju_pa4 + s*(haanju_pa5 + s*haanju_pa6)))));
		Q = haanju_one + s*(haanju_qa1 + s*(haanju_qa2 + s*(haanju_qa3 + s*(haanju_qa4 + s*(haanju_qa5 + s*haanju_qa6)))));
		if (hx >= 0)
		{
			return haanju_erx + P / Q;
		}
		else
		{
			return -haanju_erx - P / Q;
		}
	}

	if (ix >= 0x40180000)
	{
		/* inf>|x|>=6 */
		if (hx >= 0)
		{
			return haanju_one - haanju_tiny;
		}
		else
		{
			return haanju_tiny - haanju_one;
		}
	}

	x = fabs(x);
	s = haanju_one / (x * x);

	if (ix < 0x4006DB6E)
	{
		/* |x| < 1/0.35 */
		R = haanju_ra0 + s * (haanju_ra1 + s * (haanju_ra2 + s * (haanju_ra3 + s * (haanju_ra4 + s * (haanju_ra5 + s * (haanju_ra6 + s * haanju_ra7))))));
		S = haanju_one + s * (haanju_sa1 + s * (haanju_sa2 + s * (haanju_sa3 + s * (haanju_sa4 + s * (haanju_sa5 + s * (haanju_sa6 + s * (haanju_sa7 + s * haanju_sa8)))))));
	}
	else
	{
		/* |x| >= 1/0.35 */
		R = haanju_rb0 + s * (haanju_rb1 + s * (haanju_rb2 + s * (haanju_rb3 + s * (haanju_rb4 + s * (haanju_rb5 + s * haanju_rb6)))));
		S = haanju_one + s * (haanju_sb1 + s * (haanju_sb2 + s * (haanju_sb3 + s * (haanju_sb4 + s * (haanju_sb5 + s * (haanju_sb6 + s * haanju_sb7))))));
	}
	z = x;
	*(1 - n0 + (int*)&z) = 0;
	r = exp(-z * z - 0.5625) * exp((z - x) * (z + x) + R / S);
	if (hx >= 0)
	{
		return haanju_one - r / x;
	}
	else
	{
		return r / x - haanju_one;
	}
}

/************************************************************************
 Method Name: erfc
 Description:
	- complement error function
 Input Arguments:
	- any real value
 Return Values:
	- complement of error
************************************************************************/
double hj::erfc(double x)
{
	bool bFastErfc = false;

	// FAST VERSION
	if (bFastErfc) { return 1 - hj::erf(x); }

	// ACCURATE VERSION
	int n0, hx, ix;
	double R, S, P, Q, s, y, z, r;
	n0 = ((*(int*)&haanju_one) >> 29) ^ 1;
	hx = *(n0 + (int*)&x);
	ix = hx & 0x7fffffff;

	/* erfc(nan)=nan */
	/* erfc(+-inf)=0,2 */
	if (ix >= 0x7ff00000) { return (double)(((unsigned)hx >> 31) << 1) + haanju_one / x; }

	if (ix < 0x3feb0000)
	{
		/* |x|<0.84375 */
		if (ix < 0x3c700000) { return haanju_one - x; } /* |x|<2**-56 */
		z = x * x;
		r = haanju_pp0 + z * (haanju_pp1 + z * (haanju_pp2 + z * (haanju_pp3 + z * haanju_pp4)));
		s = haanju_one + z * (haanju_qq1 + z * (haanju_qq2 + z * (haanju_qq3 + z * (haanju_qq4 + z * haanju_qq5))));
		y = r / s;
		if (hx < 0x3fd00000)
		{
			/* x<1/4 */
			return haanju_one - (x + x * y);
		}
		else
		{
			r = x * y;
			r += (x - haanju_half);
			return haanju_half - r;
		}
	}

	if (ix < 0x3ff40000)
	{
		/* 0.84375 <= |x| < 1.25 */
		s = fabs(x) - haanju_one;
		P = haanju_pa0 + s * (haanju_pa1 + s * (haanju_pa2 + s * (haanju_pa3 + s * (haanju_pa4 + s * (haanju_pa5 + s * haanju_pa6)))));
		Q = haanju_one + s * (haanju_qa1 + s * (haanju_qa2 + s * (haanju_qa3 + s * (haanju_qa4 + s * (haanju_qa5 + s * haanju_qa6)))));
		if (hx >= 0)
		{
			z = haanju_one - haanju_erx;
			return z - P / Q;
		}
		else
		{
			z = haanju_erx + P / Q;
			return haanju_one + z;
		}
	}

	if (ix < 0x403c0000)
	{
		/* |x|<28 */
		x = fabs(x);
		s = haanju_one / (x * x);
		if (ix < 0x4006DB6D)
		{
			/* |x| < 1/.35 ~ 2.857143*/
			R = haanju_ra0 + s * (haanju_ra1 + s * (haanju_ra2 + s * (haanju_ra3 + s * (haanju_ra4 + s * (haanju_ra5 + s * (haanju_ra6 + s * haanju_ra7))))));
			S = haanju_one + s * (haanju_sa1 + s * (haanju_sa2 + s * (haanju_sa3 + s * (haanju_sa4 + s * (haanju_sa5 + s * (haanju_sa6 + s * (haanju_sa7 + s * haanju_sa8)))))));
		}
		else
		{
			/* |x| >= 1/.35 ~ 2.857143 */
			if (hx<0 && ix >= 0x40180000) { return haanju_two - haanju_tiny; } /* x < -6 */
			R = haanju_rb0 + s * (haanju_rb1 + s * (haanju_rb2 + s * (haanju_rb3 + s * (haanju_rb4 + s * (haanju_rb5 + s * haanju_rb6)))));
			S = haanju_one + s * (haanju_sb1 + s * (haanju_sb2 + s * (haanju_sb3 + s * (haanju_sb4 + s * (haanju_sb5 + s *(haanju_sb6 + s * haanju_sb7))))));
		}
		z = x;
		*(1 - n0 + (int*)&z) = 0;
		r = exp(-z * z - 0.5625) * exp((z - x) * (z + x) + R / S);
		if (hx > 0)
		{
			return r / x;
		}
		else
		{
			return haanju_two - r / x;
		}
	}
	else
	{
		if (hx > 0)
		{
			return haanju_tiny * haanju_tiny;
		}
		else
		{
			return haanju_two - haanju_tiny;
		}
	}
}

/************************************************************************
 Method Name: histogram
 Description:
	- generate RGB color histogram with input patch
 Input Arguments:
	- singleChannelImage: input image patch
	- numBin: the number of bins
 Return Values:
	- color histrogram contained in cv::Mat
************************************************************************/
cv::Mat hj::histogram(cv::Mat singleChannelImage, int numBin)
{
	assert(1 == singleChannelImage.channels() && numBin > 0);
	cv::Mat vecHistogram = cv::Mat::zeros(numBin, 1, CV_64FC1);
	double binSize = 256.0 / numBin;

	uchar *pData = singleChannelImage.data;
	int curHistogramBinIdx = 0;
	for (int idx = 0; idx < singleChannelImage.rows * singleChannelImage.cols; idx++, pData++)
	{
		curHistogramBinIdx = (int)std::floor((double)(*pData) / binSize);
		vecHistogram.at<double>(curHistogramBinIdx, 0)++;
	}

	return vecHistogram;
}

/************************************************************************
 Method Name: NSSD
 Description:
	- Normalized SSD between two (single channel) patches
 Input Arguments:
	- matPatch1: operand 1
	- matPatch2: operand 2
	- mask     : mask image for weighting
 Return Values:
	- double: result NSSD
************************************************************************/
double hj::NSSD(cv::Mat matPatch1, cv::Mat matPatch2, cv::InputArray mask)
{
	cv::Mat compPatch = matPatch2;
	if (matPatch1.size != matPatch2.size)
	{
		cv::resize(compPatch, compPatch, cv::Size(matPatch1.rows, matPatch1.cols));
	}

	double  ssd = 0.0;
	bool    bMasking = false;
	cv::Mat matMask = cv::Mat::ones(cv::Size(matPatch1.cols, matPatch1.rows), CV_32FC1);
	if (!mask.empty())
	{
		bMasking = true;
		matMask = mask.getMat();
	}
	//double invNumPixels = 1.0 / (matPatch1.channels() * matPatch1.rows * matPatch1.cols);
	double numPixels = 0.0;
	//for (int channelIdx = 0; channelIdx < matPatch1.channels(); channelIdx++)
	//{
	for (int rowIdx = 0; rowIdx < matPatch1.rows; rowIdx++)
	{
		for (int colIdx = 0; colIdx < matPatch1.cols; colIdx++)
		{
			if (bMasking && 0 == matMask.at<float>(rowIdx, colIdx)) { continue; }
			ssd += matMask.at<float>(rowIdx, colIdx) * abs((double)matPatch1.at<unsigned char>(rowIdx, colIdx) - (double)compPatch.at<unsigned char>(rowIdx, colIdx));
			numPixels++;
		}
	}
	//}
	ssd /= numPixels;

	return ssd;
}

/************************************************************************
 Method Name: IsLineSegmentIntersect
 Description:
	- Check line-line intersetction
 Input Arguments:
	- line1: the first line
	- line2: the second line
 Return Values:
	- wheather the lines are intersect or not
************************************************************************/
bool hj::IsLineSegmentIntersect(hj::Line3D &line1, hj::Line3D &line2)
{
	// from: "Tricks of the Windows Game Programming Gurus"
	double s1x, s1y, s2x, s2y;
	s1x = line1.second.x - line1.first.x;
	s1y = line1.second.y - line1.first.y;
	s2x = line2.second.x - line2.first.x;
	s2y = line2.second.y - line2.first.y;

	double s, t;
	s = (-s1y * (line1.first.x - line2.first.x) + s1x * (line1.first.y - line2.first.y)) / (-s2x * s1y + s1x * s2y);
	t = (s2x * (line1.first.y - line2.first.y) - s2y * (line1.first.x - line2.first.x)) / (-s2x * s1y + s1x * s2y);

	if (0 <= s && s <= 1 && 0 <= t && t <= 1) { return true; }
	return false;
}

/************************************************************************
 Method Name: Triangulation
 Description:
	- find the point which has the smallest total distance to each line
 Input Arguments:
	- line1: the first line
	- line2: the second line
 Return Values:
	- the distance between lines
************************************************************************/
double hj::Triangulation(hj::Line3D &line1, hj::Line3D &line2, hj::Point3D &midPoint3D)
{
	hj::Point3D line1Direct = line1.first - line1.second;
	hj::Point3D line2Direct = line2.first - line2.second;
	hj::Point3D lineOffset = line2.second - line1.second;

	cv::Mat matA(2, 2, CV_32FC1);
	matA.at<float>(0, 0) = (float)line1Direct.dot(line1Direct);
	matA.at<float>(0, 1) = (float)line1Direct.dot(-line2Direct);
	matA.at<float>(1, 0) = (float)line2Direct.dot(line1Direct);
	matA.at<float>(1, 1) = (float)line2Direct.dot(-line2Direct);

	cv::Mat vecB(2, 1, CV_32FC1);
	vecB.at<float>(0, 0) = (float)line1Direct.dot(lineOffset);
	vecB.at<float>(1, 0) = (float)line2Direct.dot(lineOffset);

	cv::Mat vecT = matA.inv() * vecB;

	line1Direct *= (double)vecT.at<float>(0, 0);
	line2Direct *= (double)vecT.at<float>(1, 0);
	hj::Point3D closePoint1 = line1.second + line1Direct;
	hj::Point3D closePoint2 = line2.second + line2Direct;

	midPoint3D = (closePoint1 + closePoint2) / 2;

	return (closePoint1 - closePoint2).norm_L2();
}


/************************************************************************
 Method Name: EstimateDetectionHeight
 Description:
	- Estimate a height in 3D space and (optionally) estimate top center
	  point of the box in 3D space.
 Input Arguments:
	- box       : input box
	- calibInfo : calibration information for 3D back-projection
	- z         : coordinate of z-axis of bottom center
	- location3D: (output) 3D point of top center of the input box
 Return Values:
	- height in mm unit.
************************************************************************/
double hj::EstimateBoxHeight(const hj::Rect box, hj::CCalibrationInfo &calibInfo, double z, hj::Point3D *location3D)
{
	hj::Point2D topCenter = box.topCenter(), bottomCenter = box.bottomCenter();
	hj::Point3D P11, P12, P21, P22;

	// top point
	P11.z = z;
	P12.z = P11.z + 2000;
	calibInfo.cCamModel.imageToWorld(topCenter.x, topCenter.y, P11.z, P11.x, P11.y);
	calibInfo.cCamModel.imageToWorld(topCenter.x, topCenter.y, P12.z, P12.x, P12.y);

	// bottom point
	P21.z = z;
	calibInfo.cCamModel.imageToWorld(bottomCenter.x, bottomCenter.y, P21.z, P21.x, P21.y);
	P22 = P21;
	P22.z = 2000;

	if (NULL != location3D) { *location3D = P21; }

	hj::Point3D topPoint(0.0, 0.0, 0.0);
	hj::Triangulation(hj::Line3D(P11, P12), hj::Line3D(P21, P22), topPoint);

	return (topPoint - P21).norm_L2();
}


/************************************************************************
 Method Name: PolarToEuclidean
 Description:
	- 
 Input Arguments:
	- 
 Return Values:
	- 
************************************************************************/
hj::Point2D hj::PolarToEuclidean(double _distance, double _degree)
{
	hj::Point2D pointEuclidean;
	pointEuclidean.x = _distance * std::cos((CV_PI / 180.0) * _degree);
	pointEuclidean.y = _distance * std::sin((CV_PI / 180.0) * _degree);
	return pointEuclidean;
}


/************************************************************************
 Method Name: MinDistanceBetweenPointSets
 Description:
	-
 Input Arguments:
	-
 Return Values:
	-
************************************************************************/
double  hj::MinDistanceBetweenPointSets(
	std::vector<hj::Point2D> &_vecPoints1, 
	std::vector<hj::Point2D> &_vecPoints2)
{
	double minDistance = DBL_MAX;
	double curDistance = 0.0;
	for (int i1 = 0; i1 < _vecPoints1.size(); i1++)
	{
		for (int i2 = 0; i2 < _vecPoints2.size(); i2++)
		{
			curDistance = (_vecPoints1[i1] - _vecPoints2[i2]).norm_L2();
			if (curDistance < minDistance)
			{
				minDistance = curDistance;
			}
		}
	}
	return minDistance;
}



//()()
//('')HAANJU.YOO


