/******************************************************************************
* Title        : types
* Author       : Haanju Yoo
* Initial Date : 2013.08.29 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.16)
* Description  :
*	Including basic types for overall tracking operation. Types, which are only
*	used in 3D association, are in an additional file named 'types_3D.hpp'.
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

#ifndef __HAANJU_TYPES_HPP__
#define __HAANJU_TYPES_HPP__

#include <list>
#include <deque>
#include <time.h>
#include <opencv2\imgproc\imgproc.hpp>

#include "cameraModel.h"

#define MAX_NUM_SAME_THREAD (10)
#define HJ_LIDAR_RESOLUTION (360*2)

namespace hj
{

enum DETECTION_TYPE { FULLBODY = 0, HEAD, PARTS };
enum INPUT_SOURCE   { GRABBING = 0, PILSNU, DKU, PETS };

/////////////////////////////////////////////////////////////////////////
// GEOMETRY TYPES
/////////////////////////////////////////////////////////////////////////
class Point2D
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	Point2D() : x(0), y(0) {}
	Point2D(double x, double y) : x(x), y(y) {}
	Point2D(cv::Point2f a) : x((double)a.x), y((double)a.y) {}

	// operators
	Point2D& operator=(const Point2D &a)     { x = a.x; y = a.y; return *this; }
	Point2D& operator=(const cv::Point &a)   { x = a.x; y = a.y; return *this; }
	Point2D& operator=(const cv::Point2f &a) { x = (double)a.x; y = (double)a.y; return *this; }
	Point2D& operator+=(const Point2D &a)    { x = x + a.x; y = y + a.y; return *this; }
	Point2D& operator-=(const Point2D &a)    { x = x - a.x; y = y - a.y; return *this; }
	Point2D& operator+=(const double s)      { x = x + s; y = y + s; return *this; }
	Point2D& operator-=(const double s)      { x = x - s; y = y - s; return *this; }
	Point2D& operator*=(const double s)      { x = x * s; y = y * s; return *this; }
	Point2D& operator/=(const double s)      { x = x / s; y = y / s; return *this; }
	Point2D  operator+(const Point2D &a)     { return Point2D(x + a.x, y + a.y); }
	Point2D  operator+(const double s)       { return Point2D(x + s, y + s); }
	Point2D  operator-(const Point2D &a)     { return Point2D(x - a.x, y - a.y); }
	Point2D  operator-(const double s)       { return Point2D(x - s, y - s); }
	Point2D  operator-()                     { return Point2D(-x, -y); }
	Point2D  operator*(const double s)       { return Point2D(x * s, y * s); }
	Point2D  operator/(const double s)       { return Point2D(x / s, y / s); }
	bool     operator==(const Point2D &a)    { return (x == a.x && y == a.y); }
	bool     operator==(const cv::Point &a)  { return (x == a.x && y == a.y); }

	// methods
	double  norm_L2()             { return std::sqrt(x * x + y * y); }
	double  dot(const Point2D &a) { return x * a.x + y * a.y; }
	Point2D scale(double scale)   { return Point2D(x * scale, y * scale); }
	bool    onView(const unsigned int width, const unsigned int height)
	{
		if (x < 0) { return false; }
		if (x >= (double)width) { return false; }
		if (y < 0) { return false; }
		if (y >= (double)height) { return false; }
		return true;
	}	

	// data conversion
	cv::Point cv() { return cv::Point((int)x, (int)y); }

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	// data
	double x;
	double y;
};


class Point3D
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	Point3D() : x(0), y(0), z(0) {}
	Point3D(double x, double y, double z) :x(x), y(y), z(z) {}

	// operators
	Point3D& operator=(const Point3D &a)      { x = a.x; y = a.y; z = a.z; return *this; }
	Point3D& operator=(const cv::Point3d &a)  { x = a.x; y = a.y; z = a.z; return *this; }
	Point3D& operator+=(const Point3D &a)     { x = x + a.x; y = y + a.y; z = z + a.z; return *this; }
	Point3D& operator-=(const Point3D &a)     { x = x - a.x; y = y - a.y; z = z - a.z; return *this; }
	Point3D& operator+=(const double s)       { x = x + s; y = y + s; z = z + s; return *this; }
	Point3D& operator-=(const double s)       { x = x - s; y = y - s; z = z - s; return *this; }
	Point3D& operator*=(const double s)       { x = x * s; y = y * s; z = z * s; return *this; }
	Point3D& operator/=(const double s)       { x = x / s; y = y / s; z = z / s; return *this; }
	Point3D  operator+(const Point3D &a)      { return Point3D(x + a.x, y + a.y, z + a.z); }
	Point3D  operator+(const double s)        { return Point3D(x + s, y + s, z + s); }
	Point3D  operator-(const Point3D &a)      { return Point3D(x - a.x, y - a.y, z - a.z); }
	Point3D  operator-(const double s)        { return Point3D(x - s, y - s, z - s); }
	Point3D  operator-()                      { return Point3D(-x, -y, -z); }
	Point3D  operator*(const double s)        { return Point3D(x * s, y * s, z * s); }
	Point3D  operator/(const double s)        { return Point3D(x / s, y / s, z / s); }
	bool     operator==(const Point3D &a)     { return (x == a.x && y == a.y && z == a.z); }
	bool     operator==(const cv::Point3d &a) { return (x == a.x && y == a.y && z == a.z); }

	// methods
	double norm_L2() { return std::sqrt(x * x + y * y + z * z); }
	double dot(const Point3D &a) { return x * a.x + y * a.y + z * a.z; }

	// data conversion
	cv::Point3d cv() { return cv::Point3d(x, y, z); }

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	// data
	double x;
	double y;
	double z;
};

typedef std::pair<Point2D, Point2D> Line2D;
typedef std::pair<Point3D, Point3D> Line3D;
struct FOV
{
	Point3D corner[4];
};

class Rect
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	Rect() : x(0), y(0), w(0), h(0) {}
	Rect(double x, double y, double w, double h) :x(x), y(y), w(w), h(h) {}

	// operators
	Rect& operator=(const Rect &a)      { x = a.x; y = a.y; w = a.w; h = a.h; return *this; }
	Rect& operator=(const cv::Rect &a)  { x = a.x; y = a.y; w = a.width; h = a.height; return *this; }
	bool  operator==(const Rect &a)     const { return (x == a.x && y == a.y && w == a.w && h == a.h); }
	bool  operator==(const cv::Rect &a) const { return (x == a.x && y == a.y && w == a.width && h == a.height); }
	Rect& operator*=(const double s)    { x = x * s; y = y * s; w = w * s; h = h * s; return *this; }
	Rect operator*(const double s)      const { return Rect(x * s, y * s, w * s, h * s); }

	// methods	
	Point2D bottomCenter() const { return Point2D(x + std::ceil(w / 2.0), y + h); }
	Point2D topCenter()    const { return Point2D(x + std::ceil(w / 2.0), y); }
	Point2D center()       const { return Point2D(x + std::ceil(w / 2.0), y + std::ceil(h / 2.0)); }
	Point2D reconstructionPoint() const
	{
		return this->bottomCenter();
		//switch(PSN_INPUT_TYPE)
		//{
		//case 1:
		//	return this->bottomCenter();
		//	break;
		//default:
		//	return this->center();
		//	break;
		//}
	}
	Rect cropWithSize(const double width, const double height) const
	{
		double newX = std::max(0.0, x);
		double newY = std::max(0.0, y);
		double newW = std::min(width - newX - 1, w);
		double newH = std::min(height - newY - 1, h);
		return Rect(newX, newY, newW, newH);
	}
	Rect scale(double scale) const
	{
		return Rect(x * scale, y * scale, w * scale, h * scale);
	}
	double area() { return w * h; }
	bool contain(const Point2D &a)     const { return (a.x >= x && a.x < x + w && a.y >= y && a.y < y + h); }
	bool contain(const cv::Point2f &a) const { return ((double)a.x >= x && (double)a.x < x + w && (double)a.y >= y && (double)a.y < y + h); }
	bool overlap(const Rect &a) const
	{
		return (std::max(x + w, a.x + a.w) - std::min(x, a.x) < w + a.w) && (std::max(y + h, a.y + a.h) - std::min(y, a.y) < h + a.h) ? true : false;
	}
	double distance(const Rect &a) const
	{
		Point3D descriptor1 = Point3D(x + w / 2.0, y + h / 2.0, w);
		Point3D descriptor2 = Point3D(a.x + a.w / 2.0, a.y + a.h / 2.0, a.w);
		return (descriptor1 - descriptor2).norm_L2() / std::min(w, a.w);
	}
	double overlappedArea(const Rect &a) const
	{
		double overlappedWidth = std::min(x + w, a.x + a.w) - std::max(x, a.x);
		if (0.0 >= overlappedWidth) { return 0.0; }
		double overlappedHeight = std::min(y + h, a.y + a.h) - std::max(y, a.y);
		if (0.0 >= overlappedHeight) { return 0.0; }
		return overlappedWidth * overlappedHeight;
	}

	// conversion
	cv::Rect cv() const { return cv::Rect((int)x, (int)y, (int)w, (int)h); }

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	double x;
	double y;
	double w;
	double h;
};


/////////////////////////////////////////////////////////////////////////
// OBJECTS
/////////////////////////////////////////////////////////////////////////
class CDetection
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CDetection() : box(0.0, 0.0, 0.0, 0.0), score(0.0) {}
	CDetection(const CDetection &c) { *this = c; };

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	Rect   box;
	double score;
};
typedef std::vector<CDetection> DetectionSet;


class CLidarDetection
{
public:
	CLidarDetection() : center_(Point2D(0.0, 0.0)) {};
	~CLidarDetection() {};

	void SetPoints(std::vector<Point2D> &_vecPoints) { vecPoints_ = _vecPoints; };
	void SetCenter(hj::Point2D _center) { center_ = _center; }

	std::vector<Point2D> GetPoints() { return vecPoints_; }
	hj::Point2D GetCenter() { return center_; }
private:
	Point2D center_;
	std::vector<Point2D> vecPoints_;
};
typedef std::vector<CLidarDetection> LidarDetectionSet;


class CObject2DInfo
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CObject2DInfo() : id(0), box(0.0, 0.0, 0.0, 0.0), head(0.0, 0.0, 0.0, 0.0), score(0.0) {}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	Rect         box;
	Rect         head;
	double       score;
	std::vector<cv::Point2f> prevFeatures;
	std::vector<cv::Point2f> currFeatures;
};


class CObject3DInfo
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CObject3DInfo(int numCameras) : id(0)
	{
		numCam = numCameras;
		recentPoint2Ds.resize(numCam);
		point3DBox.resize(numCam);
		rectInViews.resize(numCam, Rect(0.0, 0.0, 0.0, 0.0));
		bVisibleInViews.resize(numCam, false);
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int id;
	size_t numCam;
	std::vector<Point3D> recentPoints;
	std::vector<Point3D> curDetectionPosition;

	// values for each camera
	std::vector<std::vector<Point2D>> recentPoint2Ds;
	std::vector<std::vector<Point2D>> point3DBox;	
	std::vector<Rect> rectInViews;
	std::vector<bool> bVisibleInViews;
};


class CObjectLidarInfo
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CObjectLidarInfo() : id(0), frameIdx(0) {}
	CObjectLidarInfo(int _id) : id(_id), frameIdx(0) {}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	int id;
	int frameIdx;
	//std::vector<Point3D> recentPoints;
	Point3D location;
};


/////////////////////////////////////////////////////////////////////////
// RESULTS
/////////////////////////////////////////////////////////////////////////
class CLidarScanResult
{
public:
	CLidarScanResult()
	{
		Clear();
	}
	~CLidarScanResult() {}

	CLidarScanResult& operator=(const CLidarScanResult &a)
	{ 
		memcpy(this->arrDistances_, a.arrDistances_, sizeof(std::pair<double, double>)*HJ_LIDAR_RESOLUTION);
		size_ = a.size();
		return *this;
	}

	void   SetDistance(int _pos, double _val) 
	{
		arrDistances_[_pos].second = _val; 
		if (_pos >= size_) { size_ = (size_t)_pos + 1; }
	}
	void   SetAngle(int _pos, double _val)
	{
		arrDistances_[_pos].first  = _val;
		if (_pos >= size_) { size_ = (size_t)_pos + 1; }
	}
	double GetDistance(int _pos) { return arrDistances_[_pos].second; }
	double GetAngle(int _pos)    { return arrDistances_[_pos].first; }
	bool   Clear()
	{
		memset(arrDistances_, 0, sizeof(std::pair<double, double>) * HJ_LIDAR_RESOLUTION);
		size_ = 0;
		return true;
	}
	size_t size() const { return size_; }

private:
	std::pair<double, double> arrDistances_[HJ_LIDAR_RESOLUTION]; // angle / distance
	hj::Point2D arrCoordinates_[HJ_LIDAR_RESOLUTION];
	size_t size_;
};


class CTrackLidarResult
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CTrackLidarResult() : frameIdx(0), timeStamp(0) {}

	// operator
	CTrackLidarResult& operator=(const CTrackLidarResult &a)
	{ 
		frameIdx         = a.frameIdx; 
		timeStamp        = a.timeStamp;
		objectLidarInfos = a.objectLidarInfos;
		return *this; 
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:	
	unsigned int frameIdx;
	unsigned int timeStamp;
	std::vector<CObjectLidarInfo> objectLidarInfos;
};


// 2D tracking result at each frame
class CTrack2DResult
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CTrack2DResult() : camID(0), frameIdx(0), timeStamp(0) {}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int camID;
	unsigned int frameIdx;
	unsigned int timeStamp;
	std::vector<CObject2DInfo> object2DInfos;

	std::vector<Rect> vecDetectionRects;
	std::vector<Rect> vecTrackerRects;
	cv::Mat matMatchingCost;
};


// 3D tracking result at each frame
class CTrack3DResult
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	// constructors
	CTrack3DResult() : frameIdx(0), timeStamp(0), processingTime(0.0) {}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int frameIdx;
	unsigned int timeStamp;
	double       processingTime;
	std::vector<CObject3DInfo> object3DInfos;
};


// tracking result
class CTrackingInfo
{

};


/////////////////////////////////////////////////////////////////////////
// CALIBRATION
/////////////////////////////////////////////////////////////////////////
class CCalibrationInfo
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CCalibrationInfo() {}
	~CCalibrationInfo()
	{
		if (!matProjectionSensitivity.empty()) { matProjectionSensitivity.release(); }
		if (!matDistanceFromBoundary.empty())  { matDistanceFromBoundary.release(); }
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	unsigned int nCamIdx;
	Etiseo::CameraModel cCamModel;
	cv::Mat matProjectionSensitivity;
	cv::Mat matDistanceFromBoundary;
};

class CLidarCalibrationInfo
{
	//----------------------------------------------------------------
	// METHODS
	//----------------------------------------------------------------
public:
	CLidarCalibrationInfo()  {}
	~CLidarCalibrationInfo() {}
	bool Load(std::string _strFilePath)
	{
		float center[2];
		try
		{
			FILE *fp;			
			fopen_s(&fp, _strFilePath.c_str(), "r");
			for (int i = 0; i < 9; i++)
			{
				fscanf_s(fp, "%f,", &H_[i]);
			}
			fscanf_s(fp, "%f,%f", &center[0], &center[1]);
			fclose(fp);
		}
		catch (int e)
		{
			printf("File read error with %d\n", e);
			return false;
		}		
		center_.x = (double)center[0];
		center_.y = (double)center[1];
		
		bInit_ = true;

		return bInit_;
	}

	hj::Point2D SensorToWorld(hj::Point2D _point)
	{
		_point.x; _point.y;
		double Xw   = (double)H_[0] * _point.x + (double)H_[3] * _point.y + (double)H_[6];
		double Yw   = (double)H_[1] * _point.x + (double)H_[4] * _point.y + (double)H_[7];
		double invS = (double)H_[2] * _point.x + (double)H_[5] * _point.y + (double)H_[8];

		Xw /= invS;
		Yw /= invS;
		return hj::Point2D(Xw, Yw);
	}

	//----------------------------------------------------------------
	// VARIABLES
	//----------------------------------------------------------------
public:
	bool bInit_;
	float H_[9];	
	// H_[0], H_[3], H_[6]
	// H_[1], H_[4], H_[7]
	// H_[2], H_[5], H_[8]
	hj::Point2D center_;
};


/////////////////////////////////////////////////////////////////////////
// BUFFER RELATED
/////////////////////////////////////////////////////////////////////////
class CCircularIndex
{
public:
	CCircularIndex() : size_(0), currentIndex_(0) {}
	~CCircularIndex() {}
	// manupulation
	void setSize(int size) { size_ = size; currentIndex_ = 0; }
	void operator++(void)  { if (size_ && ++currentIndex_ >= size_) { currentIndex_ -= size_; } }
	void operator--(void)  { if (size_ && --currentIndex_ < 0) { currentIndex_ += size_; } }
	void operator++(int)   { ++*this; }
	void operator--(int)   { --*this; }
	// position access
	int size() { return size_; }
	int next()
	{
		int nextPos = currentIndex_ + (1 & size_);
		if (nextPos >= size_) { nextPos -= size_; }
		return nextPos;
	}
	int current() { return currentIndex_; }
	int previous()
	{
		int prevPos = currentIndex_ - (1 & size_);
		if (prevPos < 0) { prevPos += size_; }
		return prevPos;
	}
private:
	int size_;
	int currentIndex_;
};

class CFIFOIndicator
{
public:
	CFIFOIndicator::CFIFOIndicator(int size)
	{
		size_ = size;
		this->end();
	};
	CFIFOIndicator::~CFIFOIndicator() {};

	int current() { return pos_; }
	int next() { return ++pos_; }
	int previous() { return --pos_; }
	int end() { return pos_ = size_ - 1; }
private:
	int size_;
	int pos_;
};

class CMatFIFOBuffer
{
public:
	CMatFIFOBuffer() : bInit_(false) {}
	CMatFIFOBuffer(int bufferSize) : bInit_(true), bufferSize_(bufferSize) {}
	~CMatFIFOBuffer() { clear(); }

	bool set(int bufferSize)
	{
		if (bInit_) { this->clear(); }
		try
		{
			bufferSize_ = bufferSize;
		}
		catch (int e)
		{
			printf("An execption occured in CircularBuffer::set. Exeption number is %d.\n", e);
			return false;
		}

		return bInit_ = true;
	}

	bool clear()
	{
		if (!bInit_) { return true; }
		try
		{
			for (int bufferIdx = 0; bufferIdx < buffer_.size(); bufferIdx++)
			{
				this->remove(bufferIdx);
			}
		}
		catch (int e)
		{
			printf("An execption occured in CircularBuffer::clear. Exeption number is %d.\n", e);
			return false;
		}
		bufferSize_ = 0;
		bInit_ = false;

		return true;
	}

	bool insert(cv::Mat newMat)
	{
		//assert(!bInit_
		//	&& newMat.rows == elementSize_.height
		//	&& newMat.cols == elementSize_.width
		//	&& newMat.type() == elementType_);

		try
		{
			cv::Mat newBufferMat = newMat.clone();
			buffer_.push_back(newBufferMat);

			// circulation
			if (bufferSize_ < buffer_.size())
			{
				buffer_.pop_front();
			}
		}
		catch (int e)
		{
			printf("An execption occured in CircularBuffer::insert. Exeption number is %d.\n", e);
			return false;
		}

		return true;
	}

	bool insert_resize(cv::Mat newMat, cv::Size resizeParam)
	{
		//assert(!bInit_
		//	&& newMat.rows == elementSize_.height
		//	&& newMat.cols == elementSize_.width
		//	&& newMat.type() == elementType_);

		try
		{
			cv::Mat newBufferMat;
			cv::resize(newMat, newBufferMat, resizeParam);
			buffer_.push_back(newBufferMat);

			// circulation
			if (bufferSize_ < buffer_.size())
			{
				buffer_.pop_front();
			}
		}
		catch (int e)
		{
			printf("An execption occured in CircularBuffer::insert. Exeption number is %d.\n", e);
			return false;
		}

		return true;
	}

	bool remove(int pos)
	{
		assert(bInit_ && pos < bufferSize_);
		if (pos >= buffer_.size())
		{
			return true;
		}
		if (!buffer_[pos].empty())
		{
			buffer_[pos].release();
		}
		return true;
	}

	size_t size()
	{
		return bufferSize_;
	}

	size_t num_elements()
	{
		return buffer_.size();
	}

	cv::Mat front()
	{
		assert(bInit_);
		return buffer_.front();
	}

	cv::Mat back()
	{
		assert(bInit_);
		return buffer_.back();
	}

	cv::Mat get(int pos)
	{
		assert(bInit_ && pos < bufferSize_);
		return buffer_[pos];
	}

	int get_back_idx()
	{
		return (int)buffer_.size();
	}

	CFIFOIndicator get_indicator()
	{
		CFIFOIndicator newIndicator(bufferSize_);
		return newIndicator;
	}

	/* iterators */
	typedef std::deque<cv::Mat>::iterator iterator;
	typedef std::deque<cv::Mat>::const_iterator const_iterator;
	typedef std::deque<cv::Mat>::reverse_iterator reverse_iterator;
	typedef std::deque<cv::Mat>::const_reverse_iterator const_reverse_iterator;

	iterator begin() { return buffer_.begin(); }
	iterator end()   { return buffer_.end();  }
	const_iterator begin() const { return buffer_.begin(); }
	const_iterator end()   const { return buffer_.end(); }
	reverse_iterator rbegin() { return buffer_.rbegin(); }
	reverse_iterator rend()   { return buffer_.rend(); }
	const_reverse_iterator rbegin() const { return buffer_.rbegin(); }
	const_reverse_iterator rend()   const { return buffer_.rend(); }

private:
	bool bInit_;
	int  bufferSize_;
	std::deque<cv::Mat> buffer_;
};


/////////////////////////////////////////////////////////////////////////
// PARAMETERS
/////////////////////////////////////////////////////////////////////////
struct stViewInformation
{	
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stViewInformation()
		: nNumCameras(0)
	{};
	~stViewInformation() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int nNumCameras;
	std::vector<int> vecCamIDs;
};

struct stLidarInformation
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stLidarInformation()
		: nNumSensors(0)
	{};
	~stLidarInformation() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int nNumSensors;
	int nMotorSpeed;
	std::vector<int> vecPortNumbers;
	std::vector<int> vecSensorIDs;
};

struct stParamSocket
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamSocket()
		: strIP("192.168.0.13")
		, nPort(5007)
	{};
	~stParamSocket() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	std::string strIP;
	int nPort;
};

struct stParamFrameGrabber
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamFrameGrabber() 
		: nCamIndex(0)
		, bExistFrameRange(false)
		, bDoRealtimeOperation(false)
		, nStartFrameIndex(0)
		, nEndFrameIndex(0)
		, nInputSource(GRABBING)
		, strInputDir("") 
	{};
	~stParamFrameGrabber() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int  nCamIndex;
	bool bExistFrameRange;
	bool bDoRealtimeOperation;
	int  nStartFrameIndex;
	int  nEndFrameIndex;
	INPUT_SOURCE nInputSource;
	std::string strInputDir;
};

struct stParamLidarGrabber
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamLidarGrabber()
		: nPortNumber(0)
		, nSensorIndex (0)
		, nMotorSpeed(200)
		, bExistFrameRange(false)
		, bDoRealtimeOperation(false)
		, nStartFrameIndex(0)
		, nEndFrameIndex(0)
		, nInputSource(hj::GRABBING)
		, strInputDir("")
	{};
	~stParamLidarGrabber() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int  nPortNumber;
	int  nSensorIndex;
	int  nMotorSpeed;
	bool bExistFrameRange;
	bool bDoRealtimeOperation;
	int  nStartFrameIndex;
	int  nEndFrameIndex;
	INPUT_SOURCE nInputSource;
	std::string strInputDir;
};

struct stParamLidarDetector
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamLidarDetector()
		: dClusteringWindowSize(500)
		, dLineClusteringSize(100)
		, dMaxObjectRadius(2000)
		, cropZone(hj::Rect(0.0, 0.0, 0.0, 0.0))
	{};
	~stParamLidarDetector() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	double dClusteringWindowSize;
	double dLineClusteringSize;
	double dMaxObjectRadius;
	std::vector<CLidarCalibrationInfo*> vecStCalibInfo;
	hj::Rect cropZone;
};

struct stParamLidarTracker
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamLidarTracker()
		: dMaxMovingSpeed(1000.0)
		, nMaxTimeJump(9)
	{};
	~stParamLidarTracker() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	double dMaxMovingSpeed;
	int    nMaxTimeJump;	
};

struct stParamDetect2D
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamDetect2D()
		: nCamIndex(0)
		, strDetectorDir("")
		, strDetectionDir("")
		, pCalibrationInfo(NULL)
		, dImageRescale(1.0)
		, dImageRescaleRecover(1.0)
		, nDetectionThreshold(0)
		, nDetectionType(hj::FULLBODY)
		, dNMSOverlapRatio_1(0.0)
		, dNMSOverlapRatio_2(0.0)
	{};
	~stParamDetect2D() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int nCamIndex;
	stViewInformation stViewInfo;
	std::string strDetectorDir;
	std::string strDetectionDir;
	CCalibrationInfo *pCalibrationInfo;

	/* speed-up */
	double dImageRescale;
	double dImageRescaleRecover;

//	int    nPSN_DET_DETECTOR_TYPE; // 0:Head / 1:Full-body / 2: Full-body from Head
	//int    nPSN_DET_OPERATION_TYPE; // 0:no openmp - acfDetect / 1:openmp - acfDetect /2:TEXT
	//int    nPSN_DET_INPUT_TYPE; // 0:ETRI %d_%04d.jpg / 1:PETS2009 View_%03d\\frame_%04d.jpg / 2:VIDEO / 3: %d_%04d.jpg / others : %d_%d.jpg
	//double nPSN_DET_NMS_OVERLAP_RATIO_1;
	//double nPSN_DET_NMS_OVERLAP_RATIO_2;

	int nDetectionThreshold;
	DETECTION_TYPE nDetectionType;
	double dNMSOverlapRatio_1;
	double dNMSOverlapRatio_2;
};

struct stParamTrack2D
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamTrack2D()
		: nCamIndex(0)
		, pCalibrationInfo(NULL)
		, dImageRescale(1.0)
		, dImageRescaleRecover(1.0)
		, dDetectionMinHeight(1400.0)
		, dDetectionMaxHeight(2300.0)
		, nMaxTrackletLength(5)
		, nMinNumFeatures(4)
		, nMaxNumFeatures(100)
		, nBackTrackingLength(4)
		, dFeatureTrackWindowSizeRatio(1.0)
		, dMaxBoxDistance(1.0)
		, dMinBoxOverlapRatio(0.3)
		, dMaxBoxCenterDiffRatio(0.5)
		, dMinOpticalFlowMajorityRatio(0.5)
		, dMaxDetectionDistance(500.0)
		, dMaxHeightDifference(400.0)
		, bVisualize(false)
	{};
	~stParamTrack2D() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	int nCamIndex;
	CCalibrationInfo *pCalibrationInfo;

	/* speed-up */
	double dImageRescale;
	double dImageRescaleRecover;

	/* detection validation  */
	double dDetectionMinHeight;
	double dDetectionMaxHeight;

	/* bi-directional tracking */
	int    nMaxTrackletLength;
	int    nMinNumFeatures;
	int    nMaxNumFeatures;
	int    nBackTrackingLength;	
	double dFeatureTrackWindowSizeRatio;
	double dMaxBoxDistance;
	double dMinBoxOverlapRatio;
	double dMaxBoxCenterDiffRatio;
	double dMinOpticalFlowMajorityRatio;

	/* matching score related */
	double dMaxDetectionDistance;
	double dMaxHeightDifference;

	/* visualization for debugging */
	bool   bVisualize;	
};

struct stParamAssociate3D
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamAssociate3D()
		: nDetectionType(hj::FULLBODY)
		, bVisualize(false)
		/* optimization */
		, nProcWindowSize(10)
		, nKBestSize(50)
		, nMaxNumTracksInOptimization(1000)
		, nMaxNumTracksInUnconfirmedTree(4)
		, nMaxNumTracksInConfirmedTree(100)
		, nNumFramesForConfirmation(3)
		, bDoBranchCut(false)
		/* reconstruction */
		, nMinTrackletLength(1)
		, dMaxTrackletDistance(1000.0)		
		, bConsiderSensitivity(false)
		, dMaxSensitivityError(20.0)
		, dMinTargetProximity(500.0)
		, dDefaultHeight(1700.0)
		/* clustering */
		, dClusteringDistance(500.0)
		/* linking */
		, dMinLinkingProbability(1.0e-6)
		, nMaxTimeJump(5)
		, dCostRGBMinDistance(0.2)
		, dCostRGBCoef(100.0)
		, dCostRGBDecayCoef(0.1)
		, dCostTrackletLinkingMinDistance(1500.0)
		, dCostTrackletLinkingCoef(0.1)
		/* probability */
		, dMinConstructProbability(1.0e-3)
		, dFPRate(0.01)
		, dFNRate(0.2)
		/* enter/exit */
		, nEnterPenaltyFreeLength(2)
		, dBoundaryDistance(1000.0)
		, dEnterProbMax(1.0e-1)
		, dEnterProbDecayCoef(1.0e-3)
		, dExitProbMax(1.0e-2)		
		, dExitProbDistDecayCoef(1.0e-3)
		, dExitProbLengthDecayCoef(1.0e-1)
		, dCostEnterMax(1000.0)
		, dCostExitMax(1000.0)
		, nMaxOutpoint(3)
		/* calibration */
		, dDetectionError(4.0)
		, dCalibrationError(500.0)
		/* dynamic */
		, dKalmanProcessNoiseSigma(0.1)
		, dKalmanMeasurementNoiseSigma(0.1)
		, dKalmanPostErrorCovariance(0.1)
		, nKalmanConfidenceLevel(1)
		, dVelocityLearningRate(1.0e-1)
		, dFrameRate(6.0)
		, dMaxMovingSpeed(5000.0)
		, dMinMovingSpeed(100.0)
		/* appearance */
		, nImagePatchWidth(20)
		, nNumRGBHistogramBins(16)
		/* operation */
		, nSolverTimeLimit(0)
		, nSolverMaxIter(0)
		/* result */
		, nResultTrajectoryLength(30)
		, bShowTreeID(false)
	{		
	};
	~stParamAssociate3D() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	stViewInformation  stViewInfo;
	stLidarInformation stLidarInfo;
	std::vector<hj::CCalibrationInfo*> vecPCalibrationInfo;

	/* detection related */
	hj::DETECTION_TYPE nDetectionType;

	/* visualization */
	bool   bVisualize;

	/* optimization */
	int    nProcWindowSize;
	int    nKBestSize;
	int    nMaxNumTracksInOptimization;
	int    nMaxNumTracksInUnconfirmedTree;
	int    nMaxNumTracksInConfirmedTree;
	int    nNumFramesForConfirmation;
	bool   bDoBranchCut;

	/* reconstruction */
	int    nMinTrackletLength;
	double dMaxTrackletDistance;	
	bool   bConsiderSensitivity;
	double dMaxSensitivityError;
	double dMinTargetProximity;
	double dDefaultHeight;

	/* clustering */
	double dClusteringDistance;

	/* linking */
	double dMinLinkingProbability;
	int    nMaxTimeJump;
	double dCostRGBMinDistance;
	double dCostRGBCoef;
	double dCostRGBDecayCoef;
	double dCostTrackletLinkingMinDistance;
	double dCostTrackletLinkingCoef;

	/* probability */
	double dMinConstructProbability;
	double dFPRate;
	double dFNRate;

	/* enter/exit */
	int    nEnterPenaltyFreeLength;
	double dBoundaryDistance;
	double dEnterProbMax;
	double dEnterProbDecayCoef;
	double dExitProbMax;	
	double dExitProbDistDecayCoef;
	double dExitProbLengthDecayCoef;
	double dCostEnterMax;
	double dCostExitMax;
	int    nMaxOutpoint;
	
	/* calibration */
	double dDetectionError;
	double dCalibrationError;

	/* dynamic */
	double dKalmanProcessNoiseSigma;
	double dKalmanMeasurementNoiseSigma;
	double dKalmanPostErrorCovariance;
	int    nKalmanConfidenceLevel;
	double dVelocityLearningRate;
	double dFrameRate;
	double dMaxMovingSpeed;
	double dMinMovingSpeed;
	
	/* appearance */
	int    nImagePatchWidth;
	int    nNumRGBHistogramBins;

	/* operation */
	int    nSolverTimeLimit;
	int    nSolverMaxIter;	

	/* result */
	int    nResultTrajectoryLength;
	bool   bShowTreeID;
};

struct stParamEvaluator
{
	//------------------------------------------------
	// METHODS
	//------------------------------------------------
	stParamEvaluator()
		: strGroundTruthPath("")
		, cropZone(0.0, 0.0, 0.0, 0.0)
		, cropZoneMargin(500.0)
	{};
	~stParamEvaluator() {};

	//------------------------------------------------
	// VARIABLES
	//------------------------------------------------
	std::string strGroundTruthPath;
	hj::Rect    cropZone;
	double      cropZoneMargin;
};

}


#endif

//()()
//('')HAANJU.YOO


