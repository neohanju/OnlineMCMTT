/******************************************************************************
* Title        : CCameraModel class implementation
* Author       : (C)Silogic (initial author) / Haanju Yoo (modified for compact ver.)
* Initial Date : 2014.08.31 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.06)
* Description  : 
*	- Class for camera calibration information and its operations
*	- Originally, this program is part of the ETISEO project.
*	- See http://www.etiseo.net  http://www.silogic.fr
*	- Modified to the compact version, that only read XML (not write), by Haanju Yoo.
******************************************************************************/

#ifndef _CAMERA_MODEL_OBJECT_H_
#define _CAMERA_MODEL_OBJECT_H_

#include <string>

#define XML_TAG_CAMERA      "Camera"
#define XML_TAG_NAME        "name"

#define XML_TAG_GEOMETRY    "Geometry"
#define XML_TAG_WIDTH       "width"
#define XML_TAG_HEIGHT      "height"
#define XML_TAG_NCX         "ncx"
#define XML_TAG_NFX         "nfx"
#define XML_TAG_DX          "dx"
#define XML_TAG_DY          "dy"
#define XML_TAG_DPX         "dpx"
#define XML_TAG_DPY         "dpy"

#define XML_TAG_INTRINSIC   "Intrinsic"
#define XML_TAG_FOCAL       "focal" 
#define XML_TAG_KAPPA1      "kappa1" 
#define XML_TAG_CX          "cx" 
#define XML_TAG_CY          "cy" 
#define XML_TAG_SX          "sx" 

#define XML_TAG_EXTRINSIC   "Extrinsic"
#define XML_TAG_TX          "tx" 
#define XML_TAG_TY          "ty" 
#define XML_TAG_TZ          "tz" 
#define XML_TAG_RX          "rx" 
#define XML_TAG_RY          "ry" 
#define XML_TAG_RZ          "rz" 

namespace Etiseo {

	//!  A root class handling camera model
	class CameraModel
	{
	public:

		//! Constructor
		CameraModel();
		//! Destructor
		virtual ~CameraModel();

		/////////////////////////////////////////////////////////////////////////
		// GETTER
		/////////////////////////////////////////////////////////////////////////
		
		std::string name() { return mName; }
		const std::string& name() const { return mName; }		

		inline int    width()  const { return mImgWidth; }
		inline int    height() const { return mImgHeight; }
		inline double ncx()    const { return mNcx; }
		inline double nfx()    const { return mNfx; }
		inline double dx()     const { return mDx; }
		inline double dy()     const { return mDy; }
		inline double dpx()    const { return mDpx; }
		inline double dpy()    const { return mDpy; }
		inline double cx()     const { return mCx; }
		inline double cy()     const { return mCy; }
		inline double sx()     const { return mSx; }
		inline double focal()  const { return mFocal; }
		inline double kappa1() const { return mKappa1; }
		inline double tx()     const { return mTx; }
		inline double ty()     const { return mTy; }
		inline double tz()     const { return mTz; }
		inline double rx()     const { return mRx; }
		inline double ry()     const { return mRy; }
		inline double rz()     const { return mRz; }
		inline double cposx()  const { return mCposx; }
		inline double cposy()  const { return mCposy; }
		inline double cposz()  const { return mCposz; }

		/////////////////////////////////////////////////////////////////////////
		// SETTER
		/////////////////////////////////////////////////////////////////////////
		
		void setName(const std::string& name) { mName = name; }

		void setGeometry(int width, int height, double ncx, double nfx,	double dx, double dy, double dpx, double dpy);
		void setIntrinsic(double focal, double kappa1, double cx, double cy, double sx);
		void setExtrinsic(double tx, double ty, double tz, double rx, double ry, double rz);
		
		/////////////////////////////////////////////////////////////////////////
		// FILE INTERFACE
		/////////////////////////////////////////////////////////////////////////

		// loading from an XML
		virtual bool fromXml(std::string _filePath);

		// loading from a .dat files = output of the Tsai calibration 
		virtual bool fromDat(std::istream& is, int width, int height);

		
		/////////////////////////////////////////////////////////////////////////
		// COORDINATE MANIPULATION
		/////////////////////////////////////////////////////////////////////////

		// from image coordinate to world coordinate (2D -> 3D)
		bool imageToWorld(double Xi, double Yi, double Zw, double& Xw, double &Yw);
		// from world coordinate to image coordinate (3D -> 2D)
		bool worldToImage(double Xw, double Yw, double Zw, double& Xi, double& Yi);

		// convert from undistorted to distorted image (2D -> 2D)
		bool undistortedToDistortedImageCoord(double Xfu, double Yfu, double& Xfd, double& Yfd);
		// convert from distorted to undistorted image (2D -> 2D)
		bool distortedToUndistortedImageCoord(double Xfd, double Yfd, double& Xfu, double& Yfu);

		// from world coordinate to camera coordinate (3D -> 3D)
		bool worldToCameraCoord(double xw, double yw, double zw, double& xc, double& yc, double& zc);
		// from camera coordinate to world coordinate (3D -> 3D)
		bool cameraToWorldCoord(double xc, double yc, double zc, double& xw, double& yw, double& zw);

	protected:

		virtual void internalInit();

		// Coordinate manipulation, intermediate transformation :
		// convert from distorted to undistorted sensor plane coordinates 
		void distortedToUndistortedSensorCoord(double Xd, double Yd, double& Xu, double& Yu);
		// convert from undistorted to distorted sensor plane coordinates
		void undistortedToDistortedSensorCoord(double Xu, double Yu, double& Xd, double& Yd);

	private:

		bool			isInit;
		std::string		mName;

		// geometry
		int				mImgWidth;
		int				mImgHeight;
		double			mNcx;
		double			mNfx;
		double			mDx;
		double			mDy;
		double			mDpx;
		double			mDpy;

		// intrinsic 
		double			mFocal;
		double			mKappa1;
		double			mCx;
		double			mCy;
		double			mSx;

		// extrinsic 
		double			mTx;
		double			mTy;
		double			mTz;
		double			mRx;
		double			mRy;
		double			mRz;

		// for computation
		double			mR11;
		double			mR12;
		double			mR13;
		double			mR21;
		double			mR22;
		double			mR23;
		double			mR31;
		double			mR32;
		double			mR33;

		//camera position
		double			mCposx;
		double			mCposy;
		double			mCposz;

	};
};

#endif
