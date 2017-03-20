#include "haanju_matOperations.hpp"

namespace hj
{

/************************************************************************
 Method Name: appendRow
 Description:
	- append additional row into cv::Mat
 Input Arguments:
	- dstMat: the matrix to which the append row
	- row: appending row
 Return Values:
	- void
************************************************************************/
void appendRow(cv::Mat &dstMat, cv::Mat row)
{
	// append with empty matrix
	if (dstMat.empty())
	{
		dstMat = row.clone();
		return;
	}

	// normal appending
	if (dstMat.cols == row.cols)
	{
		cv::vconcat(dstMat, row, dstMat);
		return;
	}

	// expand dstMat
	if (dstMat.cols < row.cols)
	{
		cv::Mat newMat(dstMat.rows + 1, row.cols, dstMat.type());
		newMat(cv::Rect(0, 0, dstMat.rows, dstMat.cols)) = dstMat.clone();
		newMat.row(dstMat.rows) = row.clone();
		dstMat.release();
		dstMat = newMat.clone();
	}

	// expand row
	cv::Mat newRow(1, dstMat.cols, row.type());
	newRow(cv::Rect(0, 0, 1, row.cols)) = row.clone();
	cv::vconcat(dstMat, row, newRow);
}


/************************************************************************
 Method Name: appendCol
 Description:
	- append additional column into cv::Mat
 Input Arguments:
	- dstMat: the matrix to which the append column
	- col: appending column
 Return Values:
	- void
************************************************************************/
void appendCol(cv::Mat &dstMat, cv::Mat col)
{
	// append with empty matrix
	if (dstMat.empty())
	{
		dstMat = col.clone();
		return;
	}

	// normal appending
	if (dstMat.rows == col.rows)
	{
		cv::hconcat(dstMat, col, dstMat);
		return;
	}

	// expand dstMat
	if (dstMat.rows < col.rows)
	{
		cv::Mat newMat(col.rows, dstMat.cols + 1, dstMat.type());
		newMat(cv::Rect(0, 0, dstMat.rows, dstMat.cols)) = dstMat.clone();
		newMat.col(dstMat.cols) = col.clone();
		dstMat.release();
		dstMat = newMat.clone();
	}

	// expand row
	cv::Mat newRow(dstMat.rows, 1, col.type());
	newRow(cv::Rect(0, 0, col.rows, 1)) = col.clone();
	cv::hconcat(dstMat, col, newRow);
}

}


//()()
//('')HAANJU.YOO


