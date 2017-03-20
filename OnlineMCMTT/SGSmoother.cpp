#include "SGSmoother.h"

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>

namespace hj
{
	
CSGSmoother::CSGSmoother(void)
	: span_(SGS_DEFAULT_SPAN)
	, degree_(SGS_DEFAULT_DEGREE)
	, precomputedQsets_(NULL)
{
	Qset_.cols = 0;
	Qset_.rows = 0;
}

CSGSmoother::CSGSmoother(int span, int degree, std::vector<double> *initialData)
	: span_(span)
	, degree_(degree)
	, precomputedQsets_(NULL)
{
	Qset_.cols = 0;
	Qset_.rows = 0;
	if (NULL == initialData) { return; }
	data_.insert(data_.begin(), initialData->begin(), initialData->end());
}

CSGSmoother::~CSGSmoother(void)
{
	data_.clear();
	smoothedData_.clear();
}

int CSGSmoother::Reset(std::vector<double> &queueNewData)
{
	data_.clear();
	smoothedData_.clear();
	return Insert(queueNewData);
}

int CSGSmoother::Insert(double newData)
{
	data_.push_back(newData);
	return Smoothing();
}

int CSGSmoother::Insert(std::vector<double> &queueNewData)
{
	data_.insert(data_.end(), queueNewData.begin(), queueNewData.end());
	return Smoothing();
}

int CSGSmoother::ReplaceBack(double replaceData)
{
	data_.back() = replaceData;
	smoothedData_.pop_back();
	return Smoothing();
}

void CSGSmoother::PopBack(void)
{
	data_.pop_back();
	smoothedData_.pop_back();
	size_ = smoothedData_.size();
}

void CSGSmoother::SetSmoother(std::deque<double> &data, std::deque<double> &smoothedData, int span, int degree)
{
	data_ = data;
	smoothedData_ = smoothedData;
	span_ = span;
	degree_ = degree;
	size_ = smoothedData_.size();
}

double CSGSmoother::GetResult(int pos)
{
	assert(pos < smoothedData_.size());
	return smoothedData_[pos];
}

std::vector<double> CSGSmoother::GetResults(int startPos, int endPos)
{
	assert(endPos < smoothedData_.size() && startPos <= endPos);
	int numPoints = endPos - startPos + 1;
	std::vector<double> results;
	results.reserve(numPoints);
	for (int pos = startPos; pos <= endPos; pos++)
	{
		results.push_back(smoothedData_[pos]);
	}
	return results;
}

void CSGSmoother::GetSmoother(std::deque<double> &data, std::deque<double> &smoothedData, int &span, int &degree)
{
	data = data_;
	smoothedData = smoothedData_;
	span = span_;
	degree = degree_;
}

bool CSGSmoother::UpdateQ(int windowSize)
{
	if (Qset_.rows == windowSize) { return false; }
	if (NULL != precomputedQsets_ && precomputedQsets_->size() >= windowSize)
	{
		Qset_ = (*precomputedQsets_)[windowSize - 1];
	}
	else
	{
		Qset_ = CalculateQ(windowSize, degree_);
	}
	return true;
}

Qset CSGSmoother::CalculateQ(int windowSize, int degree)
{
	Qset resutQ_;
	int halfWindowSize = (windowSize - 1) / 2;

	// Qmid
	resutQ_.Qmid.clear();
	resutQ_.Qmid.resize(windowSize, 1.0 / (double)windowSize);

	// direction => row
	int pos = 0;
	std::vector<double> V(windowSize * (degree + 1), 1.0);
	for (int order = 1; order <= degree; order++)
	{
		for (int timeIdx = -halfWindowSize, pos = order; timeIdx <= halfWindowSize; timeIdx++, pos += degree + 1)
		{
			V[pos] = std::pow((double)timeIdx, (double)order);
		}
	}

	// find Q (recuded QR factorization, we don't need R)
	resutQ_.Q.clear();
	resutQ_.rows = windowSize;
	resutQ_.cols = degree + 1;
	resutQ_.Q.resize(resutQ_.rows* resutQ_.cols, 0.0);

	int preColPos = 0;
	double columnNorm = 0.0;
	for (int colIdx = 0; colIdx < resutQ_.cols; colIdx++)
	{
		std::vector<double> projectionsMagnitude(std::max(0, colIdx), 0.0);
		for (int rowIdx = 0, pos = colIdx; rowIdx < resutQ_.rows; rowIdx++, pos += resutQ_.cols)
		{
			resutQ_.Q[pos] = V[pos];
			// for Q(preCol)' * V
			for (int preColIdx = 0, preColPos = pos - colIdx; preColIdx < colIdx; preColIdx++, preColPos++)
			{
				projectionsMagnitude[preColIdx] += resutQ_.Q[preColPos] * V[pos];
			}
		}

		// orthogonalization
		columnNorm = 0.0;
		for (int rowIdx = 0, pos = colIdx; rowIdx < resutQ_.rows; rowIdx++, pos += resutQ_.cols)
		{
			preColPos = pos - colIdx;
			for (int preColIdx = 0; preColIdx < colIdx; preColIdx++, preColPos++)
			{
				resutQ_.Q[pos] -= projectionsMagnitude[preColIdx] * resutQ_.Q[preColPos];
			}
			columnNorm += resutQ_.Q[pos] * resutQ_.Q[pos];
		}
		columnNorm = std::sqrt(columnNorm);

		// normalization
		for (int rowIdx = 0, pos = colIdx; rowIdx < resutQ_.rows; rowIdx++, pos += resutQ_.cols)
		{
			resutQ_.Q[pos] /= columnNorm;
		}
	}

	// Qbegin: q(1:hf,:)*q'
	resutQ_.Qbegin.clear();
	resutQ_.Qbegin.resize(halfWindowSize * resutQ_.rows, 0.0);

	// Qend: q((hf+2):end,:)*q'
	resutQ_.Qend.clear();
	resutQ_.Qend.resize(halfWindowSize * resutQ_.rows, 0.0);

	int posFrontHalf = 0;
	int posBackHalf = (halfWindowSize + 1) * resutQ_.cols;
	for (int rowIdx = 0, pos = 0; rowIdx < halfWindowSize; rowIdx++)
	{
		int posQ = 0;
		for (int colIdx = 0; colIdx < resutQ_.rows; colIdx++, pos++)
		{
			for (int elementIdx = 0; elementIdx < resutQ_.cols; elementIdx++, posQ++)
			{
				resutQ_.Qbegin[pos] += resutQ_.Q[posQ] * resutQ_.Q[posFrontHalf + elementIdx];
				resutQ_.Qend[pos] += resutQ_.Q[posQ] * resutQ_.Q[posBackHalf + elementIdx];
			}
		}
		posFrontHalf += resutQ_.cols;
		posBackHalf += resutQ_.cols;
	}

	return resutQ_;
}

int CSGSmoother::Smoothing(void)
{
	int refreshPos = 0;
	int numData = (int)data_.size();
	int windowSize = std::min(span_, numData);
	windowSize -= (windowSize + 1) % 2; // will subtract 1 if frame is even.
	if (windowSize <= degree_) // bypass
	{
		refreshPos = (int)smoothedData_.size();
		for (int idx = (int)smoothedData_.size(); idx < data_.size(); idx++)
		{
			smoothedData_.push_back(data_[idx]);
		}
		return refreshPos;
	}

	// smoothing
	int halfWindowSize = (windowSize - 1) / 2;
	int midStartPos = 0;
	bool bEntireUpdate = UpdateQ(windowSize);
	std::vector<double> smoothedMid;
	if (bEntireUpdate)
	{
		// begin
		smoothedData_.clear();
		smoothedData_.resize(halfWindowSize, 0.0);
		for (int pos = 0, smoothDataPos = 0; smoothDataPos < halfWindowSize; smoothDataPos++)
		{
			for (int colIdx = 0; colIdx < windowSize; colIdx++, pos++)
			{
				smoothedData_[smoothDataPos] += Qset_.Qbegin[pos] * data_[colIdx];
			}
		}

		// middle
		smoothedMid = Filter(Qset_.Qmid, data_);
		midStartPos = windowSize - 1;
	}
	else
	{
		// middle
		refreshPos = (int)smoothedData_.size() - halfWindowSize;
		smoothedMid = Filter(Qset_.Qmid, data_, (int)smoothedData_.size());
		smoothedData_.erase(smoothedData_.begin() + refreshPos, smoothedData_.end());
	}

	// middle
	smoothedData_.insert(smoothedData_.end(), smoothedMid.begin() + midStartPos, smoothedMid.end());

	// end
	for (int pos = 0, smoothDataPos = (int)data_.size() - halfWindowSize; smoothDataPos < data_.size(); smoothDataPos++)
	{
		double curSmoothedData = 0.0;
		for (int colIdx = (int)data_.size() - windowSize; colIdx < (int)data_.size(); colIdx++, pos++)
		{
			curSmoothedData += Qset_.Qend[pos] * data_[colIdx];
		}
		smoothedData_.push_back(curSmoothedData);
	}
	size_ = smoothedData_.size();

	return refreshPos;
}

std::vector<double> CSGSmoother::Filter(std::vector<double> &coefficients, std::deque<double> &data, int startPos)
{
	std::vector<double> results(data.size() - startPos, 0.0);
	for (int resultPos = 0; resultPos < data.size() - startPos; resultPos++)
	{
		int dataPos = startPos + resultPos;
		for (int coeffPos = 0; coeffPos < coefficients.size() && 0 <= dataPos; coeffPos++, dataPos--)
		{
			results[resultPos] += coefficients[coeffPos] * data[dataPos];
		}
	}
	return results;
}

}

//()()
//('')HAANJU.YOO

