#include <numeric>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "HungarianMethod.h"

namespace hj
{

bool IsOne(float value) { return 1.0f == value ? true : false; }

inline void CHungarianMethod::MatSum(std::vector<std::vector<float>> &argInputMat, std::vector<float> &vecDest, unsigned int nDirection)
{
	unsigned int numCols = (unsigned int)argInputMat[0].size();
	unsigned int numRows = (unsigned int)argInputMat.size();
	vecDest.clear();

	switch (nDirection)
	{
	case 1: // column sum
		vecDest.resize(numCols, 0);
		for (unsigned int colIdx = 0; colIdx < numCols; colIdx++)
		{
			for (unsigned int rowIdx = 0; rowIdx < numRows; rowIdx++)
			{
				vecDest[colIdx] += argInputMat[rowIdx][colIdx];
			}
		}
		break;
	case 2: // row sum
		vecDest.resize(numRows, 0);
		for (unsigned int rowIdx = 0; rowIdx < numRows; rowIdx++)
		{
			for (unsigned int colIdx = 0; colIdx < numCols; colIdx++)
			{
				vecDest[rowIdx] += argInputMat[rowIdx][colIdx];
			}
		}
		break;
	default:
		// do nothing
		break;
	}
}

CHungarianMethod::CHungarianMethod()
{
	numRows_ = 0;
	numCols_ = 0;
	bInit_ = false;
}


CHungarianMethod::~CHungarianMethod()
{
}

void CHungarianMethod::Initialize(std::vector<float> costArray, unsigned int numRows, unsigned int numCols) {
	if (!(numRows * numCols)) { return; }
	if (costArray.size() != numRows * numCols) { assert("dimension mismatch!!"); }

	numRows_ = numRows;
	numCols_ = numCols;
	matCost_.clear(); matCost_.resize(numRows_, std::vector<float>(numCols_, 0));
	//matCost_ = cv::Mat(numRows_, numCols_, CV_32F);	
	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++) {
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++) {
			matCost_[rowIdx][colIdx] = costArray[rowIdx*numCols + colIdx];
			if (_isnanf(matCost_[rowIdx][colIdx])) {
				assert("illegal cost!! 1.#IND is contained");
				return;
			}
		}
	}
	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	matchInfo_.costMatrix = matCost_;
	matchInfo_.totalCost = 0;
	bInit_ = true;
}

void CHungarianMethod::Initialize(std::vector<std::vector<float>> costMatrix)
{
	if (!costMatrix.size())
	{
		return;
	}

	numRows_ = (unsigned int)costMatrix.size();
	numCols_ = (unsigned int)costMatrix[0].size();
	matCost_.clear(); matCost_.resize(numRows_, std::vector<float>(numCols_, 0));
	//matCost_ = cv::Mat(numRows_, numCols_, CV_32F);

	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
		{
			matCost_[rowIdx][colIdx] = costMatrix[rowIdx][colIdx];
			if (_isnanf(matCost_[rowIdx][colIdx]))
			{
				assert("illegal cost!! 1.#IND is contained");
				return;
			}
		}
	}

	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	matchInfo_.costMatrix = matCost_;
	matchInfo_.totalCost = 0;

	bInit_ = true;
}

void CHungarianMethod::Initialize(float *costArray, unsigned int numRows, unsigned int numCols)
{
	if (!(numRows * numCols))
	{
		return;
	}

	numRows_ = numRows;
	numCols_ = numCols;
	matCost_.clear(); matCost_.resize(numRows_, std::vector<float>(numCols_, 0));
	//matCost_ = cv::Mat(numRows_, numCols_, CV_32F);

	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
		{
			matCost_[rowIdx][colIdx] = costArray[rowIdx*numCols + colIdx];
			if (_isnanf(matCost_[rowIdx][colIdx]))
			{
				assert("illegal cost!! 1.#IND is contained");
				return;
			}
		}
	}

	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	matchInfo_.costMatrix = matCost_;
	matchInfo_.totalCost = 0;

	bInit_ = true;
}

void CHungarianMethod::Initialize(float **cost2DArray, unsigned int numRows, unsigned int numCols)
{
	if (!(numRows * numCols))
	{
		return;
	}

	numRows_ = numRows;
	numCols_ = numCols;
	matCost_.clear(); matCost_.resize(numRows_, std::vector<float>(numCols_, 0));
	//matCost_ = cv::Mat(numRows_, numCols_, CV_32F);

	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
		{
			matCost_[rowIdx][colIdx] = cost2DArray[rowIdx][colIdx];
			if (_isnanf(matCost_[rowIdx][colIdx]))
			{
				assert("illegal cost!! 1.#IND is contained");
				return;
			}
		}
	}

	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	matchInfo_.costMatrix = matCost_;
	matchInfo_.totalCost = 0;

	bInit_ = true;
}

void CHungarianMethod::Finalize(void)
{
	if (!bInit_)
	{
		return;
	}

	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	for (unsigned int rowIdx = 0; rowIdx < matchInfo_.costMatrix.size(); rowIdx++)
	{
		matchInfo_.costMatrix[rowIdx].clear();
		matchInfo_.matchMatrix[rowIdx].clear();
	}
	matchInfo_.costMatrix.clear();
	matchInfo_.matchMatrix.clear();
	matchInfo_.matchCosts.clear();
	matchInfo_.totalCost = 0.0f;

	bInit_ = false;
}

stMatchInfo* CHungarianMethod::Match()
{
	if (!bInit_)
	{
		return &matchInfo_;
	}

	// handling infinity
	this->CostMatrixPreprocessing();

	// initialize variables
	matMatching_.clear(); matMatching_.resize(numRows_, std::vector<float>(numCols_, 0));

	// find the number in each column and row that are connected
	std::vector<unsigned int> numX;
	std::vector<unsigned int> numY;
	numX.clear(); numX.resize(numRows_, 0);
	numY.clear(); numY.resize(numCols_, 0);
	for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
	{
		for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
		{
			//unsigned int addVal = (unsigned int)_finitef(matCost_[rowIdx][colIdx]);
			unsigned int addVal = (unsigned int)matIsFinite_[rowIdx][colIdx];
			numX[rowIdx] += addVal;
			numY[colIdx] += addVal;
		}
	}

	// find the columns(vertices) and rows(vertices) that are isolated
	std::vector<unsigned int> xCon;
	std::vector<unsigned int> yCon;
	for (unsigned int idx = 0; idx < numRows_; idx++) { if (numX[idx] > 0) { xCon.push_back(idx); } }
	for (unsigned int idx = 0; idx < numCols_; idx++) { if (numY[idx] > 0) { yCon.push_back(idx); } }

	// assemble condensed performance matrix
	nSquareSize_ = std::max(numCols_, numRows_);
	matPCond_.clear(); matPCond_.resize(nSquareSize_, std::vector<float>(nSquareSize_, 0.0f));
	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
		{
			matPCond_[rowIdx][colIdx] = matCost_[rowIdx][colIdx];
		}
	}

	// ensure that a perfect matching exists
	// calculate a form of the Edge matrix
	std::vector<std::vector<float>> matEdge(nSquareSize_, std::vector<float>(nSquareSize_, std::numeric_limits<float>::infinity()));
	float fPMaxVal = 0;
	for (unsigned int rowIdx = 0; rowIdx < nSquareSize_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < nSquareSize_; colIdx++)
		{
			if (_finitef(matPCond_[rowIdx][colIdx]))
			{
				matEdge[rowIdx][colIdx] = 0.0f;
				if (matPCond_[rowIdx][colIdx] > fPMaxVal)
				{
					fPMaxVal = matPCond_[rowIdx][colIdx];
				}
			}
		}
	}

	// find the deficiency(CNUM) in the edge matrix
	unsigned int cNum = this->minLineCover(matEdge);

	// project additional vertices and edges so that a perfect matching exists
	nSquareSize_ += cNum;
	matPCond_.clear(); matPCond_.resize(nSquareSize_, std::vector<float>(nSquareSize_, fPMaxVal));
	for (unsigned int rowIdx = 0; rowIdx < xCon.size(); rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < yCon.size(); colIdx++)
		{
			matPCond_[rowIdx][colIdx] = matCost_[xCon[rowIdx]][yCon[colIdx]];
		}
	}


	///////////////////////////////////////////////////////
	// MAIN PROGRAM: CONTROLS WHICH STEP IS EXECUTED
	///////////////////////////////////////////////////////
	std::vector<std::vector<float>> matM;
	std::vector<float> rCov, cCov, zR, zC;

	bool bExitFlag = false;
	nStepNum_ = 1;
	while (!bExitFlag)
	{
		switch (nStepNum_)
		{
		case 1:
			this->step1();
			break;
		case 2:
			this->step2(rCov, cCov, matM, matPCond_);
			break;
		case 3:
			this->step3(cCov, matM);
			break;
		case 4:
			this->step4(matPCond_, rCov, cCov, matM, zR, zC);
			break;
		case 5:
			this->step5(matM, zR, zC, rCov, cCov);
			break;
		case 6:
			this->step6(matPCond_, rCov, cCov);
			break;
		case 7:
			bExitFlag = true;
			break;
		default:
			// do nothing
			assert("Wrong case!");
			break;
		}
	}

	// remove all the virtual satellites and targets and uncondense
	// the matching to the size the original performance matrix

	// restore infinity
	this->MatchResultPostProcessing();

	// extract mathing information together
	matchInfo_.rows.clear();
	matchInfo_.cols.clear();
	matchInfo_.matchCosts.clear();
	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		for (unsigned int colIdx = 0; colIdx < numCols_; colIdx++)
		{
			matMatching_[rowIdx][colIdx] = matM[rowIdx][colIdx];
			if (1 == matMatching_[rowIdx][colIdx] && _finitef(this->matCost_[rowIdx][colIdx]))
			{
				matchInfo_.rows.push_back(rowIdx);
				matchInfo_.cols.push_back(colIdx);
				matchInfo_.matchCosts.push_back(this->matCost_[rowIdx][colIdx]);
			}
		}
	}
	matchInfo_.totalCost = std::accumulate(matchInfo_.matchCosts.begin(), matchInfo_.matchCosts.end(), 0.0f);
	matchInfo_.matchMatrix = matMatching_;

	return &matchInfo_;
}


void CHungarianMethod::PrintCost(void)
{
	if (!bInit_)
	{
		return;
	}

	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		std::cout << '[';
		for (unsigned int colIdx = 0; colIdx < numCols_ - 1; colIdx++)
		{
			std::cout << matCost_[rowIdx][colIdx] << '\t';
		}
		std::cout << matCost_[rowIdx][numCols_ - 1];
		std::cout << ']' << std::endl;
	}
}


void CHungarianMethod::PrintMatch(void)
{
	if (!bInit_)
	{
		return;
	}

	for (unsigned int rowIdx = 0; rowIdx < numRows_; rowIdx++)
	{
		std::cout << '[';
		for (unsigned int colIdx = 0; colIdx < numCols_ - 1; colIdx++)
		{
			std::cout << matMatching_[rowIdx][colIdx] << '\t';
		}
		std::cout << matMatching_[rowIdx][numCols_ - 1];
		std::cout << ']' << std::endl;
	}
}

/**************************************************************************
*   STEP 1: Find the smallest number of zeros in each row
*           and subtract that minimum from its row
**************************************************************************/
void CHungarianMethod::step1(void) {
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step1" << std::endl;
#endif
	for (unsigned int rowIdx = 0; rowIdx < nSquareSize_; rowIdx++) {
		// find min value
		float minValue = *std::min_element(matPCond_[rowIdx].begin(), matPCond_[rowIdx].end());
		// subtract with min value
		if (0.0f == minValue || std::numeric_limits<float>::infinity() == minValue) { continue; }
		for (unsigned int colIdx = 0; colIdx < nSquareSize_; colIdx++) {
			matPCond_[rowIdx][colIdx] -= minValue;
		}
	}
	nStepNum_ = 2;
}


/**************************************************************************
*   STEP 2: Find a zero in P_cond. If there are no starred zeros in its
*           column or row start the zero. Repeat for each zero
**************************************************************************/
void CHungarianMethod::step2(
	std::vector<float> &rCov,
	std::vector<float> &cCov,
	std::vector<std::vector<float>> &matM,
	std::vector<std::vector<float>> &matPCond)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step2" << std::endl;
#endif
	if (0 == matPCond.size()) { return; }
	unsigned int nSizeMatrix = (unsigned int)matPCond[0].size();

	// define variables	
	rCov.clear(); rCov.resize(nSizeMatrix, 0.0f);
	cCov.clear(); cCov.resize(nSizeMatrix, 0.0f);
	matM.clear(); matM.resize(nSizeMatrix, std::vector<float>(nSizeMatrix, 0.0f));

	for (unsigned int rowIdx = 0; rowIdx < nSizeMatrix; rowIdx++) {
		for (unsigned int colIdx = 0; colIdx < nSizeMatrix; colIdx++) {
			if (0.0f == matPCond[rowIdx][colIdx] && 0.0f == rCov[rowIdx] && 0.0f == cCov[colIdx]) {
				matM[rowIdx][colIdx] = 1.0f;
				rCov[rowIdx] = 1.0f;
				cCov[colIdx] = 1.0f;
			}
		}
	}

	// re-initialize the cover vectors
	rCov.clear(); rCov.resize(nSizeMatrix, 0.0f);
	cCov.clear(); cCov.resize(nSizeMatrix, 0.0f);

	nStepNum_ = 3;
}


/**************************************************************************
*   STEP 3: Cover each column with a starred zero. If all the columns are
*           covered then the matching is maximum
**************************************************************************/
void CHungarianMethod::step3(
	std::vector<float> &cCov,
	std::vector<std::vector<float>> &matM)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step3" << std::endl;
#endif
	if (matM.size() == 0) { assert("M matrix is empty"); }

	MatSum(matM, cCov, 1);
	if (matM[0].size() == std::accumulate(cCov.begin(), cCov.end(), 0.0f)) {
		nStepNum_ = 7;
	}
	else {
		nStepNum_ = 4;
	}
}


/**************************************************************************
*   STEP 4: Find a noncovered zero and prime it.  If there is no starred
*           zero in the row containing this primed zero, Go to Step 5.
*           Otherwise, cover this row and uncover the column containing
*           the starred zero. Continue in this manner until there are no
*           uncovered zeros left. Save the smallest uncovered value and
*           Go to Step 6.
**************************************************************************/
void CHungarianMethod::step4(
	std::vector<std::vector<float>> &matPCond,
	std::vector<float> &rCov,
	std::vector<float> &cCov,
	std::vector<std::vector<float>> &matM,
	std::vector<float> &zR,
	std::vector<float> &zC)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step4" << std::endl;
#endif
	unsigned int nSizeMatrix = (unsigned int)matPCond[0].size();
	bool zFlag = true;
	while (zFlag) {
		// find the first uncovered zero
		int row = -1, col = -1;
		bool bExitFlag = false;
		for (unsigned int rowIdx = 0; rowIdx < nSizeMatrix; rowIdx++) {
			for (unsigned int colIdx = 0; colIdx < nSizeMatrix; colIdx++) {
				if (0.0f == matPCond[rowIdx][colIdx] && 0.0f == rCov[rowIdx] && 0.0f == cCov[colIdx]) {
					row = (int)rowIdx;
					col = (int)colIdx;
					bExitFlag = true;
					break;
				}
			}
			if (bExitFlag) { break; }
		}

		if (-1 == row) {
			// if there are no uncovered zeros go to step 6
			nStepNum_ = 6;
			zFlag = false;
			zR.clear(); zR.resize(1, 0.0f);
			zC.clear(); zC.resize(1, 0.0f);
		}
		else {
			// prime the uncovered zero
			matM[row][col] = 2.0f;

			// if there is a starred zero in that row, cover the row and uncover the column containing the zero
			if (0 < std::count_if(matM[row].begin(), matM[row].end(), IsOne)) {
				rCov[row] = 1.0f;
				for (unsigned int colIdx = 0; colIdx < nSizeMatrix; colIdx++) {
					if (1.0f == matM[row][colIdx]) { cCov[colIdx] = 0.0f; }
				}
			}
			else {
				nStepNum_ = 5;
				zFlag = false;
				zR.clear(); zR.resize(1, (float)row);
				zC.clear(); zC.resize(1, (float)col);
			}
		}
	}
}


/**************************************************************************
* STEP 5: Construct a series of alternating primed and starred zeros as
*         follows.  Let Z0 represent the uncovered primed zero found in Step 4.
*         Let Z1 denote the starred zero in the column of Z0 (if any).
*         Let Z2 denote the primed zero in the row of Z1 (there will always
*         be one).  Continue until the series terminates at a primed zero
*         that has no starred zero in its column.  Unstar each starred
*         zero of the series, star each primed zero of the series, erase
*         all primes and uncover every line in the matrix.  Return to Step 3.
**************************************************************************/
void CHungarianMethod::step5(
	std::vector<std::vector<float>> &matM,
	std::vector<float> &zR,
	std::vector<float> &zC,
	std::vector<float> &rCov,
	std::vector<float> &cCov)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step5" << std::endl;
#endif
	unsigned int nMatrixSize = (unsigned int)matM[0].size();
	bool zFlag = true;

	unsigned int row = 0;
	while (zFlag) {
		// find the index number of the starred zero in the column
		int rIndex = -1;
		for (unsigned int rowIdx = 0; rowIdx < nMatrixSize; rowIdx++) {
			if (1.0f == matM[rowIdx][(unsigned int)zC[row]]) {
				rIndex = (int)rowIdx;
				break;
			}
		}

		if (rIndex >= 0) {
			// save the starred zero
			row++;
			// save the row of the starred zero
			zR.push_back((float)rIndex);
			// the column of the starred zero is the same as the column of the primed zero
			zC.push_back(zC[row - 1]);
		}
		else {
			zFlag = false;
		}

		// continue if there is a starred zero in the column of the primed zero
		if (zFlag) {
			// find the column of the primed zero in the last starred zeros row
			unsigned int cIndex = 0;
			for (unsigned int colIdx = 0; colIdx < nMatrixSize; colIdx++) {
				if (2.0f == matM[(unsigned int)zR[row]][colIdx]) {
					cIndex = colIdx;
					break;
				}
			}
			row++;
			zR.push_back(zR[row - 1]);
			zC.push_back((float)cIndex);
		}
	}

	// UNSTAR all the starred zeros in the path and STAR all primed zeros
	for (unsigned int rowIdx = 0; rowIdx < zR.size(); rowIdx++) {
		if (1.0f == matM[(unsigned int)zR[rowIdx]][(unsigned int)zC[rowIdx]]) {
			matM[(unsigned int)zR[rowIdx]][(unsigned int)zC[rowIdx]] = 0.0f;
		}
		else {
			matM[(unsigned int)zR[rowIdx]][(unsigned int)zC[rowIdx]] = 1.0f;
		}
	}

	// clear the covers
	rCov.clear(); rCov.resize(nMatrixSize, 0.0f);
	cCov.clear(); cCov.resize(nMatrixSize, 0.0f);

	// remove all the primes
	for (unsigned int rowIdx = 0; rowIdx < nMatrixSize; rowIdx++) {
		for (unsigned int colIdx = 0; colIdx < nMatrixSize; colIdx++) {
			if (2 == matM[rowIdx][colIdx]) {
				matM[rowIdx][colIdx] = 0.0f;
			}
		}
	}

	nStepNum_ = 3;
}


/**************************************************************************
* STEP 6: Add the minimum uncovered value to every element of each covered
*         row, and subtract it from every element of each uncovered column.
*         Return to Step 4 without altering any stars, primes, or covered lines.
**************************************************************************/
void CHungarianMethod::step6(
	std::vector<std::vector<float>> &matPCond,
	std::vector<float> &rCov,
	std::vector<float> &cCov)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do step6" << std::endl;
#endif
	float fMinVal = std::numeric_limits<float>::infinity();
	for (unsigned int rowIdx = 0; rowIdx < nSquareSize_; rowIdx++)
	{
		if (0 != rCov[rowIdx]) { continue; }
		for (unsigned int colIdx = 0; colIdx < nSquareSize_; colIdx++)
		{
			if (0 != cCov[colIdx]) { continue; }
			if (matPCond[rowIdx][colIdx] < fMinVal)
			{
				fMinVal = matPCond[rowIdx][colIdx];
			}
		}
	}

	float rowBiasVal = 0.0f, colBiasVal = 0.0f;
	for (unsigned int rowIdx = 0; rowIdx < nSquareSize_; rowIdx++)
	{
		rowBiasVal = 0.0f;
		if (1 == rCov[rowIdx]) { rowBiasVal = fMinVal; }
		for (unsigned int colIdx = 0; colIdx < nSquareSize_; colIdx++)
		{
			colBiasVal = 0.0f;
			if (0 == cCov[colIdx]) { colBiasVal = -fMinVal; }
			matPCond[rowIdx][colIdx] += rowBiasVal + colBiasVal;
		}
	}

	nStepNum_ = 4;
}

unsigned int CHungarianMethod::minLineCover(std::vector<std::vector<float>> &argMatrix)
{
#ifdef PSN_HUNGARIAN_DEBUG_DISPLAY
	std::cout << "do minLineCover" << std::endl;
#endif
	unsigned int nMatrixSize = (unsigned int)argMatrix[0].size();
	if (0 == nMatrixSize)
	{
		assert("minLineCover: empty argMatrix");
		return 0;
	}

	std::vector<std::vector<float>> matM;
	std::vector<float> rCov;
	std::vector<float> cCov;
	std::vector<float> zR;
	std::vector<float> zC;
	unsigned int cNum = nMatrixSize;

	// step 2
	this->step2(rCov, cCov, matM, argMatrix);

	// step 3
	this->step3(cCov, matM);

	// step 4
	this->step4(argMatrix, rCov, cCov, matM, zR, zC);

	// calculate the deficiency
	cNum -= (unsigned int)(std::accumulate(rCov.begin(), rCov.end(), 0.0f) + std::accumulate(cCov.begin(), cCov.end(), 0.0f));

	return cNum;
}

void CHungarianMethod::CostMatrixPreprocessing(void)
{
	float finiteCostSum = 0.0f;
	matIsFinite_.resize(matCost_.size());
	for (int rowIdx = 0; rowIdx < matCost_.size(); rowIdx++)
	{
		matIsFinite_[rowIdx].resize(matCost_[rowIdx].size(), false);
		for (int colIdx = 0; colIdx < matCost_[rowIdx].size(); colIdx++)
		{
			if (!_finitef(matCost_[rowIdx][colIdx])) { continue; }
			matIsFinite_[rowIdx][colIdx] = true;
			finiteCostSum += matCost_[rowIdx][colIdx];
		}
	}
	fInfiniteReplacementCost_ = FLT_MAX - finiteCostSum;

	for (int rowIdx = 0; rowIdx < matCost_.size(); rowIdx++) {

		for (int colIdx = 0; colIdx < matCost_[rowIdx].size(); colIdx++)
		{
			if (matIsFinite_[rowIdx][colIdx]) { continue; }
			matCost_[rowIdx][colIdx] = fInfiniteReplacementCost_;
		}
	}
}

void CHungarianMethod::MatchResultPostProcessing(void)
{
	for (int rowIdx = 0; rowIdx < matCost_.size(); rowIdx++)
	{
		for (int colIdx = 0; colIdx < matCost_[rowIdx].size(); colIdx++)
		{
			if (matCost_[rowIdx][colIdx] != fInfiniteReplacementCost_) { continue; }
			matCost_[rowIdx][colIdx] = std::numeric_limits<float>::infinity();
		}
	}
}

}


//()()
//('')HAANJU.YOO

