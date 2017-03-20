/******************************************************************************
* File name: GraphSolver.cpp
******************************************************************************/

#include "GraphSolver.h"
#include <algorithm>
#include <time.h>

// TODO: 티끌이라도 더 빠르게 하려면, delete(v) -> insert(u) 연계 시, OM 진입하는 OC는 u의 neighbor만 있음을 이용
//#define HJ_GRAPH_TYPE_MWISP
#define HJ_GRAPH_SOLVE_MONITORING

// comparator
bool HJGraphComparatorVertexWeightDescend(const hj::GraphVertex *vertex1, const hj::GraphVertex *vertex2)
{
	return vertex1->weight_ > vertex2->weight_;
}

bool HJGraphComparatorVertexDegreeDescend(const hj::GraphVertex *vertex1, const hj::GraphVertex *vertex2)
{
	return vertex1->queueNeighbor_.size() > vertex2->queueNeighbor_.size();
}

bool HJGraphComparatorVertexIDAscend(const hj::GraphVertex *vertex1, const hj::GraphVertex *vertex2)
{
	return vertex1->id_ < vertex2->id_;
}

bool HJSolutionComparatorSumWeightsDescend(const std::pair<hj::VertexSet, double> &solution1, const std::pair<hj::VertexSet, double> &solution2)
{
	return solution1.second > solution2.second;
}

namespace hj {

/////////////////////////////////////////////////////////////////////////
// Mat2D
/////////////////////////////////////////////////////////////////////////
Mat2D::Mat2D(void)
{
	rows_ = 0;
	cols_ = 0;
}

Mat2D::Mat2D(int rows, int cols)
{
	rows_ = rows;
	cols_ = cols;
	data_.resize(rows_);
	for (int rowIdx = 0; rowIdx < data_.size(); rowIdx++)
	{
		data_[rowIdx].resize(cols, 0.0);
	}
}

Mat2D::~Mat2D(void)
{
	this->clear();
}

void Mat2D::clear(void)
{
	for (int rowIdx = 0; rowIdx < this->data_.size(); rowIdx++)
	{
		this->data_[rowIdx].clear();
	}
	this->data_.clear();
	rows_ = 0;
	cols_ = 0;
}

void Mat2D::resize(int rows, int cols)
{
	this->clear();
	rows_ = rows;
	cols_ = cols;
	data_.resize(rows);
	for (int rowIdx = 0; rowIdx < rows_; rowIdx++)
	{
		data_[rowIdx].resize(cols_, 0.0);
	}
}

double& Mat2D::at(int row, int col)
{
	return data_[row][col];
}

//Mat2D Mat2D::operator=(const Mat2D &x)
//{
//	Mat2D newMat(x.rows(), x.cols());
//	for (int rowIdx = 0; rowIdx < x.rows_; rowIdx++)
//	{
//		for (int colIdx = 0; colIdx < x.cols(); colIdx++)
//		{
//			newMat.at(rowIdx, colIdx) = x.data_[rowIdx][colIdx];
//		}
//	}
//	return newMat;
//}

/////////////////////////////////////////////////////////////////////////
// GraphVertex
/////////////////////////////////////////////////////////////////////////
GraphVertex::GraphVertex(void)
{
	Initialize(0, 0.0);
}

GraphVertex::GraphVertex(size_t nID)
{
	Initialize(nID, 0.0);
}

GraphVertex::GraphVertex(size_t nID, double weight)
{
	Initialize(nID, weight);
}

void GraphVertex::CountUp(void)
{
	bInSolution_ = true;
	bInOM_ = false;
	countNeighborInSolution_++;
	for (size_t neighborIdx = 0; neighborIdx < queueNeighbor_.size(); neighborIdx++)
	{
		queueNeighbor_[neighborIdx]->countNeighborInSolution_++;
	}
}

void GraphVertex::CountDown(void)
{
	bInSolution_ = false;
	countNeighborInSolution_--;
	for (size_t neighborIdx = 0; neighborIdx < queueNeighbor_.size(); neighborIdx++)
	{
		queueNeighbor_[neighborIdx]->countNeighborInSolution_--;
	}
}

void GraphVertex::Initialize(size_t nID, double weight)
{
	id_ = nID;
	bValid_ = true;
	weight_ = weight;

	// BLS related
	bInOM_ = false;
	bInSolution_ = false;
	tabuStamp_ = 0;
	countNeighborInSolution_ = 0;;
}

/////////////////////////////////////////////////////////////////////////
// Graph
/////////////////////////////////////////////////////////////////////////
/* [[NOTICE]]: This graph data class has lazy update policy at deletion */
Graph::Graph(void)
	: nNumVertex_(0)
	, nNumEdge_(0)
	, nNewID_(0)
{
}

Graph::Graph(size_t nNumVertex)
	: nNumVertex_(nNumVertex)
	, nNumEdge_(0)
	, nNewID_(nNumVertex)
{
	for (size_t vertexIdx = 0; vertexIdx < nNumVertex; vertexIdx++)
	{
		listVertices_.push_back(GraphVertex(vertexIdx));
		vecPtVertices_.push_back(&listVertices_.back());
	}
}

Graph::~Graph(void)
{
}

bool Graph::Clear()
{
	try
	{
		listVertices_.clear();
		vecPtVertices_.clear();
		nNumVertex_ = 0;
		nNumEdge_ = 0;
		nNewID_ = 0;
		return true;
	}
	catch (int nError)
	{
		printf("[ERROR] Error is occured at Clear!!: %d\n", nError);
	}
	return false;
}

bool Graph::TopologyModified()
{
	return bTopologyModified_;
}

size_t Graph::maxDegree(void)
{
	size_t maxDegree = 0;
	for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
		vertexIter != listVertices_.end();
		vertexIter++)
	{
		maxDegree = std::max(maxDegree, vertexIter->degree());
	}
	return maxDegree;
}

size_t Graph::minDegree(void)
{
	size_t minDegree = 0;
	for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
		vertexIter != listVertices_.end();
		vertexIter++)
	{
		minDegree = std::min(minDegree, vertexIter->degree());
	}
	return minDegree;
}

double Graph::AverageVertexDegree(void)
{
	double averageDegree = 0;
	for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
		vertexIter != listVertices_.end();
		vertexIter++)
	{
		averageDegree += (double)vertexIter->degree();
	}
	averageDegree /= (double)listVertices_.size();
	return averageDegree;
}

GraphVertex* Graph::AddVertex(double weight)
{
	listVertices_.push_back(GraphVertex(nNewID_++, weight));
	vecPtVertices_.push_back(&listVertices_.back());
	nNumVertex_++;
	bTopologyModified_ = true;
	return vecPtVertices_.back();
}

VertexSet Graph::AddVertices(size_t numVertex)
{
	VertexSet queueAddedVertices(numVertex);
	for (size_t vIdx = 0; vIdx < numVertex; vIdx++)
	{
		listVertices_.push_back(GraphVertex(nNewID_++, 0.0));
		vecPtVertices_.push_back(&listVertices_.back());
		nNumVertex_++;
		queueAddedVertices[vIdx] = vecPtVertices_.back();
	}
	bTopologyModified_ = true;
	return queueAddedVertices;
}


bool Graph::DeleteVertex(GraphVertex* vertex)
{
	vertex->bValid_ = false;
	return true;
}

bool Graph::AddEdge(GraphVertex* vertex1, GraphVertex* vertex2)
{
	try
	{
		if (vertex1->bValid_ && vertex2->bValid_)
		{
			vertex1->queueNeighbor_.push_back(vertex2);
			vertex2->queueNeighbor_.push_back(vertex1);
			bTopologyModified_ = true;
			nNumEdge_++;
			return true;
		}
	}
	catch (int nError)
	{
		printf("[ERROR] Error is occured at AddEdge!!: %d\n", nError);
	}
	return false;
}

bool Graph::Update(void)
{
	try
	{
		bTopologyModified_ = true;

		// delete invalid verteces from the neighbor list		
		for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
			vertexIter != listVertices_.end();
			vertexIter++)
		{
			std::deque<GraphVertex*> newNeighbors;
			for (size_t neighborIdx = 0; neighborIdx < (*vertexIter).queueNeighbor_.size(); neighborIdx++)
			{
				if ((*vertexIter).queueNeighbor_[neighborIdx]->bValid_)
				{
					newNeighbors.push_back((*vertexIter).queueNeighbor_[neighborIdx]);
				}
			}
			(*vertexIter).queueNeighbor_ = newNeighbors;
		}

		// delete invalid verteces
		int vertexIdx = 0;
		for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
			vertexIter != listVertices_.end();
			/* do in the loop */)
		{
			if ((*vertexIter).bValid_)
			{
				vertexIter++;
				vertexIdx++;
				continue;
			}
			listVertices_.erase(vertexIter);
			vecPtVertices_.erase(vecPtVertices_.begin() + vertexIdx);
			nNumVertex_--;
		}
	}
	catch (int nError)
	{
		printf("[ERROR] Error is occured at Update!!: %d\n", nError);
		return false;
	}

	return true;
}

VertexSet Graph::GetAllVerteces(void)
{
	return vecPtVertices_;
}

VertexSet Graph::GetNeighbors(GraphVertex* vertex)
{
	return vertex->queueNeighbor_;
}

bool Graph::SetWeight(GraphVertex* vertex, double weight)
{
	try
	{
		if (vertex->bValid_)
		{
			vertex->weight_ = weight;
			return true;
		}
	}
	catch (int nError)
	{
		printf("[ERROR] Error is occured at SetWeight!!: %d\n", nError);
	}
	return false;
}

double Graph::GetWeight(GraphVertex* vertex)
{
	try
	{
		if (vertex->bValid_) { return vertex->weight_; }
	}
	catch (int nError)
	{
		printf("[ERROR] Error is occured at GetWeight!!: %d\n", nError);
	}
	return 0.0;
}

void Graph::ClearVertexCounters(void)
{
	for (std::list<GraphVertex>::iterator vertexIter = listVertices_.begin();
		vertexIter != listVertices_.end();
		vertexIter++)
	{
		vertexIter->countNeighborInSolution_ = 0;
	}
}

/////////////////////////////////////////////////////////////////////////
// CGraphSolver CLASS MANAGEMENT
/////////////////////////////////////////////////////////////////////////
CGraphSolver::CGraphSolver(void)
	: pGraph_(NULL)
	, bInit_(false)
	, bHasInitialSolution_(false)
	, bTimeoutSet_(false)
	, bQP_(false)
{
}

CGraphSolver::~CGraphSolver(void)
{
}

void CGraphSolver::Initialize(Graph* pGraph, HJ_GRAPH_SOLVER_TYPE solverType)
{
	if (bInit_) { Clear(); }

	pGraph_ = pGraph;
	bHasInitialSolution_ = false;
	solverType_ = solverType;

	// solver dependent initialization
	switch (solverType_)
	{
	case HJ_GRAPH_SOLVER_BLS:
		BLS_C.clear();
		BLS_PA.clear();
		BLS_OC.clear();
		BLS_OM.clear();
		BLS_nIter = 0;
		break;
	case HJ_GRAPH_SOLVER_BLS4QP:
		BLS_C.clear();
		BLS_PA.clear();
		BLS_OC.clear();
		BLS_OM.clear();
		BLS_nIter = 0;
		break;
	default:
		break;
	}

	bInit_ = true;
}

void CGraphSolver::Initialize(Graph* pGraph, HJ_GRAPH_SOLVER_TYPE solverType, Mat2D &matQ)
{
	Initialize(pGraph, solverType);
	bQP_ = true;
	matQ_ = matQ;
}

void CGraphSolver::Finalize(void)
{
	Clear();
}

void CGraphSolver::Clear(void)
{
	pGraph_ = NULL;

	// solver dependent finalization
	switch (solverType_)
	{
	case HJ_GRAPH_SOLVER_BLS:
		BLS_C.clear();
		BLS_PA.clear();
		BLS_OC.clear();
		BLS_OM.clear();
		break;
	case HJ_GRAPH_SOLVER_BLS4QP:
		BLS_C.clear();
		BLS_PA.clear();
		BLS_OC.clear();
		BLS_OM.clear();
		bQP_ = false;
		matQ_.clear();
		break;
	default:
		break;
	}

	// clear flag
	bInit_ = false;
	bHasInitialSolution_ = false;

	// clear result
	stResult_.solvingTime = 0.0;
	for (int sIdx = 0; sIdx < stResult_.vecSolutions.size(); sIdx++)
	{
		stResult_.vecSolutions[sIdx].first.clear();
	}
	stResult_.vecSolutions.clear();
}

void CGraphSolver::SetInitialSolutions(std::deque<VertexSet> &initialSolutions)
{
	if (0 == initialSolutions.size() || !bInit_) { return; }

	switch (solverType_)
	{
	case HJ_GRAPH_SOLVER_BLS:
		bHasInitialSolution_ = BLS_SetInitialSolutions(initialSolutions);
		break;
	default:
		break;
	}
}


/////////////////////////////////////////////////////////////////////////
// GRAPH SOLVING
/////////////////////////////////////////////////////////////////////////
stGraphSolvingResult* CGraphSolver::Solve(int timelimit)
{
	if (!bInit_ || NULL == pGraph_) { return NULL; }

	time_t timerStart = clock();
	//std::srand((int)std::time(NULL));

	switch (solverType_)
	{
	case HJ_GRAPH_SOLVER_BLS:
		RunBLS(timelimit);
		break;
	case HJ_GRAPH_SOLVER_BLS4QP:
		if (bQP_) { RunBLS(timelimit); }
		break;
	default:
		break;
	}

	time_t timerEnd = clock();
	stResult_.solvingTime = (double)(timerEnd - timerStart) / CLOCKS_PER_SEC;

	return &(stResult_);
}

//#define BLS_NUM_ITER      (100000)
#define BLS_P0            (.75)
#define BLS_T             (10)
#define BLS_PHI           (7)
#define BLS_MIN_ITERATION (200)
#define BLS_MAX_ITERATION (2000)
void CGraphSolver::RunBLS(int timelimit)
{
	// solving information related
	time_t timer_start = clock();
	bTimeoutSet_ = timelimit > 0;
	//std::srand((int)std::time(NULL));
	stResult_.bestIteration = 0;
	stResult_.maximumSearchInterval = 0;

	// parameters
	double L0 = 0.01 * (double)pGraph_->Size();
	double Lmax = 0.10 * (double)pGraph_->Size();
	double alphaS = 0.8;
	double alphaR = 0.8;

	// set iteration number
#ifdef HJ_GRAPH_TYPE_MWISP
	size_t maxIter = pGraph_->Size() * 112;
#else
	size_t maxIter = pGraph_->NumEdge() * 10;
#endif
	maxIter = std::min(std::max((size_t)BLS_MIN_ITERATION, maxIter), (size_t)BLS_MAX_ITERATION);
	stResult_.iterationNumber = (int)maxIter;

	/////////////////////////////////////////////////////////////////////////////
	// INITIAL SOLUTION
	/////////////////////////////////////////////////////////////////////////////
	double L = 0.0;
	double fc = 0.0;
	double fbest = 0.0;
	size_t w = 0;
	VertexSet Cbest, Cp;

	if (bHasInitialSolution_)
	{
		fc = bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C);
	}
	else
	{
		fc = BLS_GenerateInitialSolution();
	}
	Cbest = BLS_C;
	fbest = fc;
	Cp = BLS_C;
	w = 0;

	// save initial solution
	BLS_InsertSolution(Cbest, fbest);

#ifdef HJ_GRAPH_SOLVE_MONITORING
	printf("\n");
#endif
	/////////////////////////////////////////////////////////////////////////////
	// MAIN LOOP
	/////////////////////////////////////////////////////////////////////////////
	for (BLS_nIter = 0; BLS_nIter < maxIter; )
	{
		// timeout
		if (bTimeoutSet_ && (int)(clock() - timer_start) >= timelimit) { break; }

		//------------------------------------------------------
		// LOCAL SEARCH
		//------------------------------------------------------
		double fIncrement = 0.0;
		while (BLS_nIter < maxIter)
		{
			fIncrement = BLS_BestLocalMove();
			if (0.0 >= fIncrement) { break; }
			fc += fIncrement;
			BLS_nIter++;

			// DEBUG
#ifdef HJ_GRAPH_SOLVE_MONITORING
			printf("\r\tSolving : %zd/%zd...", BLS_nIter + 1, maxIter);
#endif
		}

		//------------------------------------------------------
		// SOLUTION CHECK
		//------------------------------------------------------
		fc = bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C); // to handle precision error		
		if (fc > fbest)
		{
			Cbest = BLS_C;
			fbest = fc;
			w = 0;

			// solving info			
			stResult_.maximumSearchInterval = std::max(stResult_.maximumSearchInterval, (int)BLS_nIter - stResult_.bestIteration);
			stResult_.bestIteration = (int)BLS_nIter;
		}
		else
		{
			w++;
		}

		//------------------------------------------------------
		// PERTURBATION
		//------------------------------------------------------
		// sort for comparison
		std::sort(BLS_C.begin(), BLS_C.end(), HJGraphComparatorVertexIDAscend);
		if (w > BLS_T)
		{
			L = Lmax;
			w = 0;
		}
		else if (BLS_C == Cp)
		{
			L++;
		}
		else
		{
			BLS_InsertSolution(BLS_C, fc);
			L = L0;
		}
		Cp = BLS_C;
		fc = BLS_Perturbation(L, w, alphaR, alphaS);
	}

	/////////////////////////////////////////////////////////////////////////////
	// TERMINATION
	/////////////////////////////////////////////////////////////////////////////
	// sort the solution (descending order of weight sum)
	std::sort(stResult_.vecSolutions.begin(), stResult_.vecSolutions.end(), HJSolutionComparatorSumWeightsDescend);
	// sort the best solution with vertex ID
	for (int sIdx = 0; sIdx < stResult_.vecSolutions.size(); sIdx++)
	{
		std::sort(stResult_.vecSolutions[sIdx].first.begin(), stResult_.vecSolutions[sIdx].first.end(), HJGraphComparatorVertexIDAscend);	// by track ID
	}

	// result packing
	time_t timer_end = clock();
	stResult_.solvingTime = (double)(timer_end - timer_start) / CLOCKS_PER_SEC;
	stResult_.numEdges = pGraph_->NumEdge();
	stResult_.numVertices = pGraph_->Size();
	stResult_.maximumDegree = pGraph_->maxDegree();
	stResult_.averageDegree = pGraph_->AverageVertexDegree();
}

void CGraphSolver::RunILS(void)
{
}

void CGraphSolver::RunAMTS(void)
{
}

void CGraphSolver::RunMCGA(void)
{
}

/////////////////////////////////////////////////////////////////////////
// MISCELLANEOUS
/////////////////////////////////////////////////////////////////////////
#define Graph_SOLUTION_DUPLICATION_RESOLUTION (1.0E-5)
bool CGraphSolver::CheckSolutionExistance(const VertexSet &vertexSet, double SumWeights)
{
	for (size_t solutionIdx = 0; solutionIdx < stResult_.vecSolutions.size(); solutionIdx++)
	{
		if (std::abs(SumWeights - stResult_.vecSolutions[solutionIdx].second) > Graph_SOLUTION_DUPLICATION_RESOLUTION)
		{
			continue;
		}
		if (vertexSet == stResult_.vecSolutions[solutionIdx].first)
		{
			return true;
		}
	}
	return false;
}

inline double CGraphSolver::SumWeights(VertexSet &vertexSet)
{
	double SumWeights = 0.0;
	for (VertexSet::iterator vertexIter = vertexSet.begin();
		vertexIter != vertexSet.end();
		vertexIter++)
	{
		SumWeights += (*vertexIter)->weight_;
	}
	return SumWeights;
}

inline double CGraphSolver::SumWeightsQP(VertexSet &vertexSet)
{
	double SumWeights = 0.0;
	for (int vertexIdx1 = 0; vertexIdx1 < vertexSet.size(); vertexIdx1++)
	{
		SumWeights += SumWeightsQP(vertexSet, vertexSet[vertexIdx1]);
	}
	return SumWeights;
}

inline double CGraphSolver::SumWeightsQP(VertexSet &vertexSet, GraphVertex *targetVertex)
{
	// unary
	double fWeightSum = targetVertex->weight_;
	// pair-wise
	for (int vIdx = 0; vIdx < vertexSet.size(); vIdx++)
	{
		fWeightSum += matQ_.at((int)targetVertex->id_, (int)vertexSet[vIdx]->id_);
	}
	return fWeightSum;
}

/************************************************************************
Method Name: PrintLog
Description:
-
Input Arguments:
-
Return Values:
-
************************************************************************/
void CGraphSolver::PrintLog(size_t type)
{
	printf("[iteration:%04zd] ", BLS_nIter);
	switch (type)
	{
	case 0:
		printf("local search");
		break;
	case 1:
		printf("directed perturbation");
		break;
	case 2:
		printf("random perturbation");
		break;
	default:
		break;
	}
	printf("\n");

	// print C
	printf("C(%d)={", (int)BLS_C.size());
	for (size_t vertexIdx = 0; vertexIdx < BLS_C.size(); vertexIdx++)
	{
		if (0 < vertexIdx) { printf(","); }
		printf("%zd", BLS_C[vertexIdx]->id_);
	}
	printf("}\n");

	// print PA
	printf("PA(%d)={", (int)BLS_PA.size());
	for (size_t vertexIdx = 0; vertexIdx < BLS_PA.size(); vertexIdx++)
	{
		if (0 < vertexIdx) { printf(","); }
		printf("%zd", BLS_PA[vertexIdx]->id_);
	}
	printf("}\n");

	// print OM
	printf("OM(%d)={", (int)BLS_OM.size());
	for (size_t pairIdx = 0; pairIdx < BLS_OM.size(); pairIdx++)
	{
		if (0 < pairIdx) { printf(","); }
		printf("(%zd,%zd)", BLS_OM[pairIdx].first->id_, BLS_OM[pairIdx].second->id_);
	}
	printf("}\n");

	// print OC
	printf("OC(%d)={", (int)BLS_OC.size());
	for (size_t vertexIdx = 0; vertexIdx < BLS_OC.size(); vertexIdx++)
	{
		if (0 < vertexIdx) { printf(","); }
		printf("%zd", BLS_OC[vertexIdx]->id_);
	}
	printf("}\n");

	printf("\n");
}

/////////////////////////////////////////////////////////////////////////
// BLS RELATED
/////////////////////////////////////////////////////////////////////////

/************************************************************************
Method Name: BLS_SetInitialSolution
Description:
- set initial solution for BLS solver. it also check the validity of
input solution. When a solution is valid_, do local search for the best
starting position of the solution search. At the final, the best initial
solution is set to the initial position of perturbation.
Input Arguments:
- initialSolution: candidate initial solution
Return Values:
- bool: return whether the input solution is valid or not
************************************************************************/
bool CGraphSolver::BLS_SetInitialSolutions(std::deque<VertexSet> &initialSolutions)
{
	std::vector<VertexSet> BLS_C_init;  BLS_C_init.reserve(initialSolutions.size());               // solution vertices
	std::vector<VertexSet> BLS_PA_init; BLS_PA_init.reserve(initialSolutions.size());              // not in C, but connected to all v in C
	std::vector<VertexSet> BLS_OC_init; BLS_OC_init.reserve(initialSolutions.size());              // not in C
	std::vector<std::deque<VertexPair>> BLS_OM_init(initialSolutions.size()); // <v in OC, u in C>	: except u, v is connected to all vertices in C
	BLS_OM_init.reserve(initialSolutions.size());
	double maximumWeight = DBL_MIN;
	size_t maximumSolutionIdx = 0;

	for (int sIdx = 0; sIdx < initialSolutions.size(); sIdx++)
	{
		VertexSet initialSolution = initialSolutions[sIdx];
		BLS_C.clear();
		BLS_PA.clear();
		BLS_OC.clear();
		BLS_OM.clear();

		// C
		for (VertexSet::iterator vertexIter = initialSolution.begin();
			vertexIter != initialSolution.end();
			vertexIter++)
		{
			BLS_C.push_back(*vertexIter);
			(*vertexIter)->CountUp();
		}

		//------------------------------------------------------
		// SOLUTION VALIDATION
		//------------------------------------------------------			
#ifdef HJ_GRAPH_TYPE_MWISP
		size_t sizeCheckInC = 0;
#else	
		size_t sizeCheckInC = this->BLS_C.size();
#endif
		for (VertexSet::iterator vertexIter = initialSolution.begin();
			vertexIter != initialSolution.end();
			vertexIter++)
		{
			if (sizeCheckInC == (*vertexIter)->countNeighborInSolution_) { continue; }

			// illegal solution
			BLS_C.clear();
			pGraph_->ClearVertexCounters();
			return false;
		}

		//------------------------------------------------------
		// CONSTRUCTION FOR FINDING LOCAL OPTIMUM
		//------------------------------------------------------		
		// PA and OC
#ifdef HJ_GRAPH_TYPE_MWISP
		size_t sizeOM = 1;
		size_t sizePA = 0;
#else
		size_t sizeOM = this->BLS_C.size() - 1;
		size_t sizePA = this->BLS_C.size();
#endif
		VertexSet allVertex = pGraph_->GetAllVerteces();
		for (VertexSet::iterator vertexIter = allVertex.begin();
			vertexIter != allVertex.end();
			vertexIter++)
		{
			// OC
			if ((*vertexIter)->bInSolution_) { continue; }
			this->BLS_OC.push_back(*vertexIter);

			// PA
			if (sizePA == (*vertexIter)->countNeighborInSolution_)
			{
				this->BLS_PA.push_back(*vertexIter);
				continue;
			}

			//OM
			if (sizeOM != (*vertexIter)->countNeighborInSolution_)
			{
				continue;
			}

#ifdef HJ_GRAPH_TYPE_MWISP
			// find pair
			for (VertexSet::iterator neighborIter = (*vertexIter)->queueNeighbor_.begin();
				neighborIter != (*vertexIter)->queueNeighbor_.end();
				neighborIter++)
			{
				if ((*neighborIter)->bInSolution_)
				{
					this->BLS_OM.push_back(std::make_pair(*vertexIter, *neighborIter));
					(*vertexIter)->bInOM_ = true;
					break;
				}
			}
#else
			// find another vertex
			VertexSet::iterator findIter;
			for (VertexSet::iterator vertexInCIter = this->BLS_C.begin();
				vertexInCIter != this->BLS_C.end();
				vertexInCIter++)
			{
				findIter = std::find((*vertexIter)->queueNeighbor_.begin(), (*vertexIter)->queueNeighbor_.end(), *vertexInCIter);
				if ((*vertexIter)->queueNeighbor_.end() != findIter) { continue; }
				this->BLS_OM.push_back(std::make_pair(*vertexIter, *vertexInCIter));
				(*vertexIter)->bInOM_ = true;
				break;
			}
#endif
		}

		//------------------------------------------------------
		// FIND LOCAL OPTIMUM
		//------------------------------------------------------
		BLS_BestLocalMove();
		double curWeights = SumWeights(BLS_C);
		if (BLS_InsertSolution(BLS_C, curWeights))
		{
			if (curWeights > maximumWeight)
			{
				maximumWeight = curWeights;
				maximumSolutionIdx = BLS_C_init.size();
			}

			BLS_C_init.push_back(BLS_C);
			BLS_PA_init.push_back(BLS_PA);
			BLS_OC_init.push_back(BLS_OC);
			BLS_OM_init.push_back(BLS_OM);
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	// PICK THE MAXIMUM SCORED SOLUTION
	/////////////////////////////////////////////////////////////////////////////
	BLS_C = BLS_C_init[maximumSolutionIdx];
	BLS_PA = BLS_PA_init[maximumSolutionIdx];
	BLS_OC = BLS_OC_init[maximumSolutionIdx];
	BLS_OM = BLS_OM_init[maximumSolutionIdx];

	return true;
}

/************************************************************************
Method Name: BLS_GenerateInitialSolution
Description:
- find initial solution of clique
Input Arguments:
- bGreedy: select the method of generating initial solution
Return Values:
- double: objective value of solution set
************************************************************************/
bool CGraphSolver::BLS_InsertSolution(const VertexSet &solution, const double solutionScore)
{
	if (0 == solution.size() || 0.0 > solutionScore || CheckSolutionExistance(solution, solutionScore))
	{
		return false;
	}
	stResult_.vecSolutions.push_back(std::make_pair(solution, solutionScore));
	return true;
}

/************************************************************************
Method Name: BLS_GenerateInitialSolution
Description:
- find initial solution of clique
Input Arguments:
- bGreedy: select the method of generating initial solution
Return Values:
- double: objective value of solution set
************************************************************************/
double CGraphSolver::BLS_GenerateInitialSolution(bool bGreedy)
{
	BLS_C.clear();	// solution set
	BLS_PA.clear();	// outer vertex fully connected with C
	BLS_OM.clear();	// exclusive vertex pair set
	BLS_OC.clear();	// complement of the clique

	VertexSet allVertex = pGraph_->GetAllVerteces();
	for (size_t insertIdx = 0; insertIdx < allVertex.size(); insertIdx++)
	{
		allVertex[insertIdx]->tabuStamp_ = 0;
		allVertex[insertIdx]->countNeighborInSolution_ = 0;
	}

	// construct C	
	if (bGreedy)
		std::sort(allVertex.begin(), allVertex.end(), HJGraphComparatorVertexWeightDescend);
	else
		std::random_shuffle(allVertex.begin(), allVertex.end());

	for (VertexSet::iterator vertexIter = allVertex.begin();
		vertexIter != allVertex.end();
		vertexIter++)
	{
		bool bInC = false;
#ifdef HJ_GRAPH_TYPE_MWISP
		if (0 == (*vertexIter)->countNeighborInSolution_)
#else
		if (BLS_C.size() == (*vertexIter)->countNeighborInSolution_)
#endif
		{
			if ((*vertexIter)->weight_ < 0)
			{
				if (bQP_)
				{
					double QPWeightSum = 0.0;
					for (VertexSet::iterator vertexInCIter = BLS_C.begin();
						vertexInCIter != BLS_C.end();
						vertexInCIter++)
					{
						QPWeightSum += matQ_.at((int)(*vertexIter)->id_, (int)(*vertexInCIter)->id_);
					}
					if (QPWeightSum + (*vertexIter)->weight_ >= 0.0) { bInC = true; }
				}
			}
			else
			{
				bInC = true;
			}
		}
		if (bInC)
		{
			BLS_C.push_back(*vertexIter);
			(*vertexIter)->CountUp();
		}
		else
		{
			BLS_OC.push_back(*vertexIter);
		}
	}

	// sort by track ID
	std::sort(BLS_C.begin(), BLS_C.end(), HJGraphComparatorVertexIDAscend);

	// construct OM
	for (VertexSet::iterator vertexIter = BLS_OC.begin();
		vertexIter != BLS_OC.end();
		vertexIter++)
	{
#ifdef HJ_GRAPH_TYPE_MWISP
		if (1 != (*vertexIter)->countNeighborInSolution_)
		{
			continue;
		}

		// find pair
		for (VertexSet::iterator neighborIter = (*vertexIter)->queueNeighbor_.begin();
			neighborIter != (*vertexIter)->queueNeighbor_.end();
			neighborIter++)
		{
			if ((*neighborIter)->bInSolution_)
			{
				this->BLS_OM.push_back(std::make_pair(*vertexIter, *neighborIter));
				(*vertexIter)->bInOM_ = true;
				break;
			}
		}
#else
		if (BLS_C.size() - 1 != (*vertexIter)->countNeighborInSolution_) { continue; }

		// find another vertex
		VertexSet::iterator findIter;
		for (VertexSet::iterator vertexInCIter = BLS_C.begin();
			vertexInCIter != BLS_C.end();
			vertexInCIter++)
		{
			findIter = std::find((*vertexIter)->queueNeighbor_.begin(), (*vertexIter)->queueNeighbor_.end(), *vertexInCIter);
			if ((*vertexIter)->queueNeighbor_.end() != findIter) { continue; }
			BLS_OM.push_back(std::make_pair(*vertexIter, *vertexInCIter));
			(*vertexIter)->bInOM_ = true;
			break;
		}
#endif
	}

	return bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C);
}

/************************************************************************
Method Name: BLS_BestLocalMove
Description:
- do the best local move
Input Arguments:
- none
Return Values:
- double: increment in objective value of solution set
************************************************************************/
double CGraphSolver::BLS_BestLocalMove(void)
{
	double maxObjectiveIncrement = 0.0;
	double curObjectiveValue = 0.0;
	GraphVertex *vertexInsert = NULL;
	GraphVertex *vertexRemove = NULL;

	/////////////////////////////////////////////////////////////////////////////
	// SEARCH
	/////////////////////////////////////////////////////////////////////////////
	if (0 < BLS_PA.size())
	{
		for (int paIdx = 0; paIdx < BLS_PA.size(); paIdx++)
		{
			// unary
			curObjectiveValue = BLS_PA[paIdx]->weight_;

			// pair-wise
			if (bQP_)
			{
				for (int cIdx = 0; cIdx < BLS_C.size(); cIdx++)
				{
					curObjectiveValue += matQ_.at((int)BLS_PA[paIdx]->id_, (int)BLS_C[cIdx]->id_);
				}
			}

			if (curObjectiveValue < maxObjectiveIncrement) { continue; }
			maxObjectiveIncrement = curObjectiveValue;
			vertexInsert = BLS_PA[paIdx];
		}
	}

	// movement in OM
	for (size_t pairIdx = 0; pairIdx < BLS_OM.size(); pairIdx++)
	{
		if (bQP_)
		{
			curObjectiveValue = SumWeightsQP(BLS_C, BLS_OM[pairIdx].first) - SumWeightsQP(BLS_C, BLS_OM[pairIdx].second);
		}
		else
		{
			curObjectiveValue = BLS_OM[pairIdx].first->weight_ - BLS_OM[pairIdx].second->weight_;
		}
		if (curObjectiveValue < maxObjectiveIncrement) { continue; }
		maxObjectiveIncrement = curObjectiveValue;
		vertexInsert = BLS_OM[pairIdx].first;
		vertexRemove = BLS_OM[pairIdx].second;
	}

	if (NULL == vertexInsert) { return 0.0; }

	// remove ('remove -> insert' is more effective than 'insert -> remove')
	BLS_VertexRemove(vertexRemove);

	// insert
	BLS_VertexInsert(vertexInsert);

	return maxObjectiveIncrement;
}

/************************************************************************
Method Name: BLS_Perturbation
Description:
- select perturbation condition and do perturb
Input Arguments:
- L: perturbation strength
- w: number of consecutive non-improving local optima visited
- alphaR: coefficient
- alphaS: coefficient
Return Values:
- double: objective value of solution set
************************************************************************/
const double BLS_wCheck = -BLS_T * std::log(BLS_P0);
double CGraphSolver::BLS_Perturbation(double L, size_t w, double alphaR, double alphaS)
{
	VertexSet newC;
	if (0 == w) { return BLS_PerturbRandom(L, alphaS); }

	double P = (double)w < BLS_wCheck ? std::exp(-(double)w / (double)BLS_T) : BLS_P0;
	double unitRand = (double)rand() / (double)RAND_MAX;
	if (P >= unitRand) { return BLS_PerturbDirected(L); }

	return BLS_PerturbRandom(L, alphaR);
}

/************************************************************************
Method Name: BLS_PerturbDirected
Description:
- function for directed perturbation
Input Arguments:
- L: perturbation strength
Return Values:
- double: objective value of solution set
************************************************************************/
double CGraphSolver::BLS_PerturbDirected(double L)
{
	double fVal = 0.0;
	GraphVertex *vertexInsert = NULL;
	GraphVertex *vertexRemove = NULL;
	for (size_t LIdx = 0; LIdx < L; LIdx++)
	{
		/////////////////////////////////////////////////////////////////////////////
		// UPDATE M AND SELECT MOVE
		/////////////////////////////////////////////////////////////////////////////
		fVal = bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C);
		std::vector<std::pair<BLS_MOVE_TYPE, size_t>> queueMove;
		queueMove.reserve(BLS_PA.size() + BLS_OM.size() + BLS_C.size());
		for (size_t vertexIdx = 0; vertexIdx < BLS_PA.size(); vertexIdx++)
		{
			// insert from PA
			if (BLS_PA[vertexIdx]->tabuStamp_ <= BLS_nIter)
			{
				queueMove.push_back(std::make_pair(M1, vertexIdx));
			}
		}
		for (size_t vertexIdx = 0; vertexIdx < BLS_OM.size(); vertexIdx++)
		{
			// switch
			if (BLS_OM[vertexIdx].first->tabuStamp_ <= BLS_nIter)
			{
				queueMove.push_back(std::make_pair(M2, vertexIdx));
			}
		}
		for (size_t vertexIdx = 0; vertexIdx < BLS_C.size(); vertexIdx++)
		{
			// remove from solution set			
			queueMove.push_back(std::make_pair(M3, vertexIdx));
		}

		if (0 == queueMove.size())
		{
			BLS_nIter++;
			break;
		}

		// random selection
		size_t selectedMoveIdx = (size_t)((double)(queueMove.size() - 1) * ((double)rand() / (double)RAND_MAX));

		/////////////////////////////////////////////////////////////////////////////
		// MOVE
		/////////////////////////////////////////////////////////////////////////////		
		if (M1 == queueMove[selectedMoveIdx].first)
		{
			vertexInsert = BLS_PA[queueMove[selectedMoveIdx].second];
			vertexRemove = NULL;
		}
		else if (M2 == queueMove[selectedMoveIdx].first)
		{
			vertexInsert = BLS_OM[queueMove[selectedMoveIdx].second].first;
			vertexRemove = BLS_OM[queueMove[selectedMoveIdx].second].second;
		}
		else
		{
			vertexInsert = NULL;
			vertexRemove = BLS_C[queueMove[selectedMoveIdx].second];
		}

		// remove
		fVal += BLS_VertexRemove(vertexRemove);

		// insert
		fVal += BLS_VertexInsert(vertexInsert);

		// increase iteration number
		BLS_nIter++;
	}

	return fVal;
}

/************************************************************************
Method Name: BLS_PerturbRandom
Description:
- function for random perturbation
Input Arguments:
- L: perturbation strength
- alpha: coefficient
Return Values:
- double: objective value of solution set
************************************************************************/
double CGraphSolver::BLS_PerturbRandom(double L, double alpha)
{
	// calculate the objective value of current solution
	double fVal = bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C);
	GraphVertex *vertexInsert = NULL;
	for (size_t LIdx = 0; LIdx < L; LIdx++)
	{
		/////////////////////////////////////////////////////////////////////////////
		// UPDATE M AND SELECT MOVE
		/////////////////////////////////////////////////////////////////////////////
		std::vector<size_t> queueMove;
		queueMove.reserve(BLS_OC.size());

		for (size_t vertexIdx = 0; vertexIdx < BLS_OC.size(); vertexIdx++)
		{
			if (BLS_OC[vertexIdx]->tabuStamp_ <= BLS_nIter)
			{
				queueMove.push_back(vertexIdx);
				continue;
			}

			double neighborSumWeightsInC = 0.0;
			for (VertexSet::iterator vertexIter = BLS_OC[vertexIdx]->queueNeighbor_.begin();
				vertexIter != BLS_OC[vertexIdx]->queueNeighbor_.end();
				vertexIter++)
			{
				if ((*vertexIter)->bInSolution_)
				{
					neighborSumWeightsInC += (*vertexIter)->weight_;
					if (bQP_)
					{
						neighborSumWeightsInC += matQ_.at((int)BLS_OC[vertexIdx]->id_, (int)(*vertexIter)->id_);
					}
				}
			}
			if (neighborSumWeightsInC >= alpha * fVal) { queueMove.push_back(vertexIdx); }
		}

		if (0 == queueMove.size())
		{
			BLS_nIter++;
			break;
		}

		// random selection
		size_t selectedMoveIdx = (size_t)((double)(queueMove.size() - 1) * ((double)rand() / (double)RAND_MAX));

		/////////////////////////////////////////////////////////////////////////////
		// MOVE
		/////////////////////////////////////////////////////////////////////////////
		fVal = BLS_VertexInsertM4(selectedMoveIdx);

		// increase iteration number
		BLS_nIter++;
	}

	return fVal;
}

/************************************************************************
Method Name: BLS_VertexRemove
Description:
- [NOTICE] Always be used with BLS_VertexInsert (for management of OM)
- remove vertex from C and repair PA, OM and OC
Input Arguments:
- vertexRemove: target vertex
Return Values:
- double: increment in objective value
************************************************************************/
double CGraphSolver::BLS_VertexRemove(GraphVertex *vertexRemove)
{
	if (NULL == vertexRemove) { return 0.0; }
	double fIncrements = bQP_ ? -SumWeightsQP(BLS_C, vertexRemove) : -vertexRemove->weight_;

	// remove from C
	VertexSet::iterator removeIter = std::find(BLS_C.begin(), BLS_C.end(), vertexRemove);
	BLS_C.erase(removeIter);
	vertexRemove->CountDown();

	// insert to PA
	BLS_PA.push_back(vertexRemove);

	// insert to OC and update tabu_list
	BLS_OC.push_back(vertexRemove);
	vertexRemove->tabuStamp_ = BLS_GetMoveProhibitionNumber(BLS_nIter);

	if (0 == BLS_C.size())
	{
		BLS_PA = BLS_OC;
		BLS_OM.clear();
		return fIncrements;
	}

	// update OM (-> OM and PA)
	std::deque<VertexPair> newOM;
	for (std::deque<VertexPair>::iterator pairIter = BLS_OM.begin();
		pairIter != BLS_OM.end();
		pairIter++)
	{
		if (vertexRemove != (*pairIter).second)
		{
			newOM.push_back(*pairIter);
			continue;
		}
		BLS_PA.push_back((*pairIter).first);
		(*pairIter).first->bInOM_ = false;
	}

#ifdef HJ_GRAPH_TYPE_MWISP
	size_t sizeCheckOM = 1;
#else
	size_t sizeCheckOM = this->BLS_C.size() - 1;
#endif	

	// add to OM (<- non OM)
	for (VertexSet::iterator vertexIter = BLS_OC.begin();
		vertexIter != BLS_OC.end();
		vertexIter++)
	{
		if ((*vertexIter)->bInOM_ || (*vertexIter)->countNeighborInSolution_ != sizeCheckOM) { continue; }
		(*vertexIter)->bInOM_ = true;

#ifdef HJ_GRAPH_TYPE_MWISP
		for (VertexSet::iterator neighborIter = (*vertexIter)->queueNeighbor_.begin();
			neighborIter != (*vertexIter)->queueNeighbor_.end();
			neighborIter++)
		{
			if ((*neighborIter)->bInSolution_)
			{
				newOM.push_back(std::make_pair(*vertexIter, *neighborIter));
				break;
			}
		}
#else
		VertexSet::iterator findIter;
		for (VertexSet::iterator solutionIter = BLS_C.begin();
			solutionIter != BLS_C.end();
			solutionIter++)
		{
			findIter = std::find((*vertexIter)->queueNeighbor_.begin(), (*vertexIter)->queueNeighbor_.end(), *solutionIter);
			if ((*vertexIter)->queueNeighbor_.end() == findIter)
			{
				newOM.push_back(std::make_pair(*vertexIter, *solutionIter));
				break;
			}
		}
#endif
	}

	BLS_OM = newOM;
	newOM.clear();
	return fIncrements;
}

/************************************************************************
Method Name: BLS_VertexInsert
Description:
- insert vertex to C by M1 or M2 and repair PA, OM and OC
Input Arguments:
- vertexInsert: target vertex
Return Values:
- double: increment in objective value
************************************************************************/
double CGraphSolver::BLS_VertexInsert(GraphVertex *vertexInsert)
{
	if (NULL == vertexInsert) { return 0.0; }

	// remove from OC
	VertexSet::iterator removeIter = std::find(BLS_OC.begin(), BLS_OC.end(), vertexInsert);
	BLS_OC.erase(removeIter);

	// insert to C
	BLS_C.push_back(vertexInsert);
	vertexInsert->CountUp();

#ifdef HJ_GRAPH_TYPE_MWISP
	size_t sizeCheckPA = 0;
	size_t sizeCheckOM = 1;
#else
	size_t sizeCheckPA = BLS_C.size();
	size_t sizeCheckOM = BLS_C.size() - 1;
#endif	

	// update OM
	std::deque<VertexPair> newOM;
	for (std::deque<VertexPair>::iterator pairIter = BLS_OM.begin();
		pairIter != BLS_OM.end();
		pairIter++)
	{
		if (sizeCheckOM == (*pairIter).first->countNeighborInSolution_)
		{
			newOM.push_back(*pairIter);
		}
		else
		{
			(*pairIter).first->bInOM_ = false;
		}
	}
	BLS_OM = newOM;
	newOM.clear();

	// update PA
	VertexSet newPA;
	for (VertexSet::iterator vertexIter = BLS_PA.begin();
		vertexIter != BLS_PA.end();
		vertexIter++)
	{
		if (vertexInsert == *vertexIter) { continue; }

		if (sizeCheckPA == (*vertexIter)->countNeighborInSolution_)
		{
			newPA.push_back(*vertexIter);
		}
		else
		{
			BLS_OM.push_back(std::make_pair(*vertexIter, vertexInsert));
			(*vertexIter)->bInOM_ = true;
		}
	}
	BLS_PA = newPA;
	newPA.clear();

	double fIncrement = bQP_ ? SumWeightsQP(BLS_C, vertexInsert) : vertexInsert->weight_;

	return fIncrement;
}

/************************************************************************
Method Name: BLS_VertexInsertM4
Description:
- insert vertex to C by M4(from OC) and repair PA, OM and OC
Input Arguments:
- vertexIdxInOC: index of target vertex in OC
Return Values:
- double: objective value of solution set
************************************************************************/
double CGraphSolver::BLS_VertexInsertM4(size_t vertexIdxInOC)
{
	/////////////////////////////////////////////////////////////////////////////
	// REPAIR C
	/////////////////////////////////////////////////////////////////////////////
	GraphVertex *vertexInsert = BLS_OC[vertexIdxInOC];
	for (VertexSet::iterator vertexIter = vertexInsert->queueNeighbor_.begin();
		vertexIter != vertexInsert->queueNeighbor_.end();
		vertexIter++)
	{
		(*vertexIter)->bInSolution_ = false;	// toggle for update C
	}

	// update C
	VertexSet newC;
	for (VertexSet::iterator vertexIter = BLS_C.begin();
		vertexIter != BLS_C.end();
		vertexIter++)
	{
#ifndef HJ_GRAPH_TYPE_MWISP
		(*vertexIter)->bInSolution_ = !(*vertexIter)->bInSolution_;
#endif
		if ((*vertexIter)->bInSolution_)
		{
			newC.push_back(*vertexIter);
			continue;
		}
		(*vertexIter)->CountDown();
		this->BLS_OC.push_back(*vertexIter);
	}
	newC.push_back(vertexInsert);
	vertexInsert->CountUp();
	BLS_C = newC;

	// remove from OC
	BLS_OC.erase(BLS_OC.begin() + vertexIdxInOC);

#ifdef HJ_GRAPH_TYPE_MWISP
	size_t sizeCheckPA = 0;
	size_t sizeCheckOM = 1;
#else
	size_t sizeCheckPA = BLS_C.size();
	size_t sizeCheckOM = BLS_C.size() - 1;
#endif

	/////////////////////////////////////////////////////////////////////////////
	// REPAIR PA
	/////////////////////////////////////////////////////////////////////////////
	BLS_PA.clear();
	VertexSet OMCandidate;
	for (VertexSet::iterator vertexIter = BLS_OC.begin();
		vertexIter != BLS_OC.end();
		vertexIter++)
	{
		if (sizeCheckPA == (*vertexIter)->countNeighborInSolution_
			&& !(*vertexIter)->bInSolution_)
		{
			BLS_PA.push_back(*vertexIter);
		}
		else if (sizeCheckOM == (*vertexIter)->countNeighborInSolution_
			&& !(*vertexIter)->bInOM_)
		{
			OMCandidate.push_back(*vertexIter);
		}
	}


	/////////////////////////////////////////////////////////////////////////////
	// REPAIR OM
	/////////////////////////////////////////////////////////////////////////////
	std::deque<VertexPair> newOM;
	for (std::deque<VertexPair>::iterator pairIter = BLS_OM.begin();
		pairIter != BLS_OM.end();
		pairIter++)
	{
		if (sizeCheckOM != (*pairIter).first->countNeighborInSolution_)
		{
			(*pairIter).first->bInOM_ = false;
			continue;
		}

		if ((*pairIter).second->bInSolution_)
		{
			newOM.push_back(*pairIter);
		}
		else
		{
			newOM.push_back(std::make_pair((*pairIter).first, vertexInsert));
		}
	}

	// from OM candidate
	VertexSet::iterator findIter;
	for (VertexSet::iterator vertexIter = OMCandidate.begin();
		vertexIter != OMCandidate.end();
		vertexIter++)
	{
#ifdef HJ_GRAPH_TYPE_MWISP
		for (VertexSet::iterator neighborIter = (*vertexIter)->queueNeighbor_.begin();
			neighborIter != (*vertexIter)->queueNeighbor_.end();
			neighborIter++)
		{
			if ((*neighborIter)->bInSolution_)
			{
				newOM.push_back(std::make_pair(*vertexIter, *neighborIter));
				(*vertexIter)->bInOM_ = true;
				break;
			}
		}
#else
		for (VertexSet::iterator solutionVertexIter = BLS_C.begin();
			solutionVertexIter != BLS_C.end();
			solutionVertexIter++)
		{
			findIter = std::find((*vertexIter)->queueNeighbor_.begin(), (*vertexIter)->queueNeighbor_.end(), *solutionVertexIter);
			if ((*vertexIter)->queueNeighbor_.end() != findIter)
			{
				// found
				continue;
			}
			newOM.push_back(std::make_pair(*vertexIter, *solutionVertexIter));
			(*vertexIter)->bInOM_ = true;
			break;
		}
#endif
	}
	BLS_OM = newOM;

	return bQP_ ? SumWeightsQP(BLS_C) : SumWeights(BLS_C);
}

/************************************************************************
 Method Name: BLS_GetMoveProhibitionNumber
 Description:
	- get the ending iteration number of move prohibition
 Input Arguments:
	- sizeOM: size of OM set
	- nIter: current iteration
 Return Values:
	- size_t: ending iteration number of move prohibition
************************************************************************/
size_t CGraphSolver::BLS_GetMoveProhibitionNumber(size_t nIter)
{
	return nIter + BLS_PHI + (size_t)((double)BLS_OM.size() * ((double)rand() / (double)RAND_MAX));
}

/////////////////////////////////////////////////////////////////////////
// BLS4QP RELATED
/////////////////////////////////////////////////////////////////////////

} // end of hj namespace

//()()
//('')HAANJU.YOO
