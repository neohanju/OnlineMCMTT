/******************************************************************************
* Title        : CGraphSolver class implementation
* Author       : Haanju Yoo
* Initial Date : 2014.08.31 (ver. 0.9)
* Version Num. : 1.0 (since 2016.09.06)
* Description  : Class for defining and solving graphs
******************************************************************************/

#pragma once

#include <vector>
#include <queue>
#include <list>

typedef enum {
	HJ_GRAPH_SOLVER_BLS = 0,
	HJ_GRAPH_SOLVER_BLS4QP,
	HJ_GRAPH_SOLVER_ILS,
	HJ_GRAPH_SOLVER_AMTS,
	HJ_GRAPH_SOLVER_MCGS
} HJ_GRAPH_SOLVER_TYPE;

namespace hj {

	class Mat2D
	{
	public:
		Mat2D(void);
		Mat2D(int rows, int cols);
		~Mat2D(void);
		void clear(void);
		void resize(int rows, int cols);
		double& at(int row, int col);
		int rows(void) const { return rows_; };
		int cols(void) const { return cols_; };
		//Mat2D operator=(const Mat2D &x);

	private:
		std::vector<std::vector<double>> data_;
		int rows_;
		int cols_;
	};

	/////////////////////////////////////////////////////////////////////////
	// TYPEDEFS
	/////////////////////////////////////////////////////////////////////////
	class GraphVertex
	{
		//------------------------------------------------------
		// METHODS
		//------------------------------------------------------
	public:
		// constructor
		GraphVertex();
		GraphVertex(size_t nID);
		GraphVertex(size_t nID, double weight);

		// interface
		void CountUp(void);
		void CountDown(void);
		size_t degree() const { return queueNeighbor_.size(); }

	private:
		// initialization
		void Initialize(size_t nID, double weight);

		//------------------------------------------------------
		// VARIABLES
		//------------------------------------------------------
	public:
		size_t id_;
		bool   bValid_;
		double weight_;
		std::deque<GraphVertex*> queueNeighbor_;

		// BLS related
		bool bInOM_;
		bool bInSolution_;
		size_t tabuStamp_;
		size_t countNeighborInSolution_;
	};
	typedef std::deque<GraphVertex*> VertexSet;

	class GraphEdge
	{
		// constructors
		GraphEdge(GraphVertex *argVertex1, GraphVertex *argVertex2) : id_(0), bValid_(true), vertex1_(argVertex1), vertex2_(argVertex2) {}

	public:
		size_t id_;
		bool   bValid_;
		GraphVertex *vertex1_;
		GraphVertex *vertex2_;
	};

	class Graph
	{
		//------------------------------------------------------
		// METHODS
		//------------------------------------------------------
	public:
		Graph(void);
		Graph(size_t nNumVertex);
		~Graph(void);

		bool   Clear(void);
		bool   TopologyModified(void);
		size_t Size(void) const { return nNumVertex_; }
		size_t NumEdge(void) const { return nNumEdge_; }
		size_t maxDegree(void);
		size_t minDegree(void);
		double AverageVertexDegree(void);

		GraphVertex* AddVertex(double weight);
		VertexSet    AddVertices(size_t numVertex);
		bool DeleteVertex(GraphVertex* vertex);
		bool AddEdge(GraphVertex* vertex1, GraphVertex* vertex2);
		bool Update(void);

		VertexSet GetAllVerteces(void);
		VertexSet GetNeighbors(GraphVertex* vertex);
		bool SetWeight(GraphVertex* vertex, double weight);
		double GetWeight(GraphVertex* vertex);
		void ClearVertexCounters(void);

		//------------------------------------------------------
		// VARIABLES
		//------------------------------------------------------
	private:
		bool bTopologyModified_;
		size_t nNumVertex_;
		size_t nNumEdge_;
		size_t nNewID_;
		std::list<GraphVertex> listVertices_;
		VertexSet vecPtVertices_;
	};

	struct stGraphSolvingResult
	{
		std::vector<std::pair<VertexSet, double>> vecSolutions;
		double solvingTime;
		int iterationNumber;
		int bestIteration;
		int maximumSearchInterval;
		size_t numEdges;
		size_t numVertices;
		size_t maximumDegree;
		double averageDegree;
	};
	typedef std::pair<GraphVertex*, GraphVertex*> VertexPair;
	typedef enum { M1 = 0, M2, M3, M4, A } BLS_MOVE_TYPE;


	/////////////////////////////////////////////////////////////////////////
	// GRAPH SOLVER CLASS DECLARATION
	/////////////////////////////////////////////////////////////////////////
	class CGraphSolver
	{
		//------------------------------------------------------
		// METHODS
		//------------------------------------------------------
	public:
		CGraphSolver(void);
		~CGraphSolver(void);

		void Initialize(Graph* pGraph, HJ_GRAPH_SOLVER_TYPE solverType);
		void Initialize(Graph* pGraph, HJ_GRAPH_SOLVER_TYPE solverType, Mat2D &matQ);
		void Finalize(void);
		void Clear(void);
		void SetInitialSolutions(std::deque<VertexSet> &initialSolution);

		stGraphSolvingResult* Solve(int timelimit = 0);
		stGraphSolvingResult GetResult(void) { return this->stResult_; };

	private:
		void RunBLS(int timelimit = 0);
		void RunILS(void);
		void RunAMTS(void);
		void RunMCGA(void);

		// miscellaneous
		bool   CheckSolutionExistance(const VertexSet &vertexSet, double SumWeights);
		double SumWeights(VertexSet &vertexSet);
		double SumWeightsQP(VertexSet &vertexSet);
		double SumWeightsQP(VertexSet &vertexSet, GraphVertex *targetVertex);
		void   PrintLog(size_t type);

		// BLS related
		bool   BLS_SetInitialSolutions(std::deque<VertexSet> &initialSolution);
		bool   BLS_InsertSolution(const VertexSet &solution, const double solutionScore);
		double BLS_GenerateInitialSolution(bool bGreedy = true);
		double BLS_BestLocalMove(void);
		double BLS_Perturbation(double L, size_t w, double alphaR, double alphaS);
		double BLS_PerturbDirected(double L);
		double BLS_PerturbRandom(double L, double alpha);
		double BLS_VertexRemove(GraphVertex *vertexRemove);
		double BLS_VertexInsert(GraphVertex *vertexInsert);
		double BLS_VertexInsertM4(size_t vertexIdxInOC);
		size_t BLS_GetMoveProhibitionNumber(size_t nIter);

		//------------------------------------------------------
		// VARIABLES
		//------------------------------------------------------
	private:
		Graph *pGraph_;
		bool  bInit_;
		bool  bHasInitialSolution_;
		bool  bTimeoutSet_;
		HJ_GRAPH_SOLVER_TYPE solverType_;
		stGraphSolvingResult stResult_;

		// BLS related
		size_t    BLS_nIter;
		VertexSet BLS_C;               // solution vertices
		VertexSet BLS_PA;              // not in C, but connected to all v in C
		VertexSet BLS_OC;              // not in C
		std::deque<VertexPair> BLS_OM; // <v in OC, u in C>	: except u, v is connected to all vertices in C

									   // QP related
		bool  bQP_;
		Mat2D matQ_; // NOTICE: diagonals must be zero
	};
}

//()()
//('')HAANJU.YOO
