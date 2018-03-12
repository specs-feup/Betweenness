#include "Betweenness.h"


Betweenness::Betweenness(WeightedDirectedGraph &graph)
{
	this->graph = &graph;
}

BetweennessResult Betweenness::Calculate(int startVertex, int endVertex)
{
	cout << "start " << startVertex << ", end " << endVertex << endl;
	int vertices = graph->GetVertices();
	int edges = graph->GetEdges();
	double *betweenness = new double[vertices];
	double *edgeBetweenness = new double[edges];
	for (size_t i = 0; i < vertices; i++)
	{
		betweenness[i] = 0;
	}
	for (size_t i = 0; i < edges; i++)
	{
		edgeBetweenness[i] = 0;
	}

	priority_queue<KeyValuePair, vector<KeyValuePair>, greater<KeyValuePair>> Q;
	stack<int> S;
	double *dist = new double[vertices];
	list<int> *Pred;	//list of predecessors is faster
	double *sp = new double[vertices];
	double *delta = new double[vertices];
	bool *isInStack = new bool[vertices];

	for (int s = startVertex; s < endVertex; s++)
	{
		//Single source shortest-paths problem
		{
			//Initialization
			{
				Pred = new list<int>[vertices];
				for (int i = 0; i < vertices; i++)
				{
					dist[i] = DBL_MAX;
					sp[i] = 0;
					isInStack[i] = false;
				}
				dist[s] = 0;
				sp[s] = 1;
				Q.push(KeyValuePair(dist[s], s));
			}

			while (!Q.empty())
			{
				KeyValuePair pair = Q.top();
				int v = pair.Value;
				Q.pop();

				//this part of code should replace updating in priority queue which is expensive?
				if (isInStack[v] == false)
				{
					isInStack[v] = true;
					S.push(v);


					vector<Edge> adjEdges = graph->GetAdjacentVertices(v);	//vector of edges is faster than list
					for (auto &adjEdge : adjEdges)
					{
						double edgeDistance = adjEdge.GetWeight();

						int w = adjEdge.GetOutput();

						//path discovery, shorter path found?
						if (dist[w] > dist[v] + edgeDistance)
						{
							dist[w] = dist[v] + edgeDistance;
							Q.push(KeyValuePair(dist[w], w));
							sp[w] = 0;
							Pred[w].clear();
						}

						//path counting, edge(v,w) on a shortest path?
						if (Utils::AreEqual(dist[w], dist[v] + edgeDistance))
						{
							//sp[w] = sp[w] + edgeDistance*sp[v];
							sp[w] = sp[w] + sp[v];
							Pred[w].push_back(adjEdge.GetId());//if stores edge ids instead of predecessors vertices, so multiplicity will be possible
						}
					}
				}
			}
		}

		//Accumulation
		{
			double c = 0;

			for (int v = 0; v < vertices; v++)
			{
				delta[v] = 0;
			}

			while (S.size() != 0)
			{
				int w = S.top();
				S.pop();

				for (auto edgeId : Pred[w])//go through edges on shortest paths, so we would know to which edge we should add edge betweenness value and to identify source vertex
				{
					int v = graph->GetEdge(edgeId)->GetInput();//predecessor of w on shortest path from s
					c = (sp[v] / sp[w] * (1 + delta[w]));//extention about source vertex importance and destination vertex importance
					
					edgeBetweenness[edgeId] = edgeBetweenness[edgeId] + c;	//save edge betweenness
					delta[v] = delta[v] + c;
				}
				if (w != s)
				{
					betweenness[w] = betweenness[w] + delta[w];
				}
			}
		}

		delete[] Pred;
	}

	delete[] dist;
	delete[] sp;
	delete[] delta;
	delete[] isInStack;
	return BetweennessResult(betweenness, edgeBetweenness);
}

BetweennessResult Betweenness::CalculateOpenMP(int startVertex, int endVertex, int threads)
{
	//Each thread is using serial version
	int vertices = graph->GetVertices();
	int edges = graph->GetEdges();
	double *betweenness = new double[vertices];
	double *edgeBetweenness = new double[edges];

	omp_lock_t *lockArray = new omp_lock_t[vertices + edges];//locks for vertices and edges betweenness arrays
	for (size_t i = 0; i < vertices; i++)
	{
		betweenness[i] = 0;
		omp_init_lock(&lockArray[i]);
	}
	for (size_t i = 0; i < edges; i++)
	{
		edgeBetweenness[i] = 0;
		omp_init_lock(&lockArray[i + vertices]);
	}

#pragma omp parallel shared(betweenness, lockArray) num_threads(threads)
	{
		int vert = endVertex - startVertex;
		int tid = omp_get_thread_num();

		int size = (int)ceill((double)vert / threads);
		int start = tid * size;
		int end = start + size;
		if (tid == threads - 1)
		{
			end = vert;
		}
		//for low number of vertices if start is higher than max vertex id then this method will not modify betweenness
		BetweennessResult result = Calculate(startVertex + start, startVertex + end);
		double *subBetweenness = result.VertexBetweenness;
		double *subEdgeBetweenness = result.EdgeBetweenness;

		for (size_t i = 0; i < vertices; i++) {
			omp_set_lock(&lockArray[i]);
			betweenness[i] = betweenness[i] + subBetweenness[i];
			omp_unset_lock(&lockArray[i]);
		}

		for (size_t i = 0; i < edges; i++)
		{
			omp_set_lock(&lockArray[i + vertices]);
			edgeBetweenness[i] = edgeBetweenness[i] + subEdgeBetweenness[i];
			omp_unset_lock(&lockArray[i + vertices]);
		}

		delete[] subBetweenness;
		delete[] subEdgeBetweenness;
	}

	for (size_t i = 0; i < vertices + edges; i++)
	{
		omp_destroy_lock(&lockArray[i]);
	}
	delete[] lockArray;

	return BetweennessResult(betweenness, edgeBetweenness);
}

BetweennessResult Betweenness::CalculateMpi(int startVertex, int endVertex, int threads, int chunkSize)
{
	//MASTER
	//while work not done
	//	receive request from slave
	//	if work done
	//		receive result
	//	else
	//		send work or if no more work then send termination message

	//SLAVE
	//while
	//	send request to master
	//	receive work
	//	if no more work
	//		terminate
	//	else
	//		calculate result
	//		send message to master to receive result
	//		send result to master

	int vertices = graph->GetVertices();
	int edges = graph->GetEdges();
	int size;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	enum MessageTag { wantWork, work, result };
	int jobStart = startVertex;
	int jobEnd = endVertex;
	int endsSended = 0;

	int workThreads = threads;

	if (rank == 0)	//MASTER
	{
		double *betweenness = new double[vertices];
		double *edgeBetweenness = new double[edges];
		for (size_t i = 0; i < vertices; i++)
		{
			betweenness[i] = 0;
		}
		for (size_t i = 0; i < edges; i++)
		{
			edgeBetweenness[i] = 0;
		}

		if (size == 1)
		{
			if (workThreads == 1)
			{
				BetweennessResult result = Calculate(jobStart, jobEnd);
				betweenness = result.VertexBetweenness;
				edgeBetweenness = result.EdgeBetweenness;
			}
			else
			{
				BetweennessResult result = CalculateOpenMP(jobStart, jobEnd, workThreads);
				betweenness = result.VertexBetweenness;
				edgeBetweenness = result.EdgeBetweenness;
			}
		}
		else
		{
			while (true)
			{
				//1 slave wants work, 0 slave wants send result
				int wantWork;
				MPI_Request recvRequest;
				MPI_Irecv(&wantWork, 1, MPI_INT, MPI_ANY_SOURCE, MessageTag::wantWork, MPI_COMM_WORLD, &recvRequest);
				MPI_Status status;
				MPI_Wait(&recvRequest, &status);
				int slaveRank = status.MPI_SOURCE;
				int flag = status.MPI_TAG;

				if (wantWork == 1)
				{
					int work[2];

					if (jobStart < jobEnd)
					{
						work[0] = jobStart;	//from which vertex, -1 if no more work
						if (jobStart + chunkSize <= jobEnd)
						{
							work[1] = jobStart + chunkSize;	//to which vertex
						}
						else
						{
							//last chunk of work sended
							work[1] = jobStart + jobEnd % chunkSize;
						}
					}
					else
					{
						work[0] = -1;
						//ended messages counter
						endsSended++;
					}


					MPI_Request sendRequest;
					MPI_Isend(&work, 2, MPI_INT, slaveRank, MessageTag::work, MPI_COMM_WORLD, &sendRequest);	//no need to wait for receive
					MPI_Wait(&sendRequest, &status);
					jobStart += chunkSize;

					if (endsSended == size - 1)
					{
						cout << "master ended" << endl;
						break;
					}
				}
				else
				{
					double *subBetweenness = new double[vertices + edges];
					MPI_Irecv(subBetweenness, vertices + edges, MPI_DOUBLE, slaveRank, MessageTag::result, MPI_COMM_WORLD, &recvRequest);
					MPI_Wait(&recvRequest, &status);

					//save subBetweenness, maybe in new thread
					for (size_t i = 0; i < vertices; i++)
					{
						betweenness[i] = betweenness[i] + subBetweenness[i];
					}

					for (size_t i = 0; i < edges; i++)
					{
						edgeBetweenness[i] = edgeBetweenness[i] + subBetweenness[i + vertices];
					}

					delete[] subBetweenness;
				}
			}
		}

		return BetweennessResult(betweenness, edgeBetweenness);
	}
	else	//SLAVE
	{
		int wantWork;
		while (true)
		{
			wantWork = 1;
			MPI_Request sendRequest;
			MPI_Isend(&wantWork, 1, MPI_INT, 0, MessageTag::wantWork, MPI_COMM_WORLD, &sendRequest);
			MPI_Status status;
			MPI_Wait(&sendRequest, &status);

			int work[2];
			MPI_Request recvRequest;
			MPI_Irecv(&work, 2, MPI_INT, 0, MessageTag::work, MPI_COMM_WORLD, &recvRequest);
			MPI_Wait(&recvRequest, &status);

			if (work[0] == -1)
			{
				//no more work
				break;
			}
			else
			{
				double *results = new double[vertices + edges];
				double *subBetweenness;
				double *subEdgeBetweenness;
				if (workThreads == 1)
				{
					BetweennessResult result = Calculate(work[0], work[1]);
					subBetweenness = result.VertexBetweenness;
					subEdgeBetweenness = result.EdgeBetweenness;
				}
				else
				{
					BetweennessResult result = CalculateOpenMP(work[0], work[1], workThreads);
					subBetweenness = result.VertexBetweenness;
					subEdgeBetweenness = result.EdgeBetweenness;
				}

				//merge two results
				for (size_t i = 0; i < vertices; i++)
				{
					results[i] = subBetweenness[i];
				}

				for (size_t i = 0; i < edges; i++)
				{
					results[i + vertices] = subEdgeBetweenness[i];
				}


				wantWork = 0;
				MPI_Isend(&wantWork, 1, MPI_INT, 0, MessageTag::wantWork, MPI_COMM_WORLD, &sendRequest);
				MPI_Wait(&sendRequest, &status);
				MPI_Isend(results, vertices + edges, MPI_DOUBLE, 0, MessageTag::result, MPI_COMM_WORLD, &sendRequest);
				MPI_Wait(&sendRequest, &status);
				delete[] results;
				delete[] subBetweenness;
				delete[] subEdgeBetweenness;
			}
		}
	}

	return BetweennessResult(nullptr, nullptr);
}



Betweenness::~Betweenness()
{
}
