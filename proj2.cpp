#include <iostream>
#include <list>
#include <queue>
#include <algorithm>
#include <utility>
using namespace  std;

/*************************
 *    ASA - LEIC-TP      *
 * 		Grupo 32		 *
 * Andre Fonseca 84698	 *
 * Leonor Loureiro 84736 *
 *************************/
#define AERO 1
#define ROAD 2



typedef  pair<int, int> nrConnections;

/**
 * Class representation of an Edge
 * An Edge or Arc is a connection between 2 vertices.
 */
class Edge{
	private:
    	int predecessor;	// vertex from which edge originates from
    	int successor;		// vertex to which the edge is going to
    	int weight;         // weight of edge
		int type;
	public:
		/**
		 * Default consctructor
		 * @param pre the predecessor of the edge
		 * @param suc the successor of the edge
		 */
        Edge(int pre, int suc, int wei, int typ);
        int getPredecessor() const;
        int getSucessor() const;
        int getWeight() const;
		int getType() const;
};

struct cheapestEdge
{
    inline bool operator() (const Edge& struct1, const Edge& struct2)
    {
		if(struct1.getWeight() == struct2.getWeight()){
			return (struct1.getType()==ROAD);
		}
        return (struct1.getWeight() < struct2.getWeight());
    }
};

/**
 * Structure used to model relations between 2 cities
 * This implementation of a directional graph is made through an adjacency list
 */
class Graph{
	private:
		int nrOfVertices;         // number of vertexs
		vector<Edge> edges;  		 // vector of edges
	public:

		/**
		 * Default Constructor
		 * @param v the size(number of vertex) of the graph
		 * @param cost the cost of building a aeroport
		 */
		Graph(int v);

		/**
		 * Defines the cost of building an aeroport
		 * @param v the vertex
		 * @param cost the cost of building a aeroport
		 */
		void setVertexCost(int v, int cost);

		/**
		 * Inserts and edge to the graph
		 * @param e Edge to be added to the graph
		 */
		void insertEdge(Edge e);

		/**
		 * Based on Kruskal's algorithm for finding the minimum
		 * spanning tree.
		 */
		pair<int,nrConnections> kruskalMST();
};

class DisjointSets{
	private:
		int *parent, *rank;
		int nrElements;
	public:
		DisjointSets(int n);
		/**
		 * Finds the parent of a node
		 * @param u the node of which we want the parent
		 */
		int find(int u);
		void merge(int u, int v);
};

DisjointSets::DisjointSets(int n){
	nrElements = n;
	parent = new int[n];
	rank = new int[n];
	for(int i = 0; i < n; i++){
		parent[i] = i;
		rank[i] = 0;
	}
}

int DisjointSets::find(int u){
	if(u != parent[u]){
		parent[u] = find(parent[u]);
	}
	return parent[u];
}

void DisjointSets::merge(int u, int v){
	u = find(u);
	v =  find(v);

	//Make the tree of smaller rank a
	//subtree of the other tree
	if (rank[u] > rank[v]){
		parent[v] = u;
	}
	else{
		parent[u] = v;
	}

	if(rank[u] == rank[v]){
		rank[v]++;
	}
}

Edge::Edge(int pre, int suc, int wei, int typ){
	predecessor = pre;
	successor = suc;
	weight = wei;
	type = typ;
}

int Edge::getPredecessor() const{
	return predecessor;
}

int Edge::getSucessor() const{
	return successor;
}

int Edge::getWeight() const{
	return weight;
}

int Edge::getType() const{
	return type;
}

Graph::Graph(int v){
	nrOfVertices = v;
}


void Graph::insertEdge(Edge e){
	edges.push_back(e);
}

pair<int,nrConnections> Graph::kruskalMST(){

	//Weight of the past of miminum cost
	int mstWeight = 0;
	int nrRoadsMST = 0;
	int nrAerosMST = 0;

	//Orders the edges by weight
	sort(edges.begin(), edges.end(), cheapestEdge());

	DisjointSets sets(nrOfVertices);

	for(vector<Edge>::iterator i = edges.begin(); i != edges.end(); i++){
		int u = (*i).getPredecessor();
		int v = (*i).getSucessor();

		int set_u = sets.find(u);
		int set_v = sets.find(v);

		// Check if u and v belong to the same set
		// and the edge selected creates a cycle
		if(set_u != set_v){

			if((*i).getType() == ROAD){
				nrRoadsMST++;
			}
			else{
				nrAerosMST++;
			}

			mstWeight += (*i).getWeight();

			sets.merge(u,v);
		}
	}


	return make_pair(mstWeight,make_pair(nrAerosMST,nrRoadsMST));
}



/*
 * Input shape:
 * int: number of cities
 * int: number of potential airports to build
 *
 * for every potential airport:
 * 	int: airport's number
 * 	int: airport's cost
 *
 * int: number with potential roads to build
 *
 * for every potential road
 * 	int: city1
 * 	int: city2
 * 	int: road's building cost
 */
int inputProcess(){

	int potentialAirports;
	int potentialRoads;
	int numberOfCities;
	int airportVertex;
	int buildingCost;
	int airportID;
	int city1;
	int city2;
	int i;

	std::cin >> numberOfCities;

	int airports[numberOfCities];
	int roadsInCity[numberOfCities];

	Graph roadGraph (numberOfCities);
	Graph fullGraph (numberOfCities+1);

//initialize all cities with no airports
	for( i = 0; i < numberOfCities; i++)
		airports[i] = -1;



// airports
	std::cin >> potentialAirports;
	airportVertex = numberOfCities;

	for(i = 0; i < potentialAirports; i++){
		std::cin >> airportID >> buildingCost;

		Edge airportEdge(airportID-1, airportVertex, buildingCost, AERO);
		fullGraph.insertEdge(airportEdge);
		airports[airportID] = buildingCost;
	}

// Roads
	std::cin >> potentialRoads;

	for(i = 0; i < potentialRoads; i++){
		std::cin >> city1 >> city2 >> buildingCost;

		Edge roadEdge(city1-1, city2-1, buildingCost, ROAD);
		fullGraph.insertEdge(roadEdge);
		roadGraph.insertEdge(roadEdge);

		roadsInCity[city1]++;
		roadsInCity[city2]++;
	}

	for(i = 0; i < numberOfCities; i++){
		if(roadsInCity[i] == 0){
			if(airports[i] == -1){
				std::cout << "Insuficiente" << '\n';
			}
			//only do the complete Graph
		}
	}

// Runs MST on each graph in order to choose which is cheapest
	pair<int,nrConnections> roadsOutput = roadGraph.kruskalMST();
	pair<int,nrConnections> fullOutput = fullGraph.kruskalMST();

// Verifies if the output was insuficient or not and prints it
	if (roadsOutput.first > fullOutput.first || roadsOutput.second.second < numberOfCities -1){
		if(fullOutput.second.first + fullOutput.second.second < numberOfCities){
			std::cout << "Insuficiente" << '\n';
		}else{
			std::cout << fullOutput.first << '\n';
			std::cout << fullOutput.second.first <<' '<< fullOutput.second.second << '\n';
		}
	}else{
		std::cout << roadsOutput.first << '\n';
		std::cout << roadsOutput.second.first <<' '<< roadsOutput.second.second << '\n';
	}


	return 0;
}



int main(int argc, char const *argv[])
{
	//int *aeroCost;          // cost of every aeroport
	// Graph g(9);
	// g.insertEdge(Edge(0,1,4,0));
	// g.insertEdge(Edge(0,7,8,0));
	// g.insertEdge(Edge(1,2,8,0));
	// g.insertEdge(Edge(1,7,11,0));
	// g.insertEdge(Edge(2,3,7,0));
	// g.insertEdge(Edge(2,8,2,0));
	// g.insertEdge(Edge(2,5,4,0));
	// g.insertEdge(Edge(3,4,9,0));
	// g.insertEdge(Edge(3,5,14,0));
	// g.insertEdge(Edge(4,5,10,0));
	// g.insertEdge(Edge(5,6,2,0));
	// g.insertEdge(Edge(6,7,1,0));
	// g.insertEdge(Edge(6,8,6,0));
	// g.insertEdge(Edge(7,8,7,0));
	// cout << g.kruskalMST().first <<endl;
	//
	inputProcess();
	return 0;
}
