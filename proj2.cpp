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
#define AERO 1 // used to reference the Airspace connection/edge
#define ROAD 2 // used to reference a road connection/edge




/**
 * Class representation of an Edge
 * An Edge or Arc is a connection between 2 vertices.
 */
class Edge{
	private:
    	int predecessor;	// vertex from which edge originates from
    	int successor;		// vertex to which the edge is going to
    	int weight;         // weight of edge
		int type;			// type: Airport connection or road
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


/**
 * Comparison between 2 edges
 * @return True if Edge1 is Smaller or Lesser than Edge 2
 */
bool sortEdges(const Edge& edge1, const Edge& edge2)
{
	if(edge1.getWeight() == edge2.getWeight()){
		if(edge1.getType() == edge2.getType()){
			return false;
		}
		return (edge1.getType()==ROAD);
	}
    return (edge1.getWeight() < edge2.getWeight());
}



/**
 * Structure used to model relations between severeral Cities
 * This implementation of a undirected graph is made through it's Edges
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
		pair<int,pair<int,int> > kruskalModified();
};


/**
 * Disjoint Sets CLass
 * 2 sets are disjoints if there is no common element between them.
 * This class models basic union, find and merge functions between disjoint sets.
 */
class DisjointSets{
	private:
		int *parent; 			//the parent or representative of a set
		int *rank;				//rank (number of merges) of a given set
		int nrElements;			//number of elements in a given set
	public:

		/**
		 * Default constructor of DisjointSets class
		 * It receives the number of Disjoint sets with which the object starts
		 */
		DisjointSets(int size);

		/**
		 * Deconstructor
		 * destroys alocated memory for the dejoinSets
		 */
		~DisjointSets();

		/**
		 * Finds the parent of a node
		 * @param u the node of which we want the parent
		 */
		int find(int u);

		/**
		 * Merges two disjoint sets into one
		 * @param u disjoint set one
		 * @param v disjoint set two
		 */
		void merge(int u, int v);
};





DisjointSets::DisjointSets(int size){
	nrElements = size;
	parent = new int[size+1];
	rank = new int[size+1];
	for(int i = 0; i < size+1; i++){
		parent[i] = i;
		rank[i] = 0;
	}
}

DisjointSets::~DisjointSets(){
	delete [] parent;
	delete [] rank;
}

int DisjointSets::find(int u){
	if(u != parent[u]){
		parent[u] = find( parent[u] );
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

int Edge::getPredecessor() const{ return predecessor; }

int Edge::getSucessor() const{ return successor; }

int Edge::getWeight() const{ return weight; }

int Edge::getType() const{ return type; }






Graph::Graph(int v){
	nrOfVertices = v;
}

void Graph::insertEdge(Edge e){
	edges.push_back(e);
}

pair<int,pair<int,int> > Graph::kruskalModified(){

	//Weight of the past of miminum cost
	int mstWeight = 0;
	int nrRoadsMST = 0;
	int nrAerosMST = 0;

	//Orders the edges by weight
	sort(edges.begin(), edges.end(), sortEdges);

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
			}else{
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
 *
 *
 * Finds the number of airports and roads needed to be built between given cities
 * in the cheapest possible way and the minimum airports required.
 */
int inputProcess(){

	int potentialAirports = 0;
	int potentialRoads = 0;
	int numberOfCities = 0;
	int airportVertex = 0;
	int buildingCost = 0;
	int airportID = 0;
	int city1 = 0;
	int city2 = 0;
	int i = 0;

	std::cin >> numberOfCities;

	int airports[numberOfCities];
	int roadsInCity[numberOfCities];

	Graph roadGraph (numberOfCities);
	Graph fullGraph (numberOfCities+1);

//initialize all cities with no airports and 0 roads
	for( i = 0; i < numberOfCities; i++){
		airports[i] = -1;
		roadsInCity[i] = 0;
	}



// airports process
	std::cin >> potentialAirports;
	airportVertex = numberOfCities;

	for(i = 0; i < potentialAirports; i++){
		std::cin >> airportID >> buildingCost;

		fullGraph.insertEdge(Edge (airportID-1, airportVertex, buildingCost, AERO));
		airports[airportID-1] = buildingCost;
	}

// Roads process
	std::cin >> potentialRoads;

	for(i = 0; i < potentialRoads; i++){
		std::cin >> city1 >> city2 >> buildingCost;

		fullGraph.insertEdge(Edge (city1-1, city2-1, buildingCost, ROAD));
		roadGraph.insertEdge(Edge (city1-1, city2-1, buildingCost, ROAD));

		roadsInCity[city1-1]++;
		roadsInCity[city2-1]++;
	}


	for(i = 0; i < numberOfCities; i++){
		if(roadsInCity[i] == 0){
			if(airports[i] == -1){
				std::cout << "Insuficiente" << '\n';
				return 0;
			}
			//only do the complete Graph
		}
	}

// Runs MST on each graph in order to choose which is cheapest
	pair<int,pair<int,int> > roadsOutput = roadGraph.kruskalModified();
	pair<int,pair<int,int> > fullOutput = fullGraph.kruskalModified();

// Verifies if the output was insuficient or not and prints it
	if (roadsOutput.first > fullOutput.first || roadsOutput.second.second != numberOfCities -1){
		if(fullOutput.second.first + fullOutput.second.second != numberOfCities){
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
	inputProcess();
	return 0;
}
