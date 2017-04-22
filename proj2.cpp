#include <iostream>
#include <list>
#include <queue>
#include <algorithm>
using namespace  std;

/*************************
 *    ASA - LEIC-TP      *
 * 		Grupo 32		 *
 * Andre Fonseca 84698	 *
 * Leonor Loureiro 84736 *
 *************************/


/**
 * Class representation of an Edge
 * An Edge or Arc is a connection between 2 vertices.
 */
class Edge{
	public:
    	int predecessor;	// vertex from which edge originates from
    	int successor;		// vertex to which the edge is going to
    	int weight;         // weight of edge
		/**
		 * Default consctructor
		 * @param pre the predecessor of the edge
		 * @param suc the successor of the edge
		 */
        Edge(int pre, int suc, int wei);
        int getPredecessor();
        int getSucessor();
        int getWeight();
};

struct cheapestEdge
{
    inline bool operator() (const Edge& struct1, const Edge& struct2)
    {
        return (struct1.weight < struct2.weight);
    }
};

/**
 * Structure used to model relations between 2 cities
 * This implementation of a directional graph is made through an adjacency list
 */
class Graph{
	private:
		int nrOfVertices;           // number of vertexs
		vector<Edge> edges;  		// vector of edges
		int *aeroCost;              // cost of every aeroport
		int *key;					// cost of the edge that led to the vertex
		int *predecessor;			// predecessor vertex

		void sortEdges();
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
		int kruskalMST();
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
	parent = new int[n+1];
	rank = new int[n+1];
	for(int i = 0; i <= n; i++){
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

Edge::Edge(int pre, int suc, int wei){
	predecessor = pre;
	successor = suc;
	weight = wei;
}

int Edge::getPredecessor(){
	return predecessor;
}

int Edge::getSucessor(){
	return successor;
}

int Edge::getWeight(){
	return weight;
}

Graph::Graph(int v){
	nrOfVertices = v;
	aeroCost = new int[v];
	for(int i = 0; i < v; i++)
		aeroCost[i] = -1;			// can't build aeroport
}

void Graph::setVertexCost(int v, int cost){
	aeroCost[v] = cost;
}

void Graph::insertEdge(Edge e){
	edges.push_back(e);
}

int Graph::kruskalMST(){

	//Weight of the past of miminum cost
	int mstWeight = 0;
	int nrEdgesMST = 0;

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

			key[(*i).getSucessor()] = (*i).getWeight(); 
			predecessor[(*i).getSucessor()] = (*i).getPredecessor();
			cout << u << " - " << v << endl;

			nrEdgesMST++;

			mstWeight += (*i).getWeight();

			sets.merge(u,v);
		}

	}
	cout << mstWeight << endl;
	return nrEdgesMST;
}

int main(int argc, char const *argv[])
{
	Graph g(9);
	g.insertEdge(Edge(0,1,4));
	g.insertEdge(Edge(0,7,8));
	g.insertEdge(Edge(1,2,8));
	g.insertEdge(Edge(1,7,11));
	g.insertEdge(Edge(2,3,7));
	g.insertEdge(Edge(2,8,2));
	g.insertEdge(Edge(2,5,4));
	g.insertEdge(Edge(3,4,9));
	g.insertEdge(Edge(3,5,14));
	g.insertEdge(Edge(4,5,10));
	g.insertEdge(Edge(5,6,2));
	g.insertEdge(Edge(6,7,1));
	g.insertEdge(Edge(6,8,6));
	g.insertEdge(Edge(7,8,7));
	cout << g.kruskalMST() <<endl;
	return 0;
}