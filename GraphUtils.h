#pragma once
#include <list>
#include <string>
#include <fstream>
#include <iostream>

class GraphUtils{
public:
   static GraphUtils* New(int vertices);

   void printAdjMatrix();

   int nvertices;				//will be the actual number of vertices in the graph
   double** adj_cost_matrix;//cost and adjacency matrix (combined)
   double** M;//matrix for determining number of spanning trees
   short** clique_adj;	//will be adjacency matrix of complement graph
   std::ofstream outputfile;//file to which results are sent
   std::ofstream MISfile;	//writes the information from each step of the MIS algorithm
   std::ofstream MISets;		//writes just the MISs and cliques

protected:
   GraphUtils(int vertices);
};