#include "GraphUtils.h"

GraphUtils* GraphUtils::New(int vertices){
   return new GraphUtils(vertices);
}

GraphUtils::GraphUtils(int vertices){
   this->outputfile = std::ofstream( "graph.output" );
   this->MISfile = std::ofstream("MISinfo.data");
   this->MISets = std::ofstream("MIS.data");
   
   this->nvertices = vertices;
   this->adj_cost_matrix = new double*[vertices]();
   this->M = new double*[vertices]();
   this->clique_adj = new short*[vertices]();

   for(int i = 0; i < vertices; i++){
      this->adj_cost_matrix[i] = new double[vertices]();
      this->M[i] = new double[vertices]();
      this->clique_adj[i] = new short[vertices]();
   }
}