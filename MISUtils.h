#pragma once

#include <iostream>
#include <chrono>

#include "GraphUtils.h"

class MISUtils{
public:
   MISUtils(GraphUtils* gu);

   bool equation_3_8(char MIS_or_clique);

   int step2(char directed, char MIS_or_clique);
   int step4(char directed, char MIS_or_clique);
   int step5(char MIS_or_clique);

   int step2_h(char directed, char MIS_or_clique);
   int step4_h(char directed, char MIS_or_clique);
   int step5_h(char MIS_or_clique);
   bool equation_3_8_h(char MIS_or_clique);

   int MIS(char directed);
   int H_MIS(char directed);
   int clique(char directed);

   void printQPlus();
   void printQMinus();
   void printS();
   void printD();
   void printAdjMatrix(int d);

private:
   void init();

   GraphUtils* gu;

   int** Qminus;
   int** Qplus;
   int** S;
   int** D;
   int* X;

   int* cost;
   int * fS;
   int fCost;
   
   int*** adj_matrix;

   int k, maximum; //state of k and max_ind(the maximum size of an MIS)

   bool eq38 = false;		//says whether or not equation 3.8 is satisfied
   bool done = false;
   bool print_MIS = false;
};