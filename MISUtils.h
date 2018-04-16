#pragma once

#include <iostream>

#include "GraphUtils.h"

class MISUtils{
public:
   MISUtils(GraphUtils* gu);

   bool equation_3_8(char MIS_or_clique);

   int step2(char directed, char MIS_or_clique);
   int step4(char directed, char MIS_or_clique);
   int step5(char MIS_or_clique);

   int MIS(char directed);
   int H_MIS(char directed);
   int clique(char directed);

private:
   void init();

   GraphUtils* gu;

   int** Qminus;
   int** Qplus;
   int** S;
   int* X;

   int k, maximum; //state of k and max_ind(the maximum size of an MIS)
   double start, finish, duration;	//will be used for determining time of calculation
   
   bool eq38 = false;		//says whether or not equation 3.8 is satisfied
   bool done = false;
   bool print_MIS = false;
};