//---------------------begin MIS/clique 
//************************************************************************
//*
//*	Title: 			MIS
//*	Date:  			July 2000
//*	Version:			1.1
//*	
//*	Project:			Development of C++ graph package
//*
//*	Author: 			Matthew N. Johnson
//*	Coordinator: 		Professor Gary Lamont
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Function description:	This function uses the algorithm as described in 
//*					Christofides, page 35, to find the maximal 
//*					independent sets.  In addition to the output that 
//*					is displayed on the screen, the state of the 
//*					information at each successive step and iteration 
//*					is written to the file 'MISinfo.data'. However, 
//*					note that this file can be rather large. For 
//*				     example, for a 13-vertex, directed graph, the file 
//*					was 23 pages in length when sent to Microsoft 
//*					Word.  This size grows almost exponentially, as a 
//*					20-vertex, directed graph created a 400+ page 
//*					file.  If you wish to view only the maximal 
//*					independent sets, they are printed to the file 
//*				     called 'MIS.data'.  This may be more convenient if 
//*					your only concern is the final answer.
//*
//*	Algorithm:			As described in Christofides, p. 35. 
//*
//*	Order of analysis:	O((3^(n/3))*gu->nvertices^2)
//*
//*	Globals:			adjacency matrix, number of vertices
//*
//*	Local Variables:		The sets of Q+, Q-, S, and X
//*
//*	Calling Modules:		Main program
//*
//************************************************************************
//this function is common to both MIS and clique 
#include <time.h>

#include "MISUtils.h"

bool MISUtils::equation_3_8(char MIS_or_clique){
   short a, b, sum=0, qminus_size=0, qplus_size=0;
   for (a=0; a<gu->nvertices; a++){
      qminus_size += Qminus[k][a];
   }
   if (qminus_size==0){  		//3.8 unsatisfied if Q- is empty
      eq38 = false;			//S)atisfied, U)nsatisfied 
   }
   else{	//Q- is not empty			//else #1
      for (a=0; a<gu->nvertices; a++){
         qplus_size+=Qplus[k][a];
      }
      if (qplus_size==0){
         eq38 = true; 	//3.8 satisfied if Q+ is empty
      }
      else{		//Q+ is not empty		//else #2
         eq38 = false;   //assume 3.8 is unsatisfied
         sum=0;
         for (a=0; a<gu->nvertices; a++){	//for #1
            if (Qminus[k][a] != 0){ //i.e. l+1 element of Q-
               for (b=0; b<gu->nvertices; b++){	//for #2
                  if(MIS_or_clique=='c'){			//c means doing clique calculation
                     if ((gu->clique_adj[a][b]!=0)&&(Qplus[k][b]!=0)){
                        sum++;
                     }
                  }
                  else{			//doing mis calculation
                     if((gu->adj_cost_matrix[a][b]!=0)&&(Qplus[k][b]!=0)){
                        sum++; // equation 3.8 not satisfied
                     }
                  }
               }		// end for #2
            }			//end if 		
         }				//end for #1
         if (sum == 0) {
            eq38 = true; 
         }
      }					//end else #2
   }				//end else #1
   if (print_MIS){
      gu->MISfile<<"step 3:";
      if (eq38){
         gu->MISfile<<"  equation 3.8 satisfied, go to step 5"<<std::endl;
      }
      else{
         gu->MISfile<<"  equation 3.8 not satisfied, go to step 4"<<std::endl;
      }
   }
   return eq38;
}

//this function is used by MIS and clique calculation
int MISUtils::step2(char directed, char MIS_or_clique){
   short a, xk = gu->nvertices;
   for (a=0; a<gu->nvertices; a++){
      if ((Qplus[k][a]!=0)&&(Qplus[k][a]<xk)){
         xk = Qplus[k][a];//find the lowest numbered vertex in Qplus
      }
      Qplus[k+1][a]=Qplus[k][a];//set Q+ (k+1) to Q+ k, same for Q- and S
      Qminus[k+1][a]=Qminus[k][a];
      S[k+1][a]=S[k][a];
   }
   X[k]=S[k+1][xk-1]= xk;//add xk to the set X of vertices and to S
   for (a=0;a<gu->nvertices;a++){
      if (MIS_or_clique=='c'){
         if (directed=='Y'){
            if(gu->clique_adj[a][(X[k])-1]==1) {
               Qminus[k+1][a]=Qplus[k+1][a]=0;	//remove correspondences of xk
            }
         }
         if(gu->clique_adj[(X[k])-1][a]==1) {
            Qminus[k+1][a]=Qplus[k+1][a]=0;	//remove correspondences of xk
         }
      }
      else{
         if (directed=='Y'){
            if (gu->adj_cost_matrix[a][(X[k]-1)]!=0){
               Qminus[k+1][a]=Qplus[k+1][a]=0;//remove correspondences of xk(to and from)
            }
         }
         if(gu->adj_cost_matrix[(X[k])-1][a]!=0) {
            Qminus[k+1][a]=Qplus[k+1][a]=0;		//remove correspondences
         }
      }
   }
   Qplus[k+1][(X[k])-1]=0; //remove xk from Q+(k+1)
   k++;						//increment k
   if (print_MIS)		//print the values of the sets, if requested
   {
      gu->MISfile<<"step 2:"<<std::endl<<"  xi"<<k<<" = "<<xk<<std::endl<<"  Q+"<<k<<":";
      for (a=0; a<gu->nvertices; a++){
         if (Qplus[k][a]!=0){
            gu->MISfile.width(3);
            gu->MISfile<<Qplus[k][a];		//the Q+ set
         }
      }
      gu->MISfile<<std::endl<<"  Q-"<<k<<":";
      for (a=0; a<gu->nvertices; a++){
         if (Qminus[k][a]!=0){
            gu->MISfile.width(3);
            gu->MISfile<<Qminus[k][a];		//the Q- set
         }
      }
      gu->MISfile<<std::endl<<"  S "<<k<<":";
      for (a=0; a<gu->nvertices; a++){
         if (S[k][a]!=0){
            gu->MISfile.width(3);
            gu->MISfile<<S[k][a];			//the S set
         }
      }
      gu->MISfile<<std::endl<<"  k = "<<k<<std::endl;
   }
   eq38 = this->equation_3_8(MIS_or_clique);
   return k;
}

//MIS and clique
int MISUtils::step5(char MIS_or_clique){
   short a;
   k--;						//decrement k
   if (print_MIS){
      gu->MISfile<<"step 5:"<<std::endl<<"  k = "<<k<<std::endl;
   }
   if (k>=0){
      Qplus[k][(X[k])-1]=0;			//remove xk from Q+ and...
      Qminus[k][(X[k])-1] = X[k];	//add it to Q-
      for (a=0; a < gu->nvertices; a++){
         S[k][a]=S[k+1][a];				//set S(k+1) to S(k)
      }
      S[k][(X[k])-1]=0;				//remove xk from S
      done = false;
      if (print_MIS){
         gu->MISfile<<"  Q+"<<k<<":";
         for (a=0; a<gu->nvertices; a++){
            if (Qplus[k][a]!=0)
            {
               gu->MISfile.width(3);
               gu->MISfile<<Qplus[k][a];		//print Q+
            }
         }
         gu->MISfile<<std::endl<<"  Q-"<<k<<":";
         for (a=0; a<gu->nvertices; a++){
            if (Qminus[k][a]!=0){
               gu->MISfile.width(3);
               gu->MISfile<<Qminus[k][a];		//print Q-
            }
         }
         gu->MISfile<<std::endl<<"  S "<<k<<":";
         for (a=0; a<gu->nvertices; a++){
            if (S[k][a]!=0)
            {
               gu->MISfile.width(3);
               gu->MISfile<<S[k][a];			//print S
            }
         }
         gu->MISfile<<std::endl;
      }
      eq38=equation_3_8(MIS_or_clique);
   }
   else{	//k<0, so algorithm is done
      done = true; 
      if (print_MIS){
         gu->MISfile<<"k<=0 and Q+ is empty, stop"<<std::endl;
      }
   }
   return k;
}

//MIS and clique 
int MISUtils::step4(char directed, char MIS_or_clique)
{
   short a, sSize=0, Qplus_size=0, Qminus_size=0;
   for (a=0; a < gu->nvertices; a++)
   {
      Qplus_size+=Qplus[k][a];	//count to see if the sets are empty
      Qminus_size+=Qminus[k][a];
   }
   if (Qplus_size!=0)	//i.e. Q+ is not empty
   {
      if(print_MIS)
      {
         gu->MISfile<<"step 4:  Q+ is not empty, go to step 2"<<std::endl;
      }
      k=step2(directed, MIS_or_clique);
   }
   else if(Qminus_size==0)//Q+ and Q- are empty
   {
      if (print_MIS)
      {
         gu->MISfile<<"step 4:  Q+ and Q- are empty"<<std::endl<<"  Maximal independent set:  [";
      }
      if (MIS_or_clique=='m')
      {
         gu->MISets<<"Maximal independent set:  [";
         std::cout<<"Maximal independent set:  [";
      }
      else
      {
         gu->MISets<<"clique:  [";
         std::cout<<"clique:  [";
      }
      for(a=0;a<gu->nvertices;a++)
      {
         if (S[k][a]!=0)//print the mis
         {
            std::cout.width(3);
            std::cout<<S[k][a];
            if (print_MIS)
            {
               gu->MISfile.width(3);
               gu->MISfile<<S[k][a];
            }
            gu->MISets.width(3);
            gu->MISets<<S[k][a];
            sSize++;
         }
      }
      std::cout<<"]"<<std::endl;
      gu->MISets<<"]"<<std::endl;
      if (print_MIS)
      {
         gu->MISfile<<"]"<<"  Go to step 5"<<std::endl;
      }
      if (sSize>maximum)
      {
         maximum=sSize;
      }
      k=step5(MIS_or_clique);
   }
   else	//Q+ empty, but Q- not
   {
      if (print_MIS)
      {
         gu->MISfile<<"step 4:  Q+ empty, but Q- not empty, go to step 5"<<std::endl;
      }
      k=step5(MIS_or_clique);
   }
   //set the values that are returned
   return k;
}

void MISUtils::updateEdgeCount() {

}

void MISUtils::init(){
   Qminus = new int*[gu->nvertices]();
   Qplus = new int*[gu->nvertices]();
   S = new int*[gu->nvertices]();
   D = new int*[gu->nvertices]();
   X = new int[gu->nvertices]();
  
   adj_matrix = new int**[gu->nvertices / 2 + 1]();

   for(int i = 0; i < gu->nvertices; i++){ // dynamically initialize 2d array
      Qminus[i] = new int[gu->nvertices]();
      Qplus[i] = new int[gu->nvertices]();
      S[i] = new int[gu->nvertices]();
      D[i] = new int[gu->nvertices]();
      adj_matrix[i] = new int*[gu->nvertices]();
      for (int j = 0; j < gu->nvertices; j++) {
         adj_matrix[i][j] = new int[gu->nvertices]();
      }
   }

   k = 0; //set state of k and max_ind(the maximum size of an MIS)
   cost = 0;
   maximum = 0; 
   eq38 = false;		//says whether or not equation 3.8 is satisfied
   done = false;
   print_MIS = false;

   gu->printAdjMatrix();
}

//MIS only
int MISUtils::MIS(char directed){
   this->init();
   char print_MIS;

   do{	//ask whether or not to print the information from each step
      std::cout<<"Do you want to print the intermediate results of the MIS calculation? (Y/N)"<<std::endl;
      std::cin>>print_MIS;
      print_MIS=toupper(print_MIS);
      if ((print_MIS!='N')&&(print_MIS!='Y')){
         std::cout<<"Illegal input, please try again"<<std::endl;
      }
   }while((print_MIS!='N')&&(print_MIS!='Y'));
   start=clock();
   for (int i=0; i<gu->nvertices; i++){
      this->Qplus[0][i]=(i+1);//set Q+0 to all vertices
   }
   if (print_MIS == 'Y'){
      this->print_MIS = true;
      gu->MISfile<<"step 1:"<<std::endl<<"  k = 0"<<std::endl<<"  Q+0:";
      for (int i=0; i<gu->nvertices; i++){
         gu->MISfile.width(3);
         gu->MISfile<<Qplus[k][i];	//print Q+
      }
      gu->MISfile<<std::endl<<"  Q-0:"<<std::endl<<"  S 0:"<<std::endl;//Q- and S always empty initially
   }
   k=step2(directed,'m');
   do{
      if (eq38){	//equation 3.8 is satisfied, go to step 5
         k=step5('m');
      }
      else{		//equation 3.8 is not satisfied, go to step 4
         k=step4(directed, 'm');
      }
   }while(!done);
   finish=clock();
   duration=(finish-start)/CLOCKS_PER_SEC;		//calculate amount of time it took
   std::cout<<"independence number = "<< maximum << std::endl  << "MIS caculation took " << duration << " seconds." << std::endl;
   if (this->print_MIS){
      gu->MISfile<<"independence number = "<<maximum<<std::endl;
   }
   gu->MISets<<"  independence number = "<<maximum<<std::endl<<"MIS calculation took "<<duration<<" seconds."<<std::endl;
   return 0;
}

void MISUtils::printQPlus() {
   gu->MISfile << std::endl << "  Q+" << k << ":";
   for (int a = 0; a<gu->nvertices; a++) {
      if (Qplus[k][a] != 0) {
         gu->MISfile.width(3);
         gu->MISfile << Qplus[k][a];		//the Q+ set
      }
   }
}

void MISUtils::printQMinus() {
   gu->MISfile << std::endl << "  Q-" << k << ":";
   for (int a = 0; a<gu->nvertices; a++) {
      if (Qminus[k][a] != 0) {
         gu->MISfile.width(3);
         gu->MISfile << Qminus[k][a];		//the Q- set
      }
   }
}

void MISUtils::printS() {
   gu->MISfile << std::endl << "  S " << k << ":";
   for (int a = 0; a<gu->nvertices; a++) {
      if (S[k][a] != 0) {
         gu->MISfile.width(3);
         gu->MISfile << S[k][a];			//the S set
      }
   }
}

void MISUtils::printD() {
   gu->MISfile << std::endl << "  D " << k << ":";
   for (int a = 0; a<gu->nvertices; a++) {
      if (D[k][a] != 0) {
         gu->MISfile.width(3);
         gu->MISfile << D[k][a];			//the S set
      }
   }
}

void MISUtils::printAdjMatrix(int d) {
   for (int i = 0; i < gu->nvertices; i++) {
      for (int j = 0; j < gu->nvertices; j++) {
         std::cout << adj_matrix[d][i][j] << " ";
      }
      std::cout << "\n";
   }
   std::cout << std::endl;
}

bool MISUtils::equation_3_8_h(char MIS_or_clique) {
   short a, b, sum = 0, qminus_size = 0, qplus_size = 0;

   for (a = 0; a<gu->nvertices; a++) {
      qminus_size += Qminus[k][a];
   }
   if (qminus_size == 0) {  		//3.8 unsatisfied if Q- is empty
      eq38 = false;			//S)atisfied, U)nsatisfied 
   }
   else {	//Q- is not empty			//else #1
      for (a = 0; a<gu->nvertices; a++) {
         qplus_size += Qplus[k][a];
      }
      if (qplus_size == 0) {
         eq38 = true; 	//3.8 satisfied if Q+ is empty
      }
      else {		//Q+ is not empty		//else #2
         eq38 = false;   //assume 3.8 is unsatisfied
         /*
         sum = 0;
         for (a = 0; a<gu->nvertices; a++) {	//for #1
            if (Qminus[k][a] != 0) { //i.e. l+1 element of Q-
               for (b = 0; b<gu->nvertices; b++) {	//for #2
                  if ((adj_matrix[k][a][b] != 0) && (Qplus[k][b] != 0)) {
                     sum++; // equation 3.8 not satisfied
                  }
               }		// end for #2
            }			//end if 		
         }				//end for #1
         if (sum == 0) {
            eq38 = true;
         }
         */
      }					//end else #2
   }				//end else #1
   if (print_MIS) {
      gu->MISfile << "step 3:";
      if (eq38) {
         gu->MISfile << "  equation 3.8 satisfied, go to step 5" << std::endl;
      }
      else {
         gu->MISfile << "  equation 3.8 not satisfied, go to step 4" << std::endl;
      }
   }
   return eq38;
}

//MIS and clique
int MISUtils::step5_h(char MIS_or_clique) {
   short a;
   k--;						//decrement k
   if (print_MIS) {
      gu->MISfile << "step 5:" << std::endl << "  k = " << k;
   }
   if (k >= 0) {
      Qplus[k][(X[k]) - 1] = 0;			//remove xk from Q+ and...
      Qminus[k][(X[k]) - 1] = X[k];	//add it to Q-
      for (a = 0; a < gu->nvertices; a++) {
         S[k][a] = S[k + 1][a];				//set S(k+1) to S(k)
         D[k + 1][a] = 0;
      }
      S[k][(X[k]) - 1] = 0;				//remove xk from S
      done = false;
      if (print_MIS) {
         printQPlus();
         printQMinus();
         printS();
         printD();
         gu->MISfile << std::endl;
      }
      eq38 = equation_3_8_h(MIS_or_clique);
   }
   else {	//k<0, so algorithm is done
      done = true;
      if (print_MIS) {
         gu->MISfile << "k<=0 and Q+ is empty, stop" << std::endl;
      }
   }
   return k;
}

//this function is used by MIS and clique calculation
int MISUtils::step2_h(char directed, char MIS_or_clique) {
   short a, xk = gu->nvertices;

   for (a = 0; a<gu->nvertices; a++) {
      D[k + 1][a] = D[k][a];

      if ((Qplus[k][a] != 0) && (Qplus[k][a]<xk)) {
         xk = Qplus[k][a];//find the lowest numbered vertex in Qplus
      }
      for (int j = 0; j < gu->nvertices; j++) {
         adj_matrix[k + 1][a][j] = adj_matrix[k][a][j];
      }
   }

   // create new adj_matrix to reflect removal of xk neighbors
   for (int i = 0; i < gu->nvertices; i++) {
      if (adj_matrix[k][xk - 1][i] != 0) {
         D[k + 1][i] = i + 1;
         for (int j = 0; j < gu->nvertices; j++) {
            adj_matrix[k + 1][j][i] = 0;
            adj_matrix[k + 1][i][j] = 0;
         }
      }
   }

   if (xk == 7) {
      printAdjMatrix(k + 1);
   }
   else {
      printAdjMatrix(k + 1);
   }

   D[k + 1][xk - 1] = xk;
   // determine smallest degree in adjacency matrix
   int min = gu->nvertices;
   for (int i = 0; i < gu->nvertices; i++) {
      int sum = 0;
      if (D[k + 1][i] == 0) {
         for (int j = 0; j < gu->nvertices; j++) {
            sum += adj_matrix[k + 1][i][j];
         }
         if (sum < min) {
            min = sum;
         }
      }      
   }

   // set Q+k to only those with the smallest degree vertices
   for (int i = 0; i < gu->nvertices; i++) {
      int sum = 0;
      for (int j = 0; j < gu->nvertices; j++) {
         sum += adj_matrix[k+1][i][j];
      }
      if (sum == min) {
         if (D[k + 1][i] == 0) {
            Qplus[k + 1][i] = i + 1;
         }
      }
   }

   for (a = 0; a<gu->nvertices; a++) {
      //set Q- and S
      //Qminus[k + 1][a] = Qminus[k][a];
      S[k + 1][a] = S[k][a];
   }

   X[k] = S[k + 1][xk - 1] = xk;//add xk to the set X of vertices and to S
   for (a = 0; a<gu->nvertices; a++) {
      if (adj_matrix[k + 1][(X[k]) - 1][a] != 0) {
         Qminus[k + 1][a] = Qplus[k + 1][a] = 0;		//remove correspondences
      }
   }
   Qplus[k + 1][(X[k]) - 1] = 0; //remove xk from Q+(k+1)
   k++;						//increment k
   if (print_MIS)		//print the values of the sets, if requested
   {
      gu->MISfile << "step 2:" << std::endl << "  xi" << k << " = " << xk;

      printQPlus();
      printQMinus();
      printS();
      printD();

      gu->MISfile << std::endl << "  k = " << k << std::endl;
   }
   eq38 = this->equation_3_8_h(MIS_or_clique);
   return k;
}

//MIS and clique 
int MISUtils::step4_h(char directed, char MIS_or_clique)
{
   short a, sSize = 0, Qplus_size = 0, Qminus_size = 0;
   for (a = 0; a < gu->nvertices; a++)
   {
      Qplus_size += Qplus[k][a];	//count to see if the sets are empty
      Qminus_size += Qminus[k][a];
   }
   if (Qplus_size != 0)	//i.e. Q+ is not empty
   {
      if (print_MIS)
      {
         gu->MISfile << "step 4:  Q+ is not empty, go to step 2" << std::endl;
      }
      k = step2_h(directed, MIS_or_clique);
   }
   else if (Qminus_size == 0)//Q+ and Q- are empty
   {
      if (print_MIS)
      {
         gu->MISfile << "step 4:  Q+ and Q- are empty" << std::endl << "  Maximal independent set:  [";
      }
      if (MIS_or_clique == 'm')
      {
         gu->MISets << "Maximal independent set:  [";
         std::cout << "Maximal independent set:  [";
      }
      else
      {
         gu->MISets << "clique:  [";
         std::cout << "clique:  [";
      }
      for (a = 0; a<gu->nvertices; a++)
      {
         if (S[k][a] != 0)//print the mis
         {
            std::cout.width(3);
            std::cout << S[k][a];
            if (print_MIS)
            {
               gu->MISfile.width(3);
               gu->MISfile << S[k][a];
            }
            gu->MISets.width(3);
            gu->MISets << S[k][a];
            sSize++;
         }
      }
      std::cout << "]" << std::endl;
      gu->MISets << "]" << std::endl;
      if (print_MIS)
      {
         gu->MISfile << "]" << "  Go to step 5" << std::endl;
      }
      if (sSize>maximum)
      {
         maximum = sSize;
      }
      k = step5_h(MIS_or_clique);
   }
   else	//Q+ empty, but Q- not
   {
      if (print_MIS)
      {
         gu->MISfile << "step 4:  Q+ empty, but Q- not empty, go to step 5" << std::endl;
      }
      k = step5_h(MIS_or_clique);
   }
   //set the values that are returned
   return k;
}


int MISUtils::H_MIS(char directed){
   this->init();
   char print_MIS;

   do {	//ask whether or not to print the information from each step
      std::cout << "Do you want to print the intermediate results of the MIS calculation? (Y/N)" << std::endl;
      std::cin >> print_MIS;
      print_MIS = toupper( print_MIS );
      if( ( print_MIS != 'N' ) && ( print_MIS != 'Y' ) ){
         std::cout << "Illegal input, please try again" << std::endl;
      }
   } while( ( print_MIS != 'N' ) && ( print_MIS != 'Y' ) );

   start = clock();
   // determine smallest degree in adjacency matrix
   int min = gu->nvertices;
   for (int i = 0; i < gu->nvertices; i++) {
      int sum = 0;
      for (int j = 0; j < gu->nvertices; j++) {
         sum += gu->adj_cost_matrix[i][j];
      }
      if (sum < min) {
         min = sum;
      }
   }

   // set Q+0 to only those with the smallest degree vertices
   int q = 0;
   for (int i = 0; i < gu->nvertices; i++) {
      int sum = 0;
      for (int j = 0; j < gu->nvertices; j++) {
         sum += gu->adj_cost_matrix[i][j];
         adj_matrix[0][i][j] = gu->adj_cost_matrix[i][j];
      }
      if (sum == min) {
         Qplus[0][i] = i + 1;
         q++;
      }
   }

   /*
   for(int i = 0; i < gu->nvertices; i++ ){
      Qplus[0][i] = ( i + 1 );//set Q+0 to all vertices
   }
   */
   if( print_MIS == 'Y' ){
      this->print_MIS = true;
      gu->MISfile << "step 1:" << std::endl << "  k = 0";
      printQPlus();
      printQMinus();
      printS();
      gu->MISfile << std::endl;
   }


   k = step2_h(directed, 'm');
   do {
      if(eq38){	//equation 3.8 is satisfied, go to step 5
         k = step5_h('m');
      }
      else{		//equation 3.8 is not satisfied, go to step 4
         k = step4_h(directed, 'm');
      }
   } while(!done);
   finish = clock();
   duration = ( finish - start ) / CLOCKS_PER_SEC;		//calculate amount of time it took
   std::cout << "independence number = " << maximum << std::endl << "MIS calculation took " << duration << " seconds." << std::endl;
   if(this->print_MIS){
      gu->MISfile << "independence number = " << maximum << std::endl;
   }
   gu->MISets << "  independence number = " << maximum << std::endl << "MIS calculation took " << duration << " seconds." << std::endl;
   return 0;
}

//*****************************************************************************
//*
//*	Title: 			clique
//*	Date:  			July 2000
//*	Version:			1.1
//*	
//*	Project:			Development of C++ graph package
//*
//*	Author: 			Matthew N. Johnson
//*	Coordinator: 		Professor Gary Lamont
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Function description:	finds the cliques of a graph by considering the 
//*					complementary graph to the one entered in the 
//*					data. 
//*
//*	Algorithm:			Uses the same algorithm as MIS (Christofides, p. 
//*					35) but uses the information of the complementary 
//*					graph.
//*
//*	Order of analysis:	O((3^(n/3))*gu->nvertices^2)
//*
//*	Globals:			adjacency matrix, number of vertices
//*
//*	Local Variables:		maximum
//*
//*	Calling Modules:		main program
//*
//*****************************************************************************
//uses equation_3_8, step2, step4 and step5 from above
int MISUtils::clique(char directed)
{
   this->init();
   for (int a=0; a<gu->nvertices; a++){//make gu->clique_adj the adjacency matrix for the complement graph
      for (int b=0; b<gu->nvertices; b++){
         if (gu->adj_cost_matrix[a][b]!=0){
            gu->clique_adj[a][b]=0;
         }
         else{
            gu->clique_adj[a][b]=1;
         }
      }
      gu->clique_adj[a][a]=1;		//a node is adjacent to itself
      Qplus[0][a]=(a+1);		//Qplus equals all vertices
   }
   k=step2(directed, 'c');
   do{
      if (eq38){
         k=step5('c');
      }
      else{
         k=step4(directed, 'c');
      }
   }while(!done);
   std::cout<<"clique number = "<<maximum<<std::endl;
   gu->MISets<<"  clique number = "<<maximum<<std::endl;
   return 0;
}   

MISUtils::MISUtils(GraphUtils* gu){
   this->gu = gu;
}
