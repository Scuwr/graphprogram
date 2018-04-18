//*********************************************************************
//**			copyright (c) 1994/1996/1999/2000					 **
//**																 **
//**		Permission to use, copy and distribute this software	 **
//**		for educational purposes without fee is hereby granted   **
//**		provided that the copyright notice and this permission   **
//**		notice appear on all copies.  Permission to modify the   **
//**		software is granted, but not the right to distribute the **
//**		modified code.  All modifications are to be distributed  **
//**		as changes to released versions by AFIT.				 **
//**																 **
//*********************************************************************
//*********************************************************************
//*
//*	Title: 		Graph Package
//*	Date:  		July 2000
//*	Version:		1.1
//*	Project:		Development of C++ graph package
//*
//*	Authors: 		Professor Gary B. Lamont, James R. Holten 
//*				(1985); Timothy S. Weber (1986); Andy Beard, 
//*				Mark Huson, and Mike Proicou (1990); Andrew E. 
//*				McNear (1992); Matthew N. Johnson(1999,2000)
//*
//*	Coordinator: 	Professor Gary Lamont
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Objective:		This is the main program that coordinates the
//*				actions of the graphing program.  It calls all
//*				functions.  
//*
//*	Parameters-	Inputs   File or terminal data describing a graph		
//*			outputs  Those outputs described in the comments that
//*				   precede the main driver routine.
//*
//*History: An original graph program was written in Pascal in 1985 in 
//*	   support of the AFIT courses in data structures and algorithms.
//*	   The Pascal code was converted to Ada in 1986 with new 
//*	   functions added.  This version represented the graph matrix 
//*	   using a 2-dimensional array (adjacency matrix).  In 1989, 
//*	   the capability was added to compute the maximal independent
//*	   subsets as outlined in the Christofides text.  The software to 
//*	   represent the graph using a linked-list, sparse matrix was 
//*	   developed in 1991.  The user interface was expanded and code 
//*	   was written to compute all spannning trees.  In 1992, several 
//*	   routines were changed to improve functionality.  The help 
//*	   screens were moved into separate file, making future editing 
//*	   much simpler.  The screen displays were improved to better 
//*	   clarify data inputs.  The storing of user-entered data into 
//*	   input files was changed to add descriptive meaning to each 
//*	   data value for ease of reconstruction of a pictorial graph 
//*	   from the output data file.  Many of the modules of the 
//*	   previous version were refined in 1994.  In 1999, the Ada 
//*	   version was rewritten and updated in C++. The function to find 
//*	   all spanning trees was changed to give the number of spanning 
//*	   trees to save space and unnecessary calculation.  
//*
//*********************************************************************
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

#include "GraphUtils.h"
#include "MISUtils.h"

using namespace std;

const int vertices = 100;				//will be the actual number of vertices in the graph

GraphUtils* gu;
//*******************************************************************
//*
//*	Title: 			graph_menus
//*	Date:  			July 2000
//*	Version:		1.1
//*	
//*	Project:		Development of C++ graph package
//*
//*	Author: 		Matthew N. Johnson
//*	Coordinator: 		Professor Gary Lamont
//*
//*	Operating system:	Unix
//*	Language:  		C++
//*
//*	Function description:	This consists of several menus that appear 
//*				throughout the program.  It includes the main
//*				title that is displayed initially, the menu that
//*				asks for type of input, the help menu (if any),
//*				and the menu displaying available functions.
//*
//*	Calling Modules:	main_program
//*
//********************************************************************
//----------------------------- begin Graph_menus
void graph_menus_title(void)
{
	cout<<endl<<"***********************************************************************"<<endl<<endl;
	cout<<"\t\tAir Force Institute of Technology"<<endl;
	cout<<"\t\t Graduate School of Engineering"<<endl;
	cout<<"\t  Department of Electrical and Computer Engineering" <<endl<<endl;
	cout<<"***********************************************************************"<<endl<<endl;
	cout<<"\t\t STANDARD GRAPH ROUTINES IN C++"<<endl<<endl;
	cout<<"\t\tWelcome to the AFIT Graph Program"<<endl<<endl;
}

void graph_menus_help(void)
{
	cout<<endl<<"***********************************************************************"<<endl;
	cout<<"\t\t\tSorry, no help currently"<<endl;
	cout<<"***********************************************************************"<<endl<<endl;
}

char graph_menus_method (void)
{
	char method = '?';
	do {
	cout<<"***********************************************************************"<<endl;
	cout<<"\tWhat type of graph input method do you desire?" <<endl;
	cout<<"\t\t\tYour options are:"<<endl;
	cout<<" F) for file Input, I) for interactive input, H) for help, X) to exit"<<endl;
	cout<<"***********************************************************************"<<endl;
	cin>>method;
	method = toupper (method);
	}while((method!='F')&&(method!='I')&&(method!='H')&&(method!='X'));
	return method;
}

char graph_menus_functions(void)
{
	char operation = '?';
	do	{
	cout<<"***********************************************************************"<<endl<<endl;
	cout<<"\t\tWhat functions are desired?"<<endl<<endl;
	cout<<"\tThe available options for directed/undirected graphs are:" <<endl<<endl;
	cout<<"\t\tS for strong components"<<endl;
	cout<<"\t\tD for depth first search"<<endl;
   cout<<"\t\tC for maximum clique"<<endl;
	cout<<"\t\tM for maximum independent set"<<endl;
   cout<<"\t\tH for maximum independent set with heuristics"<<endl;
	cout<<"\t\tT for the number of spanning trees"<<endl;
	cout<<"\t\tK for Kruskal's minimum spanning tree"<<endl;
	cout<<"\t\tP for Prim's minimum spanning tree"<<endl;
	cout<<"\t\tB for both Kruskal's and Prim's MST"<<endl;
	cout<<"\t\tA for all of the above"<<endl;
	cout<<"\t\tJ for Dijkstra's algorithm"<<endl;
	cout<<"\t\tN for entry of a new graph"<<endl;
	cout<<"\t\tX to terminate program"<<endl<<endl;
	cout<<"All outputs are stored in graph.output except MISs which are stored in "<<endl;
	cout<<"MIS.data and the information about each step of the MIS calculation which,"<<endl;
	cout<<"if desired, is stored in MISinfo.data"<<endl;
	cout<<"***********************************************************************"<<endl<<endl;
	cin>>operation;
	operation = toupper(operation);
	}while ((operation != 'S')&& (operation != 'D') && (operation != 'C') && (operation != 'M') && (operation != 'H') && (operation  != 'T') && (operation != 'K') && (operation != 'P') && (operation != 'B')&& (operation != 'A'  )&& (operation!='J')&&(operation != 'N')&& (operation != 'X'));
	return operation;
}
//------------------------end graph menus
//------------------------begin file input 
//*****************************************************************************
//*
//*	Title: 			file_input
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
//*	Function description:	This program takes file input from a source file 
//*					and enters it in the program as the graph data.  
//*					However, currently, there is no guard against a 
//*				     filename that does not exist.  If such a file name 
//*					is typed in, there will be no data for the graph.  
//*				     Running Kruskal's algorithm seems to be the safest 
//*					method to check that there is data, as other 
//*					operations may encounter infinite loops, causing 
//*					the program to hang. 
//*
//*	Calling Modules:		Main program
//*
//****************************************************************************
void print_adj_matrix(void);

char file_input(string inputfile){
	//string inputfile;	//name of the file data is taken from
	char weighted='?', directed='?', file_opened='n';//test to see if file was found
	short a, b, current_vtx, out_deg, to_vtx, arc_index;//iterators
	double arc_weight;//allow for non-integer values of weights
	
	do{		//make sure a valid file is found
		cout<<"Please enter the name of graph input file:"<<endl;
		cin>>inputfile;
		ifstream infile(inputfile.c_str());
		if (infile.fail())
		{	
			file_opened='n';
			cout<<"**Error:  Failed to open "<<inputfile<<"!!"<<endl<<endl;
		}
		else
		{
			file_opened='y';
		}
	}while(file_opened!='y');	//loop until a good file is found
	ifstream infile(inputfile.c_str());	//re-initialize once proper name is taken
	infile>>directed;		// directed?
	//skip reading up to 50 characters and to the next line, to allow comments
	infile.ignore(50,'\n');	
	infile>>weighted;		// weighted?
	infile.ignore(50,'\n');	

   int nvertices;
	infile>>nvertices; 	// number of vertices
   gu = GraphUtils::New(nvertices);

	for (a=0; a < gu->nvertices; a++)
	{	//set cost/adjacency and M matrices to 0
		for (b=0; b<gu->nvertices; b++)
		{
			gu->adj_cost_matrix[a][b]=gu->M[a][b]=0; 
		}
	}
	for (current_vtx=0; current_vtx<gu->nvertices; current_vtx++)//begin for #1, 
	{
		infile.ignore(50,'\n');	
		infile>>out_deg; //out degree of current vertex
		if (directed!='Y')
		{
			gu->M[current_vtx][current_vtx]=out_deg;
		}
		for(arc_index=1; arc_index<=out_deg; arc_index++)//begin for #2 
		{
			infile.ignore(50,'\n');	
			infile>>to_vtx;
			if (directed=='Y')
			{
				gu->adj_cost_matrix[current_vtx][to_vtx-1]=1;
				gu->M[current_vtx][to_vtx-1]--;	//decrement for a link between vertices
				gu->M[to_vtx-1][current_vtx]--;
				gu->M[to_vtx-1][to_vtx-1]++;	//increment for indegree
			}
			else	//not directed
			{
gu->adj_cost_matrix[current_vtx][to_vtx-1]=gu->adj_cost_matrix[to_vtx-1][current_vtx]=1; 
				gu->M[current_vtx][to_vtx-1]--; //decrement for a link between vertices
			}
			if(weighted=='Y')
			{
				infile.ignore(50,'\n');	
				infile>>arc_weight;
				if (directed=='Y')
				{
					gu->adj_cost_matrix[current_vtx][to_vtx-1]=arc_weight;
				}
				else
				{
			gu->adj_cost_matrix[current_vtx][to_vtx-1]=gu->adj_cost_matrix[to_vtx-1][current_vtx]=arc_weight;
				}
			}
		}		//end for #2 
	}		//end for #1
	gu->outputfile<<"Data for "<<inputfile<<":"<<endl;
	return directed;
}
//------------------------end file input
//------------------------begin interactive input 
//*****************************************************************************
//*
//*	Title: 			interactive_input
//*	Date:  			July 2000
//*	Version:			1.1
//*	
//*	Project:			Development of C++ graph package 
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Function description:	Obtains graph information via keyboard input from 
//*					the user.  This information is used throughout 
//*					the rest of the program to run other functions.
//*
//*	Calling Modules:		graph_program
//*
//*****************************************************************************
char interactive_input(void)
{
	string savefile;	//will be the name of the file to which data is saved
	char save ='?', correct_data='?', weighted='?', directed='?'; //save data or not, test to see if input was correct
	bool save_graph_data;
	short a, b, current_vtx, out_deg, to_vtx, arc_index=1;	//iterators
	double tempdata[2*vertices+1], arc_weight=1;
//to vertex and weight to every vertex + overall out-degree=2*vertices+1
	int counter;	//will count the number of useful inputs in tempdata
	do{
		cout<<"Do you want to save the graph input data? (Y/N)"<<endl;
		cin>>save;
		save = toupper(save);
		if ((save!='Y')&&(save!='N'))
		{
			cout<<"Illegal input, please try again."<<endl;
		}
	} while ((save != 'Y') && (save != 'N'));	
	save_graph_data = (save == 'Y'); //qualifies bool 
	if (save_graph_data)
	{
	cout<<"What do you want to call this data file?"<<endl;//ask user for name of file
	cin>>savefile;
	cout<< "You are about to enter graph data in file:  "<<savefile<<endl;
	}
	else//avoid trying to initialize a null string as a file name
	{
		savefile="void.data";
	}
	ofstream outFile(savefile.c_str());
	do{		//ask if directed
		cout<<"Is this a directed graph? (Y/N)"<<endl;
		cin>>directed;
		directed = toupper (directed);
		if ((directed!='Y')&&(directed!='N'))
		{
		cout<<"Illegal input, please try again."<<endl;
		}
	} while ((directed != 'Y')&& (directed != 'N')); 
	if (save_graph_data)
	{
		outFile<<directed<<"\t//directed"<<endl;	//write to file and comment
	}
	do{				//find out if the graph is weighted
		if (directed=='Y') 
		{
			cout<<"Do the arcs have weights? (Y/N)"<<endl;
		}
		else
		{
			cout<<"Do the links have weights? (Y/N)"<<endl;
		}
		cin>>weighted;
		weighted = toupper(weighted);
		if ((weighted != 'Y')&& (weighted != 'N'))
		{
			cout<<"Illegal input, please try again"<<endl;
		}
	}while ((weighted != 'Y')&& (weighted != 'N'));
	if (save_graph_data)
	{
		outFile<<weighted<<"\t//weighted"<<endl;	//write to file and comment
	}
	do{
		cout<<"How many vertices/nodes does the graph have:"<<endl;
      int nvertices;
		cin>>nvertices;
      gu = GraphUtils::New(nvertices);
		if ((gu->nvertices<=0)||(gu->nvertices>vertices))
		{
		cout<<"**Error: You must have more than 0 and less than "<<vertices+1<<" vertices"<<endl;
		}
	}while((gu->nvertices<=0)||(gu->nvertices>vertices));//check to see there is an appropriate number of nodes
	//doesn't check to make sure gu->nvertices is an integer (will crash program)
	if (save_graph_data)
	{
		outFile<<gu->nvertices<<"\t//number of vertices"<<endl;
	}
	if (directed=='Y')
	{
		cout<<"Input the arcs coming out of the ";
	}
	else
	{
		cout<<"Input the links coming out of the ";
	}
	cout<<gu->nvertices<<" vertices by their destination number."<<endl;
	for (a=0; a < gu->nvertices; a++)
	{//set cost/adjacency and M matrices to 0
		for (b=0; b<gu->nvertices; b++)
		{
			gu->adj_cost_matrix[a][b]=gu->M[a][b]=0;	
		}
	}
	for(current_vtx=0;current_vtx<gu->nvertices; current_vtx++)	//begin for #1
	{
	do{											//do #1
		for (a=0; a<=2*vertices; a++)
		{
			tempdata[a]=0;	//clear tempdata
		}
		counter=0;				//will be the number of entries to be saved
		if (directed=='Y')
		{
		cout<<"Input the outdegree of vertex "<<(current_vtx+1)<<" (include loops):"<<endl;
		}
		else
		{
			cout<<"Input the degree of vertex "<<(current_vtx+1)<<" (include loops):"<<endl;
		}
		cin>>out_deg;
		tempdata[counter]=out_deg;
		counter++;					//increment for out degree
		for (arc_index=1;arc_index<=out_deg;arc_index++) //begin for #2 
		{
			if (directed=='Y')
			{
			cout<<"arc "<<arc_index<<" out of vertex "<<(current_vtx+1)<<" goes to vertex: (if a loop, put in same vertex)"<<endl;
			}
			else
			{
 			cout<<"link "<<arc_index<<" out of vertex "<<(current_vtx+1)<<" goes to vertex: (if a loop, put in same vertex)"<<endl;
			}
			cin>>to_vtx;
			tempdata[counter]=to_vtx;
			counter++;		//increment for to vertex
			if(weighted=='Y')
			{
				cout<<"and has weight:"<<endl;
				cin>>arc_weight;
				tempdata[counter]=arc_weight;
				counter++;			//increment for weight
			}
		}		//end for #2
	    do{						//do #2
			cout<<"Was the data for vertex "<<current_vtx+1<<" correct? (Y/N)"<<endl;
			cin>>correct_data;
			correct_data=toupper(correct_data);
			if ((correct_data!='Y')&&(correct_data!='N'))
	        {
				cout<<"Illegal input, please try again."<<endl;
			}
		}while((correct_data!='Y')&&(correct_data!='N'));	//end do #2
	}while(correct_data!='Y');		//end do #1
//now input the data when it is confirmed to be correct
	if (directed!='Y')		
	{
		gu->M[current_vtx][current_vtx]=tempdata[0];	//the out degree
	}
	if (weighted=='Y')
	{
		for (a=1; a<counter; a+=2)//take every other value if weighted
		{
			if(directed=='Y')
			{
				gu->adj_cost_matrix[current_vtx][int(tempdata[a])-1]=tempdata[a+1];//a+1 is a weight
				gu->M[current_vtx][int(tempdata[a])-1]--;
				gu->M[int(tempdata[a])-1][current_vtx]--;
				gu->M[int(tempdata[a])-1][int(tempdata[a])-1]++;
			}
			else	//not directed
			{
				gu->adj_cost_matrix[current_vtx][int(tempdata[a])-1]=gu->adj_cost_matrix[int(tempdata[a])-1][current_vtx]=tempdata[a+1];
				gu->M[current_vtx][int(tempdata[a])-1]--;
			}
		}
	}
	else		//i.e. not weighted
	{
		for (a=1;a<counter; a++)
		{
			if(directed=='Y')
			{
				gu->adj_cost_matrix[current_vtx][int(tempdata[a])-1]=1;
				gu->M[current_vtx][int(tempdata[a])-1]--;
				gu->M[int(tempdata[a])-1][current_vtx]--;
				gu->M[int(tempdata[a])-1][int(tempdata[a])-1]++;
			}
			else
			{
				gu->adj_cost_matrix[current_vtx][int(tempdata[a])-1]=gu->adj_cost_matrix[int(tempdata[a])-1][current_vtx]=1;
				gu->M[current_vtx][int(tempdata[a])-1]--;
			}
		}
	}			
	if (save_graph_data)//saves and comments the graph data
	{
		outFile<<tempdata[0]<<"\t//(out) degree of vertex "<<current_vtx+1<<endl;
		for (a=1;a<counter;a++)
		{
			if (weighted=='Y')
			{
				if (a%2==0)	//every other input is a weight, starting with tempdata[2]
				{
					outFile<<tempdata[a]<<"\t//  arc/link weight"<<endl;
				}
				else//the odd-numbered inputs are the to-vertexes, starting with tempdata[1]
				{
					outFile<<tempdata[a]<<"\t//arc/link to vertex "<<tempdata[a]<<endl;
				}
			}
			else//if not weighted, every input is a to-vertex
			{
				outFile<<tempdata[a]<<"\t//arc/link to vertex "<<tempdata[a]<<endl;
			}
		}
	}
	}		//end for #1
	if (save_graph_data)
	{
		gu->outputfile<<"Data for "<<savefile<<":"<<endl;
	}
	else
	{
		gu->outputfile<<"Data for your graph:"<<endl;
	}
	return directed;
}
//--------------------------------------end interactive input

void print_adj_matrix(void)
{
	short a, b;
	gu->outputfile<<"adjacency/cost matrix:"<<endl;
	for (a=0; a<gu->nvertices; a++)
	{	
		for (b=0; b<gu->nvertices; b++)
		{
			gu->outputfile.width(4);
			gu->outputfile<<gu->adj_cost_matrix[a][b];
		}
		gu->outputfile<<endl;
	}
}
//--------------------------------------begin strong component
//************************************************************************
//*
//*	Title: 			strongComponent
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
//*	Function description:	Finds the strong components of a graph along with 
//*					the basis and contra basis of the graph.  For 
//*					nondirected graphs, it returns all vertices for
//*					strong component, basis and contrabasis.
//*
//*	Algorithm:			Strong components are determined from the Q and R 
//*					matrices.  If a vertex can reach another vertex 
//*					and be reached from that vertex, the two vertices 
//*					are in the same strong component.  This is 
//*					determined by comparing a particular spot in the 
//*					matrices, if Q[i][j]=1 and R[i][j]=1, xi and xj 
//*					are in the same strong component.  Basis is found 
//*					with the Reaching matrix.  Rows of R are compared 
//*					to find the vertex that can reach the most other 
//*				      vertices.  If a vertex can reach all the vertices, 
//*					it is the basis.  Otherwise, other vertices are 
//*					added until all the vertices can be reached.  The 
//*					contrabasis is found in the same manner, except 
//*					that the Q (reachability) matrix is used instead 
//*					of the R matrix.
//*
//*	Order of analysis:	O(gu->nvertices^3)
//*
//*	Globals:			adjacency matrix
//*
//*	Local Variables:		Reaching and reachability matrices
//*
//*	Calling Modules:		graph_program
//*
//************************************************************************
void basis_contra(int Q_R[vertices][vertices], int SC_mat[vertices][vertices], int type)//type = 1 for basis, 2 for contrabasis
{
	short a, b, c, current_vtx=0, number=0, number_reached=0, most_reached=0, contra_base=0, complement=0, complementor=0;
	list <int> contra_basis;//stores a list of the vertices currently in the (contra)basis (actual numbers)
	list <int>::iterator c_b;
	list <int> not_reached;//stores a list of vertices not reached(reachable) (indexes to vertices)
	list <int>::iterator nr;
	for (a = 0; a < gu->nvertices; a++)
	{
		number_reached = 0;
		for (b=0; b < gu->nvertices; b++)
		{
			number_reached +=Q_R[a][b];//will equal 1 if a can reach b, 0 if not
		}
		if (number_reached > most_reached)
		{
			most_reached = number_reached;  
			contra_base=(a+1);	//gives the row (i.e. vertex) reaching(reached from) the most 
					//vertices, must be part of (contra)basis
		}
	}	
	contra_basis.push_back(contra_base); //begins a list of the vertices in the (contra)basis
	for (a=0; a<gu->nvertices; a++)
	{
		if(Q_R[contra_base-1][a]==0)
		{
			not_reached.push_back(a);	//list of vertices not reached(reachable) by vertices
		}								//in the (contra)basis currently
	}
	if (not_reached.size()>0)//i.e. one vertex/strong component is not the sole (contra)basis //if #1
	{
		do{
			for (current_vtx=0; current_vtx<gu->nvertices; current_vtx++)		//for #1
			{
				for(nr=not_reached.begin();nr!=not_reached.end();nr++)//for #2
				{
		//find a vertex that does reach a not-reached(reachable) vertex
					if (Q_R[current_vtx][(*nr)]==1)
					{
						number=0;
						for (c_b=contra_basis.begin(); c_b!=contra_basis.end(); c_b++)
						{
							if((*c_b)!=current_vtx+1)//current vertex not that particular element of the basis
							{
								number++;
							}
						}
					}    //if
				//find the vertex that is not in the (contra)basis and that reaches(is reached by) the most 
				//vertices not yet found by vertices in the current (contra)basis
					if (number==contra_basis.size())//current vertex is not in the (contra)basis
					{
						b=0;
						for (c=0; c<gu->nvertices; c++)
						{
							for(c_b=contra_basis.begin();c_b!=contra_basis.end();c_b++)
							{
								if ((Q_R[current_vtx][c]==1)&&(Q_R[(*c_b)-1][c]==0))
					//current vertex reaches(reachable from) c and the members of the (contra)basis do(are) not
								{
									b++;
								}
							}
						}
						if (b>complement)
						{
							complement=b;			//highest total from above
							complementor=current_vtx;	//add this to the (contra)basis
						}
					}
				}			//end for #2
			}				// end for #1
			//once the best new candidate is found, add it and proceed
			contra_basis.push_back(complementor+1);	//store actual value
			for (a=0; a<gu->nvertices; a++)
			{
				if (Q_R[complementor][a]==1)
				{
					not_reached.remove(a);//take out vertices that are now reached(reachable)
				}
			}
			complement=0;//reset for another pass
		}while(not_reached.size()>0); //untill all the vertices are reached(reachable)
		if (type==1)	//means calculating the basis, so indicate this in print
		{
			gu->outputfile<<"  basis: (choose one from each []):"<<endl;
			cout<<"basis: (choose one from each [])"<<endl;
		}
		else				//i.e. calculating contrabasis
		{
			gu->outputfile<<"  contrabasis: (choose one from each []):"<<endl;
			cout<<"contrabasis: (choose one from each [])"<<endl;
		}
		for (c_b=contra_basis.begin(); c_b!=contra_basis.end(); c_b++)
		{//output the (contra)basis
			cout<<"[";
			gu->outputfile<<"[";
			for (a=0; a<gu->nvertices; a++)
			{
				if (SC_mat[(*c_b)-1][a]!=0)
				{
					cout.width(3);
					cout<<SC_mat[(*c_b)-1][a];
					gu->outputfile.width(3);
					gu->outputfile<<SC_mat[(*c_b)-1][a];
				}
			}
			cout<<"]";
			gu->outputfile<<"]";
		}
		cout<<endl;
		gu->outputfile<<endl;
	}			//end if #1
	else			//one vertex/strong component is the (contra)basis
	{

		if (type==1)//calculating the basis
		{
			gu->outputfile<<"  basis: (choose one from within the []):"<<endl<<"[";
			cout<<"basis: (choose one from within [])"<<endl<<"[";
		}
		else		//calculating the contrabasis
		{
			gu->outputfile<<"  contrabasis: (choose one from within the []):"<<endl<<"[";
			cout<<"contrabasis: (choose one from within the [])"<<endl<<"[";
		}
		for (a=0; a<gu->nvertices; a++)
		{
			if (SC_mat[(contra_base)-1][a]!=0)
			{
				cout.width(3);
				cout<<SC_mat[(contra_base)-1][a];
				gu->outputfile.width(3);
				gu->outputfile<<SC_mat[(contra_base)-1][a];
			}
		}
		cout<<"]"<<endl;
		gu->outputfile<<"]"<<endl;
	}
}

void strong_component(char directed)
{
	int a, b, c, no_of_elements=0, similarities=0;	//iterators
	int R[vertices][vertices]={0}, R2[vertices][vertices], Q[vertices][vertices], SC_mat[vertices][vertices]={0};
	for (a=0; a<gu->nvertices; a++)
	{
		R[a][a]=1;	//a vertex can reach itself obviously
	}
	for (a=0; a<gu->nvertices; a++)	//find correspondence up to degree 2
	{
		for (b=0; b<gu->nvertices; b++)
		{
			if (gu->adj_cost_matrix[a][b]!=0)//vertex a and b are connected...
			{
				R[a][b]=1;
				for (c=0; c<gu->nvertices; c++) 
				{
					if (gu->adj_cost_matrix[b][c]!=0)//b and c are connected...
					{
						R[a][c]=1;//therefore, a can reach c
					}
				}
			}
		}
	}
	do{					//find up to correspondence^n
		for (a=0; a<gu->nvertices; a++)		
		{
			for (b=0; b<gu->nvertices; b++)	
			{
				R2[a][b]=R[a][b];		
			}					
		}						
		for (a=0; a<gu->nvertices; a++)
		{
			for (b=0; b<gu->nvertices; b++)
			{
				if (R[a][b]==1)		//a and b are connected
				{
					for (c=0; c<gu->nvertices; c++)
					{
						if (R[b][c]==1)//b and c are connected
						{
							R[a][c]=1;	//so a can reach c
						}	//if
					}		//for
				}			//if
			}				//for 
		}					//for
		similarities=0;	
		//similarities counts number of changes from one iteration to the next,
		// if equal to gu->nvertices^2, all correspondences are found
		for (a=0; a<gu->nvertices; a++)		
		{
			for (b=0; b<gu->nvertices; b++)	
			{
				if (R2[a][b]==R[a][b])	
				{
					similarities++;			
				}				
			}					
		}						
	}while(similarities<(gu->nvertices*gu->nvertices));
	//display R matrix and 	create reachability matrix(transpose of R matrix)
	gu->outputfile<<"R matrix"<<endl;
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			gu->outputfile.width(3);
			gu->outputfile<<R[a][b];
			Q[b][a]=R[a][b];	
		}
		gu->outputfile<<endl;
	}
	//display Q matrix
	gu->outputfile<<"Q matrix"<<endl;
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			gu->outputfile.width(3);
			gu->outputfile<<Q[a][b];
		}
		gu->outputfile<<endl;
	}
	//now find Strong Component from Q and R matrices, to be output later
	for (a=0;a<gu->nvertices;a++)
	{
		for (b=0;b<gu->nvertices;b++)
		{
			if ((R[a][b]==1)&&(Q[a][b]==1))  //compare a row of R and Q
			{//a can reach b, and a can reached from b, so in the same SC
				SC_mat[a][b]=b+1; 
			}
		}
	}
	basis_contra(R, SC_mat,1);				//find the basis
	basis_contra(Q, SC_mat,2);				//find the contrabasis
	//if a and b are in the same strong component, their rows in SC_mat will be the same
	//so set one to all zeros, so that the same strong component isn't printed more than once
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=a+1; b<gu->nvertices; b++)
		{
			if((SC_mat[b][a]==SC_mat[a][a])&&(SC_mat[a][a]!=0))
			{//a and b are in the same strong component
				for (c=0; c<gu->nvertices; c++)
				{			
					SC_mat[b][c]=0;//discount correspondence of b, because
				}					//they are already counted by a
			}
		}
	}
	cout<<"strong components:"<<endl;			//output the strong components
	gu->outputfile<<"strong components:"<<endl;
	for (a=0; a<gu->nvertices; a++)
	{
		no_of_elements=0;
		for (b=0; b<gu->nvertices; b++)
		{
			no_of_elements+=SC_mat[a][b];
		}
		if (no_of_elements!=0)//already eliminated double counting above, so it gives unique bases
		{
			cout<<"[";
			gu->outputfile<<"[";
			for (c=0; c<gu->nvertices; c++)
			{
				if (SC_mat[a][c]!=0)
				{
					gu->outputfile.width(3);
					gu->outputfile<<SC_mat[a][c];
					cout.width(3);
					cout<<SC_mat[a][c];
				}
			}
			cout<<"]"<<endl;
			gu->outputfile<<"]"<<endl;
		}
	}	
}
//---------------------end strong component

//---------------------begin dfs
//************************************************************************
//*
//*	Title: 			dfs
//*
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
//*	Function description:	performs a depth first search on the vertices of 
//*				     the graph.  In the case of an undirected graph, it 
//*					will usually return the vertices in only one or a 
//*					few branches.  Directed graphs will return a new 
//*					list on each line that starts with the lowest 
//*					numbered vertex that has not yet been reached and 
//*					begins a search from that vertex. 
//*
//*	Order of analysis:	at most O(n^2)
//*
//*	Globals:			the adjacency matrix	
//*
//*	Local Variables:		a list of discovered vertices
//*
//*	Calling Modules:		main program
//*
//************************************************************************
int find_new_tree(list <int> discovered)//new tree altogether
{
	short possible_start=1, new_tree_start=-1;
	list <int>::iterator d;
	discovered.push_front(gu->nvertices+1);//in case all the elements from 1 to n are present
	discovered.sort();//gives an ordered list of the discovered vertices
	//finds the lowest numbered vertex not yet in the tree
	do{
	for (d=discovered.begin(); d!=discovered.end(); d++)
		{
			if ((possible_start!=(*d))&&(new_tree_start==-1))//second half assures it is the first one
			{
				new_tree_start=possible_start;
			}
		possible_start++;
		}
	}while(new_tree_start==-1);	//assures taking the lowest one
	cout<<"new tree:"<<endl<<new_tree_start;
	gu->outputfile<<"new tree:"<<endl<<new_tree_start;
	return new_tree_start;
}

int find_new_start(list <int> discovered)//new start in same tree
{
	short a, correspondences_discovered=0, fully_used_vertices=0, total_correspondences=0, new_start=0;
	list <int>::iterator d;
	list <int>::iterator dd;
	do{
		for (d=discovered.begin(); d!=discovered.end(); d++)
		{
			correspondences_discovered=total_correspondences=0;//reset values
 			for (a=0; a<gu->nvertices; a++)
			{
				if (gu->adj_cost_matrix[(*d)-1][a]!=0)
				{
					total_correspondences++;	//increments for # of correspondence
					for (dd=discovered.begin(); dd!=discovered.end(); dd++)
					{
						if ((a+1)==(*dd))
						{
							correspondences_discovered++;//increments if correspondence already 
						}		//discovered
					}
				}
			}
			if (correspondences_discovered==total_correspondences)	//all correspondences have been discovered
			{
				fully_used_vertices++;
			}
			else		//vertex has some correspondence not yet discovered
			{
				new_start=(*d);
			}
		}
	}while((new_start==0)&&(fully_used_vertices<discovered.size()));
	if ( fully_used_vertices == discovered.size() )	//all correspondences used, find a new tree
	{
		new_start=-1;
	}
	return new_start;
}

void dfs(short start)
{
	short a, b, number_of_corr, tempstart, dfs_mat[vertices][vertices]={0};
	char found='?';
	list <int> discovered;
	list <int>::iterator d;
	double ave=0.0, sum;
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			if (gu->adj_cost_matrix[a][b]!=0)
			{
				dfs_mat[a][b]=1;	//create simple, 1,0 matrix that will be changed
			}
		}
	}
	do{
		cout<<"At which vertex would you like to start the depth-first search?"<<endl;
		cin>>start;
		if((start<=0)||(start>gu->nvertices))
		{
			cout<<"**Error: The start must be between 1 and "<<gu->nvertices<<endl;
		}
	}while((start<=0)||(start>gu->nvertices));
	discovered.push_front(start);//indicates a vertex has been found in the search
	cout<<"Depth first search starting with vertex "<<start<<":"<<endl<<start;
	gu->outputfile<<"Depth first search starting with vertex "<<start<<":"<<endl<<start;
	do{			//do #1
		  do{	//do #2
				tempstart=start;		//save tempstart in case changed later
				for (a=0; a<gu->nvertices; a++)
				{
					if (gu->adj_cost_matrix[start-1][a]!=0)//a is a correspondence to the start
					{
						found='?';
						for (d=discovered.begin(); d!=discovered.end(); d++)
						{
							if ((a+1)==(*d))
							{
								found='y';//vertex a has been discovered
							}
						}
						if (found!='y')//hadn't yet been found
						{
							dfs_mat[start-1][a]=dfs_mat[a][start-1]=-1;
							//-1 indicates a vertex is now found, and won't be considered again
							cout<<" --> "<<a+1;
							gu->outputfile<<" --> "<<a+1;
							discovered.push_front(a+1);//indicate the vertex is discovered
							start=a+1;//start from the new vertex just found
						}
						else
						{
							dfs_mat[start-1][a]=dfs_mat[a][start-1]=-1;
						//already have that vertex, so don't consider again
						}
					}//if
					number_of_corr=sum=ave=0;
					if (tempstart!=start)//there is a correspondence to the new start
					{
						for (b=0; b<gu->nvertices; b++)
						{	
							sum+=dfs_mat[start-1][b];
						//finds the sum of the values for the correspondences of the start
							if (dfs_mat[start-1][b]!=0)
							{
								number_of_corr++;//finds the number of elements in the sum
							}
						}
						if (number_of_corr>0)
						{
							ave=sum/number_of_corr; 
						}
					}//if 
 				}	//for
		  }while((ave>-1)&&(number_of_corr!=0));//ave=-1 when all correspondences of a vertex are found
  //number_of_corr==0 when there are no correspondences to the start (end while #1)
		if ((ave==-1)||(number_of_corr==0)&&(discovered.size()<gu->nvertices))
	//either all correspondences are used(ave==-1) or there are none(number_of_corr==0)
		{
			cout<<endl;
			gu->outputfile<<endl;
			start=find_new_start(discovered);
			if (start!=-1)	//start=-1 if all correspondences of all the discovered vertices 
			{	//have been used, if so, start a new tree
				cout<<start;
				gu->outputfile<<start;
			}
			else //i.e. start==-1, so start a new tree
			{
				if (discovered.size()<gu->nvertices)
				{
					start=find_new_tree(discovered);//find a new tree, send the list of vertices already discovered
					discovered.push_front(start);	//add the new start to the list of discovered vertices
				}
			}
		}
	}while(discovered.size()<gu->nvertices);		//end while #2
	cout<<endl;
	gu->outputfile<<endl;
}
//---------------------end dfs

//----------------------------begin calculation of number of spanning trees 
//*****************************************************************************
//*	Title: 			spanning_trees
//*
//*	Date:  			July 2000
//*	Version:			1.1
//*	
//*	Project:			Development of C++ graph package
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Function description:	This function determines the number of spanning 
//*					trees in a given graph.  It does so by taking the 
//*					determinant of the matrix M, which is determined 
//*					from graph input.  M is defined as follows:  a 
//*					diagonal entry, mii, is the degree of vertex xi, if not 
//*					directed, and the in degree if directed, and an entry mij 
//*					is minus the number of parallel links between xi and xj.  
//*					It contains only n-1 rows/columns, where n is the number 
//*					of vertices. 
//*
//*	Algorithm:		Uses Gauss-Jordan elimination to find an upper 
//*					triangular matrix and then multiplies the 
//*					pivots to find the determinant.  It is assumed 
//*					that no vertex has degree(in and out together) of 
//*					zero.  This makes this abbreviated method viable 
//*					for the matrix.
//*
//*	Order of analysis:
//*
//*	Calling Modules:		Main program, kruskal, prim
//*
//*****************************************************************************
double spanning_trees(int print_option)//print if option==1, don't if option==0
//option is 0 when called from kruskal() or prim()
{
 	short a, b, h;	//all iterators
	double multiplier;	//factor for multiplying and eliminating columns
	double determinant=1.0;//determinant of the matrix
	for (a=0; a<gu->nvertices-1; a++)
	{
		for(h=1; h<gu->nvertices-1-a; h++)
		{
			multiplier=1.0;
			if (gu->M[a+h][a]!=0)		//an element below a pivot is not zero
			{
 				multiplier = ( ( -(gu->M[a+h][a]) ) / (gu->M[a][a]));
				for (b=a; b<gu->nvertices-1; b++)
 				{
					gu->M[a+h][b]+=((gu->M[a][b])*multiplier);//eliminate a ro
				}
			}
		}
	}
	for (a=0; a<gu->nvertices-1; a++)
	{
		determinant*=gu->M[a][a];
	}
	if(determinant < 1)
		determinant = 0;
	if (print_option==1)	//when called from main, not kruskal or prim or SC
	{
		cout<<"number of spanning trees = "<<determinant<<endl;
		gu->outputfile<<"number of spanning trees = "<<determinant<<endl;
	}
	return determinant;
}
//--------------------------end calculation of number of spanning trees 
//*******************************
//* Title:		remove_loops
//* Date:		May 2000
//* Version:	1.0
//*
//* Function description:   It rules out one side of a loop(that is two links between
//* two vertices) If weighted, get rid of heavier side, else arbitrarily from 
//* the higher numbered one to the lower numbered one. This is to ease in the 
//* calculation of the minimum spanning trees.
//* 
//* Calling modules: kruskal(), prim()
//********************************
double remove_loops(char directed, double tempcost_matrix[vertices][vertices])
{
	short a, b;//iterators
	if (directed=='Y')
	{
		for (a=0; a<gu->nvertices; a++)
		{
			for (b=gu->nvertices-1; b>=a; b--)//only consider the lower triangle of the matrix
			{
				if ((tempcost_matrix[a][b]!=0)&&(tempcost_matrix[b][a]!=0))
				{
					if(tempcost_matrix[a][b]>tempcost_matrix[b][a])//if unweighted, are equal
					{
						tempcost_matrix[a][b]=999999;
					}
					else
					{
						tempcost_matrix[b][a]=999999;//defaults here if unweighted
					}
				}
			}
		}
	}
	return tempcost_matrix[vertices][vertices];
}
//--------------------------begin kruskal and prim algorithms
//*****************************************************************************
//*
//*	Title: 			Prim
//*	Date:  			July 2000
//*	Version:			1.1
//*
//*	Project:			Development of C++ graph package
//*
//*	Operating system:		Unix
//*	Language:  			C++
//*
//*	Function description:	Determines the minimum spanning tree using Prim's
//*					algorithm as found in Christofides, p. 139-139.  
//*					A problem could arise in a situation where 
//*					an undirected, weighted graph had an area where 
//*					there were two parallel links between two 
//*					vertices with different weights.  In this case, 
//*					if the lower cost is entered first, it would be 
//*					erased and the second, higher cost would 
//*					overwrite it.  To avoid this, simply input the 
//*					lower cost second.
//*
//*	Order of analysis:	O(gu->nvertices^2)
//*
//*	Globals:			number of vertices (gu->nvertices) the adjacency 
//*					matrix (adj_matrix)
//*	
//*	Local Variables:		a temporary adjacency matrix (tempadj_matrix) 
//*
//*	Calling Modules:		graph program
//*
//************************************************************************
void prim(char directed)
{
	double possible;
	possible=spanning_trees(0);//test to see if there are even spanning trees(undirected only)
	if ((possible==0)&&(directed!='Y'))
	{
		cout<<"Prim: The graph is disconnected, so there are no spanning trees"<<endl;
		gu->outputfile<<"Prim: The graph is disconnected, so there are no spanning trees"<<endl;
		return;
	}
	list <int> visited;//a list of the vertices that are visited
	list <int>::iterator v;
	short low_column = 0, low_row=0;//will store the vertices to be used
	short a, b, f, g, h, nlinks=0;	//iterators and counters
	double total_cost=0, lowest_cost=999998;	//total cost so far and current lowest cost under consideration
	double beta_matrix[vertices][vertices];	//will store the costs of links under consideration
	double tempcost_matrix[vertices][vertices];//stores costs and will be changed
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			tempcost_matrix[a][b]=gu->adj_cost_matrix[a][b]; //set a temporary matrix to be changed
			beta_matrix[a][b]=1000000;	//equivalent to infinity in Prim's description
		}
	}
	tempcost_matrix[vertices][vertices]=remove_loops(directed, tempcost_matrix);
	if (directed=='Y')
	{
		cout << "Prim - include arcs:"<<endl;
		gu->outputfile << "Prim - include arcs:"<<endl; 
	}
	else
	{
		cout << "Prim - include links:"<<endl;
		gu->outputfile << "Prim - include links:"<<endl; 
	}
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			if ( (gu->adj_cost_matrix[a][b]!=0) && (gu->adj_cost_matrix[a][b]<lowest_cost) )
			{//arbitrarily start with the cheapest link
				lowest_cost = gu->adj_cost_matrix[a][b];
				low_row = a;	//a and b are the two ends of the cheapest link
				low_column = b;
			}
		}
	}
	total_cost += gu->adj_cost_matrix[low_row][low_column];
	cout<<"("<<low_row+1<<" --> "<<low_column+1<<")";
	gu->outputfile<<"("<<low_row+1<<" --> "<<low_column+1<<") ";
	visited.push_back(low_row);//add the two endpoints of the cheapest link to
	visited.push_back(low_column);//the list of visited vertices
	nlinks++;
	//find all the correspondences and costs of the start vertex (low_row)
	for (a=0; a<gu->nvertices; a++)
	{
		if (gu->adj_cost_matrix[low_row][a] != 0) //new possibilities for links
		{
			beta_matrix[low_row][a]=gu->adj_cost_matrix[low_row][a];
		}
		if (directed=='Y')
		{
			if (gu->adj_cost_matrix[a][low_row] != 0)
			{
				beta_matrix[a][low_row]=gu->adj_cost_matrix[a][low_row];
			}
		}
	}
	//indicate that the link is used
	tempcost_matrix[low_row][low_column]=tempcost_matrix[low_column][low_row]=-1;
	beta_matrix[low_row][low_column]=beta_matrix[low_column][low_row]=999999; 
	do{
		lowest_cost = 999998;
		for (f=0; f<gu->nvertices; f++)
		{
			//find the correspondences of the two endpoints of a used link
			if (tempcost_matrix[low_column][f]>0)//greater than zero if a link exists and hasn't been used
			{
			beta_matrix[low_column][f]=gu->adj_cost_matrix[low_column][f];
			}
			if (tempcost_matrix[low_row][f]>0)
			{
				beta_matrix[low_row][f]=gu->adj_cost_matrix[low_row][f];
			}
			if (directed=='Y')		//check correspondences to and from
			{
				if (tempcost_matrix[f][low_column]>0)
				{
			beta_matrix[f][low_column]=gu->adj_cost_matrix[f][low_column];
				}
				if (tempcost_matrix[f][low_row]>0)
				{
				beta_matrix[f][low_row]=gu->adj_cost_matrix[f][low_row];
				}
			}
		}
		//find the next cheapest link off the subtree, add it to the tree
		for (g=0; g<gu->nvertices; g++)
		{
			for (h=0; h<gu->nvertices; h++)
			{
				if (beta_matrix[g][h]<lowest_cost)
				{
					lowest_cost = beta_matrix[g][h];
					low_row = g;		//the two ends of the cheapest link
					low_column = h;
				}
			}
		}
		nlinks++;	//add the link
		a=0;
		b=0;
		// check to see if the endpoints have been visited to avoid circuits
		for (v=visited.begin(); v!= visited.end(); v++)
		{
			if (low_row!=(*v))
			{
				a++;
			}
			if (low_column!=(*v))
			{
				b++;
			}
		}
		if ( (a==visited.size()) || (b==visited.size()) ) 
		//at most, one of the above can be true because Prim's algorithm grows from a single subtree
		{
			if (a==visited.size())//a (the low_row) has not been visited
			{
				visited.push_back(low_row);
			}
			else//b (the low_column) has not been visited
			{
				visited.push_back(low_column);
			}
			//add the link to the tree, add its cost to the total cost, indicate
			//that the link is used with -1
			total_cost += tempcost_matrix[low_row][low_column];
			cout<<"("<<low_row+1<<" --> " <<low_column+1<<")";
			gu->outputfile<<"("<<low_row+1<<" --> " <<low_column+1<<") ";
			tempcost_matrix[low_row][low_column]=tempcost_matrix[low_column][low_row]=-1;
			beta_matrix[low_row][low_column]=beta_matrix[low_column][low_row]=999999;
		}
		else//i.e. both endpoints of the next cheapest link have already been visited
		{
			nlinks--;//indicate that the link should not be added, but was considered
			tempcost_matrix[low_row][low_column]=tempcost_matrix[low_column][low_row]=-1;
			beta_matrix[low_row][low_column]=beta_matrix[low_column][low_row]=999999;
		}
	}while(nlinks<gu->nvertices-1);
	cout<<endl<<"total cost is: "<<total_cost<<endl;
	gu->outputfile<<endl<<"total cost is: "<<total_cost<<endl;
}
//************************************************************************
//*	Title: 			kruskal
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
//*	Function description:	Determines the minimum spanning tree using 
//*					Kruskal's algorithm. It uses a system that 
//*					assigns values to each vertex according to the
//*					root of the subtree it is a part of.  It uses 
//*					this information to decide if a link can be 
//*					added.  A link can be added, if either: 1.  One 
//*					or both of the vertices are not yet included in 
//*					any subtree, or 2. The two vertices are part of 
//*					different subtrees. 
//*
//*	Order of analysis:	O(nlinks*gu->nvertices^2)
//*
//*	Globals:			number of vertices (gu->nvertices), adjacency matrix
//* 					(adj_matrix)
//*
//*	Local Variables:		a temporary adjacency matrix (tempadj_matrix) 
//*
//*	Calling Modules:		main program
//************************************************************************
void kruskal(char directed)
{
	double possible;
	possible=spanning_trees(0);//check to see if the graph is connected(undirected only)
	if ((possible==0)&&(directed!='Y'))
	{
		cout<<"Kruskal: The graph is disconnected, so there are no spanning trees"<<endl;
		gu->outputfile<<"Kruskal: The graph is disconnected, so there are no spanning trees"<<endl;
		return;
	} 
	double tempcost_matrix[vertices][vertices], total_cost=0.0;//will store the cost values
	short a, b, c, d, c_visited=0, d_visited=0, old_root=0, nlinks=0;//mostly iterators
 	list <int> visited;//a list of the vertices that have been used/considered
	list <int>::iterator vi;
 	list <double> costList;//will be an ordered list of the costs, so that they can be used lowest to highest
	list <double>::iterator cl;
	int roots[vertices];//indicates the root of the subtree to which each vertex belongs
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			tempcost_matrix[a][b]=gu->adj_cost_matrix[a][b]; //initialize all costs
		}
		roots[a]=-1;//-1 indicates a root has not yet been assigned 
	}
	tempcost_matrix[vertices][vertices]=remove_loops(directed, tempcost_matrix);
	if (directed!='Y')
	{	//get lower triangular matrix so costs are only counted once (undirected only)
		for (a=0; a<gu->nvertices; a++)
		{
			for (b=0; b<=a; b++)
			{
				tempcost_matrix[a][b]=0;
			}	
		}
	}
	for (a=0; a<gu->nvertices; a++)
	{
		for (b=0; b<gu->nvertices; b++)
		{
			if (tempcost_matrix[a][b] != 0)
			{
				costList.push_back(tempcost_matrix[a][b]);//make a list of all the costs
			}
		}
	}
	costList.sort();//put all the costs in order to be used in sequence
	if (directed=='Y')
	{
		cout << "Kruskal - include arcs:"<<endl;
		gu->outputfile << "Kruskal - include arcs:"<<endl; 
	}
	else
	{
		cout << "Kruskal - include links:"<<endl;
		gu->outputfile << "Kruskal - include links:"<<endl; 
	}
do{
for (cl=costList.begin();cl!= costList.end();cl++)	//for #1
{
	for (c=0; c<gu->nvertices; c++)				//for #2
	{
		for (d=0; d<gu->nvertices; d++)			//for #3
		{
	if((tempcost_matrix[c][d]==(*cl))&&(nlinks<gu->nvertices-1))//if #1
//the connection from c to d is the cheapest left and there aren't too many links
		{
				nlinks++;//increment the number of links
			if (nlinks==1)//if the first link, it must be counted
			{
				roots[c]=roots[d]=c;// call c the root of this tree
				visited.push_back(c);//indicate these vertices have now been visited
				visited.push_back(d);
				total_cost+=tempcost_matrix[c][d];//add the cost of this link to the total
				cout<<"("<<c+1<<" --> "<<d+1<<")";
				gu->outputfile<<"("<<c+1<<" --> "<<d+1<<") ";
				tempcost_matrix[c][d]=tempcost_matrix[c][d]=-1;//indicate that a link has been used
			}
			else//there is more than one link			//else #1
			{
				c_visited=0;
				d_visited=0;
				for (vi=visited.begin(); vi!=visited.end(); vi++)
				{
					if (c!=(*vi))
					{
						c_visited++;//increment if c isn't the value on the visited list
					}
					if (d!=(*vi))
					{
						d_visited++;//increment if d isn't the value on the visited list
					}
				}
				if((c_visited==visited.size())||(d_visited==visited.size()))//if #2
				{//either c or d (or both) has (have) not been visited
					if((c_visited==visited.size())&&(d_visited==visited.size()))
					{//starting a new subtree(neither c nor d have been visited)
						roots[c]=roots[d]=c;//make c the root of this new tree
						visited.push_back(c);	//add the vertices to the list of visited
						visited.push_back(d);
					}
					else if(c_visited==visited.size())
					{		//c has not been visited, but d has
						roots[c]=roots[d];//add d-c link to the subtree of d
						visited.push_back(c);
					}
					else if(d_visited==visited.size())
					{		//d has not been visited, but c has
						roots[d]=roots[c];//add c-d link to the subtree of c
						visited.push_back(d);
					}
					total_cost+=tempcost_matrix[c][d];//add the cost
					cout<<"("<<c+1<<" --> "<<d+1<<")";
					gu->outputfile<<"("<<c+1<<" --> "<<d+1<<") ";
					tempcost_matrix[c][d]=tempcost_matrix[c][d]=-1;//indicate the link is used
				}										//end if #2
				else if (roots[c]!=roots[d]) //i.e. c and d have both been visited
				//two possibilities:1.joining two subtrees(allowable). . .
 				{ 						
					old_root=roots[d];
					for (b=0; b<gu->nvertices; b++)
					{
						if (roots[b]==old_root)	//i.e. d is considered the root of the 
						{							//subtree to which b belongs
							roots[b]=roots[c];//make c the root of the combined subtrees
						}
					}
					total_cost+=tempcost_matrix[c][d];//add the cost of the link
					tempcost_matrix[c][d]=tempcost_matrix[d][c]=-1;//indicate that the link is used
					cout<<"("<<c+1<<" --> "<<d+1<<")";
					gu->outputfile<<"("<<c+1<<" --> "<<d+1<<") ";
				}
 				else //2. Are in the same tree, i.e. creates a circuit (not allowed)
				{
					tempcost_matrix[c][d]=tempcost_matrix[d][c]=-1; //indicate the link is used
					nlinks--;		//remove link
				}		//end else 
			}	//end else #1
				}					//end if #1
			}						//end for #3
		}							//end for #2
	}								//end for #1
	
	}while(nlinks<gu->nvertices-1);
	cout<<endl<< "total cost is: " << total_cost << endl ;
	gu->outputfile<<endl<< "total cost is: " << total_cost << endl ;
}
//--------------------------end kruskal and prim
//--------------------------begin dijkstra
//************************************************************************
//*	Title: 				dijkstra
//*	Date:  				June 2000
//*	Version:			1.0
//*	
//*	Project:			Development of C++ graph package
//*
//*	Author: 			Matthew N. Johnson
//*	Coordinator: 		Professor Gary Lamont
//*
//*	Operating system:	Unix
//*	Language:  			C++
//*
//*	Function description:	Determines the minimum path between two vertices 
//*							using Dijkstra's algorithm. It labels vertices 
//*							according to the minimum of their current cost (initially
//*							set to infinity for all except the start) and the 
//*							cost to add their link to the current tree.  Once the pre-
//*							determined end vertex is labelled, its cost is the minimum 
//*							path from the start vertex to the end vertex.
//*				
//*	Order of analysis:	O(nlinks*gu->nvertices)
//*
//*	Globals:			number of vertices (gu->nvertices), adjacency/cost matrix
//* 					(gu->adj_cost_matrix)
//*
//*	Local Variables:	a temporarily labelled list (temporary) 
//*
//*	Calling Modules:	main program
//************************************************************************
//steps two and three in Dijkstra's algorithm
int dijkstra2_3(list<int>temporary, int p_vertex, double labels[vertices], double newlabels[vertices])
{
	short a, x_star=-1;	//iterator and x_star value (to be returned)
	double min_label=10000000;		//keeps track of minimum label
	list <int>::iterator temp;
	for (temp=temporary.begin(); temp!=temporary.end(); temp++)//go through all temporarily labeled vertices
	{
		if (gu->adj_cost_matrix[p_vertex][(*temp)]!=0)	//a temporarily labeled vertex is an element of 
		{											//the correspondence of the p_vertex
			if(labels[(*temp)]>(labels[p_vertex]+gu->adj_cost_matrix[p_vertex][(*temp)]))
			{
				labels[(*temp)]=labels[p_vertex]+gu->adj_cost_matrix[p_vertex][(*temp)];
			}//update its label, if necessary
		}
		if (labels[(*temp)]<min_label)
		{
			min_label=labels[(*temp)];	//keep track of the cheapest label
			x_star=(*temp);					//and the vertex with that label
		}
	}
	for (a=0; a<gu->nvertices; a++)
	{
		newlabels[a]=labels[a];		//the returned values
	}
	return x_star;
}
	
void dijkstra(short start_vertex, short end_vertex)//steps one, four and five in Dijkstra's algorithm
{
	short a, x_star, p_vertex;
	list <int> temporary;
	double labels[vertices];
	
	do{
		cout<<"What is the start vertex for Dijkstra's algorithm?"<<endl;
		cin>>start_vertex;
		if (start_vertex<1 || start_vertex>gu->nvertices)
		{
			cout<<"**Error: start vertex must be between 1 and "<<gu->nvertices<<".  Please try again."<<endl;
		}
	}while(start_vertex<1 || start_vertex>gu->nvertices);
	do{
		cout<<"What is the end vertex?"<<endl;
		cin>>end_vertex;
		if (end_vertex<1 || end_vertex>gu->nvertices)
		{
			cout<<"**Error: end vertex must be between 1 and "<<gu->nvertices<<".  Please try again."<<endl;
		}
	}while(end_vertex<1 || end_vertex>gu->nvertices);
	if (start_vertex == end_vertex)
	{
		cout<<"Dijkstra - Minimum cost from vertex "<<start_vertex<<" to vertex "<<end_vertex<<" is 0"<<endl;
		gu->outputfile<<"Dijkstra - Minimum cost from vertex "<<start_vertex<<" to vertex "<<end_vertex<<" is 0"<<endl;
		return;
	}	
	start_vertex--;end_vertex--;//to avoid having to subtract one for all referencing
	for (a=0; a<gu->nvertices; a++)
	{
		labels[a]=10000000;	//equivalent to infinity in algorithm
		temporary.push_back(a);
	}
	labels[start_vertex]=0;		//start vertex's label is 0
	temporary.remove(start_vertex);	//equivalent to setting the label to permanent
	p_vertex=start_vertex;
	do{
		x_star=dijkstra2_3(temporary,p_vertex,labels,labels);
		if (x_star==-1)		//no correspondences on the list
		{
			cout<<"A path does not exist between the two vertices"<<endl;
			gu->outputfile<<"A path does not exist between the two vertices"<<endl;
			return;
		}
		temporary.remove(x_star);	//set x_star to permanent
		p_vertex=x_star;			//make x_star the new p
	}while(x_star!=end_vertex);
	cout<<"Dijkstra - Minimum cost from vertex "<<start_vertex+1<<" to vertex "<<end_vertex+1<<" is "<<labels[end_vertex]<<endl;
	gu->outputfile<<"Dijkstra - Minimum cost from vertex "<<start_vertex+1<<" to vertex "<<end_vertex+1<<" is "<<labels[end_vertex]<<endl;
}
//---------------------------end Dijkstra's algorithm		

//---------------------------begin main program 
//*****************************************************************************
//*
//*	Title:		Main driver Program
//*	Date:		July 2000
//*	Version:	1.1 (in C++)
//*	Name:		main
//*	Function:	This is the main driver program for the graph package of 
//*			programs.  It executes the following functions:
//*			1.	Prints the opening title screen. (graph_menus_title)
//*			2.	Displays the help screen describing the program's 
//*				capabilities. (graph_menus_help)
//*			3.	Determines the method the user will use to enter the 
//*				graph data.  (graph_menus_method)
//*			4.	Determines the operation to perform on the graph data. 
//*				(graph_menus_functions)
//*			5.	Opens the input file.  (open_input_file) 
//*			6.	Computes the depth first search of the graph.
//*			7.	Computes the strong components of the graph. 
//*			8.	Computes the basis and contrabasis of the graph.
//*			9.	Computes the maximal independent sets and cliques of 
//*					the graph.
//*			10.	Computes the number of spanning trees of the graph.
//*			11.	Computes the minimum spanning tree of the graph using 
//*					Kruskal's algorithm.
//*			12.	Computes the minimum spanning tree of the graph using 
//*					Prim's algorithm.
//*			13. Computes the minimum path between two vertices using 
//*					Dijkstra's algorithm
//*
//*	Modules called: graph_menus_title, graph_menus_help, graph_menus_method, 
//*			    input_file, interactive_input, strong_component, dfs, 
//*			    MIS, clique, spanning trees, kruskal, prim, dijkstra
//*****************************************************************************
int main(){
	graph_menus_title();//display the title
	char choice = '?', function = '?', directed='?';
	short start = 0, start_vertex = 0, end_vertex = 0;
	string inputfile = "?";
do{
	do{
		choice=graph_menus_method();	//obtain input method
		if ('F' == choice)		//file input method
		{
			directed=file_input(inputfile);
		}
		if ('I' == choice)		//interactive input method
		{
			directed=interactive_input();
		}
		if ('H' == choice)		//help menu (if one exists)
		{
			graph_menus_help();
		}
		if ('X' == choice)		//exit
		{
			return choice;
		}
	}while((choice != 'X')&&(choice!='F')&&(choice!='I'));
	
   MISUtils mu(gu);
   
   do{
		function=graph_menus_functions();	//obtain function to be performed
		if ('S' == function)			//strong component
		{
			strong_component(directed);
		}
		if ('D' == function)			//depth-first search
		{			
			dfs(start);
		}
      if( 'C' == function ) {
         mu.clique(directed );
      }
		if ('M' == function)			//MIS and clique
		{
			mu.MIS(directed);
		}
      if( 'H' == function ) {
         mu.H_MIS(directed);
      }
		if ('T' == function)			//number of spanning trees
		{
			spanning_trees(1);
		}
		if ('K' == function)			//kruskal's algorithm
		{
			 kruskal(directed);
		}
		if ('P' == function)			//prim's algorithm
		{
			prim(directed);
		}
		if ('B' == function)			//both kruskal and prim
		{
			kruskal(directed);
			prim(directed);
		}
		if ('A' == function)			//all of the above functions
		{
			strong_component(directed);
			dfs(start);
			mu.MIS(directed);
         mu.H_MIS(directed);
			mu.clique(directed);
			spanning_trees(1);
			kruskal(directed);
			prim(directed);
		}
		if ('J'==function)	//dijkstra's algorithm
		{
			dijkstra(start_vertex, end_vertex);
		}
	}while ((function !='X')&&(choice != 'X')&&(function != 'N'));
		gu->outputfile<<endl;
}while((function !='X')&&(choice != 'X')||(function == 'N'));
	return 0;
}
