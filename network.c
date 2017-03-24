/*
 *  Created on: Mar 22, 2016
 *      Author: Simon Reynders(260502362)
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

struct node{
	bool inTree;
};

struct edge{

	bool isFirst,isLast;
	int forwardSource,forwardDestination,backwardsSource,backwardsDestination,_multiplicity;
	float edgeCost;
	float edgeReliability;
	float edgeEff;
	struct edge * ptrStart;
	struct edge * ptrEnd;
	struct edge * ptrNext;
	struct edge * ptrPrv;
};

struct minSpanTree{

	//Structure for holding all elements in the minimum spanning tree


	//the parameters of the problem
	float paramaterCost;
	float paramaterReliability;

	//number of nodes in the tree at any give moment
	//max number of nodes allowed in the tree considering N
	int nodesInTree;
	int maxPossibleNodes;

	//The running cost and reliability of the edges in the tree
	float networkCost;
	float runningReliability;

	//Largest and smallest costs
	float largestCost;
	float smallestCost;

	//pointers to the edges that compose the minimum spanning tree
	struct edge * treeBegin;
	struct edge * treeEnd;
	struct edge * lastEntryFromLL;
};

void printNetworkMatrix(int ELEMENTS,struct minSpanTree * mst);

int getNumberOfElements(FILE * text);
struct edge * readFromFile(FILE * text,int ELEMENTS,struct edge * edgePtr,float *A_B,float *REQ_COST,float *REQ_RELIABILITY);
struct edge * readArrayRow(char firstRow[],int COLUMN,int ROW,int ELEMENTS,struct edge * edgePtr,struct edge * prvEdge,struct edge * startEdge,bool firstCall,bool cost);
float readValue(FILE * text,char BUFFER[],char NUMERIC_KEY[]);

struct minSpanTree * initMinTree(struct node * nodeArray,float paramCost,float paramRel,int elements);
void findPossibleMST(int seed,struct minSpanTree * mst,struct edge * edgeListStart,struct node * nodeArray,struct edge * treeElement);
struct edge * seekEdge(int seed,struct edge * edgePtr);
void efficiencyProblem(struct edge * edgePtr);
struct edge * addNode(bool initialNode,int number,struct node * nodeArray,struct minSpanTree * mst,struct edge * edgeFromGraph,struct edge * prvLeaf);

int addParallelEdge(struct edge * edgeToAdd,struct minSpanTree* mst);
void refineSolution(struct minSpanTree * mst,int problemType);
struct edge * findCheapestEdge(struct minSpanTree * mst);
struct edge * findBestEdge(struct minSpanTree * mst,struct edge * prvEdge);

float parallelReliability(struct edge * scan,int multiplicity,float edgeReliability);
float currentNetEff(struct minSpanTree * mst, struct edge * consider);
float reCalculateRel(struct edge * edgePtr,struct minSpanTree * mst);
struct edge * findEdgeMultiplicity(struct minSpanTree * mst);


int main(){


	//prompts user to provide a text file
	FILE * text;
	char toOpen[30];
	char quit[]="exit";
	printf("Please enter file you want processed or type exit to quit: \n");
	gets(toOpen);
	if(!strcmp(toOpen,quit)){
		exit(1);
	}
	text=fopen(toOpen,"r");
	if(text==NULL){
		printf("Text file is either missing or has wrong name.\n");
		exit(0);
	}

	//gets the first thing in the text file: N
	int elements=getNumberOfElements(text);

	//init. an array of N elements which will represent the nodes in the minimum spanning tree(MST)
	struct node nodeArray[elements];
	int i;
	for(i=0;i<elements;i++){
		//set all the nodes to not be in the tree
		nodeArray[i].inTree=false;
	}
	float A_B=0.0,REQ_COST=0.0, REQ_RELIABILITY=0.0;

	//we allocate an edge object; this will be the first edge in the MST;
	//we map the ptr elements for it and then read information from the file
	//the information is put into many edge structures and they are connected
	//to form a linked list.
	struct edge * edgePtr=malloc(sizeof(struct edge));
	edgePtr->ptrStart=edgePtr;
	edgePtr->ptrPrv=edgePtr;
	edgePtr->isFirst=true;
	readFromFile(text,elements,edgePtr,&A_B,&REQ_COST,&REQ_RELIABILITY);

	//done with the file so clost it;
	fclose(text);

	//get the efficiency of each edge (reliability/cost);
	efficiencyProblem(edgePtr);

	//init. the MST
	struct minSpanTree * mst=initMinTree(nodeArray,REQ_COST,REQ_RELIABILITY,elements);

	//find a MST
	findPossibleMST(0,mst,edgePtr,nodeArray,NULL);

	//checks if the MST leaves us with a negative or positive cost budget
	if(mst->networkCost>mst->paramaterCost){
		printf("No minimum spanning tree with cost less than cost cap exists. Program will now terminate\n");
		exit(0);
	}

	//prints the MST
	printNetworkMatrix(elements,mst);
	printf("Network cost: %3.f \t Reliability of network: %3.3f\n\n\n",mst->networkCost,mst->runningReliability);

	//refines the MST; adds parallel edges
	refineSolution(mst,A_B);

	//print the refined network
	printNetworkMatrix(elements,mst);
	printf("Network cost: %3.f \t Reliability of network: %3.3f\n",mst->networkCost,mst->runningReliability);

	//free all that memory
	struct edge * toFree=mst->treeEnd;
	struct edge * toFree2;
	while(toFree!=NULL){
		toFree2=toFree->ptrPrv;
		free(toFree);
		toFree=toFree2;
	}
	toFree=edgePtr->ptrEnd;
	while(toFree->isLast){
		toFree2=toFree->ptrPrv;
		free(toFree);
		toFree=toFree2;
	}
	free(toFree);
	free(toFree2);
	free(edgePtr);
	free(nodeArray);

return(0);

}


void refineSolution(struct minSpanTree * mst,int problemType){
	//Depending on the problem type, we stop adding edges when certain
	//conditions are met.
	/*
	 * If we are looking to make the most reliable network while remaining
	 * under the cost cap, problemType=0, otherwise problemType=1.
	 * The terminating conditions for the loops reflect these problem types.
	 * If problemt type 0, we want to stop as soon as we reach a certain reliability
	 * parameter. If looking for a solution to problem type 1, we continue looking
	 * for the best edge to add until we can no longer add without breaching the
	 * network cost cap.
	 *
	 * Note that if addParallelEdge in either case returns a 0, then we have reached
	 * some parameter cap and can no longer continue adding edges in parallel.
	 */
	int k=1;
	if(problemType<1){
		while(mst->runningReliability<mst->paramaterReliability && k==1){
			k=addParallelEdge(findCheapestEdge(mst),mst);
		}
	}
	else{
		while(mst->networkCost<mst->paramaterCost && k==1){
			k=addParallelEdge(findCheapestEdge(mst),mst);
		}
	}
}

struct edge * findCheapestEdge(struct minSpanTree * mst){

	/*
	 * We are looking for the next best edge to add here.
	 * We have two place holders: one for the best item seen yet
	 * and one for the scan that changes with each edge in the linked list.
	 * We initialize them to the first item in the tree.
	 * Calculate the cost per percentage reliability gain of making an edge parallel
	 * and compare each edge with these values. The function parallelReliability finds
	 * the reliability of that edge when its multiplicity is increased and thus made parallel.
	 */
		struct edge * scan=mst->treeBegin;
		struct edge * bestPresentEdge=scan;
		float bestPresentRel=parallelReliability(scan,scan->_multiplicity,scan->edgeReliability);
		float bestPresentEff=bestPresentEdge->edgeCost/((bestPresentRel-bestPresentEdge->edgeReliability)/bestPresentEdge->edgeReliability);
		while(scan!=NULL){
			//continue through the tree until we hit a NULL character

			float nextRel=parallelReliability(scan,scan->_multiplicity,scan->edgeReliability);
			float nextEff=scan->edgeCost/((nextRel-scan->edgeReliability)/scan->edgeReliability);
			if(bestPresentEff>nextEff && scan->_multiplicity<3){
				//if we find a better element and its multiplicity is not yet 3, take it
				//else move on.
				bestPresentRel=nextRel;
				bestPresentEff=nextEff;
				bestPresentEdge=scan;
			}
			scan=scan->ptrNext;
		}
		//This is here to catch the case that no edge in the previous if statement met the
		//condition. In this case, the bestPresentEdge is set to the first edge in the tree
		//regardless of cost or multiplicity. This is here to make sure that edge would not get
		//returned, but a NULL character would. For example, if all edges were added in parallel
		//three times, the bestPresentEdge would be set to an edge with multiplicty 3 and we cannot
		//have that, so we must set it to a NULL character.
		if(bestPresentEdge->edgeCost+mst->networkCost<mst->paramaterCost && !(bestPresentEdge->_multiplicity>2)){
			bestPresentEdge->edgeReliability=bestPresentRel;
			bestPresentEdge->_multiplicity++;
			return(bestPresentEdge);
		}
		else{
			return(NULL);
		}
}
int addParallelEdge(struct edge * edgeToAdd,struct minSpanTree* mst){

	/*The best edge found in the findCheapestEdge function must be added to the
	 * MST. To do this we need to update the reliability of the tree to reflect the
	 * change in reliability made by adding an edge in parallel. Additionally, the cost
	 * must be updated.
	 * A check is in place to reflect the fact that findCheapestEdge may have returned a
	 * NULL character. If this is so, this function returns a 0, otherwise it ran to completion
	 * and returns a 1.
	 */

	if(edgeToAdd==NULL){
		return(0);
	}
	struct edge * step=mst->treeEnd;
	mst->runningReliability=1;
	mst->networkCost=mst->networkCost+edgeToAdd->edgeCost;
	while(step!=NULL){
		mst->runningReliability=mst->runningReliability*step->edgeReliability;
		step=step->ptrPrv;
	}
	return(1);
}

float parallelReliability(struct edge * scan,int multiplicity,float edgeReliability){

	//a function to calculate reliability of parallel edges. If the edge input is and edge
	//with multiplicity one, it finds the reliability if that edge had multiplicity 2.

	multiplicity++;
	int i;
	float pi=1;
	for(i=0;i<multiplicity;i++){
		pi=pi*(1-edgeReliability);
	}
	return(1-pi);
}


struct minSpanTree * initMinTree(struct node * nodeArray,float paramCost,float paramRel,int elements){

	/*This initializes the MST. Nothing special.
	 * This is more here for cosmetic issues. Didn't like
	 * having all of these things in main();
	 */


	struct minSpanTree * minTree=malloc(sizeof(struct minSpanTree));
	minTree->nodesInTree=0;
	minTree->paramaterCost=paramCost;
	minTree->paramaterReliability=paramRel;
	minTree->networkCost=0;
	minTree->runningReliability=1;
	minTree->maxPossibleNodes=elements;
	return(minTree);
}


struct edge * addNode(bool initialNode,int number,struct node * nodeArray,struct minSpanTree * mst,struct edge * edgeFromGraph,struct edge * prvLeaf){

	/*This function does a lot of work.
	 * From the function findPossibleMST, it will ad a node to the MST.
	 * If findPossibleMST is called for the first time, there are no
	 * edges in the MST yet: we must add it. In that case, initialNode will
	 * be true and the code to get the ball rolling on the MST happens below:
	 * -A node gets added to the list according to Prim's algo. ( I choose node 0)
	 * -The ptr to the beginning of the tree is set to the edge I have allocated here.
	 * -->NOTE, THAT EDGE IS THE FIRST AND SUBSEQUENT EDGE ELEMENT IN THE MST.
	 * -The edge is then returned, as all we have done is included a node in the tree
	 * -That edge will have its date fields filled one we have found a suitable edge from node 0
	 * -nodeArray[number] represents a node getting add to the tree
	 *
	 * If initialNode is false:
	 * -The tree has already been initialized and we have found a suitable edge to add to the tree
	 * -The information is to be saved in an edge given to the tree through the parameter prvLeaf
	 * -prvLeaf was made in the previous call to addNode. In fact, if this is the first edge we are
	 * adding, then initialNode was true the first time this function was called and we are filling
	 * the pointer to the edge structure it returned.
	 * -All the required fields in the MST are changed accordingly
	 * -The edge elements and pointers are set accordingly
	 * -Another edge element was created upon the calling of addNode; this way we can set
	 * prvLeaf's pointers to the next edge item in the MST here. We set prvLeaf->next to point to the
	 * edge *ptr we have just made. This was we make our linked list. We can easily enough do the reverse
	 * and set ptr->previous to point to prvLeaf. Thus we have a doubly linked list.
	 *
	 */

	struct edge * ptr=malloc(sizeof(struct edge));
	if(initialNode){
		mst->largestCost=0;
		mst->smallestCost=mst->paramaterCost;
		mst->nodesInTree++;
		mst->treeBegin=ptr;
		ptr->ptrPrv=NULL;
		nodeArray[number].inTree=true;
		return(ptr);
	}
	else{
		mst->nodesInTree++;
		mst->networkCost=mst->networkCost+edgeFromGraph->edgeCost;
		mst->runningReliability=mst->runningReliability*edgeFromGraph->edgeReliability;
		mst->treeEnd=prvLeaf;

		if(edgeFromGraph->edgeCost>mst->largestCost){
			mst->largestCost=edgeFromGraph->edgeCost;
		}
		if(edgeFromGraph->edgeCost<mst->smallestCost){
			mst->smallestCost=edgeFromGraph->edgeCost;
		}
		prvLeaf->forwardDestination=edgeFromGraph->forwardDestination;
		prvLeaf->_multiplicity=edgeFromGraph->_multiplicity;
		prvLeaf->forwardSource=edgeFromGraph->forwardSource;
		prvLeaf->backwardsDestination=edgeFromGraph->backwardsDestination;
		prvLeaf->backwardsSource=edgeFromGraph->backwardsSource;
		prvLeaf->edgeCost=edgeFromGraph->edgeCost;
		prvLeaf->edgeEff=edgeFromGraph->edgeEff;
		prvLeaf->edgeReliability=edgeFromGraph->edgeReliability;
		prvLeaf->ptrNext=ptr;
		prvLeaf->ptrStart=mst->treeBegin;
		ptr->ptrPrv=prvLeaf;
		ptr->ptrNext=NULL;
		nodeArray[edgeFromGraph->forwardDestination].inTree=true;
		nodeArray[edgeFromGraph->forwardSource].inTree=true;


		//HERE we free the memory of the edge we are adding to the MST
		//I realize we could just move the pointers around so the edge from
		//the graph just becomes linked to the items in the MST but I couldn't
		//get it to work and it was late.
		struct edge * temp;
		temp=edgeFromGraph->ptrPrv;
		temp->ptrNext=edgeFromGraph->ptrNext;
		temp=edgeFromGraph->ptrNext;
		edgeFromGraph->ptrPrv=temp;

		return(ptr);

	}
}


void findPossibleMST(int seed,struct minSpanTree * mst,struct edge * edgeListStart,struct node * nodeArray,struct edge * emptyElement){

	/*
	 * Finds a minimum spanning tree. First there are no elements in the minimum spanning tree so
	 * we call addNode which adds a node for us and returns an edge element for us to fill with
	 * the edge we will add when we find a suitable edge to connect in the network.
	 * That edge element will be passed to the next call of this function.
	 */

		if(mst->nodesInTree==0){
			struct edge * ptr=malloc(sizeof(struct edge));
			ptr=addNode(true,seed,nodeArray,mst,NULL,NULL);
			findPossibleMST(seed,mst,edgeListStart,nodeArray,ptr);
		}

		/*
		 * Here, addNode has been called once to add a node to the tree. Now we must find an edge
		 * that goes from the node in the tree to a node not in the tree which meets the criteria
		 * to add. This criteria is found in the if statement below and loops while there are items
		 * in the linked list of edges.
		 * The condition for the if statement is rather long but can be thought of as follows:
		 * If we say node 0 is in the tree and node 1 is not, we may want to add that edge, but
		 * as we scan through the linked list of edges, we may come accross an edge that originates
		 * from 1 and ends in 0. However as this is an undirected graph, an edge from 1 to 0 is the same
		 * as an edge from 0 to 1, so we must consider it.
		 * The mathematic criteria for adding an edge is if its reliability to cost ratio is the highest
		 * we can find and its cost is also low.
		 */

		if(mst->nodesInTree<mst->maxPossibleNodes){
			//look through LL to find cheapest edge, with a source in MST, with destination not in MST
			struct edge * index=edgeListStart;
			struct edge * presentBest=edgeListStart;
			while(!index->isLast){
				if(((nodeArray[index->forwardSource].inTree && !nodeArray[index->forwardDestination].inTree)||(nodeArray[index->backwardsSource].inTree && !nodeArray[index->backwardsDestination].inTree))  && index->edgeEff>presentBest->edgeEff && index->edgeCost<presentBest->edgeCost){
					presentBest=index;
				}
				index=index->ptrNext;
			}
			struct edge * ptr=malloc(sizeof(struct edge));
			ptr=addNode(false,seed,nodeArray,mst,presentBest,emptyElement);
			findPossibleMST(0,mst,edgeListStart,nodeArray,ptr);
		}


}

void efficiencyProblem(struct edge * edgePtr){

	//sets the edges edgeEff value to edge's reliability/edge's cost

	while(!edgePtr->isLast){
		edgePtr->edgeEff=edgePtr->edgeReliability/edgePtr->edgeCost;
		edgePtr=edgePtr->ptrNext;
	}

}



struct edge * readFromFile(FILE * text,int ELEMENTS,struct edge * edgePtr,float *A_B,float *REQ_COST,float *REQ_RELIABILITY){

	//This simply reads information from the file.

	//Bools are used to determine if reading the cost matrix or reliability matrix
	bool R_COST=false;
	bool R_REQ=false;

	//input buffer
	char BUFFER[256];
	int COLUMN=1;
	int ROW=0;

	//numeric and alphabetical keys for text files
	char NUMERIC_KEY[]="0123456789";
	char ALPHA_KEY[]="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	//points to the last element of the linked list
	struct edge * lastPos;

	while(EOF!=fscanf(text,"%s",BUFFER)){

		char * pointerToItem;

		//if the buffer has anything in it from the alpha key...
		if((strpbrk(BUFFER,ALPHA_KEY)!=NULL)){

			//if buffer has the element "C=["...
			if((pointerToItem=strstr(BUFFER,"C=["))){

				//we are looking at the cost matrix so set the r_cost true
				//readArrayRow into the edge structures
				R_COST=true;
				lastPos=readArrayRow(pointerToItem+3,COLUMN,ROW,ELEMENTS,edgePtr,edgePtr->ptrPrv,edgePtr->ptrStart,true,R_COST);//reads the rowString element by element and places in the cost matrix
				ROW++;
			}
			if((pointerToItem=strstr(BUFFER,"R=["))){

				//same as above

				ROW=0;
				COLUMN=1;
				R_REQ=true;
				lastPos=readArrayRow(pointerToItem+3,COLUMN,ROW,ELEMENTS,edgePtr,NULL,NULL,true,R_COST);
				ROW++;
			}
			//if the buffer has the "a_b"...
			if((strstr(BUFFER,"a_b"))!=NULL){
				*A_B=readValue(text,BUFFER,NUMERIC_KEY);
			}
			//if the buffer has "req_reliability"...
			if((pointerToItem=strstr(BUFFER,"Req_Reliability"))){
				*REQ_RELIABILITY=readValue(text,BUFFER,NUMERIC_KEY);
			}
			//if the buffer has "req_cost"...
			if((pointerToItem=strstr(BUFFER,"Req_Cost"))){
				*REQ_COST=readValue(text,BUFFER,NUMERIC_KEY);
			}
		}
		else{
			//otherwise, we are still reading elements from the cost or reliability matrices

			//if we are reading from the cost matrix, r_cost will be true and r_req will be false
			if(R_COST && (!R_REQ)){
				COLUMN++;
				//continues reading in to edge elements
				//we pass it the last accessed edge, pointed to by lasPos
				//at the end of readArrayRow, we get a new edge structure pointed to by lastPos
				lastPos=readArrayRow(BUFFER,COLUMN,ROW,ELEMENTS,lastPos,lastPos->ptrPrv,lastPos->ptrStart,true,R_COST);
				if(ROW==ELEMENTS-2){
					//if we have read all there is to read of the cost matrix...
					//set r_cost to false because we are no longer reading cost matrices and set the lastPos structure to be the
					//terminating link in the linked list of edges.
					R_COST=false;
					edgePtr->ptrEnd=lastPos;
					lastPos->ptrNext=lastPos->ptrStart;
					edgePtr->ptrEnd->isLast=true;
				}
				ROW++;
			}
			//same as above but for the reliability matrix
			if((!R_COST) && R_REQ){
				COLUMN++;
				lastPos=readArrayRow(BUFFER,COLUMN,ROW,ELEMENTS,lastPos,lastPos->ptrPrv,lastPos->ptrStart,true,R_COST);
				if(ROW==ELEMENTS-2){
					R_REQ=false;
				}
				ROW++;
			}
		}
		//clear the buffer each time
		memset(BUFFER,'\0',256);
	}
	//return the pointer that points to the first edge in the linked list of edges
	return(edgePtr);
}

float readValue(FILE * text,char BUFFER[],char NUMERIC_KEY[]){
	//reads single values
	char * pointerToItem;
	while((strpbrk(BUFFER,NUMERIC_KEY)==NULL)){
		fscanf(text,"%s",BUFFER);
	}
	pointerToItem=strpbrk(BUFFER,NUMERIC_KEY);
	float number=strtof(pointerToItem,NULL);
	return(number);
}

struct edge * readArrayRow(char firstRow[],int COLUMN,int ROW,int ELEMENTS,struct edge * edgePtr,struct edge * prvEdge,struct edge * startEdge,bool firstCall, bool cost){
	//reads a row of the cost/reliability matrices

	//if we are reading from the cost matrix...
	if(cost){
		char * VALUE_POINTER;

		//break the string into tokens delimited by commas
		VALUE_POINTER=strtok(firstRow,",");
		struct edge * temp=edgePtr;
		struct edge * prv=prvEdge;
		do{
			//create a new edge structure for the next inevitable edge
			//set the edge pointed to by temp information:
			/*-source
			 * -destination
			 * -cost
			 * -reliability
			 * -multiplicity
			 * -the start pointer which points to the first edge in the edge's linked list
			 * -the next pointer should point to the next edge element which will be newEdge
			 *-the previous pointer should point to the previous edge element which will be pointed to by pointer prv
			*/
			struct edge * newEdge=malloc(sizeof(struct edge));
			float number=strtof(VALUE_POINTER,NULL);
			newEdge->isFirst=false;
			newEdge->isLast=false;
			temp->isLast=false;
			temp->edgeCost=number;
			temp->forwardSource=ROW;
			temp->forwardDestination=COLUMN;
			temp->backwardsDestination=ROW;
			temp->backwardsSource=COLUMN;
			temp->_multiplicity=1;
			temp->ptrStart=startEdge;
			temp->ptrNext=newEdge;
			temp->ptrPrv=prv;

			prv=temp;
			temp=newEdge;

			COLUMN++;
			VALUE_POINTER=strtok(NULL,",");
		}while(VALUE_POINTER!=NULL);

		//if we have reached the end of row of that matrix, return temp
		//this will be stored as lastPos in readFromFile();
		temp->ptrPrv=prv;
		temp->ptrStart=prv->ptrStart;
		return(temp);
	}
	else{

		//if we come here, we are reading from the reliability matrix
		char * VALUE_POINTER;
		VALUE_POINTER=strtok(firstRow,",");

		do{
			float number=strtof(VALUE_POINTER,NULL);
			edgePtr->edgeReliability=number;
			edgePtr=edgePtr->ptrNext;
			COLUMN++;
			VALUE_POINTER=strtok(NULL,",");
		}while(VALUE_POINTER!=NULL);
		return(edgePtr);
	}
}

int getNumberOfElements(FILE * text){

	//gets the number of elements N
	char N_VALUE[5];
	int ELEMENTS;
	if(text==NULL){
		perror("Error opening file");
		exit(0);
	}
	fscanf(text,"%s",N_VALUE);
	sscanf(N_VALUE, "%*c%*c%d",&ELEMENTS);
	return(ELEMENTS);
}

void printNetworkMatrix(int ELEMENTS,struct minSpanTree * mst){

	//just prints the network matrix

	int i,j;

	int output[ELEMENTS][ELEMENTS];
	for(i=0;i<ELEMENTS;i++){
		for(j=0;j<ELEMENTS;j++){
			output[i][j]=0;
		}
	}
	struct edge * thing;
	thing=mst->treeEnd;
	while(thing!=NULL){
		output[thing->forwardSource][thing->forwardDestination]=thing->_multiplicity;
		output[thing->forwardDestination][thing->forwardSource]=thing->_multiplicity;
		thing=thing->ptrPrv;
	}

	for(i=0;i<ELEMENTS;i++){
		printf("%3.1i-[",i);
		for(j=0;j<ELEMENTS;j++){
			if(j<ELEMENTS-1){
				printf("%i,",output[i][j]);
			}
			else{
				printf("%i]",output[i][j]);
			}
		}
		printf("\n");
	}
}

