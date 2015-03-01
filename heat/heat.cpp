#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <ctime>
#include <cstdio>
#define np 8

double* init()
{
	double* A=(double*)calloc(sizeof(double),np*np);
	int i;
	for( i=0; i<np*np; i++ )
	 {
		  A[i]= i;
		  if(i>np) A[i]=i*2;
	 }
	return A;
}
inline int map(int row,int col,int sizex)
{
	return (row*sizex)+col;
}
void send(int ox,int oy,int yBlockSize,int xBlockSize, int sizex,int sizey,int i,double* A,int xSplit,int ySplit)
{
	int x=ox*yBlockSize;
	int y=oy*xBlockSize;
	int j;
	//borders (process 0,0)
	if(ox==0 && oy==0)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x+j,y,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
		
	}
	//process (0,xBlockSize-1) top right
	if(ox==0 && oy==xSplit-1)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x+j,y-1,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//top row, sans first and last (covered in the 2 conditions above)
	if(ox==0)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x+j,y-1,sizex)],xBlockSize+2,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//left border, process (ySplit-1,0) bottom left
	if(oy==0 && ox==ySplit-1)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x-1+j,y,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//left border sans bottom left
	if(oy==0)
	{
		for(j=0;j<yBlockSize+2;j++)
		{MPI_Send(&A[map(x-1+j,y,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//bottom row, process (ySplit-1,xSplit-1)
	if(ox == ySplit-1 && oy==xSplit-1)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x-1+j,y-1,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//bottom row wizout bottom right
	if(ox==ySplit-1)
	{
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Send(&A[map(x-1+j,y-1,sizex)],xBlockSize+2,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//right border, sans top n bottom extremes (already covered above)
	if(oy==xSplit-1)
	{
		for(j=0;j<yBlockSize+2;j++)
		{MPI_Send(&A[map(x-1+j,y-1,sizex)],xBlockSize+1,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
	}
	//the rest of the inner processes grid
	for(j=0;j<yBlockSize+2;j++)
		{MPI_Send(&A[map(x-1+j,y-1,sizex)],xBlockSize+2,MPI_DOUBLE,i,j,MPI_COMM_WORLD);}return;
}

double* receive(int ox,int oy,int yBlockSize,int xBlockSize, int sizex,int sizey,int xSplit,int ySplit,double* localTmp)
{int j;
//borders (process 0,0)
	if(ox==0 && oy==0)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}

	//process (0,xBlockSize-1) top right
	if(ox==0 && oy==xSplit-1)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//top row, sans first and last (covered in the 2 conditions above)
	if(ox==0)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+2),xBlockSize+2,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//left border, process (ySplit-1,0) bottom left
	if(oy==0 && ox==ySplit-1)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//left border sans bottom left
	if(oy==0)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+2));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+2));
		for(j=0;j<yBlockSize+2;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//bottom row, process (ySplit-1,xSplit-1) bottom right
	if(ox == ySplit-1 && oy==xSplit-1)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//bottom row wizout bottom right
	if(ox==ySplit-1)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+1));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+1));
		for(j=0;j<yBlockSize+1;j++)
		{MPI_Recv(localA+j*(xBlockSize+2),xBlockSize+2,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//right border, sans top n bottom extremes (already covered above)
	if(oy==xSplit-1)
	{
		double* localA=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+2));
		localTmp=(double*)calloc(sizeof(double),(xBlockSize+1)*(yBlockSize+2));
		for(j=0;j<yBlockSize+2;j++)
		{MPI_Recv(localA+j*(xBlockSize+1),xBlockSize+1,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;
	}
	//the rest of the inner processes grid
	double* localA=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+2));
	localTmp=(double*)calloc(sizeof(double),(xBlockSize+2)*(yBlockSize+2));
	for(j=0;j<yBlockSize+2;j++)
		{MPI_Recv(localA+j*(xBlockSize+2),xBlockSize+2,MPI_DOUBLE,0,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);}return localA;


}

double operate(int ox,int oy,int yBlockSize,int xBlockSize,double* localA,double* localTmp,int xSplit,int ySplit)
{int i,j;double sum=0;
	if(ox==0 && oy==0)//top left
	{
			for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize;j++)
			{
				
				localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
				sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			}
		}return sum;
	}

	//process (0,xBlockSize-1) top right
	if(ox==0 && oy==xSplit-1)
	{
		for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize;j++)
			{
			localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
				sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			}
		}
	return sum;
	}
	//top row, sans first and last (covered in the 2 conditions above)
	if(ox==0)
	{
	
		for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize+1;j++)
			{
				localTmp[i*(xBlockSize+2)+j]=0.25*
					(localA[i*(xBlockSize+2)+j-1]+ //left
					localA[i*(xBlockSize+2)+j+1]+//right
					localA[(i-1)*(xBlockSize+2)+j]+//top
					localA[(i+1)*(xBlockSize+2)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j])*(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j]);
			}
		}
		return sum;
	}
	//left border, process (ySplit-1,0) bottom left
	if(oy==0 && ox==ySplit-1)
	{
	
		for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize;j++)
			{
				localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			}
		}
	return sum;
	}

	//left border sans bottom left
	if(oy==0)
	{
	for(i=1;i<yBlockSize+1;i++)
		{
			for(j=1;j<xBlockSize+1;j++)
			{
				localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			}
		}
	
	return sum;}

	//bottom row, process (ySplit-1,xSplit-1) bottom right
	if(ox == ySplit-1 && oy==xSplit-1)
	{
	for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize;j++)
			{
				localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			}
		}
	return sum;}

	//bottom row wizout bottom right
	if(ox==ySplit-1)
	{
		for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize+1;j++)
			{
				localTmp[i*(xBlockSize+2)+j]=0.25*
					(localA[i*(xBlockSize+2)+j-1]+ //left
					localA[i*(xBlockSize+2)+j+1]+//right
					localA[(i-1)*(xBlockSize+2)+j]+//top
					localA[(i+1)*(xBlockSize+2)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j])*(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j]);
			
			}
		}
	
	return sum;}

	//right border, sans top n bottom extremes (already covered above)
	if(oy==xSplit-1)
	{
	
		for(i=1;i<yBlockSize;i++)
		{
			for(j=1;j<xBlockSize;j++)
			{
				localTmp[i*(xBlockSize+1)+j]=0.25*
					(localA[i*(xBlockSize+1)+j-1]+ //left
					localA[i*(xBlockSize+1)+j+1]+//right
					localA[(i-1)*(xBlockSize+1)+j]+//top
					localA[(i+1)*(xBlockSize+1)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j])*(localTmp[i*(xBlockSize+1)+j] - localA [i*(xBlockSize+1)+j]);
			
			}
		}

	return sum;}

	//the rest of the inner processes grid
		for(i=1;i<yBlockSize+1;i++)
		{
			for(j=1;j<xBlockSize+1;j++)
			{
				localTmp[i*(xBlockSize+2)+j]=0.25*
					(localA[i*(xBlockSize+2)+j-1]+ //left
					localA[i*(xBlockSize+2)+j+1]+//right
					localA[(i-1)*(xBlockSize+2)+j]+//top
					localA[(i+1)*(xBlockSize+2)+j]);//bottom
					sum+=(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j])*(localTmp[i*(xBlockSize+2)+j] - localA [i*(xBlockSize+2)+j]);
			}
		}
		return sum;

}
void jacobi(double* A, double* tmp, int sizex, int sizey,int xSplit,int ySplit,int argc,char* argv[])
{
	int rank,threads;int coords[2];
	int i=0;int j=0;
	MPI_Request send_request,recv_request;
	//Init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &threads);
	//x and y BlockSize are the size of each dimension of a single block that each process should loop on, this is useful bcoz a block can be identified by the index of its first element
	int xBlockSize=sizex/xSplit;
	int yBlockSize=sizey/ySplit;
    //create virtual topolgy 
	MPI_Comm vt;
	int dim[2], cycle[2], reorder;
	dim[0]=xSplit; 
	dim[1]=ySplit;
	cycle[0]=false; cycle[1]=false;//acyclic
	reorder=false;//same ordering of ranks
	MPI_Cart_create(MPI_COMM_WORLD,2,dim,cycle,reorder,&vt);
	//get coords from rank
	MPI_Cart_coords(vt,rank,2,coords);
	//printf("my rank is %d and coords are %d,%d\n",rank,coords[0],coords[1]);

	if(rank==0)
	{
		//send all chunks of A to all processes,
		//printf("total of %d threads\n",threads);
		for(i=1;i<threads;i++)
		{
				MPI_Cart_coords(vt,i,2,coords);
				//sending is done row by row, because rowns and columns are not sequential in memory
				send(coords[0],coords[1],yBlockSize,xBlockSize,sizex,sizey,i,A,xSplit,ySplit);
			
		}
		MPI_Cart_coords(vt,0,2,coords);
		//operate
	
		//operate(coords[0],coords[1],yBlockSize,xBlockSize,localA,localTmp,xSplit,ySplit);

		//recieve results chunks from processees
	}
	else
	{ 
		double* localTmp=(double*)calloc(sizeof(double),1);
		double* localA=receive(coords[0],coords[1],yBlockSize,xBlockSize,sizex,sizey,xSplit,ySplit,localTmp);
		//printf("my rank is %d my coord are %dx%d the first item is %f second row item is %f\n",rank,coords[0],coords[1],localA[0],localA[xBlockSize+1]);
		/*
		printf("rank %d\n", rank);
		for(int i=0;i<yBlockSize;i++)
		{
			for(int j=0;j<xBlockSize+1;j++)
			{
				printf("%f ",localA[j+i*(xBlockSize+1)]);
			}
			printf("\n");
		}*/
		
		//operate on it
		double sum=operate(coords[0],coords[1],yBlockSize,xBlockSize,localA,localTmp,xSplit,ySplit);
	

		//send it to 0
	}



	

	MPI_Finalize();
}

int main(int argc, char* argv[])
{
	double* A=init();
	double* tmp=(double*)calloc(sizeof(double),np*np);
	static const char filename[] = "dist.txt";
    FILE *file = fopen ( filename, "r" );
	char line[4];
	fgets(line, sizeof line,file);
	int xSplit=line[0]-'0';
	int ySplit=line[2]-'0';
	jacobi(A,tmp,np,np,xSplit,ySplit,argc,argv); 
}