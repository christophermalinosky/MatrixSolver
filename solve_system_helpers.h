#ifndef SOLVE_SYSTEM_HELPERS
#define SOLVE_SYSTEM_HELPERS

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

//Defines a matrix struct.  The matrix has m rows and n columns.
typedef struct matrix{
  int m;        //Number of rows in the matrix
  int n;        //Number of columns in the matrix

  double* data; //This is a pointer to the array data.  This stores the data as a linear array, so
                //accsible as data[k].  The data is stored in row-major order, so the array is laid out as
                //a_11 a_12 a_13 ...  a_1N a_21 a_22 a_23 ... a_2N ... a_NN. 

  double** entry; //points the the data of the matrix, but represents a 2-D array.  So, entry[i][j] gives the entry of the matrix
                  //in row i, column j.  Incidentally, entry[i] is the address of entry[i][0].
                  //Also, keep in mind that everything here is 0-indexed. 
  
}matrix;


//Correctly set upt the entry array in a matrix if its m,n,and data entries are already set up.
//Not hard.  Just tedious. 
void setup_entry_array(matrix* A){
  A->entry=malloc(sizeof(double*) * A->m);
  int i;
  for (i=0;i<A->m;i++){
    A->entry[i]=A->data+i*A->n;
  }
}

//Read a matrix out of a plain text file.
//text file must contain the dimenstions of the matrix in the first line,
//and then the matrix data thereafter. 

matrix read_matrix(char* filename){

  FILE* fp=fopen(filename,"r");

  if (fp==NULL){
    fprintf(stderr,"Error opening %s for reading. \n",filename);
    exit(1);
  }
  matrix A;

  if (fscanf(fp,"%d,%d",&(A.m),&(A.n))!=2){
    fprintf(stderr,"Error reading line 1 of %s.\n",filename);
    exit(1);
  }

  int i;
  int j;
  A.data=malloc(sizeof(double)*A.m*A.n);
  for (i=0;i<A.m;i++){
    if (fscanf(fp,"%lf",A.data+A.n*i)!=1){
      fprintf(stderr,"Error reading %s.\n",filename);
      exit(1);
    }
    for (j=1;j<A.n;j++){
      if (fscanf(fp,",%lf",A.data+A.n*i+j)!=1){
	fprintf(stderr,"Error reading %s.\n",filename);
	exit(1);
      }
    }
  }

  setup_entry_array(&A);
  
  fclose(fp);
  return A;


}


//Initialize a new matrix of size m x n.  The data is will be allocated for you, and the entry array all set up.
//But, the data in there is pure garbage. 
matrix new_matrix(int m, int n){

  matrix A;
  A.m=m;
  A.n=n;
  A.data=calloc(m*n,sizeof(double));

  setup_entry_array(&A);
  
  return A;
}

//Free a matrix after it is allocated.
void free_matrix(matrix A){
  free(A.data);
  free(A.entry);
}

//Print a matrix to file.  
void fprintmatrix(FILE* fp,matrix A,char delim){

  int i;
  int j;
  for (i=0;i<A.m;i++){
    for (j=0;j<A.n;j++){
      fprintf(fp,"%10lf",A.data[i*A.n+j]);
      if (j<A.n-1){
	fprintf(fp,"%c",delim);
      }
    }
    fprintf(fp,"\n");
  }
    
}


//Write the matrix A to the filename.
void write_matrix(char* filename, matrix A){

  FILE* fp=fopen(filename,"w");
  if (fp==NULL){
    fprintf(stderr,"Error opening %s for writing.\n",filename);
    exit(1);
  }

  fprintf(fp,"%d,%d\n",A.m,A.n);
  fprintmatrix(fp,A,',');
  
  fclose(fp); 


}

//Print a matrix to the command prompt. 
void printmatrix(matrix A){
  fprintmatrix(stdout,A,' ');
}

//Create a random matrix of size m x n with entries between min and max. 
matrix random_matrix(int m,int n,int min, int max){

  matrix A=new_matrix(m,n);
  int i;
  int j;

  srand(time(NULL));
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      A.data[i*n+j]=rand()%(max-min)+min;
    }
  }

  
  return A;
  
}

//returns a copy of the matrix A
matrix copy_matrix(matrix A){
  matrix B=new_matrix(A.m,A.n);
  memcpy(B.data,A.data,A.m*A.n*sizeof(double));
  return B;
}


//Calculates AB and returns it. 
matrix matrix_mult(matrix A, matrix B){
  if (A.n != B.m){
    fprintf(stderr,"Error, A is %d x %d while B is %d x %d.  Matrix dimensions must be consistent.\n",A.m,A.n,B.m,B.n);
    exit(1);
  }

  matrix C=new_matrix(A.m,B.n);
  int i,j,k;

  for (i=0;i<C.m;i++){
    for (k=0;k<A.n;k++){
      for (j=0;j<C.n;j++){
	C.entry[i][j]+=A.entry[i][k]*B.entry[k][j];
      }
    }
  }

  return C;
}

//Computes Ax-b.  If this is nearly the zero vector then x is close to the solution to Ax=b. 
matrix residual_matrix(matrix A, matrix b, matrix x){

  matrix R=matrix_mult(A,x);
  int i;
  for (i=0;i<R.n;i++){
    R.entry[i][0]-=b.entry[i][0];
  }
  return R;
}

//Computes ||Ax-b||  If this is close to 0 then x is very close to the solution to Ax=b (hopefully). 
double residual(matrix A, matrix b, matrix  x){
  matrix R=residual_matrix(A,b,x);
  int i;
  double r=0;
  for (i=0;i<R.n;i++){
    r+=R.entry[i][0]*R.entry[i][0];
  }
  r=pow(r,0.5);
  free_matrix(R);
  return r;
}
#endif


  
