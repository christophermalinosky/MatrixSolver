#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include "solve_system_helpers.h"
#include "nbodyutils.h"
#include<mpi.h>

MPI_Datatype MPI_Vector;

MPI_Datatype create_MPI_vector(){
  vector v;
  MPI_Datatype dtype[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  int blocklen[3]={1,1,1};
  MPI_Aint disp[3]={(long int)&(v.x)-(long int)&v,(long int)&(v.y)-(long int)&v,(long int)&(v.z)-(long int)&v};
  MPI_Datatype vector_type;
  MPI_Type_create_struct(3, blocklen, disp, dtype, &vector_type);
  MPI_Type_commit(&vector_type);
  return vector_type;
}

//print a buffer of length n. 
void printbuffer(double* buffer, int n){

  int i;
  printf("[");
  for (i=0;i<n;i++){
    printf("%f", buffer[i]);
    if(i < n - 1){
    printf(" ");
    }
  }
  printf("]\n");

  
}

void main(int argc,char** argv){
  MPI_Init(&argc,&argv);
  MPI_Vector = create_MPI_vector();
  int rank;
  int nproc;
  int i;
  int j;
  MPI_Status status;
  double div;
  int done = 1;

  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int lastRow = -1;
    int k;
    if(rank == 0){
        matrix A=read_matrix(argv[1]);
        matrix b=read_matrix(argv[2]);
        matrix x=new_matrix(b.m,1);

        if (argc<3){
            printf("Usage: solve_system A b x\n");
            printf("Here A is the file storing the matrix A\n");
            printf("b is the file storing the vector b");
            printf("and x is the output file");
            exit(1);
        } else if(A.m < nproc){
            printf("The amount of processors is more than the number of rows. Please enter a lower amount.\n");
            exit(1);
        } else if(A.m != b.m){
            printf("The vector must have an equal number of rows as the matrix\n");
            exit(1);
        }
        
        MPI_Bcast(&A.n,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&A.m,1,MPI_INT,0,MPI_COMM_WORLD);

        int numCol = A.n;
        int numRow = A.m;

        int index = 0;
        int totalValues = (numRow/nproc) + 1;
        double assignment[totalValues][numCol];
        double assignA[totalValues][1];
        int j = 0;
        while(j < A.m){
            if(j%nproc != 0){
                MPI_Send(A.entry[j], A.n, MPI_DOUBLE, j%nproc, 1, MPI_COMM_WORLD);
                MPI_Send(b.entry[j],1,MPI_DOUBLE, j%nproc, 1, MPI_COMM_WORLD);
            } else {
                for(k=0;k<numCol;k++){
                    assignment[index][k] = A.entry[j][k];
                }
                assignA[index][0] = b.entry[j][0];
                index++;
            }              
            j++;
        }
        for(j = 1; j < nproc; j++){
            MPI_Send(A.entry[0], A.n, MPI_DOUBLE, j, 2, MPI_COMM_WORLD);
        }
        
        int start = 1;
        double* recvbuf = malloc(numCol*sizeof(double));
        recvbuf = A.entry[0];
        double temp = assignA[0][0];
        
        index = 0;
        int assignIndex = 0;        
        while(done == 1){
            if(start == 0){
                MPI_Recv(recvbuf, numCol, MPI_DOUBLE, nproc-1, 0, MPI_COMM_WORLD,&status);
                MPI_Recv(&lastRow, 1, MPI_INT, nproc-1, 0, MPI_COMM_WORLD,&status);
                MPI_Recv(&temp,1,MPI_DOUBLE, nproc-1, 0, MPI_COMM_WORLD,&status);
            }

            if(lastRow < numCol - 2){
                MPI_Send(recvbuf, numCol,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                int sendRow = lastRow + 1;
                MPI_Send(&sendRow,1,MPI_INT,(rank+1)%nproc,0,MPI_COMM_WORLD);
                MPI_Send(&temp,1,MPI_DOUBLE, (rank+1)%nproc, 0, MPI_COMM_WORLD);
            }

            if(start == 0){
                int i;
                div = assignment[assignIndex][index]/recvbuf[index];
                for(i = 0; i < numCol; i++){
                    assignment[assignIndex][i] = assignment[assignIndex][i] - recvbuf[i] * (div);
                }
                assignA[assignIndex][0] = assignA[assignIndex][0] - temp * div;

                index++;
                if(index == (assignIndex * nproc) + rank){
                    if(lastRow < numCol - 2){
                        MPI_Send(assignment[assignIndex], numCol,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                        int sendRow = lastRow + 1;
                        MPI_Send(&sendRow,1,MPI_INT,(rank+1)%nproc,0,MPI_COMM_WORLD);
                        MPI_Send(assignA[assignIndex],1,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                    }
                    div = assignment[assignIndex][index];
                    for(i = index; i < numCol; i++){
                        assignment[assignIndex][i] = assignment[assignIndex][i]/div;
                    }
                    assignA[assignIndex][0] = assignA[assignIndex][0]/div;
                    index = 0;
                    assignIndex++;
                }
            } else {
                div = assignment[0][0];
                for(i = 0; i < numCol; i++){
                        assignment[0][i] = assignment[0][i]/div;
                }
                assignA[0][0] = assignA[0][0]/div;
                index = 0;
                assignIndex = 1;
                start = 0;
            }
            if((assignIndex * nproc) + rank >= numRow){
                assignIndex--;
                index = numRow - 1;
                if((assignIndex * nproc) + rank == (numRow - 1)){
                    MPI_Send(&assignA[assignIndex],1,MPI_DOUBLE,nproc-1,0,MPI_COMM_WORLD);
                    assignIndex--;
                }
                while(1){
                    if(assignIndex < 0){
                        write_matrix(argv[3],x);
                        done = 0;
                        break;
                    }
                    MPI_Recv(&temp,1,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD,&status);
                    
                    x.entry[index][0] = temp;

                    if(assignIndex != 0){
                        MPI_Send(&temp,1,MPI_DOUBLE,nproc-1,0,MPI_COMM_WORLD);
                    }

                    assignA[assignIndex][0] = assignA[assignIndex][0] - temp * assignment[assignIndex][index];
                    index--;
                    if(index == (assignIndex * nproc) + rank){
                        if(assignIndex != 0){
                            MPI_Send(assignA[assignIndex],1,MPI_DOUBLE,nproc-1,0,MPI_COMM_WORLD);
                        }
                        x.entry[index][0] = assignA[assignIndex][0];
                        index = numRow - 1;
                        assignIndex--;
                    }
                }
            }
        }
    } else {
        int numCol;
        int numRow;

            MPI_Bcast(&numCol,1,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast(&numRow,1,MPI_INT,0,MPI_COMM_WORLD);

            int totalValues = (numRow/nproc) + 1;
            double assignment[totalValues][numCol];
            double assignA[totalValues][1];

            int index = 0;
            double* recvbuf = malloc(numCol*sizeof(double));
            while(1) {
                MPI_Recv(recvbuf, numCol, MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
               
                if(status.MPI_TAG == 2) {
                    break;
                }
                MPI_Recv(assignA[index], 1, MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);
                for(j=0;j<numCol;j++){
                    assignment[index][j] = recvbuf[j];
                }
                index++;
            }

            index = 0;
            int assignIndex = 0;
            double temp;

            while(done == 1){
                MPI_Recv(recvbuf, numCol, MPI_DOUBLE, (rank-1)%nproc, 0, MPI_COMM_WORLD,&status);
                MPI_Recv(&lastRow, 1, MPI_INT, (rank-1)%nproc, 0, MPI_COMM_WORLD,&status);
                MPI_Recv(&temp, 1, MPI_DOUBLE, (rank-1)%nproc, 0, MPI_COMM_WORLD,&status);
                
                if(lastRow < numCol - 2){
                    MPI_Send(recvbuf, numCol,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                    int sendRow = lastRow + 1;
                    MPI_Send(&sendRow,1,MPI_INT,(rank+1)%nproc,0,MPI_COMM_WORLD);
                    MPI_Send(&temp,1,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                }

                int i;
                div = assignment[assignIndex][index]/recvbuf[index];
                for(i = 0; i < numCol; i++){
                    assignment[assignIndex][i] = assignment[assignIndex][i] - recvbuf[i] * (div);
                }
                assignA[assignIndex][0] = assignA[assignIndex][0] - temp * div;
                index++;
                if(index == (assignIndex * nproc) + rank){
                    
                    if(lastRow < numCol - 2){
                        MPI_Send(assignment[assignIndex], numCol,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                        int sendRow = lastRow + 1;
                        MPI_Send(&sendRow,1,MPI_INT,(rank+1)%nproc,0,MPI_COMM_WORLD);
                        MPI_Send(assignA[assignIndex],1,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD);
                    }

                    div = assignment[assignIndex][index];
                    for(i = index; i < numCol; i++){
                        assignment[assignIndex][i] = assignment[assignIndex][i]/div;
                    }
                    assignA[assignIndex][0] = assignA[assignIndex][0]/div;
                    index = 0;
                    assignIndex++;
                    int lol = (assignIndex * nproc) + rank;
                    
                }
                if((assignIndex * nproc) + rank >= numRow){                   
                    assignIndex--;
                    index = numRow - 1;
                    if((assignIndex * nproc) + rank == (numRow - 1)){
                        MPI_Send(&assignA[assignIndex],1,MPI_DOUBLE,(rank-1)%nproc,0,MPI_COMM_WORLD);
                        assignIndex--;
                    }
                    while(1){
                        if(assignIndex < 0){
                            done = 0;
                            break;
                        }
                        MPI_Recv(&temp,1,MPI_DOUBLE,(rank+1)%nproc,0,MPI_COMM_WORLD,&status);\
                        MPI_Send(&temp,1,MPI_DOUBLE,(rank-1)%nproc,0,MPI_COMM_WORLD);

                        assignA[assignIndex][0] = assignA[assignIndex][0] - temp * assignment[assignIndex][index];
                        index--;
                        if(index == (assignIndex * nproc) + rank){
                            MPI_Send(assignA[assignIndex],1,MPI_DOUBLE,(rank-1)%nproc,0,MPI_COMM_WORLD);
                            index = numRow - 1;
                            assignIndex--;
                        }
                    }
                }
            }
        }
    MPI_Finalize();
}
