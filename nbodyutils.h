#ifndef __NBODYUTILS__
#define __NBODYUTILS__

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include "vector.h"
#include<mpi.h>



typedef struct{
  int numsteps;     //number of steps to take
  double endtime;   //ending time of simulation
  double stepsize;  //stepsize
  double theta;
  int numprocs;     //number of processors for this simulation
  int rank;         //my rank in the world of processors
  int (*body_assmts)[2];  //which processors are assigned to which bodies?
  int *global_displs;
  int  *body_count;  //how many bodies are assigned to each processor?
  int method;                //ODE solver method
}sim_opts;


typedef struct{
  int N;    //how many bodies there are
  double G; //Gravitational constant for units
  vector* X0; //vector array of initial positions
  vector* V0; //vector array of initial velocities
  vector* X;  //vector array [n*times] of all positions at all times
  double* M;  //array of masses
  double* times; //array of times in simulation
  int numsteps; //number of steps taken
  int* assigned;
}nbody_dataset;

void print_usage(){
  printf("Example usage:\n");
  printf("nbody input_file output_file <--nsteps=x> <--endtime=x> <--stepsize=x> <options>\n\n");
  printf("two of nsteps, endtime, and stepsize must be supplied.\n\n");
  printf("endtime        the ending time of the simulation, taking the initial state as time 0.\n");
  printf("stepsize       the length of time between steps\n");
  printf("nsteps         the number of timesteps to take\n");
  printf("\n");
  printf("Time units may be specified in earth (d)ays, (y)ears, or (h)ours by appending d,y,or h to the time unit.  If no time unit is supplied, the default is seconds.\n");
  printf("We assume that the input velocity is in something/second\n\n");
  printf("For example endtime=1y means 1 earth year, endtime=1d means 1 earth day, endtime=1h means 1 earth hour\n\n\n");
  printf("Time units may be applied to stepsize as well.\n");
  printf("Options are:\n");
  printf("--method=<RK4|Euler|ImEuler>  Default is Euler.\n");
   
}


void print_options(sim_opts o){
  
  printf("Number of timesteps: %d\n",o.numsteps);
  printf("Endtime: %lf (s)\n",o.endtime);
  printf("Stepsize: %lf (s)\n",o.stepsize);
  printf("Theta: %lf\n",o.theta);
  printf("Solution method: ");
  switch(o.method){
  case 1:
    printf("RK4");
    break;
  case 2:
    printf("Euler");
    break;
  case 3:
    printf("Improved Euler");
    break;
  }
  printf("\n");
  
}



void read_sim_opts(int argc,char** argv, sim_opts* s){
  int end_time_specified=0;
  double T=0;
  int n_steps_specified=0;
  int stepsize_specified=0;

  float dummy;
    double dummylf;
  s->method=2;
  #ifndef serial
  MPI_Comm_size(MPI_COMM_WORLD,&(s->numprocs));
  MPI_Comm_rank(MPI_COMM_WORLD,&(s->rank));
  #endif 
  
  
  int i;

  for (i=3;i<argc;i++){
    if (sscanf(argv[i],"--endtime=%lf",&T)==1){
      switch (argv[i][strlen(argv[i])-1]){
      case 'y':
	    T*=3.15569e7;
    	break;
      case 'd':
	    T*=86400;
    	break;
      case 'h':
	    T*=3600;
    	break;
      }
      s->endtime=T;
      end_time_specified=1;
    }
    else if (sscanf(argv[i],"--nsteps=%f",&(dummy))==1){
      s->numsteps=(int)dummy;
      n_steps_specified=1;
    }
    else if (sscanf(argv[i],"--theta=%lf",&(dummylf))==1){
      s->theta=dummylf;
    }
    else if (sscanf(argv[i],"--stepsize=%lf",&(T))==1){
      stepsize_specified=1;
      switch(argv[i][strlen(argv[i])-1]){
      case 'y':
    	T*=3.15569e7;
	    break;
      case 'd':
	    T*=86400;
	    break;
      case 'h':
	    T*=3600;
	    break;
      }
      s->stepsize=T;
    }
    else if(strcmp(argv[i],"--method=RK4")==0){
      s->method=1;
    }
    else if (strcmp(argv[i],"--method=Euler")==0){
      s->method=2;
    }
    else if (strcmp(argv[i],"--method=ImEuler")==0){
      s->method=3;
    }
    else{
      print_usage();
      fprintf(stderr,"Error parsing option %s\n",argv[i]);
      exit(1);
    }
  }

  if (end_time_specified && n_steps_specified && stepsize_specified){
    fprintf(stderr,"You can't specify stepsize, number of steps, and end time.\n");
    print_usage();
    exit(1);
  }
  else if (end_time_specified && n_steps_specified){
    s->stepsize=s->endtime/s->numsteps;
  }
  else if (end_time_specified && stepsize_specified){
    s->numsteps=s->endtime/s->stepsize;
  }
  else if (stepsize_specified && n_steps_specified){
    s->endtime=s->stepsize*s->numsteps;
  }
  else{
    fprintf(stderr,"You must specify two of end time, number of steps, and step size.\n");
    print_usage();
    exit(1);
  }


}

double dist(vector v, vector w){
  return sqrt(vsnorm(minus(w,v)));

}

void free_nbody_dataset(nbody_dataset d){
  free(d.X0);
  free(d.V0);
  free(d.X);
  free(d.M);
  free(d.times);
}



//Load data from file into a dataset. 
void load_data(char* filename, nbody_dataset* d){
  FILE* fp=fopen(filename,"r");
  if (fp==NULL){
    fprintf(stderr,"Error opening %s\n",filename);
    exit(1);
  }
  if (fscanf(fp,"%d ",&(d->N))!=1){
    fprintf(stderr,"Error reading first line of %s to get number of bodies. Sorry!\n",filename);
    exit(1);
  }
  d->X0=malloc(sizeof(vector)*(d->N));
  d->V0=malloc(sizeof(vector)*(d->N));
  d->M=malloc(sizeof(double)*(d->N));
  d->assigned=malloc(sizeof(int)*(d->N));
  d->G=6.674e-11;
  vector* x0=d->X0;
  vector* v0=d->V0;
  double* M=d->M;
  int numread;
  int i=0;
  char line[10000];
  while (fgets(line,10000,fp)!=NULL){
    if (line[0]!='#'){
      numread=sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&(x0[i].x),&(x0[i].y),&(x0[i].z),&(v0[i].x),&(v0[i].y),&(v0[i].z),&(M[i]));
      if (numread!=7){
	fprintf(stderr,"Error reading line %d,%s,of %s.\n",i,line,filename);
	exit(1);
      }
      i++;
    }
  }

  if (i!=d->N){
    fprintf(stderr,"Error.  %s said there were %d bodies, but read %d.\n",filename,d->N,i);
    exit(1);
  }
}


//Write data from dataset to file. 
void write_data(char* filename,nbody_dataset* d ){

  FILE* fp=fopen(filename,"w+");
  if (fp==NULL){
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(1);
  }
  fprintf(fp,"v1.0\n");
  fprintf(fp,"nbody=%d\n",d->N);
  vector (*x)[d->N] = (vector (*)[d->N])d->X;
 
  int ti; int i;
  for (ti=0;ti<=d->numsteps;ti++){
    for (i=0;i<d->N;i++){
      fprintf(fp,"%lf,%lf,%lf\n",x[ti][i].x,x[ti][i].y,x[ti][i].z);
    }
  }

  fclose(fp);
}



#endif
