#ifndef VECTORH
#define VECTORH

typedef struct vector{

  double x;
  double y;
  double z;
}vector;

vector zerovec(){
  vector z;
  z.x=0;
  z.y=0;
  z.z=0;
  return z;
}

vector plus( vector v1, vector v2){
  vector s;
  s.x=v1.x+v2.x;
  s.y=v1.y+v2.y;
  s.z=v1.z+v2.z;
  return s;
}

vector minus(vector v1, vector v2){
  vector s;
  s.x=v1.x-v2.x;
  s.y=v1.y-v2.y;
  s.z=v1.z-v2.z;
  return s;
}

vector mult(vector v1, double alpha){
  vector s;
  s.x=alpha*v1.x;
  s.y=alpha*v1.y;
  s.z=alpha*v1.z;
  return s;
}

double vsnorm(vector v){
  return v.x*v.x+v.y*v.y+v.z*v.z;

}

double vnorm(vector v){
  return sqrt(vsnorm(v));
}

void fprintvec(FILE* fp, char* fmt, vector v){
  fprintf(fp,fmt,v.x,v.y,v.z);
}

void printvec(vector v){
  fprintvec(stdout,"(%lf,%lf,%lf)",v);
}

#endif
