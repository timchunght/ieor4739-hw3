#ifndef POWER

#define POWER

typedef struct powerbag{
  int n;
  double *vector;
  double *newvector;
  double *matrix;
  double topeigvalue;
  int ID;
}powerbag;

void PWRpoweralg_new(powerbag *pbag);
void PWRfreespace(powerbag **ppbag);
void PWRfree(double **pvector);
int PWRreadnload_new(char *filename, int ID, powerbag **ppbag);

int PWRallocatespace(int n, double **pvector, double **pnewvector, double **pmatrix);
int PWRreadnload(char *filename, int *pn, double **pvector, double **pnewvector, double **pmatrix);

void PWRpoweriteration(int ID, int k, 
		    int n, double *vector, double *newvector, double *matrix,
		       double *peigvalue, double *perror);

void PWRpoweralg(int n, double *vector, double *newvector, double *matrix,
		 double *peigvalue);

void PWRshowvector(int n, double *vector);

void PWRcompute_error(int n, double *perror, double *newvector, double *vector);
#endif

/*class powerunit{
 public:
  powerunit(char *inputname);
  ~powerunit();
  int readnload(char *file);
  int allocatespace();
  void iterate();
  void showvector(){ for(int j = 0; j < n; j++) printf("%g ",vector[j]);
  printf("\n");}
  char *getname(){return name;}
 private:
  void releasespace();
  int n;
  double *matrix;
  double *vector;
  double *newvector;
  char *name;
};
*/
