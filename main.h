#include "stdint.h"
#define __STDC_FORMAT_MACROS
#include "inttypes.h"

#ifndef PRId64
#define PRId64 "ld"
#endif
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif


size_t getFilesize(const char* filename) {
    struct stat st;
    stat(filename, &st);
    return st.st_size;
}

inline bool fexists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int n;
bigint ntimestep,natoms;
int size_one,nchunk,triclinic;
double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;
int boundary[3][2];
//char boundstr[9];

char *lineD = NULL;
char *lineS = NULL;

int VB_tic;

int *vel_Jout = NULL;
int *vel_Jcoll = NULL;
int *vel_Jevap = NULL;
int *vel_Jcond = NULL;

double *indexID = NULL;
int *indices = NULL;
int *JevapL_oldlist = NULL;
int *JevapL_newlist = NULL;
int *JcondL_oldlist = NULL;
int *JcondL_newlist = NULL;
int *JcollL_oldlist = NULL;
int *JcollL_newlist = NULL;
int *JoutL_oldlist = NULL;
int *JoutL_newlist = NULL;

int data,vflag;
//int *Jold = NULL;
//int *Jnew = NULL;

int *JL = NULL;

FILE *screen = stdout;
FILE *error = stderr;
FILE *fpdatvel = NULL;
FILE *fpposition = NULL;

int NATOMS,SBYTES;

int imass,imol,iid,itype,ix,iy,iz,ivx,ivy,ivz,ike,ipe,is1,is2,is3,is4,is5,is6;
int STRflag,output_every_flag,nflag;
int start,stop,N,output,step_size,threads,iV,output_every;
double LB_start,LB_stop,iLB,VB_start,VB_stop,iVB,vmin,vmax;
//double *coord = NULL;

char *STR = NULL;

int maxbuf = 0;
double *buf = NULL;

const double Na = 6.02214086e23; // 1/mol
const double kb = 1.38064852e-23; // J/K
const double R = Na * kb;
const double Pi = 3.141592653589793;

