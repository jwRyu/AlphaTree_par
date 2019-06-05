#ifndef _COMMON_H
#define _COMMON_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include </usr/local/bin/cfitsio-3.47/fitsio.h>
#include <math.h>
//#include <FreeImage.h>

#define PI 3.14159265358979323846
#define false 0
#define true  1
#define MAXTHREADS 128
#define CONNECTIVITY  6

/*********************************************/
/*#define USEFLOATPOINT 1
#define USEFLOATPOINT_ONLYPOSITIVE 0
typedef float greyval_t; //intensities
typedef long pixel_t; //coordinates
*/
/********************************************/

/*********************************************/
#define USEFLOATPOINT 1
#define USEFLOATPOINT_ONLYPOSITIVE 0
//#define USEFLOATPOINT_SORT 1
//#define USEFLOATPOINT_ONLYPOSITIVE_SORT 1
typedef float greyval_t; //intensities
typedef long pixel_t; //coordinates

/********************************************/

#define NOTPROCESSED 0
#define PROCESSED 1

typedef unsigned char ubyte;
typedef short bool;
int MULFACTOR;

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
//#define LWB(self, nthreads) (width*(((self)*depth)/nthreads))
//#define UPB(self, nthreads) (width*(((self+1)*depth)/nthreads))

//#define LWB(self, nthreads) ((width*height)*(((self)*depth)/nthreads))
//#define UPB(self, nthreads) ((width*height)*(((self+1)*depth)/nthreads))

#define LWB(self, nthreads) (size2D*(((self)*depth)/nthreads)) // careful: it must hold that depth >= nthreads!
#define UPB(self, nthreads) (size2D*(((self+1)*depth)/nthreads))

#define bottom (-1)

//#define M(px) (mapPixelToSorted[px])

float floatpointMUL;
int nthreads;
pthread_t threadID[MAXTHREADS];
int nthreadsRef; // number of threads of the refinement, equals to the number of quantized grey levels numQTZLEVELS

long width, height, depth, size, size2D;  /* precondition: width <= size/nthreads */

greyval_t *gval;
greyval_t *outRef;
int *gval_qu;

int numQTZLEVELS;
int inputQTZLEVELS;
double lambda;
pixel_t minpos;

short bitsPerPixel;

typedef struct MaxNode
{
    pixel_t parent;
    //short gvalqu;
    long Area;
} MaxNode;

MaxNode *node_qu;
MaxNode *node_qus;
bool *reached_qu;

// Used for parallel counting sort
pixel_t *SORTED;
pixel_t *SORTEDRS[2];
unsigned int NUMBUCKETS;

// Used for the creation of the quantized image
pixel_t *pxStartPosition;
pixel_t *pxEndPosition;

// Used for the refinement with Najman Couprie
typedef struct Node
{
	//greyval_t Level;
	greyval_t filter;
    pixel_t parent;
    long Area;
} Node;

Node *node_ref;
pixel_t *zpar;

/***** Timings ****/
clock_t start;
struct tms tstruct;
float musec;
long tickspersec;

clock_t start_sort, end_sort;
clock_t start_quimg, end_quimg;
clock_t start_qutree, end_qutree;
clock_t start_ref, end_ref;

#endif
