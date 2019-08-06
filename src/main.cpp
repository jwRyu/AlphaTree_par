// #include <cstdio>
#include <ctime>
//#include <opencv2/opencv.hpp>
//#include <filesystem>
#include <iostream>
//#include <ctime>
//#include <chrono>
//#include <fstream>
#include <stdio.h>
#include <string.h>
#include <fitsio.h>
#include "defines.h"
#include "AlphaTree.h"
//#include <opencv.h>

using namespace std;
//using namespace cv;

#define DEBUG 0

#define OUTPUT_FNAME "./AlphaTree.dat"
#define OUTIMG_FNAME "./outimg.jpg"
#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/OrientedAlphaTree_test"

#define REPEAT 10	//program repetition for accurate runtime measuring

void RandomizedHDRimage(uint64* hdrimg, uint8* ldrimg, int64 imgsize)
{
	uint64 pix;

	for (int64 i = 0; i < imgsize; i++)
	{
		pix = ((uint64)ldrimg[i]) << 56;
		//pix = ((uint64)(rand() & 0xff) << 56); //tmp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		pix |= ((uint64)(rand() & 0xff) << 48);
		pix |= ((uint64)(rand() & 0xff) << 40);
		pix |= ((uint64)(rand() & 0xff) << 32);
		pix |= ((uint64)(rand() & 0xff) << 24);
		pix |= ((uint64)(rand() & 0xff) << 16);
		pix |= ((uint64)(rand() & 0xff) << 8);
		pix |= ((uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
}

int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, bnaxis, check = 1;
    long npixels = 1, firstpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
    int64 imgsize;
    double *img, *bpix, value;
//    uint8 *img8;
//    uint16 *img16;
//    uint32 *img32;
//    uint64 *img64;
    int image2=1;

    if (argc != 5) {
      printf("Usage: imarith image1 { image2 | value } oper outimage \n");
      printf("\n");
      printf("Perform 'image1 oper image2' or 'image1 oper value'\n");
      printf("creating a new output image.  Supported arithmetic\n");
      printf("operators are add, sub, mul, div (first character required\n");
      printf("\n");
      printf("Examples: \n");
      printf("  imarith in1.fits in2.fits a out.fits - add the 2 files\n");
      printf("  imarith in1.fits 1000.0 mul out.fits - mult in1 by 1000\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }
    fits_open_file(&bfptr, argv[2], READONLY, &status);
    if (status) {
      value = atof(argv[2]);
      if (value == 0.0) {
	printf("Error: second argument is neither an image name"
	       " nor a valid numerical value.\n");
	return(status);
      }
      image2 = 0;
      status = 0;
    }

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    if (image2) fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_size(afptr, 3, anaxes, &status);
    if (image2) fits_get_img_size(bfptr, 3, bnaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }
         /* check that the input 2 images have the same size */
    else if ( image2 && ( anaxes[0] != bnaxes[0] ||
			  anaxes[1] != bnaxes[1] ||
			  anaxes[2] != bnaxes[2] ) ) {
       printf("Error: input images don't have same size\n");
       check = 0;
    }
/*
    if      (*argv[3] == 'a' || *argv[3] == 'A')
      op = 1;
    else if (*argv[3] == 's' || *argv[3] == 'S')
      op = 2;
    else if (*argv[3] == 'm' || *argv[3] == 'M')
      op = 3;
    else if (*argv[3] == 'd' || *argv[3] == 'D')
      op = 4;
    else {
      printf("Error: unknown arithmetic operator\n");
      check = 0;
    }
*/
    /* create the new empty output file if the above checks are OK */
    if (check && !fits_create_file(&outfptr, argv[4], &status) )
    {
      /* copy all the header keywords from first image to new output file */
      fits_copy_header(afptr, outfptr, &status);

      npixels = anaxes[0];  /* no. of pixels to read in each row */
      if(anaxis == 1)
        imgsize = anaxes[0];  /* no. of pixels to read in each row */
      else if(anaxis == 2)
        imgsize = anaxes[0] * anaxes[1];
      else
        imgsize = anaxes[0] * anaxes[1] * anaxes[2];

      img = (double *) malloc(imgsize * sizeof(double)); /* mem for 1 row */
      if (image2) bpix = (double *) malloc(npixels * sizeof(double));

      if (img == NULL || (image2 && bpix == NULL)) {
        printf("Memory allocation error\n");
        return(1);
      }


      cout << firstpix[0] << " " << firstpix[1] << " " << firstpix[2] << endl;

      //fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, img,
        //                NULL, &status);

      fits_read_pix(afptr, TDOUBLE, firstpix, imgsize, NULL, img,
                          NULL, &status);
//Alpha-Tree Filtering

//      int datamax;
//      afptr->readKey("DATAMAX", &datamax);
//      cout << "datamax: " << afptr << endl;

      AlphaTreeWrapper *tree;

/*
      img8 = (uint8*)malloc(imgsize * sizeof(uint8));
      img16 = (uint16*)malloc(imgsize * sizeof(uint16));
      img32 = (uint32*)malloc(imgsize * sizeof(uint32));
      img64 = (uint64*)malloc(imgsize * sizeof(uint64));
      for(int64 ii = 0;ii < imgsize; ii++)
      {
        img8[ii] =  (uint8)img[ii];
        img16[ii] = (uint16)img[ii];
        img32[ii] = (uint32)img[ii];
        img64[ii] = (uint64)img[ii];
      }
*/

      cout << anaxes[1] << " x " << anaxes[0] << " x " << anaxes[2] << endl;

      tree = (AlphaTreeWrapper*)malloc(sizeof(AlphaTreeWrapper));
      tree->BuildAlphaTree(img, (int)anaxes[1], (int)anaxes[0],(int)anaxes[2], (int)4, (int)0);
      //void BuildAlphaTree(Pixel *img, Imgidx height_in, Imgidx width_in, Imgidx channel_in, Imgidx connectivity_in, int8 listsz_idx)
      //tree->clear();
      free(tree);

      fits_write_pix(outfptr, TDOUBLE, firstpix, imgsize,
                       img, &status); /* write new values to output image */

      fits_close_file(outfptr, &status);
      free(img);
      if (image2) free(bpix);
//      free(img8);
//      free(img16);
//      free(img32);
//      free(img64);

    }

    fits_close_file(afptr, &status);
    if (image2) fits_close_file(bfptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
