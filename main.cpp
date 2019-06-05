#include <cstdio>
#include <ctime>
//#include <opencv2/opencv.hpp>
//#include <filesystem>
#include <iostream>
//#include <ctime>
//#include <chrono>
//#include <fstream>
//#include "defines.h"
//#include "allocator.h"
//#include "HQueue.hpp"
//#include "AlphaTree.hpp"
#include <fitsio.h>

using namespace std;
//using namespace cv;

#define DEBUG 0

#define OUTPUT_FNAME "./AlphaTree.dat"
#define OUTIMG_FNAME "./outimg.jpg"
#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/OrientedAlphaTree_test"

#define REPEAT 10	//program repetition for accurate runtime measuring

int main(int argc, char **argv)
{
//fits_init_cfitsio();
	return 0;
}

/*
int main(int argc, char **argv)
{
	AlphaTree tree;
//	Pixel *outimg;
	uint32 width, height, channel;
	uint32 cnt = 0;
	ofstream f;
	ifstream fcheck;
	uint32 contidx;
	char in;

	contidx = 0;
	//	f.open("C:/Users/jwryu/RUG/2018/AlphaTree/AlphaTree_grey_Exp.dat", std::ofstream::app);
	fcheck.open(OUTPUT_FNAME);
	if (fcheck.good())
	{
		cout << "Output file \"" << OUTPUT_FNAME << "\" already exists. Overwrite? (y/n/a)";
		//cin >> in;
		in = 'y';
		if (in == 'a')
		{
			f.open(OUTPUT_FNAME, std::ofstream::app);
			cout << "Start from : ";
			cin >> contidx;
		}
		else if (in == 'y')
			f.open(OUTPUT_FNAME);
		else
			exit(-1);
	}
	else
		f.open(OUTPUT_FNAME);
	//f.open("C:/Users/jwryu/RUG/2018/AlphaTree/AlphaTree_RGB_Exp.dat", std::ofstream::app);

	cnt = 0;
	std::string path = INPUTIMAGE_DIR;
	for (auto & p : std::experimental::filesystem::directory_iterator(path))
	{
		if (++cnt < contidx)
		{
			continue;
		}
		cv::String str1(p.path().string().c_str());
		cv::Mat cvimg = imread(str1, cv::IMREAD_ANYCOLOR);
		cv::Mat outimg(cvimg);

		height = cvimg.rows;
		width = cvimg.cols;
		channel = cvimg.channels();

		cout << cnt << ": " << str1 << ' ' << height << 'x' << width << endl;

		if (channel != 1)
		{
			cout << "input should be a 1-ch image" << endl;
			getc(stdin);
			exit(-1);
		}

		//img = (Pixel*)malloc(height * width * sizeof(Pixel));
		//imreshape((uint8*)img, cvimg.data, height, width);

		double runtime, minruntime;
		for (int testrep = 0; testrep < REPEAT; testrep++)
		{
			//memuse = max_memuse = 0;
			auto wcts = std::chrono::system_clock::now();

			tree.BuildAlphaTree((Pixel*)cvimg.data, height, width, channel, 4);

			//outimg = new Pixel[width * height];

			//outimg = new Pixel[width * height];
			//tree.AlphaFilter(outimg, 10, 30);
			//cv::String outfname(OUTIMG_FNAME);
			//imwrite(outfname,

			std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);

			//tree.AlphaFilter((Pixel*)outimg.data, 30, 100);

//			namedWindow("image", WINDOW_AUTOSIZE);
	//		imshow("image", outimg);
		//	waitKey(0);

			runtime = wctduration.count();
			minruntime = testrep == 0 ? runtime : min(runtime, minruntime);

			if (testrep < (REPEAT - 1))
				tree.clear();
		}
		//f << p.path().string().c_str() << '\t' << height << '\t' << width << '\t' << max_memuse << '\t' << nrmsd << '\t' << tree->maxSize << '\t' << tree->curSize << '\t' << minruntime << endl;

		cout << "Time Elapsed: " << minruntime << endl;

		cin >> in;

		//free(outimg);
		tree.clear();
		cvimg.release();
		str1.clear();
	}

	cin >> in;

	return 0;
}
*/
