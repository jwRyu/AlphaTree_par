#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "AlphaTree.h"
#include "HQueue.h"
#include "allocator.h"
#include "Trie.h"
using namespace std;

#define OUTPUT_FNAME "C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/tmp.dat"
//#define OUTPUT_FNAME "C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/tmptmptmp.dat"

#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Grey"
#define INPUTIMAGE_DIR_COLOUR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Colour" //colour images are used after rgb2grey conversion
#define REPEAT 8
#define RUN_TSE_ONLY 0

#define DEBUG 0

#if DEBUG
void* buf;
uint64 bufsize;
void save_buf(void* src, uint64 size)
{
	memcpy(buf, src, size);
	bufsize = size;
}
uint8 isChanged(void *src)
{
	uint64 i;
	for (i = 0; i < bufsize; i++)
	{
		if (((uint8*)buf)[i] != ((uint8*)src)[i])
			return 1;
	}
	return 0;
}
#endif

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

void RandomizedHDRimage(uint32* hdrimg, uint8* ldrimg, int64 imgsize)
{
	uint32 pix;

	for (int64 i = 0; i < imgsize; i++)
	{
		pix = ((uint64)ldrimg[i]) << 24;
		pix |= ((uint64)(rand() & 0xff) << 16);
		pix |= ((uint64)(rand() & 0xff) << 8);
		pix |= ((uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
}

/*
void DeleteAlphaTree(AlphaTree* tree)
{
	Free(tree->parentAry);
	Free(tree->node);
	Free(tree);
}*/


void sort_h1(double* h1_arr, int obs)
{
	int *hist, *h, *h1, i, j, bitfield, shamt, hsum;
	double *tmp, *p, *end;
	int64_t mask = 0xffff, num, offset1, offset2, offset3;

	offset1 = 65536;
	offset2 = 65536 << 1;
	offset3 = offset1 + offset2;

	hist = new int[65536 << 2];
	tmp = new double[obs];

	for (i = 0; i < 65536 << 2; i++)
		hist[i] = 0;

	end = h1_arr + obs;
	for (p = h1_arr; p < end; p++)
	{
		num = *((int64_t*)(p));

		hist[num & mask]++;
		hist[offset1 + ((num >> 16) & mask)]++;
		hist[offset2 + ((num >> 32) & mask)]++;
		hist[offset3 + ((num >> 48) & mask)]++;
	}

	shamt = 0;
	for (bitfield = 0; bitfield < 4; bitfield++)
	{
		hsum = 0;
		h = hist;
		h1 = h + 65536;
		while (h != h1)
		{
			hsum += *h;
			*(h++) = hsum;
		}	
		
		for (p = h1_arr + obs - 1; p >= h1_arr; p--)
		{
			num = *((int64_t*)p);
			j = --hist[(num >> shamt) & mask];
			tmp[j] = *p;
		}
		p = tmp; tmp = h1_arr; h1_arr = p;

		hist += 65536;
		shamt += 16;
	}
	hist -= 65536 << 2;

	delete[] hist;
	delete[] tmp;
}

int main(int argc, char **argv)
{
	//	AlphaTree<int32, uint8> *tree;
	AlphaTree<int32, uint64> *tree;
	int32 width, height, channel;
	int32 cnt = 0;
	ofstream f;
	ifstream fcheck;
	char in;
	int32 i, contidx;
	std::string path;
	double time_elapsed = 0;
	double pixels_processed = 0;
//	uint8 testimg[4 * 4] = { 4, 4, 2, 0, 4, 1, 1, 0, 0, 3, 0, 0, 2, 2, 0, 5 };
	uint8 testimg[4 * 4] = { 30, 37, 26, 22, 6, 4, 7, 30, 36, 26, 6, 24, 8, 35, 35, 10 };
	uint64 testimg_4[4 * 4] = { 30, 37, 26, 22, 6, 4, 7, 30, 36, 26, 6, 24, 8, 35, 35, 10 };
	uint64 testimg_2[2 * 2] = { 30, 37, 6, 4};
	uint64 testimg_64[30 * 30] = { 512372, 231570, 507352, 114529, 185524, 1029, 788748, 658566, 720989, 717518, 847815, 499829, 372402, 293441, 523032, 6885, 944424, 612075, 211972, 98664, 985184, 223966, 875069, 870690, 794572, 828932, 40616, 23714, 551095, 284978, 180980, 624162, 130240, 502427, 839004, 552390, 2419, 591446, 193930, 170408, 42749, 339308, 292055, 106215, 405547, 201640, 491947, 216391, 963379, 995142, 313421, 837912, 355320, 276629, 313683, 949994, 347141, 818517, 643771, 150945, 751058, 644008, 77438, 228652, 402256, 417640, 660188, 380027, 384192, 623806, 111673, 242377, 298844, 580963, 521363, 986232, 566030, 396037, 600030, 269684, 383579, 30574, 959072, 656253, 594689, 549265, 289945, 181701, 429003, 813174, 109687, 835014, 203804, 884826, 152083, 57397, 352371, 429217, 12456, 721680, 110819, 4208, 379917, 797415, 911521, 848010, 969867, 981779, 771739, 344064, 439920, 614758, 239988, 652002, 980276, 98068, 70337, 719684, 840124, 797871, 142600, 951506, 303380, 457111, 726044, 856057, 956805, 759245, 677583, 518320, 200533, 737326, 936941, 240212, 99109, 391675, 903776, 944183, 433441, 697801, 576217, 337286, 159255, 297869, 117066, 886324, 237202, 560821, 895522, 79939, 457421, 642679, 542074, 3604, 679924, 286368, 896661, 455284, 253800, 83340, 342031, 73046, 679402, 572373, 849384, 240757, 6171, 28047, 260317, 303353, 80343, 506673, 236394, 265251, 721115, 609949, 865870, 940533, 932071, 429963, 423205, 735073, 717459, 480795, 81597, 394811, 552894, 920556, 105266, 667005, 131107, 566078, 753721, 386717, 35119, 451125, 620795, 21653, 149019, 768958, 260535, 274886, 876468, 62227, 752622, 564005, 164778, 138698, 429974, 450127, 457621, 504981, 992288, 351882, 227532, 212444, 976079, 193617, 257823, 191434, 789981, 925565, 591306, 872606, 644159, 711525, 957635, 624795, 466053, 371883, 632277, 519206, 463405, 998605, 955666, 376048, 200511, 156420, 543584, 239046, 816402, 44140, 828118, 231225, 102844, 453861, 467956, 497364, 459681, 650319, 42111, 939082, 959199, 117586, 603530, 445072, 26834, 919493, 517857, 700304, 942987, 223117, 63517, 999855, 297006, 201399, 145341, 82582, 763359, 437579, 978885, 366189, 688665, 18319, 540951, 345744, 791992, 608677, 924691, 879300, 271002, 935990, 79742, 46066, 114383, 961335, 13077, 654115, 515510, 635536, 976374, 928251, 22129, 24473, 357723, 169491, 649894, 950009, 173273, 440208, 141083, 372417, 321450, 933184, 383836, 378762, 154988, 470233, 128949, 399304, 275216, 174115, 245794, 4183, 324145, 174509, 507646, 958183, 394333, 133578, 786511, 791629, 603448, 397289, 571904, 656913, 25519, 798155, 631138, 977056, 374975, 872631, 550722, 319430, 551543, 768007, 495988, 907524, 710120, 562530, 141129, 904231, 873523, 509902, 410430, 89477, 515071, 898753, 577574, 423812, 947961, 521120, 321264, 978141, 219495, 343890, 644692, 159594, 109577, 350564, 487773, 469794, 780675, 790148, 20872, 584694, 630937, 840802, 83437, 211799, 447690, 351300, 444478, 899715, 397982, 642861, 983318, 857453, 918034, 494949, 96559, 50368, 991216, 276280, 62973, 518685, 758472, 720370, 933273, 207145, 917415, 976433, 713534, 899121, 328567, 141484, 329444, 958762, 275800, 942434, 287416, 865192, 398220, 168026, 544938, 755772, 852202, 574631, 553030, 197874, 922874, 168656, 562797, 285194, 578283, 212858, 745139, 68497, 32818, 447800, 38022, 797817, 569395, 55018, 892841, 299712, 40334, 419080, 912035, 704005, 199409, 326426, 151546, 400600, 298367, 540384, 290738, 494808, 5009, 897743, 864798, 333117, 231654, 431441, 285709, 875256, 86638, 529725, 489292, 200297, 241610, 759031, 642165, 474890, 391833, 343723, 553915, 853089, 40728, 177784, 205019, 407155, 501672, 276024, 252898, 874633, 741321, 269500, 272004, 614931, 987822, 841396, 175153, 604012, 44771, 147014, 578649, 977794, 378193, 80427, 70921, 894669, 494510, 145243, 183324, 355631, 340861, 242659, 984289, 658682, 602408, 571522, 138276, 67983, 982462, 11477, 673492, 287280, 734347, 143178, 774901, 337271, 627491, 848732, 494567, 566939, 458613, 842166, 688456, 435834, 940562, 804601, 214, 399437, 953389, 795289, 704801, 33689, 920555, 314576, 151888, 80329, 722412, 408761, 621423, 662313, 444186, 55356, 831350, 437324, 736680, 69881, 244056, 427000, 983550, 2336, 901983, 663096, 828820, 416684, 282188, 930019, 214180, 43017, 753115, 241571, 522313, 758579, 197631, 561104, 985856, 391867, 374495, 977826, 330204, 46303, 871230, 231994, 966820, 475227, 323159, 758031, 698300, 51096, 306316, 407888, 437207, 985993, 592038, 414900, 942027, 785976, 767530, 150116, 810067, 657576, 945998, 163612, 505392, 140647, 95391, 99385, 112695, 48263, 734299, 749698, 895652, 722469, 735474, 637601, 344983, 74970, 920756, 86755, 824045, 844747, 443462, 98211, 483359, 422347, 22663, 276663, 850780, 26000, 402576, 10476, 807937, 970675, 439388, 974260, 472433, 827009, 694320, 692338, 351467, 258254, 101770, 18000, 145977, 233157, 717440, 64347, 802197, 188238, 381495, 974275, 293923, 211707, 528545, 164562, 897072, 49786, 350994, 655754, 175753, 8337, 470558, 896950, 294036, 775019, 568713, 867412, 470418, 758447, 261201, 821781, 892143, 924250, 825323, 470961, 881053, 545625, 96173, 712220, 756329, 164572, 41726, 785644, 690313, 165348, 89029, 475581, 436181, 455864, 893462, 567943, 731150, 665886, 897428, 760791, 946549, 193876, 464040, 151916, 231113, 555254, 135972, 19197, 613464, 932550, 482415, 583746, 564024, 578486, 493012, 969918, 553664, 460458, 225152, 962347, 566652, 263221, 654569, 63758, 68873, 276515, 942375, 28192, 587706, 716388, 517098, 349063, 308446, 898703, 465519, 217946, 934879, 994474, 711156, 846958, 614468, 858882, 469812, 763810, 658069, 818323, 523238, 738812, 720197, 97000, 842136, 770360, 901885, 777257, 513054, 914601, 814696, 727415, 881255, 892021, 992135, 127660, 53555, 776356, 326106, 332211, 623090, 373137, 871459, 92306, 102230, 27690, 966432, 725733, 590170, 597762, 439352, 419503, 984429, 43839, 768514, 987931, 308133, 539807, 696961, 433372, 335832, 98115, 851772, 904723, 476853, 23979, 134916, 40521, 725112, 84936, 647875, 945882, 269510, 976053, 174409, 234376, 467788, 252644, 528440, 713349, 192228, 333233, 466196, 21853, 805362, 770227, 433988, 651638, 194823, 911330, 470925, 48975, 766538, 331590, 292466, 439739, 586809, 517453, 469390, 570608, 491521, 878745, 586438, 362363, 903683, 691011, 628274, 934650, 763202, 122775, 250809, 561681, 830394, 152554, 811203, 516172, 252619, 132128, 448962, 96760, 472520, 373424, 703213, 198327, 946773, 792429, 929898, 828695, 51009, 236405, 383533, 937271, 304245, 368942, 793021, 829505, 564852, 649240, 133806, 384030, 531429, 262668, 77871, 29989, 530918, 542377, 26098, 74567, 622905, 674351, 938200, 763473, 60429, 743386, 329380, 76455, 663208, 842707, 711873, 682691, 745116, 505843, 759577, 126780, 795718, 516698, 882751, 601023, 41245, 103055, 402182, 249103, 823832, 947291, 221279, 706746, 572234, 501493, 530092, 787210, 229250, 736093, 805979, 801241, 308984, 968616, 128552, 978872, 461424, 826626, 584787, 862188, 592648, 382901, 69311, 197085, 237864, 163481, 131983, 980450, 228813, 656742, 965051, 705796 };
	int8 listsize[8] = { 2, 4, 6, 8, 10, 16, 32 , 48};
	
	srand(time(NULL));
	contidx = 0;
	//	f.open("C:/Users/jwryu/RUG/2018/AlphaTree/AlphaTree_grey_Exp.dat", std::ofstream::app);
	fcheck.open(OUTPUT_FNAME);
	if (fcheck.good())
	{
		cout << "Output file \"" << OUTPUT_FNAME << "\" already exists. Overwrite? (y/n/a)";
		in = 'a';
		cin >> in;
		if (in == 'a')
		{
			f.open(OUTPUT_FNAME, std::ofstream::app);
			cout << "Start from : ";
			contidx = 617;
			cin >> contidx;
		}
		else if (in == 'y')
			f.open(OUTPUT_FNAME);
		else
			exit(-1);
	}
	else
		f.open(OUTPUT_FNAME);

	cnt = 0;
#if RUN_TSE_ONLY
	if (mem_scheme > 0)
		break;
#endif
	for (i = 1; i < 2; i++) // grey, colour loop
	{
		if (i == 0)
			path = INPUTIMAGE_DIR;
		else
			path = INPUTIMAGE_DIR_COLOUR;

		for (auto & p : std::experimental::filesystem::directory_iterator(path))
		{
			cnt++;
			if (cnt < contidx)// || cnt % 10 != 0)
			{
				cout << cnt << ": " << p << endl;
				continue;
			}
			cv::String str1(p.path().string().c_str());
			cv::Mat cvimg;
			if (i == 0)
				cvimg = imread(str1, cv::IMREAD_GRAYSCALE);
			else
			{
				cvimg = imread(str1, cv::IMREAD_COLOR);
				cv::cvtColor(cvimg, cvimg, CV_BGR2GRAY);
			}

			/*
			cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE);// Create a window for display.
			cv::imshow("Display window", cvimg);                   // Show our image inside it.
			cv::waitKey(0);
			getc(stdin);
			*/

			height = cvimg.rows;
			width = cvimg.cols;
			channel = cvimg.channels();

			cout << cnt << ": " << str1 << ' ' << height << 'x' << width << endl;

			uint64 *hdrimg = new uint64[width * height * channel];
			RandomizedHDRimage(hdrimg, cvimg.data, width * height);

			if (channel != 1)
			{
				cout << "input should be a greyscale image" << endl;
				getc(stdin);
				exit(-1);
			}
			int8 listsz[] = { 10, 10, 10 };
			for(int8 listsz_idx = 0;listsz_idx < 2;listsz_idx++)
			{
				double runtime, minruntime;
				for (int testrep = 0; testrep < REPEAT; testrep++)
				{
					memuse = max_memuse = 0;

					tree = (AlphaTree<int32, uint64>*)Malloc(sizeof(AlphaTree<int32, uint64>));
					auto wcts = std::chrono::system_clock::now();

					//tree = (AlphaTree<int32, uint64>*)Malloc(sizeof(AlphaTree<int32, uint64>));
					tree->BuildAlphaTree(hdrimg, height, width, channel, 4, listsz_idx);
					//tree->BuildAlphaTree(testimg_2, 2, 2, channel, 4, listsz_idx);
					//tree->BuildAlphaTree(testimg_64, 30, 30, channel, 4, listsz_idx);

					//tree = (AlphaTree<int32, uint8>*)Malloc(sizeof(AlphaTree<int32, uint8>));
					//tree->BuildAlphaTree(cvimg.data, height, width, channel, 4);
					//tree->BuildAlphaTree(testimg, 4, 4, channel, 4);
					std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
					runtime = wctduration.count();
					minruntime = testrep == 0 ? runtime : min(runtime, minruntime);

					if (testrep < (REPEAT - 1))
						tree->clear();// tree->clear();
					//std::cout << testrep << endl;
				}
				f << p.path().string().c_str() << '\t' << height << '\t' << width << '\t' << max_memuse << '\t' << tree->nrmsd << '\t' << tree->maxSize << '\t' << tree->curSize << '\t' << minruntime << endl;

				pixels_processed = width * height;
				time_elapsed = minruntime;
				cout << "Time Elapsed: " << minruntime << " # Nodes: " << tree->curSize << " mean processing speed(Mpix/s): " << pixels_processed / (time_elapsed * 1000000) << endl;
				tree->clear();
			}

			cvimg.release();
			str1.clear();
			delete[] hdrimg;	
			//break;//tmp
		}
		//break;//tmp
	}

	f.close();
	Free(tree);
	return 0;
}