#pragma once
#pragma once
#include<iostream>
#include<math.h>
#include "defines.h"
#include "HQueue.h"
#include "Trie.h"
#include "HybridQueue.h"
#include <omp.h>

#define DELAYED_NODE_ALLOC		1
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define NODE_CANDIDATE		0xfffffffe

#define dimg_idx_v(pidx) ((pidx)<<1)
#define dimg_idx_h(pidx) ((pidx)<<1)+1

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))


using namespace std;
/*
#define dimg_idx_v(pidx) ((pidx)<<1)
#define dimg_idx_h(pidx) ((pidx)<<1)+1

#define dimg_idx_0(pidx) ((pidx)<<2)
#define dimg_idx_1(pidx) ((pidx)<<2)+1
#define dimg_idx_2(pidx) ((pidx)<<2)+2
#define dimg_idx_3(pidx) ((pidx)<<2)+3

#define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
#define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
#define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
#define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))
*/

#define A		1.3901
#define SIGMA	-2.1989
#define B		-0.1906
#define M		0.05


#define	IMGIDX_32BITS		0
#define	IMGIDX_64BITS		1
#define	PIXEL_8BIT			0
#define	PIXEL_16BIT			1
#define	PIXEL_32BIT			2
#define	PIXEL_64BIT			3
#define	PIXEL_FLOAT			4
#define	PIXEL_DOUBLE		5

//Memory allocation reallocation schemes
/*
#define TSE 0
#define MAXIMUM 1
#define LINEAR 2
#define EXP 3
int mem_scheme = -1;
double size_init[4] = { -1, 1, .2, .15 };
double size_mul[4] = { 1, 1, 1, 2 };
double size_add[4] = { .05, 0, 0.15, 0 };
*/

//tmptmptmptmptmptmtmp
int disp = 0;


template<class Imgidx, class Pixel>
class AlphaNode
{
public:
	Imgidx area;
	Pixel alpha;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	Imgidx parentidx;

	Imgidx rootidx;
	Pixel filter_val;

	inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	{
		this->area = area_in;
		this->alpha = level;
		this->sumPix = sumPix_in;
		this->minPix = minPix_in;
		this->maxPix = maxPix_in;
		this->parentidx = -1;
	}
	inline void add(AlphaNode* q)
	{
		if(disp) printf("Entering add_this\n");
		if(disp) this->print(this);
		this->area += q->area;
		this->sumPix += q->sumPix;
		this->minPix = min(this->minPix, q->minPix);
		this->maxPix = max(this->maxPix, q->maxPix);
		if(disp) printf("Exiting add_this\n");
		if(disp) this->print(this);
		if(disp) getchar();
	}
	inline void add(Pixel pix_val)
	{
		this->area++;
		this->sumPix += (double)pix_val;
		this->minPix = min(this->minPix, pix_val);
		this->maxPix = max(this->maxPix, pix_val);
	}
	inline void copy(AlphaNode* q)
	{
		this->area = q->area;
		this->sumPix = q->sumPix;
		this->minPix = q->minPix;
		this->maxPix = q->maxPix;
	}
	inline void connect_to_parent(AlphaNode* pPar, Imgidx iPar)
	{
		if(disp) printf("conn2par_enter\n");
		this->parentidx = iPar;
		pPar->add(this);
	}

	void print(AlphaNode* node)
	{
		printf("Node idx: %lld\nparent: %lld\narea: %lld\nsumpix: %lf\nminpix: %lld\nmaxpix: %lld\n\n"
		,(uint64)(this - node)
		,(uint64)this->parentidx
		,(uint64)this->area
		,(double)this->sumPix
		,(uint64)this->minPix
		,(uint64)this->maxPix);
	}
};

template<class Imgidx, class Pixel>
class RankItem
{
public:
	Pixel alpha;
	Imgidx dimgidx;
	Imgidx p, q;

	inline void operator=(const RankItem& item)
	{
		this->alpha = item.alpha;
		this->dimgidx = item.dimgidx;
		this->p = item.p;
		this->q = item.q;
	}
};


template<class Imgidx, class Pixel>
class ATree
{
	inline Pixel abs_diff(Pixel p, Pixel q)
	{
		if (p > q)
		return p - q;
		else
		return q - p;
	}
	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>*& rankinfo, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;

		if (connectivity == 4)
		nredges = (width - 1) * height + width * (height - 1);
		else
		nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);

		//rankinfo = (RankItem1<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem1<Imgidx, Pixel>));
		tmp = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc((sizeof(Pixel) >> 1) * 65536 * sizeof(Imgidx));

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 2;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  -
			//   -  p  3
			//   0  1  2
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				dimgidx++; //skip 0
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width - 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + width + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					rankinfo[contidx].p = imgidx;
					rankinfo[contidx].q = imgidx + 1;
					rankinfo[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width - 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + width;
				rankinfo[contidx++].dimgidx = dimgidx;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				rankinfo[contidx].alpha = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				rankinfo[contidx].p = imgidx;
				rankinfo[contidx].q = imgidx + 1;
				rankinfo[contidx++].dimgidx = dimgidx++;
				imgidx++;
			}
		}

		for (i = 0; i < (Imgidx)((sizeof(Pixel) >> 1) * 65536); i++)
		hist[i] = 0;
		for (i = 0; i < nredges; i++)
		{
			hidx = rankinfo[i].alpha;
			hist[hidx & mask]++;

			shamt = (Pixel)16;
			h_offset = (Pixel)65536;
			for (nbyte = 2; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
			{
				hidx = hidx >> shamt;
				hist[h_offset + (hidx & mask)]++;
				h_offset += (Pixel)65536;
			}
		}

		h = hist;
		shamt = 0;
		for (nbyte = 0; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
		{
			hsum = 0;
			for (i = 0; i < (Imgidx)(hist_size); i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = (rankinfo[i].alpha >> shamt) & mask;
				j = --h[hidx];
				tmp[j] = rankinfo[i]; //slow!
			}
			r = rankinfo; rankinfo = tmp; tmp = r; //swap p,q
			shamt += 16;
			h += hist_size;
		}

		for (i = 0; i < nredges; i++)
		{
			//rank_to_alpha[i] = rankinfo[i].alpha;
			rank[rankinfo[i].dimgidx] = i;
		}

		//std::cout << "in_func2: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

		//Free(rankinfo);
		Free(tmp);
		Free(hist);
	}
	void compute_dimg(Pixel* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + stride_w] - (int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + stride_w] - (int64)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  -
			//   -  p  3
			//   0  1  2
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				dimgidx++; //skip 0
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				dhist[dimg[dimgidx++]]++;
				imgidx++;
				for (j = 1; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width + 1] - (int64)img[imgidx]));//2
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width - 1] - (int64)img[imgidx]));//0
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + width] - (int64)img[imgidx]));//1
				dhist[dimg[dimgidx]]++;
				dimgidx += 3;//skip 2,3
				imgidx++;
			}

			//bottom
			for (j = 0; j < width - 1; j++)
			{
				dimgidx += 3; //skip 0,1,2
				dimg[dimgidx] = (Pixel)(abs((int64)img[imgidx + 1] - (int64)img[imgidx]));//3
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
	}
	void init_isAvailable(uint8* isAvailable)
	{
		int32 i, j, k;
		int32 imgsize = width * height;

		if (connectivity == 4)
		{

			//		    Neighbour Index
			// 			-      3      -
			// 			2    pixel    1
			// 			-      0      -
			//
			//			Neighbour indices to bit field
			//			7 6 5 4 3 2 1 0
			//         MSB			 LSB
			//			0: Neighbour pixel not available (corner of Image or partition)
			//			1: available
			for (i = 0; i < ((imgsize + 1) >> 1); i++)
			isAvailable[i] = 0xff;

			set_field(isAvailable, 0, 0x3);
			set_field(isAvailable, width - 1, 0x5);
			set_field(isAvailable, width * (height - 1), 0xa);
			set_field(isAvailable, width * height - 1, 0xc);


			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				set_field(isAvailable, i, 0x7);
				set_field(isAvailable, j, 0xe);
				j++;
			}

			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				set_field(isAvailable, j, 0xb);
				set_field(isAvailable, k, 0xd);
				j += width;
				k += width;
			}
		}
		else
		{
			//		    Neighbour Index
			// 			6      5      4
			// 			7    pixel    3
			// 			0      1      2
			//
			//			Neighbour indices to bit field
			//			7 6 5 4 3 2 1 0
			//         MSB			 LSB
			//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
			//			1: available

			//initialize to all available
			for (i = 0; i < imgsize; i++)
			isAvailable[i] = 0xff;

			//four corners
			isAvailable[0] = 0x0e;
			isAvailable[width - 1] = 0x83;
			isAvailable[width*(height - 1)] = 0x38;
			isAvailable[width * height - 1] = 0xe0;

			//top and bottom row
			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				isAvailable[i] = 0x8f;
				isAvailable[j] = 0xf8;
				j++;
			}

			//leftest and rightest column
			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				isAvailable[j] = 0x3e;
				isAvailable[k] = 0xe3;
				j += width;
				k += width;
			}
		}
	}

	void init_isAvailable_par(uint8* isAvailable, int npartitions)
	{
		int32 i, j, k;
		Imgidx imgsize = width * height;
		Imgidx wstride = width / npartitions;
		Imgidx hstride = height / npartitions;

		init_isAvailable(isAvailable);

		if (connectivity == 4)
		{
			//hor partitions
			j = (hstride - 1) * width;
			for (i = 0; i < npartitions - 1; i++)
			{
				k = j + width;
				//set_field(isAvailable, j, 0xa);
				//set_field(isAvailable, j + width, 0x3);
				for (; j < k; j++)
				{
					set_field(isAvailable, j, 0xe);
					set_field(isAvailable, j + width, 0x7);
				}
				//set_field(isAvailable, j, 0xc);
				//set_field(isAvailable, j + width, 0x5);

				j += (hstride - 1) * width;
			}

			//ver partitions

			for (i = 0; i < npartitions - 1; i++)
			{
				j = (i + 1) * wstride - 1;
				//k = j + height;
				//set_field(isAvailable, j, 0x5);
				//set_field(isAvailable, j + 1, 0x3);
				for (; j < imgsize; j += width)
				{
					set_field(isAvailable, j, 0xd);
					set_field(isAvailable, j + 1, 0xb);
				}
				//set_field(isAvailable, j, 0xc);
				//set_field(isAvailable, j + width, 0xa);
			}

			k = 0;
			for(i = 0;i < height;i++)
			{
				for(j = 0; j< width; j++)
				{
					uint16 p = get_field(isAvailable, k++);
					std::cout << p << '\t';
				}
				std::cout << std::endl;
			}

			/*		//partition crossings - no need
			k = (hstride - 1) * width;
			for(i = 0;i < npartitions - 1;i++)
			{
			k--;
			for(j = 0;j < npartitions - 1;j++)
			{
			k += wstride;
			set_field(isAvailable, k, 0xc);
			set_field(isAvailable, k + 1, 0xa);
			set_field(isAvailable, k + width, 0x5);
			set_field(isAvailable, k + width + 1, 0x3);
		}
		k += (hstride - 1) * width;
	}
	*/
}
else
{
	//		    Neighbour Index
	// 			6      5      4
	// 			7    pixel    3
	// 			0      1      2
	//
	//			Neighbour indices to bit field
	//			7 6 5 4 3 2 1 0
	//         MSB			 LSB
	//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
	//			1: available

	//initialize to all available
	for (i = 0; i < imgsize; i++)
	isAvailable[i] = 0xff;

	//four corners
	isAvailable[0] = 0x0e;
	isAvailable[width - 1] = 0x83;
	isAvailable[width*(height - 1)] = 0x38;
	isAvailable[width * height - 1] = 0xe0;

	//top and bottom row
	j = width * (height - 1) + 1;
	for (i = 1; i < width - 1; i++)
	{
		isAvailable[i] = 0x8f;
		isAvailable[j] = 0xf8;
		j++;
	}

	//leftest and rightest column
	j = width;
	k = (width << 1) - 1;
	for (i = 1; i < height - 1; i++)
	{
		isAvailable[j] = 0x3e;
		isAvailable[k] = 0xe3;
		j += width;
		k += width;
	}
}
}

inline uint8 is_available(uint8 isAvailable, uint8 iNeighbour)
{
	//return	(((isAvailable[idx >> 1] >> ((idx & 1) << 2)) & 0x0f) >> iNeighbour) & 1;
	return	(isAvailable >> iNeighbour) & 1;
}

inline void set_field(uint8* arr, Imgidx idx, uint8 in)
{
	uint8 shamt = (idx & 1) << 2;
	arr[idx >> 1] &= (in << shamt) | ((uint8)(0x0f) << (4 - shamt));
}

inline uint8 get_field(uint8* arr, Imgidx idx)
{
	return (arr[idx >> 1] >> ((idx & 1) << 2)) & 0x0f;
}

inline void push_neighbour(HQueue<Imgidx> *hqueue, Imgidx* levelroot, Pixel* dimg, Imgidx idx, Imgidx dimgidx)
{
	Pixel dissim = dimg[dimgidx];
	hqueue->push(idx, dissim);
	if (levelroot[dissim] == NULL_LEVELROOT)
	#if DELAYED_NODE_ALLOC
	levelroot[dissim] = NODE_CANDIDATE;
	#else
	levelroot[dissim] = NewAlphaNode1(dissim);
	#endif
}

inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
{
	AlphaNode<Imgidx, Pixel>* pNode;
	pNode = node + iNode;
	parentAry[pidx] = iNode;
	if (pNode->area) //possibly unnecessary branch..
	pNode->add(pix_val);
	else
	pNode->set(1, level, (double)pix_val, pix_val, pix_val);
}


#if DELAYED_NODE_ALLOC
inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx *levelroot, int32 level)
{
	AlphaNode<Imgidx, Pixel>* pNode;
	Imgidx iNode = levelroot[level];
	if (iNode == NODE_CANDIDATE)
	{
		iNode = NewAlphaNode();
		levelroot[level] = iNode;
		parentAry[pidx] = iNode;
		pNode = node + iNode;

		pNode->set(1, level, (double)pix_val, pix_val, pix_val);
	}
	else
	{
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		pNode->add(pix_val);
	}
}
#else
inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode)
{
	AlphaNode<Imgidx, Pixel> *pNode = &node[iNode];
	parentAry[pidx] = iNode;
	pNode->add(pix_val);
}
#endif


inline void connectNode2Node(Imgidx iChild, Imgidx iPar, Pixel level)
{
	AlphaNode<Imgidx, Pixel> *pPar, *pChild;
	pChild = node + iChild;
	pPar = node + iPar;
	pChild->parentidx = iPar;
	if (pPar->area)
	pPar->add(pChild);
	else
	{
		pPar->alpha = level;
		pPar->copy(pChild);
	}
}
#if DELAYED_NODE_ALLOC
inline void connectNode2Node(Imgidx* levelroot, Imgidx iChild, int32 level)
{
	AlphaNode<Imgidx, Pixel> *pPar, *pChild;
	Imgidx iPar = levelroot[level];
	if (iPar == NODE_CANDIDATE)
	{
		iPar = NewAlphaNode();
		levelroot[level] = iPar;
		pPar = node + iPar;
		pChild = node + iChild;
		pChild->parentidx = iPar;
		pPar->alpha = level;
		pPar->copy(pChild);
	}
	else
	{
		pPar = node + iPar;
		pChild = node + iChild;
		pChild->parentidx = iPar;
		pPar->add(pChild);
	}
}
#else
inline void connectNode2Node(AlphaNode<Imgidx, Pixel>* pPar, Imgidx iPar, AlphaNode<Imgidx, Pixel>* pNode)
{
	pNode->parentidx = iPar;
	pPar->add(pNode);
}
#endif
#if DELAYED_NODE_ALLOC
inline Imgidx NewAlphaNode()
{
	//		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

	if (curSize == maxSize)
	{
		std::cout << "Reallocating...\n";
		maxSize = min(2 * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

		node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		//			pNew = node + curSize;
	}
	return curSize++;
}
#else
inline Imgidx NewAlphaNode(Pixel level) //Fix it later - no need to initialize
{
	AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

	if (curSize == maxSize)
	{
		std::cout << "Reallocating...\n";
		maxSize = min(height * width, maxSize + (Imgidx)(height * width * 0.1));

		node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		pNew = node + curSize;
	}
	pNew->alpha = level;
	pNew->minPix = (uint8)-1;
	pNew->maxPix = 0;
	pNew->sumPix = 0.0;
	//		pNew->parentidx = 0;
	pNew->area = 0;

	return curSize++;
}
#endif
inline uint8 is_visited(uint8* isVisited, Imgidx p)
{
	return (isVisited[p >> 3] >> (p & 7)) & 1;
}
inline void visit(uint8* isVisited, Imgidx p)
{
	isVisited[p >> 3] = isVisited[p >> 3] | (1 << (p & 7));
}

void Flood_HQueue(Pixel* img)
{
	Imgidx imgsize, dimgsize, nredges, x0;
	int32 numlevels, max_level, current_level, next_level;
	HQueue<Imgidx>* hqueue;
	Imgidx *dhist;
	Pixel *dimg;
	Imgidx iChild, *levelroot;
	uint8 *isVisited, *isAvailable, isAv;;
	Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

	imgsize = width * height;
	nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
	dimgsize = (connectivity >> 1) * width * height;
	numlevels = 1 << (8 * sizeof(uint8));

	dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
	dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
	levelroot = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	isVisited = (uint8*)Malloc((size_t)((imgsize + 7) >> 3));
	if (connectivity == 4)
	isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
	else
	isAvailable = (uint8*)Malloc((size_t)(imgsize));
	for (p = 0; p < numlevels; p++)
	levelroot[p] = NULL_LEVELROOT;
	memset(dhist, 0, (size_t)numlevels * sizeof(int32));
	memset(isVisited, 0, (size_t)((imgsize + 7) >> 3));
	init_isAvailable(isAvailable);

	max_level = (uint8)(numlevels - 1);

	compute_dimg(dimg, dhist, img);

	dhist[max_level]++;
	hqueue = new HQueue<Imgidx>(nredges + 1, dhist, numlevels);
	curSize = 0;

	/////////////////////////////////////////
	//tree size estimation (TSE)
	nrmsd = 0;
	for (p = 0; p < numlevels; p++)
	nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
	nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
	maxSize = min(2 * imgsize, (int32)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + M)));
	//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	/////////////////////////////////////////
	//	printf("NRMSD: %f\tEst. NTS: %f\tEst. Tree size: %d\n", nrmsd, ((A * exp(SIGMA * nrmsd) + B) + M), tree->maxSize);
	//maxSize = (int32)(2 * imgsize * size_init[mem_scheme]);
	Free(dhist);

	parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
	node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

	#if DELAYED_NODE_ALLOC
	//		levelroot[max_level + 1] = NODE_CANDIDATE;
	levelroot[max_level + 1] = NewAlphaNode();
	AlphaNode<Imgidx, Pixel> *pNode = node + levelroot[max_level+1];
	pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
	//(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
	//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

	#else
	levelroot[max_level + 1] = NewAlphaNode((uint8)max_level);
	node[levelroot[max_level + 1]].parentidx = levelroot[max_level + 1];
	#endif

	current_level = max_level;
	x0 = imgsize >> 1;
	hqueue->push(x0, current_level);

	iChild = levelroot[max_level + 1];
	while (current_level <= max_level)
	{
		while (hqueue->min_level <= current_level)
		{
			p = hqueue->pop();
			if (is_visited(isVisited, p))
			{
				hqueue->find_min_level();
				continue;
			}
			visit(isVisited, p);
			#if !HQUEUE_COST_AMORTIZE
			hqueue->find_min_level();
			#endif
			if (connectivity == 4)
			{
				isAv = get_field(isAvailable, p);
				q = p << 1;

				if (is_available(isAv, 0) && !is_visited(isVisited, p + width))		push_neighbour(hqueue, levelroot, dimg, p + width, q);
				if (is_available(isAv, 1) && !is_visited(isVisited, p + 1))				push_neighbour(hqueue, levelroot, dimg, p + 1, q + 1);
				if (is_available(isAv, 2) && !is_visited(isVisited, p - 1))				push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
				if (is_available(isAv, 3) && !is_visited(isVisited, p - width))		push_neighbour(hqueue, levelroot, dimg, p - width, q - (width << 1));
			}
			else
			{
				isAv = isAvailable[p];
				q = p << 2;
				if (is_available(isAv, 0) && !is_visited(isVisited, p + wstride1))		push_neighbour(hqueue, levelroot, dimg, p + wstride1, q);
				if (is_available(isAv, 1) && !is_visited(isVisited, p + width))				push_neighbour(hqueue, levelroot, dimg, p + width, q + 1);
				if (is_available(isAv, 2) && !is_visited(isVisited, p + wstride0))		push_neighbour(hqueue, levelroot, dimg, p + wstride0, q + 2);
				if (is_available(isAv, 3) && !is_visited(isVisited, p + 1))						push_neighbour(hqueue, levelroot, dimg, p + 1, q + 3);
				if (is_available(isAv, 4) && !is_visited(isVisited, p - wstride1))		push_neighbour(hqueue, levelroot, dimg, p - wstride1, q - wstride_d + 4);
				if (is_available(isAv, 5) && !is_visited(isVisited, p - width))				push_neighbour(hqueue, levelroot, dimg, p - width, q - wstride_d + 1);
				if (is_available(isAv, 6) && !is_visited(isVisited, p - wstride0))		push_neighbour(hqueue, levelroot, dimg, p - wstride0, q - wstride_d - 2);
				if (is_available(isAv, 7) && !is_visited(isVisited, p - 1))						push_neighbour(hqueue, levelroot, dimg, p - 1, q - 1);
			}

			if (current_level > hqueue->min_level)
			current_level = hqueue->min_level;
			#if HQUEUE_COST_AMORTIZE
			else
			hqueue->find_min_level();
			#endif

			#if DELAYED_NODE_ALLOC
			connectPix2Node(p, img[p], levelroot, current_level);
			#else
			connectPix2Node(p, img[p], node + levelroot[current_level], levelroot[current_level]);
			#endif

		}
		//std::cout << "checking redundancy..." << std::endl;
		//std::cout << "node[iChild].parentidx: " << node[iChild].parentidx << std::endl;
		//std::cout << "levelroot[current_level]: " << levelroot[current_level] << std::endl;
		if (node[iChild].parentidx == levelroot[current_level] &&
			node[levelroot[current_level]].area == node[iChild].area)
			{
				levelroot[current_level] = iChild;
				#if DELAYED_NODE_ALLOC
				curSize--;
				#endif
			}

			next_level = current_level + 1;
			while (next_level <= max_level && (levelroot[next_level] == NULL_LEVELROOT))
			next_level++;

			#if DELAYED_NODE_ALLOC
			connectNode2Node(levelroot, levelroot[current_level], next_level);
			#else
			connectNode2Node(node + levelroot[next_level], levelroot[next_level], node + levelroot[current_level]);
			#endif


			iChild = levelroot[current_level];
			levelroot[current_level] = NULL_LEVELROOT;
			current_level = next_level;

		}
		node[iChild].parentidx = iChild;
		rootidx = iChild;
		curSize--;

		for (p = 0; p < imgsize; p++)
		{
			if (node[parentAry[p]].alpha)//Singleton 0-CC
			{
				x0 = NewAlphaNode();
				(&node[x0])->set(1, 0, (double)img[p], img[p], img[p]);
				node[x0].parentidx = parentAry[p];
				parentAry[p] = x0;
			}
		}

		delete hqueue;
		Free(dimg);
		Free(levelroot);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_Trie(Pixel* img)
	{
		Imgidx imgsize, dimgsize, nredges, x0;
		Imgidx current_rank, next_rank;
		Trie<Imgidx, trieidx> *queue;
		RankItem<Imgidx, Pixel>* rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx dimgidx;
		Imgidx *rank, top_rank;
		int8 incidence, shamt, mask;
		//		Pixel *rank2alpha;
		Imgidx iChild, *levelroot;
		uint8 *isVisited, *isVisited_edges, *isAvailable, isAv;
		Imgidx trietop, p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;
		Imgidx incidence_map[4];

		if (connectivity == 4)
		{
			incidence_map[0] = width;
			incidence_map[1] = 1;
			shamt = 1;
			mask = 1;
		}
		else
		{
			incidence_map[0] = width - 1;
			incidence_map[1] = width;
			incidence_map[2] = width + 1;
			incidence_map[3] = 1;
			shamt = 2;
			mask = 3;
		}

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;

		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize)));
		isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
		if (connectivity == 4)
		isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; //
		memset(isVisited, 0, (size_t)((imgsize)));
		memset(isVisited_edges, 0, (size_t)((nredges)));
		init_isAvailable(isAvailable);

		compute_dimg(rank, rankinfo, img);
		queue = new Trie<Imgidx, trieidx>(nredges);
		//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

		//trie = new Trie<Imgidx, trieidx>(nredges << 1);
		//trie->push(nredges, 0);


		// 		//tree size estimation (TSE)
		//
		// 		nrmsd = 0;
		// 		for (p = 0; p < numlevels; p++)
		// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		// 		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
		// 		//maxSize = imgsize;
		// 		Free(dhist);
		maxSize = imgsize + nredges;

		parentAry = (int32*)Malloc((size_t)imgsize * sizeof(int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;


		for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
		pNode->set(1, 0, (double)img[p], img[p], img[p]);

		maxpixval = ~(1 << ((sizeof(Pixel) << 3) - 1));
		for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
		pNode->set(0, pRank->alpha, 0.0, maxpixval, 0);

		isVisited[0] = 1;
		if (connectivity == 4)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
			// 			trie->push(rank[0], 0);
			// 			trie->push(rank[1], 0);
		}
		else
		{
			// 			trie->push(rank[1], 0);
			// 			trie->push(rank[2], 0);
			// 			trie->push(rank[3], 0);
			queue->push(rank[1]);
			queue->push(rank[2]);
			queue->push(rank[3]);

		}

		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank);
		isVisited_edges[current_rank] = 1;
		//connectPix2Node(parentAry, 0, img[0], current_rank, rankinfo[current_rank].alpha);
		//x0 = imgsize >> 1;
		iChild = current_rank;

		while (1)//(current_rank <= nredges)
		{
			while (1)//((trie->top() >> 1) <= current_rank)
			{
				//trietop = queue->top();		//remove tmp variables later if possible
				//incidence = queue->top_incidence();	//0 is outgoing, 1 is incoming
				top_rank = queue->top();	//remove tmp variables later if possible
				//dimgidx = rankinfo[top_rank].dimgidx;
				//p = (dimgidx >> shamt) + (incidence_map[dimgidx & mask] & (incidence - 1)); //current pixel idx
				if (isVisited_edges[top_rank])
				break;
				pRank = rankinfo + top_rank;
				if (isVisited[pRank->p])
				{
					if (isVisited[pRank->q])
					break;
					p = pRank->q;
				}
				else
				p = pRank->p;

				isVisited[p] = 1;
				#if !HQUEUE_COST_AMORTIZE
				//find_min_level();
				#endif
				if (connectivity == 4)
				{
					isAv = get_field(isAvailable, p);
					q = p << 1;

					if (is_available(isAv, 0) && !isVisited[p + width])
					queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])
					queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])
					queue->push(rank[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])
					queue->push(rank[q - (width << 1)]);
				}
				else
				{
					isAv = isAvailable[p];
					q = p << 2;

					if (is_available(isAv, 0) && !isVisited[p + wstride1])
					queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + width])
					queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + wstride0])
					queue->push(rank[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p + 1])
					queue->push(rank[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - wstride1])
					queue->push(rank[q - wstride_d + 4]);
					if (is_available(isAv, 5) && !isVisited[p - width])
					queue->push(rank[q - wstride_d + 1]);
					if (is_available(isAv, 6) && !isVisited[p - wstride0])
					queue->push(rank[q - wstride_d - 2]);
					if (is_available(isAv, 7) && !isVisited[p - 1])
					queue->push(rank[q - 1]);
				}

				next_rank = queue->top();
				node[p].connect_to_parent(&node_in[next_rank], next_rank);

				if (current_rank == next_rank)
				break;
				current_rank = next_rank;
			}

			//Redundant node removal
			// 			if (node[iChild].parentidx == levelroot[current_rank] &&
			// 				node[levelroot[current_rank]].area == node[iChild].area)
			// 			{
			// 				levelroot[current_rank] = iChild;
			// #if DELAYED_NODE_ALLOC
			// 				curSize--;
			// 				//memset((uint8*)(node + curSize), 0, sizeof(AlphaNode));
			// #endif
			// 			}

			queue->pop();
			next_rank = queue->top();

			//Redundant node removal
			if (node_in[iChild].parentidx == current_rank &&
				node_in[iChild].area == node_in[current_rank].area)
					current_rank = iChild;

				node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
				//connectNode2Node(current_rank, next_rank, rankinfo[next_rank].alpha);
				// #if DELAYED_NODE_ALLOC
				// 			connectNode2Node(current_rank, next_rank, rankinfo[current_rank].alpha);
				// #else
				// 			connectNode2Node(node + levelroot[next_rank], levelroot[next_rank], node + levelroot[current_rank]);
				//#endif
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
				next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				node_in[next_rank].parentidx = next_rank;
				break;
			}

			iChild = current_rank;
			current_rank = next_rank;
		}


		for (pNode = node; pNode < node + maxSize; pNode++)
			pNode->parentidx += imgsize;
		node_in[next_rank].parentidx = -1;

		// 		curSize = 0;
		// 		for (p = 0; p < nredges; p++)
		// 		{
		// 			if (node[p].parentidx > 0)
		// 				curSize++;
		// 		}

		// 		q = nredges;
		// 		for (p = 0; p < imgsize; p++)
		// 		{
		// 			if (node[parentAry[p]].alpha)//Singleton 0-CC
		// 			{
		// 				(&node[q])->set(1, 0, (double)img[p], img[p], img[p]);
		// 				node[q++].parentidx = parentAry[p];
		// 				parentAry[p] = q;
		// 				curSize++;
		// 			}
		// 		}


		//		delete hqueue;
		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
		//		Free(levelroot);
		Free(isVisited);
		Free(isVisited_edges);
		Free(isAvailable);
	}

		void Flood_HybridQueue(Pixel* img, int8 listsize)
		{
			Imgidx imgsize, dimgsize, nredges;
			Imgidx current_rank, next_rank;
			HybridQueue<Imgidx, trieidx> *queue;
			RankItem<Imgidx, Pixel> *rankinfo, *pRank;
			AlphaNode<Imgidx, Pixel> *pNode;
			Pixel maxpixval;
			Imgidx *rank, top_rank;
			int8 nbits;

			Imgidx iChild;
			uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
			Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

			imgsize = width * height;
			nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
			dimgsize = (connectivity >> 1) * width * height;

			rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
			rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
			//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
			//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
			isVisited = (uint8*)Malloc((size_t)((imgsize)));
			//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
			if (connectivity == 4)
			isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
			else
			isAvailable = (uint8*)Malloc((size_t)(imgsize));

			//	levelroot[p] = NULL_LEVELROOT;
			//levelroot[nredges] = NODE_CANDIDATE; //
			memset(isVisited, 0, (size_t)((imgsize)));
			//		memset(isVisited_edges, 0, (size_t)((nredges)));
			init_isAvailable(isAvailable);

			//	std::cout << "before:" << rankinfo->alpha << std::endl;
			//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
			compute_dimg(rank, rankinfo, img);
			//std::cout << "after:" << rankinfo->alpha << std::endl;
			//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

			queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);
			//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

			//trie = new Trie<Imgidx, trieidx>(nredges << 1);
			//trie->push(nredges, 0);



			// 		//tree size estimation (TSE)
			//
			// 		nrmsd = 0;
			// 		for (p = 0; p < numlevels; p++)
			// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
			// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
			// 		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
			// 		//maxSize = imgsize;
			// 		Free(dhist);
			maxSize = imgsize + nredges;
			num_node = maxSize;
			num_node_in = nredges;

			parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
			node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			node_in = node + imgsize;


			for (pNode = node, p = 0; pNode < node_in; pNode++, p++)
			pNode->set(1, 0, (double)img[p], img[p], img[p]);

			nbits = ((sizeof(Pixel) << 3) - 1);
			maxpixval = ~(1 << nbits);
			for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0);
			//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

			isVisited[0] = 1;
			if (connectivity == 4)
			{
				queue->push(rank[0]);
				queue->push(rank[1]);
			}
			else
			{
				queue->push(rank[1]);
				queue->push(rank[2]);
				queue->push(rank[3]);
			}

			current_rank = queue->top();
			node[0].connect_to_parent(&node_in[current_rank], current_rank);
			//isVisited_edges[current_rank] = 1;
			iChild = current_rank;

			while (1)//(current_rank <= nredges)
			{
				while (1)//((trie->top() >> 1) <= current_rank)
				{
					top_rank = queue->top();	//remove tmp variables later if possible
					//	if (isVisited_edges[top_rank])
					//		break;
					pRank = rankinfo + top_rank;
					if (isVisited[pRank->p])
					{
						if (isVisited[pRank->q])
						break;
						p = pRank->q;
					}
					else
					p = pRank->p;

					isVisited[p] = 1;
					#if !HQUEUE_COST_AMORTIZE
					//find_min_level();
					#endif
					if (connectivity == 4)
					{
						isAv = get_field(isAvailable, p);
						q = p << 1;

						if (is_available(isAv, 0) && !isVisited[p + width])
						queue->push(rank[q]);
						if (is_available(isAv, 1) && !isVisited[p + 1])
						queue->push(rank[q + 1]);
						if (is_available(isAv, 2) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
						if (is_available(isAv, 3) && !isVisited[p - width])
						queue->push(rank[q - (width << 1)]);
					}
					else
					{
						isAv = isAvailable[p];
						q = p << 2;

						if (is_available(isAv, 0) && !isVisited[p + wstride1])
						queue->push(rank[q]);
						if (is_available(isAv, 1) && !isVisited[p + width])
						queue->push(rank[q + 1]);
						if (is_available(isAv, 2) && !isVisited[p + wstride0])
						queue->push(rank[q + 2]);
						if (is_available(isAv, 3) && !isVisited[p + 1])
						queue->push(rank[q + 3]);
						if (is_available(isAv, 4) && !isVisited[p - wstride1])
						queue->push(rank[q - wstride_d + 4]);
						if (is_available(isAv, 5) && !isVisited[p - width])
						queue->push(rank[q - wstride_d + 1]);
						if (is_available(isAv, 6) && !isVisited[p - wstride0])
						queue->push(rank[q - wstride_d - 2]);
						if (is_available(isAv, 7) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
					}

					next_rank = queue->top();
					/*				if(p == 27579)
					{
					node[p].print(node);
					node[5773602].print(node);
					getchar();

					disp = 0;
				}
				*/
				node[p].connect_to_parent(&node_in[next_rank], next_rank);

				/*				disp = 0;
				if(p == 27579)
				{
				node[p].print(node);
				node_in[node[p].parentidx].print(node);
				getchar();
				}
				*/
				if (current_rank == next_rank)
				break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			//Redundant node removal
			if (node_in[iChild].parentidx == current_rank &&
				node_in[iChild].area == node_in[current_rank].area)
				current_rank = iChild;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
					next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;

				break;
			}

			iChild = current_rank;
			current_rank = next_rank;
		}

		for (pNode = node; pNode < node + maxSize; pNode++)
		pNode->parentidx += imgsize;
		node_in[next_rank].parentidx = -1;

		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
		//		Free(levelroot);
		Free(isVisited);
		//Free(isVisited_edges);
		Free(isAvailable);
	}

	void par_Flood_HybridQueue(Pixel* img, int8 listsize)
	{
		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank, next_rank;
		HybridQueue<Imgidx, trieidx> *queue, **queues;
		RankItem<Imgidx, Pixel> *rankinfo, *pRank;
		AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank, top_rank;
		int8 nbits;

		Imgidx iChild;
		uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
		Imgidx p, q, wstride_d = width << 2, wstride0 = width + 1, wstride1 = width - 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;

		rankinfo = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		//rank2alpha = (Pixel*)Malloc((size_t)nredges * sizeof(Pixel));
		//levelroot = (Imgidx*)Malloc((Imgidx)(nredges + 1) * sizeof(Imgidx));
		isVisited = (uint8*)Malloc((size_t)((imgsize)));
		//isVisited_edges = (uint8*)Malloc((size_t)((nredges)));
		if (connectivity == 4)
		isAvailable = (uint8*)Malloc((size_t)((imgsize + 1) >> 1));
		else
		isAvailable = (uint8*)Malloc((size_t)(imgsize));

		//	levelroot[p] = NULL_LEVELROOT;
		//levelroot[nredges] = NODE_CANDIDATE; //
		memset(isVisited, 0, (size_t)((imgsize)));
		//		memset(isVisited_edges, 0, (size_t)((nredges)));

		int8 npartition = 2;

		init_isAvailable_par(isAvailable, npartition); //yogi

		//	int8 num_thread = omp_get_num_threads();

		int64 blksz_x = (width - (width % npartition)) / npartition;
		int64 blksz_y = (height - (height % npartition)) / npartition;
		int64 blksz_xn = blksz_x + (width % npartition);
		int64 blksz_yn = blksz_y + (height % npartition);
		int64 blksz = blksz_x * blksz_y;
		int64 blksz_lastcol = blksz_xn * blksz_y;
		int64 blksz_lastrow = blksz_x * blksz_yn;
		int64 blksz_last = blksz_xn * blksz_yn;

		queues = new HybridQueue<Imgidx, trieidx>*[npartition * npartition];//(nredges, listsize);
		p = 0;
		for(int x = 0;x < npartition; x++)
		{
			for(int y = 0;y < npartition - 1; y++)
			queues[p++] = new HybridQueue<Imgidx, trieidx>(blksz, listsize);
			queues[p++] = new HybridQueue<Imgidx, trieidx>(blksz_lastrow, listsize);
		}
		for(int y = 0;y < npartition - 1; y++)
		queues[p++] = new HybridQueue<Imgidx, trieidx>(blksz_lastcol, listsize);
		queues[p++] = new HybridQueue<Imgidx, trieidx>(blksz_last, listsize);


		//		#pragma omp parallel default(shared)
		{
			//	std::cout << "before:" << rankinfo->alpha << std::endl;
			//		std::cout << "before: " << rankinfo << " " << rankinfo[0].alpha << std::endl;
			compute_dimg(rank, rankinfo, img);
			//std::cout << "after:" << rankinfo->alpha << std::endl;
			//		std::cout << "after: " << rankinfo << " " << rankinfo[0].alpha << std::endl;

			queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);
			//		queue = new HybridQueue<Imgidx, trieidx>(nredges, listsize);

			//trie = new Trie<Imgidx, trieidx>(nredges << 1);
			//trie->push(nredges, 0);



			// 		//tree size estimation (TSE)
			//
			// 		nrmsd = 0;
			// 		for (p = 0; p < numlevels; p++)
			// 			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
			// 		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
			// 		maxSize = min(imgsize, (Imgidx)(imgsize * A * (exp(SIGMA * nrmsd) + B + M)));
			// 		//maxSize = imgsize;
			// 		Free(dhist);
			maxSize = imgsize + nredges;
			num_node = maxSize;
			num_node_in = nredges;

			parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(int32));
			node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			node_in = node + imgsize;

			std::cout << "initiating nodes" << std::endl;
			//		pNode = node;pNode++
			//		#pragma omp parallel for
			for (p = 0; p < imgsize;p++)
			{
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
				//node[p].print(node);
			}


			pNode = node_in;
			nbits = ((sizeof(Pixel) << 3) - 1);
			maxpixval = ~(1 << nbits);
			for (pRank = rankinfo; pNode < node + maxSize; pNode++, pRank++)
			pNode->set(0, pRank->alpha, 0.0, maxpixval, 0);
			//inline void set(Imgidx area_in, Pixel level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)

			isVisited[0] = 1;
			if (connectivity == 4)
			{
				queue->push(rank[0]);
				queue->push(rank[1]);
			}
			else
			{
				queue->push(rank[1]);
				queue->push(rank[2]);
				queue->push(rank[3]);
			}

			current_rank = queue->top();
			node[0].connect_to_parent(&node_in[current_rank], current_rank);
			//isVisited_edges[current_rank] = 1;
			iChild = current_rank;

			while (1)//(current_rank <= nredges)
			{
				while (1)//trie->top() >> 1) <= current_rank)
				{
					top_rank = queue->top();	//remove tmp variables later if possible
					//	if (isVisited_edges[top_rank])
					//		break;
					pRank = rankinfo + top_rank;
					if (isVisited[pRank->p])
					{
						if (isVisited[pRank->q])
						break;
						p = pRank->q;
					}
					else
					p = pRank->p;

					isVisited[p] = 1;
					#if !HQUEUE_COST_AMORTIZE
					//find_min_level();
					#endif
					if (connectivity == 4)
					{
						isAv = get_field(isAvailable, p);
						q = p << 1;

						if (is_available(isAv, 0) && !isVisited[p + width])
						queue->push(rank[q]);
						if (is_available(isAv, 1) && !isVisited[p + 1])
						queue->push(rank[q + 1]);
						if (is_available(isAv, 2) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
						if (is_available(isAv, 3) && !isVisited[p - width])
						queue->push(rank[q - (width << 1)]);
					}
					else
					{
						isAv = isAvailable[p];
						q = p << 2;

						if (is_available(isAv, 0) && !isVisited[p + wstride1])
						queue->push(rank[q]);
						if (is_available(isAv, 1) && !isVisited[p + width])
						queue->push(rank[q + 1]);
						if (is_available(isAv, 2) && !isVisited[p + wstride0])
						queue->push(rank[q + 2]);
						if (is_available(isAv, 3) && !isVisited[p + 1])
						queue->push(rank[q + 3]);
						if (is_available(isAv, 4) && !isVisited[p - wstride1])
						queue->push(rank[q - wstride_d + 4]);
						if (is_available(isAv, 5) && !isVisited[p - width])
						queue->push(rank[q - wstride_d + 1]);
						if (is_available(isAv, 6) && !isVisited[p - wstride0])
						queue->push(rank[q - wstride_d - 2]);
						if (is_available(isAv, 7) && !isVisited[p - 1])
						queue->push(rank[q - 1]);
					}

					next_rank = queue->top();
					/*				if(p == 27579)
					{
					node[p].print(node);
					node[5773602].print(node);
					getchar();

					disp = 0;
				}
				*/
				node[p].connect_to_parent(&node_in[next_rank], next_rank);

				/*				disp = 0;
				if(p == 27579)
				{
				node[p].print(node);
				node_in[node[p].parentidx].print(node);
				getchar();
			}
			*/
			if (current_rank == next_rank)
			break;
			current_rank = next_rank;
		}

		queue->pop();
		next_rank = queue->top();

		//Redundant node removal
		if (node_in[iChild].parentidx == current_rank &&
			node_in[iChild].area == node_in[current_rank].area)
				current_rank = iChild;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank);
			if (node_in[next_rank].area == imgsize)
			{
				if (node_in[current_rank].area == imgsize)
				next_rank = current_rank;
				//				node_in[next_rank].parentidx = next_rank;
				// node_in[next_rank].parentidx = -1;

				break;
			}

			iChild = current_rank;
			current_rank = next_rank;
		}

		for (pNode = node; pNode < node + maxSize; pNode++)
		pNode->parentidx += imgsize;
		node_in[next_rank].parentidx = -1;

		delete queue;
		Free(rank);
		Free(rankinfo);
		//Free(rank2alpha);
		//		Free(levelroot);
		Free(isVisited);
		//Free(isVisited_edges);
		Free(isAvailable);
	}

}

public:
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Imgidx, Pixel> *node, *node_in;
	Imgidx num_node, num_node_in;
	Imgidx rootidx;
	Imgidx* parentAry;
	double nrmsd;

	ATree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~ATree() { if (node) { Free(node); Free(parentAry); } }
	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }

	void BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, int connectivity_in, int par)
	{
		this->height = (Imgidx)height_in;
		this->width = (Imgidx)width_in;
		this->channel = (Imgidx)channel_in;
		this->connectivity = (Imgidx)connectivity_in;
		curSize = 0;
		if (connectivity != 4 && connectivity != 8)
		{
			//std::cout << "connectivity should be 4 or 8\n" << std::endl;
			return;
		}
		if (sizeof(Pixel) == 1)
		{
			//tmp
			//std::cout << "buildalphatree - hqueue" << std::endl;
			Flood_HQueue(img);
		}
		else
		{
			if (par == 0)
			//				Flood_Trie1(img);
			Flood_HybridQueue(img, LISTSIZE_DEFAULT);
			else
			par_Flood_HybridQueue(img, LISTSIZE_DEFAULT);
		}
	}

	void AreaFilter(Pixel *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		//	Pixel val;
		AlphaNode<Imgidx, Pixel> *pNode;

		for(idx_lim = 0;idx_lim < num_node_in && node_in[idx_lim].alpha < alpha;idx_lim++)
		;

		for(i = 0;i < num_node;i++)
		{
			node[i].rootidx = 0;
			node[i].filter_val = 0;
		}

		imgsize = height * width;
		iarea = (Imgidx)(area * (double)imgsize);
		iarea = min(imgsize, max(0, iarea));
		//val = 1;
		for(i = 0;i < imgsize; i++)
		{
			pNode = &node[i];
			while(pNode->parentidx != -1 &&
				pNode->alpha < alpha
				//&& pNode->alpha < alpha
			)
			{
				//					if(val == 9260)
				{
					//						pNode->print(node);
				}
				pNode = &node[pNode->parentidx];
			}
			// 				if(pNode->filter_val == 0 && val < 200)
			//					pNode->filter_val = val++;

			//	printf("%lld\n",(uint64)val);
			// 				outimg[i] = ((double)pNode->area / (double)imgsize) * 255;
			outimg[i] = min(255,(Pixel)pNode->area/200.0);
		}
		//	printf("%lld\n",(uint64)val);
	}

	void AreaFilter(double *outimg, double area, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		Imgidx iarea;
		//	Pixel val;
		AlphaNode<Imgidx, Pixel> *pNode;

		for(idx_lim = 0;idx_lim < num_node_in && node_in[idx_lim].alpha < alpha;idx_lim++)
		;

		for(i = 0;i < num_node;i++)
		{
			node[i].rootidx = 0;
			node[i].filter_val = 0;
		}

		imgsize = height * width;
		iarea = (Imgidx)(area * (double)imgsize);
		iarea = min(imgsize, max(0, iarea));
		//val = 1;
		for(i = 0;i < imgsize; i++)
		{
			pNode = &node[i];
			while(pNode->parentidx != -1 &&
				pNode->alpha < alpha
				//&& pNode->alpha < alpha
			)
			{
				//					if(val == 9260)
				{
					//						pNode->print(node);
				}
				pNode = &node[pNode->parentidx];
			}
			// 				if(pNode->filter_val == 0 && val < 200)
			//					pNode->filter_val = val++;

			//	printf("%lld\n",(uint64)val);
			// 				outimg[i] = ((double)pNode->area / (double)imgsize) * 255;
			outimg[i] = min(255,(double)pNode->area/200.0);
		}
		//	printf("%lld\n",(uint64)val);
	}

	void print_tree()
	{
		for(int i = 0;i < maxSize; i++)
		{
			node[i].print(node);
		}
	}

};


class AlphaTree
{
	void *tree;
	int imgidx, pixel;
public:
	AlphaTree() : tree(0) {}
	~AlphaTree() { if (tree) { Free(tree);} }
	inline void clear()
	{
		/*		if(tree)
		{
		switch(imgidx)
		{
		case (IMGIDX_32BIT):
		case (IMGIDX_64BIT):
	}
}
tree->\
}
*/
}

void BuildAlphaTree(uint8 *img, int height, int width, int channel, int connectivity, int par)
{
	pixel = PIXEL_8BIT;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		//tmp
		//std::cout << "Entered wrapper build" << std::endl;

		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint8>));
		((ATree<int32,uint8>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint8>));
		((ATree<int64,uint8>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
}

void BuildAlphaTree(uint16 *img, int height, int width, int channel, int connectivity, int par)
{
	pixel = PIXEL_16BIT;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint16>));
		((ATree<int32,uint16>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint16>));
		((ATree<int64,uint16>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
}
void BuildAlphaTree(uint32 *img, int height, int width, int channel, int connectivity, int par)
{
	pixel = PIXEL_32BIT;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint32>));
		((ATree<int32,uint32>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint32>));
		((ATree<int64,uint32>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
}
void BuildAlphaTree(uint64 *img, int height, int width, int channel, int connectivity, int par)
{
	pixel = PIXEL_64BIT;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint64>));
		((ATree<int32,uint64>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint64>));
		((ATree<int64,uint64>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
	}
}


void BuildAlphaTree(float *img, int height, int width, int channel, int connectivity, int par)
{
	pixel = PIXEL_FLOAT;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint32>));
		((ATree<int32,uint32>*)tree)->BuildAlphaTree((uint32*)img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint32>));
		((ATree<int64,uint32>*)tree)->BuildAlphaTree((uint32*)img, height, width, channel, connectivity, par);
	}
}


void BuildAlphaTree(double *img, int height, int width, int channel, int connectivity, int par)
{

	//		printf("Entered BuildAlphaTree\n");
	pixel = PIXEL_DOUBLE;
	if((int64)height * (int64)width < (int64)0xefffffff)
	{
		imgidx = IMGIDX_32BITS;
		tree = Malloc(sizeof(ATree<int32,uint64>));
		((ATree<int32,uint64>*)tree)->BuildAlphaTree((uint64*)img, height, width, channel, connectivity, par);
	}
	else
	{
		imgidx = IMGIDX_64BITS;
		tree = Malloc(sizeof(ATree<int64,uint64>));
		((ATree<int64,uint64>*)tree)->BuildAlphaTree((uint64*)img, height, width, channel, connectivity, par);
	}
}

void AreaFilter(double *outimg, double area, double alpha)
{
	if(imgidx == IMGIDX_32BITS)
	{
		if(pixel == PIXEL_8BIT)
		((ATree<int32,uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
		else if(pixel == PIXEL_16BIT)
		((ATree<int32,uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
		else if(pixel == PIXEL_32BIT)
		((ATree<int32,uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
		else if(pixel == PIXEL_64BIT)
		((ATree<int32,uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
		else if(pixel == PIXEL_FLOAT)
		((ATree<int32,uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
		else
		((ATree<int32,uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);

	}
	else
	{
		if(pixel == PIXEL_8BIT)
		((ATree<int64,uint8>*)tree)->AreaFilter((uint8*)outimg, area, alpha);
		else if(pixel == PIXEL_16BIT)
		((ATree<int64,uint16>*)tree)->AreaFilter((uint16*)outimg, area, alpha);
		else if(pixel == PIXEL_32BIT)
		((ATree<int64,uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
		else if(pixel == PIXEL_64BIT)
		((ATree<int64,uint64>*)tree)->AreaFilter((uint64*)outimg, area, alpha);
		else if(pixel == PIXEL_FLOAT)
		((ATree<int64,uint32>*)tree)->AreaFilter((uint32*)outimg, area, alpha);
		else
		((ATree<int64,uint64>*)tree)->AreaFilter((double*)outimg, area, alpha);
	}

}

void print_tree()
{
	if(imgidx == IMGIDX_32BITS)
	{
		if(pixel == PIXEL_8BIT)
		((ATree<int32,uint8>*)tree)->print_tree();
		else if(pixel == PIXEL_16BIT)
		((ATree<int32,uint16>*)tree)->print_tree();
		else if(pixel == PIXEL_32BIT)
		((ATree<int32,uint32>*)tree)->print_tree();
		else if(pixel == PIXEL_64BIT)
		((ATree<int32,uint64>*)tree)->print_tree();
		else if(pixel == PIXEL_FLOAT)
		((ATree<int32,uint32>*)tree)->print_tree();
		else
		((ATree<int32,uint64>*)tree)->print_tree();
	}
	else
	{
		if(pixel == PIXEL_8BIT)
		((ATree<int64,uint8>*)tree)->print_tree();
		else if(pixel == PIXEL_16BIT)
		((ATree<int64,uint16>*)tree)->print_tree();
		else if(pixel == PIXEL_32BIT)
		((ATree<int64,uint32>*)tree)->print_tree();
		else if(pixel == PIXEL_64BIT)
		((ATree<int64,uint64>*)tree)->print_tree();
		else if(pixel == PIXEL_FLOAT)
		((ATree<int64,uint32>*)tree)->print_tree();
		else
		((ATree<int64,uint64>*)tree)->print_tree();
	}
}
/*
void BuildAlphaTree(float *img, int height, int width, int channel, int connectivity, int par)
{
pixel = PIXEL_8BIT;
if((int64)height * (int64)width < (int64)0xefffffff)
{
imgidx = IMGIDX_32BITS;
tree = Malloc(sizeof(ATree<int32,float>));
((ATree<int32,float>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
}
else
{
imgidx = IMGIDX_64BITS;
tree = Malloc(sizeof(ATree<int64,float>));
((ATree<int64,float>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
}
}
void BuildAlphaTree(double *img, int height, int width, int channel, int connectivity, int par)
{
pixel = PIXEL_8BIT;
if((int64)height * (int64)width < (int64)0xefffffff)
{
imgidx = IMGIDX_32BITS;
tree = Malloc(sizeof(ATree<int32,double>));
((ATree<int32,double>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
}
else
{
imgidx = IMGIDX_64BITS;
tree = Malloc(sizeof(ATree<int64,double>));
((ATree<int64,double>*)tree)->BuildAlphaTree(img, height, width, channel, connectivity, par);
}
}
*/
};
