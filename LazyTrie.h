#pragma once
#include "allocator.h"
// 
// //tmptmptmptmpt tmp
// #include <fstream>
// using namespace std;

#define LAZYTRIE_DEBUG 0


template<class Imgidx, class Trieidx>
class LazyTrie
{
	Imgidx minidx;
	Trieidx **trie, *trie0;
	//Imgidx *levelsize;
	Imgidx popdist, popdist_thr;
	Imgidx triesize, *levelsize, mask_field, nbits;
	int8 shamt, numlevels;
	//delayed non-leaf node push

#if TRIE_DEBUG
	Imgidx cursize;
#endif


	// 	//tmptmptmp
	// 	ofstream f;


		//int32 curSize;//tmp

	void initLazyTrie(Imgidx triesize)
	{

#if TRIE_DEBUG
		cursize = 0;
#endif
		//curSize = 0;//tmp
		Imgidx size;
		shamt = 2;
		for (int8 nbyte = sizeof(Trieidx); nbyte; nbyte >>= 1)
			shamt++;
		nbits = 1 << shamt;
		mask_field = nbits - 1;
		numlevels = 1;
		for (size = (triesize + 1) >> shamt; size; size >>= shamt)
			numlevels++;

		trie = (Trieidx**)Malloc(sizeof(Trieidx*) * numlevels);
		levelsize = (Imgidx*)Malloc(sizeof(Imgidx) * numlevels);
		size = triesize + 1;
		for (int8 i = 0; i < numlevels; i++)
		{
			size >>= shamt;
			levelsize[i] = size + 1;
			trie[i] = (Trieidx*)Malloc(sizeof(Trieidx) * (size + 1));
			for (Imgidx j = 0; j < (size + 1); j++)
				trie[i][j] = 0;
			//levelsize[i] = lvlsz;
		}
		trie0 = trie[0];
		//push(triesize);
		//curSize = 1;
		//curSize = 0;
		//tmp!
// 		f.open("C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/trie0rrr.dat",std::ofstream::out);
// 		f << triesize << endl;
	}
public:
	LazyTrie(Imgidx triesize, Imgidx popdist_thr)
		: popdist(0), popdist_thr(popdist_thr), minidx(triesize)
	{
		initLazyTrie(triesize);
	}
	LazyTrie(Imgidx triesize)
		: popdist(0), minidx(triesize)
	{
		popdist_thr = 1200;
		initLazyTrie(triesize);
	}
	~LazyTrie()
	{
		for (int8 i = 0; i < numlevels; i++)
			Free(trie[i]);
		Free(trie);
		Free(levelsize);

		// 		f.close();//tmptmp
	}

	inline Imgidx top() { return minidx; }
	inline Imgidx min_rank() { return minidx >> 1; }
	inline Imgidx min_incidence() { return minidx & 1; }
	inline void push(Imgidx in, int8 incidence)
	{
		//curSize++; //tmp
		//tmp
/*		f << '0' << '\n' << in << endl;*/
		push((in << 1) + incidence);
	}
	inline void push(Imgidx in)
	{
		Imgidx n = in, s_in, shamt1;
		Trieidx *p;

		//		curSize++; //tmp
				//tmp
		/*		f << '0' << '\n' << in << endl;*/

				//n = (in << 1) + incidence;
		s_in = n >> shamt;

		if (n < minidx)
			minidx = n;

		p = trie0 + s_in;
		*p |= ((Trieidx)1 << (n & mask_field));

		//lazytrie
		popdist >>= 1;


// 		n = s_in;
// 		s_in >>= shamt;
// 		for (int8 i = 1; i < numlevels; i++)
// 		{
// 			p = trie[i] + s_in;
// 			shamt1 = n & mask_field;
// 			//if (((*p) >> shamt1) & 1)
// 				//break;
// 			*p |= ((Trieidx)1 << shamt1);
// 			n = s_in;
// 			s_in >>= shamt;
// 		}

#if TRIE_DEBUG
		cursize++;
#endif
	}

	inline void buildtrie()
	{
		Imgidx size, lvlsz, lvlsz1;
		Trieidx *trie_i, *trie_i1, *p, *q, *r, s, *t;
		int8 i, k;
		Imgidx j, l, m;

		for (int8 i = 0; i < numlevels - 1; i++)
		{
			m = 0;
			for (j = 0; j < levelsize[i + 1] - 1; j++)
			{
				trie[i + 1][j] = 0;
				p = &trie[i + 1][j];
				for (l = 0; l < nbits; l++)
					*p |= ((Trieidx)(trie[i][m++] != 0) << l);
			}
			p = &trie[i + 1][j];
			*p = 0;
			for (l = 0; m < levelsize[i]; l++)
				*p |= ((Trieidx)(trie[i][m++] != 0) << l);

// 			for (j = 0; j < levelsize[i]; j++)
// 			{
// 				if (trie[i][j])
// 				{
// 					p = trie[i + 1] + (j >> shamt);
// 		//			s = *p;
// 					*p |= ((Trieidx)1 << (j & mask_field));
// //					if (s != *p)
// 	//					std::cout << "huphup" << std::endl;
// 				}
// 			}
		}

		return;

// 		trie_i1 = trie0;
// 		lvlsz1 = levelsize[0];
// 		for (int8 i = 0; i < numlevels - 1; i = k)
// 		{
// 			k = i + 1;
// 			trie_i = trie_i1;
// 			trie_i1 = trie[k];
// 			lvlsz = lvlsz1;
// 			lvlsz1 = levelsize[k];
// 			q = trie_i - 1;
// 			j = (lvlsz & (~mask_field)) - 1;
// 			p = trie_i + j;
// 			t = trie[k] + (j >> shamt);
// 			while (p != q)
// 			{
// 				s = *p-- != 0;
// 				for (r = p - nbits; p != r; p--)
// 					s = (s << 1) | (*p != 0);
// 				if(*t != s)//tmp
// 					*t-- = s;
// 			}
// 			q = trie_i + j;
// 			p = trie_i + lvlsz - 1;
// 			t = trie_i1 + lvlsz1 - 1;
// 			s = *p-- != 0;
// 			while (p != q)
// 				s = (s << 1) | (*p-- != 0);
// 			if (*t != s)//tmp
// 				*t-- = s;
// 	//		*t = s;
// 		}
	}

	inline void pop()
	{
		Imgidx s_idx = minidx >> shamt, shamt1;
		Trieidx *p, tmp;
		int8 lvl;

		//curSize--;//tmp
// 		//tmp
// 		f << '1' << '\n' << (minidx>>1) << endl;

		if (1)//(popdist > popdist_thr)
		{
			buildtrie();
			shamt1 = minidx & mask_field;
			p = trie0 + s_idx;
			*p &= (~((Trieidx)1 << shamt1++));
			for (lvl = 0; !*p;)
			{
				minidx = s_idx;
				s_idx >>= shamt;
				shamt1 = minidx & mask_field;
				p = trie[++lvl] + s_idx;
				*p &= (~((Trieidx)1 << shamt1));
			}
			for (tmp = *p >> shamt1; !(tmp & 1); tmp >>= 1)
				shamt1++;
			minidx = (minidx & ~mask_field) | shamt1;
			while (lvl)
			{
				tmp = trie[--lvl][minidx];
				p = &tmp;
				s_idx = minidx;
				minidx <<= shamt;
				for (shamt1 = 0; !(tmp & 1); shamt1++)
					tmp >>= 1;
				minidx |= shamt1;
			}
		}
		else
		{

		}
#if TRIE_DEBUG
		cursize--;
#endif
	}
};

