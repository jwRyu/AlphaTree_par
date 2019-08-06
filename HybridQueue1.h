#pragma once
#include "allocator.h"
#include "defines.h"
#include "LazyTrie.h"

//#include <iostream>

//#define LISTSIZE_DEFAULT 12

#define TRACK_QUEUEING	0

//tmp
#if TRACK_QUEUEING
#include <fstream>
using namespace std;
#endif

template<class Imgidx>
struct MinList1
{
	Imgidx idx;
	MinList1 *next;
};

template<class Imgidx, class Trieidx>//, class Qidx>
class HybridQueue1
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	LazyTrie<Imgidx, Trieidx> *trie;
	int8 *queue;
	Imgidx minidx_queue;
	int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	int8 shamt, nbit;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		// 
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));	
		//queue = (int8*)Malloc((size + 1) * sizeof(int8));
		//trie = (Trie<Imgidx, int64>*)Malloc(size * sizeof(Trie<Imgidx, int64>*));
		trie = new LazyTrie<Imgidx, int64>(size);
		list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
		list[0] = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/Hybrid_queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue1(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue1(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}
	inline Imgidx top() { return list[0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		int16 i;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		// 
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (idx < trie->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list]);

				for (i = curSize_list - 1; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		trie->push(idx);
	}
	inline void pop()
	{
		int8 i;
		Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << head->idx << endl;
#endif

		if (curSize_list == 0)
		{
			list[0] = trie->top();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << head->idx << endl;
#endif
	}
	inline void pop_queue()
	{
		trie->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	// 	int8 checklist()//tmp
	// 	{
	// 		MinList<Imgidx> *p;
	// 		if (head)
	// 		{
	// 			for (p = head; p; p = p->next)
	// 			{
	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
	// 					return 1;
	// 			}
	// 			if (tail->next)
	// 				return 1;
	// 		}
	// 		int n = 0;
	// 		for (int i = 0; i < maxSize_list; i++)
	// 		{
	// 			if (list[i].idx != -1)
	// 				n++;
	// 		}
	// 		if (n != curSize_list)
	// 			return 1;
	// 		return 0;
	// 	}

	~HybridQueue1()
	{
		delete trie;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};
