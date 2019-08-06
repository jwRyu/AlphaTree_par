#pragma once

#include "defines.h"
#include "allocator.h"

template <class Imgidx>
class HQueue
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
public:
	int64 qsize;
	int16 min_level;
	HQueue(uint64 qsize_in, Imgidx *dhist, int16 numlevels)
	{
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void push(Imgidx pidx, int16 level)
	{
		min_level = min(level, min_level);
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
	}

	inline Imgidx pop()
	{
		return queue[--cur[min_level]];
	}

	inline void find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};
/*

struct neighbouridx {
	//  -  3  -
	//  2  p  1
	//  -  0  -
	uint8 neighbour;
	uint32 pidx;
};
template <>
class HQueue <neighbouridx>
{
	neighbouridx *queue;
	uint64 *bottom, *cur;
	uint64 qsize;
	uint64 min_level, max_level;
public:
	HQueue(uint64 qsize, uint64 *dhist, uint32 dhistsize, uint8 neighbours)
	{
		uint64 nn = neighbours >> 1;
		int shamt;
		queue = (neighbouridx*)Malloc((size_t)qsize * nn * sizeof(neighbouridx));
		bottom = (uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(uint64));
		cur = (uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(uint64));

		this->qsize = qsize;
		min_level = max_level = dhistsize;

		for (shamt = -1; nn; nn >>= 1)
			shamt++;

		uint64 sum_hist = 0;
		for (uint64 i = 0; i < dhistsize; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i] << shamt;
		}
		bottom[dhistsize] = 0;
		cur[dhistsize] = 1;
}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void hqueue_push(uint64 pidx, uint64 level)
	{
		min_level = min(level, min_level);
#if DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
	}

	inline T hqueue_pop()
	{
		return queue[--cur[min_level]];
	}

	inline void hqueue_find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};
inline void hqueue_push(HQueue<neighidx>* hqueue, uint32 idx, uint8 neighbor, uint32 level)
{
	hqueue->min_level = min(level, hqueue->min_level);
#if DEBUG
	assert(level < hqueue->max_level);
	assert(hqueue->cur[level] < hqueue->qsize);
#endif
	hqueue->queue[hqueue->cur[level]].pidx = idx;
	hqueue->queue[hqueue->cur[level]++].neighbour = neighbor;
}
*/

/*
void hqueue_new(HQueue<uint32>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize);
void hqueue_new(HQueue<neighidx>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize, uint8 neighbours);
*/
