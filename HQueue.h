#pragma once

#include "defines.h"
#include "allocator.h"

template <typename T>
struct HQueue
{
	T *queue;
	uint32 *bottom, *cur;
	uint64 qsize;
	uint32 min_level, max_level;
};

struct neighidx {
	//  -  3  -
	//  2  p  1
	//  -  0  -
	uint8 neighbour;
	uint32 pidx;
};

void hqueue_new(HQueue<uint32>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize);
void hqueue_new(HQueue<neighidx>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize, uint8 neighbours);

template <typename T>
void hqueue_free(HQueue<T>* hqueue)
{
	Free(hqueue->queue);
	Free(hqueue->bottom);
	Free(hqueue->cur);
	Free(hqueue);
}


inline void hqueue_push(HQueue<uint32>* hqueue, uint32 pidx, uint32 level)
{
	hqueue->min_level = min(level, hqueue->min_level);
#if DEBUG
	assert(level < hqueue->max_level);
	assert(hqueue->cur[level] < hqueue->qsize);
#endif
	hqueue->queue[hqueue->cur[level]++] = pidx;
}


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

template <typename T>
inline T hqueue_pop(HQueue<T>* hqueue)
{
	return hqueue->queue[--hqueue->cur[hqueue->min_level]];
}

template <typename T>
inline void hqueue_find_min_level(HQueue<T>* hqueue)
{
	while (hqueue->bottom[hqueue->min_level] == hqueue->cur[hqueue->min_level])
		hqueue->min_level++;
}