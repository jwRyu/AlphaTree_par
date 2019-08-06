/*
void hqueue_new(HQueue<uint32>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize)
{
	uint32 i;
	(*hqueue) = (HQueue<uint32>*)Malloc(sizeof(HQueue<uint32>));
	(*hqueue)->queue = (uint32*)Malloc((size_t)qsize * sizeof(uint32));
	(*hqueue)->bottom = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));
	(*hqueue)->cur = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));

	(*hqueue)->qsize = qsize;
	(*hqueue)->min_level = (*hqueue)->max_level = dhistsize;

	uint32 sum_hist = 0;
	for (i = 0; i < dhistsize; i++)
	{
		(*hqueue)->bottom[i] = (*hqueue)->cur[i] = sum_hist;
		sum_hist += dhist[i];
	}
	(*hqueue)->bottom[dhistsize] = 0;
	(*hqueue)->cur[dhistsize] = 1;
}

void hqueue_new(HQueue<neighidx>** hqueue, uint64 qsize, uint32 *dhist, uint32 dhistsize, uint8 neighbours)
{
	uint32 i, nn = neighbours >> 1;
	int shamt;
	(*hqueue) = (HQueue<neighidx>*)Malloc(sizeof(HQueue<neighidx>));
	(*hqueue)->queue = (neighidx*)Malloc((size_t)qsize * nn * sizeof(neighidx));
	(*hqueue)->bottom = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));
	(*hqueue)->cur = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));

	(*hqueue)->qsize = qsize;
	(*hqueue)->min_level = (*hqueue)->max_level = dhistsize;

	for (shamt = -1; nn; nn >>= 1)
		shamt++;

	uint32 sum_hist = 0;
	for (i = 0; i < dhistsize; i++)
	{
		(*hqueue)->bottom[i] = (*hqueue)->cur[i] = sum_hist;
		sum_hist += dhist[i] << shamt;
	}
	(*hqueue)->bottom[dhistsize] = 0;
	(*hqueue)->cur[dhistsize] = 1;
}

HQueue* hqueue_new(uint64 qsize, uint32 *dhist, uint32 dhistsize)
{
	uint32 i;
	HQueue* hqueue = (HQueue*)Malloc(sizeof(HQueue));
	hqueue->queue = (uint32*)Malloc((size_t)qsize * sizeof(uint32));
	hqueue->bottom = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));
	hqueue->cur = (uint32*)Malloc((size_t)(dhistsize + 1) * sizeof(uint32));

	hqueue->qsize = qsize;
	hqueue->min_level = hqueue->max_level = dhistsize;

	int sum_hist = 0;
	for (i = 0; i < dhistsize; i++)
	{
		hqueue->bottom[i] = hqueue->cur[i] = sum_hist;
		sum_hist += dhist[i];
	}
	hqueue->bottom[dhistsize] = 0;
	hqueue->cur[dhistsize] = 1;

	return hqueue;
}

void hqueue_free(HQueue* hqueue)
{
	Free(hqueue->queue);
	Free(hqueue->bottom);
	Free(hqueue->cur);
	Free(hqueue);
}
*/
