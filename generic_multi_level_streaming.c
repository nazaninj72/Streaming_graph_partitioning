#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
/*------ struct definitions -----*/
struct graph {
	int n, m, fmt;
	uint64_t* vtx;
	uint32_t* adjncy;
	uint32_t* vwgt;
	uint32_t* adjwgt;
	int id;
};

struct Bins {
	int *Bin;
	uint32_t *binwgt;
};
struct Bins* Bins(int s, int k) {
	struct Bins* bins = malloc(sizeof(struct Bins));
	bins->Bin = (int*)calloc(s + 1, sizeof(int));
	bins->binwgt = (uint32_t*)calloc(k + 1, sizeof(uint32_t));
	return bins;
}
struct graph* graph_coarse(int K) {
	struct graph* Gprime = malloc(sizeof(struct graph));
	Gprime->n = K;
	Gprime->vtx = (uint64_t*)calloc(K + 2, sizeof(uint64_t));
	Gprime->vwgt = malloc(sizeof(uint32_t) * (K + 2));

	return Gprime;
}
struct graph* graph_new(int n, int m, int fmt) {
	struct graph* G = malloc(sizeof(struct graph));
	G->n = n;
	G->m = m;
	G->vtx = (uint64_t*)calloc(n + 2, sizeof(uint64_t));
	G->adjncy = malloc(sizeof(uint32_t) * (2 * m + 2));
	G->adjwgt = malloc(sizeof(uint32_t) * (2 * m + 2));
	G->vwgt = (uint32_t*)calloc(n + 2, sizeof(uint32_t));
	G->fmt = fmt;
	return G;
}


/*
 * Function:  showgraph
 * --------------------
 * displays the graph information
 *
 *
 *
 *  graph*G: graph
 *
 *
 *
 */
void showgraph(struct graph *G) {
	int i = 0;
	printf("\n \n The VTX array is : \n");
	for (i = 0; i < G->n + 2; i++)
	{
		printf("vtx[%d] = %u ", i, G->vtx[i]);
	}
	printf("\n \n The vwgt array is : \n");
	for (i = 0; i < G->n + 1; i++)
	{
		printf("vwgt[%d] = %u ", i, G->vwgt[i]);
	}
	printf("\n \n The adjncy array is : \n");
	printf("G->vtx[G->n+1] : %d \n ", G->vtx[G->n + 1]);
	for (i = 0; i < G->vtx[G->n + 1]; i++)
	{
		printf("adjncy[%d] = %u ", i, G->adjncy[i]);
	}
	printf("\n \n The adjwgt array is : \n");
	for (i = 0; i < G->vtx[G->n + 1]; i++)
	{
		printf("adjwgt[%d] = %u ", i, G->adjwgt[i]);
	}
}


/*
 * Function:  readgraphfromfile
 * --------------------
 * read graph information which is in csr format into the graph structure
 *
 *
 * fname[]: name of the input file
 *  returns: graph information
 *
 *
 */
struct graph* readgraphfromfile(char fname[])
{
	FILE *ifp;
	char *mode = "r";
	int n;
	int m;
	int fmt;
	char ch = 'a';
	ifp = fopen(fname, "r");
	if (ifp == NULL)
	{
		printf("\nFailed to open file.\n");
		exit(1);
	}

	fscanf(ifp, "%d %d %d%c", &n, &m, &fmt, &ch);

	struct graph *G = graph_new(n, m, fmt);

	char *token = malloc(1024);
	uint64_t Vindx = 1;
	uint32_t tindx = 0;
	uint64_t edgeindx = 1;

	G->vtx[0] = 0;
	G->vwgt[0] = 0;
	G->adjncy[0] = 0;
	G->adjwgt[0] = 0;
	G->n = n;

	int i;
	int j;
	for (j = 1; j <= G->n + 1; j++) {
		G->vwgt[j] = 1;
	}
	while (Vindx <= G->n) {

		fscanf(ifp, "%c", &ch);
		if (ch == '\n')
		{
			G->vtx[Vindx++] = edgeindx;

		}
		else
		{
			fseek(ifp, -sizeof(char), SEEK_CUR);
			fscanf(ifp, "%s%c", token, &ch);
			uint32_t val = atoi(token);

			if (tindx % 2 == 0)
			{

				if (tindx == 0)//it means that this is my first value in the line
				{
					G->vtx[Vindx] = edgeindx;
				}
				G->adjncy[edgeindx] = val;
			}
			else
			{
				G->adjwgt[edgeindx++] = val;
			}
			if (ch == '\n') {

				tindx = 0;

				Vindx++;

				G->vtx[Vindx] = edgeindx;

			}
			else
				tindx++;

		}
	}

	G->vtx[Vindx] = edgeindx;
	fclose(ifp);
	return G;
}

/*
 * Function:  LDGpartitioning
 * --------------------
 * partitions the graph using linear deterministic algorithm, it assigns vertex v to partition i such that it will have most number of neighbors
 *
 *
 *  K: number of partitions
 *  graph*G : graph that needs to be partitioned
 *
 *	capacity: capacity of the bins
 *  returns: bin information of the graph
 *
 *
 */
struct Bins* LDGpartitioning(int K, struct graph *G, int capacity)
{
	float edgecutsum = 0.0;
	float imbalancesum = 0.0;
	int i;
	int j;
	int u;
	uint32_t v;
	int b;
	int bmax;
	uint32_t r;
	int kstart;
	int kend;
	int k;
	int cutcount = 0;
	int temp;
	int maxwgt = 0;
	double epsilon = 0.000001;
	//initialize bin ,bin affinity and bin weight to zero
	struct Bins *binha = Bins(G->n, K);
	uint32_t *BinAff;
	BinAff = (uint32_t*)calloc(K + 1, sizeof(uint32_t));
	int *Blist;
	Blist = (uint32_t*)calloc(K + 1, sizeof(uint32_t));

	int bi;
	//random number generation
	uint64_t* rvertices;
	rvertices = malloc(sizeof(uint64_t) * (G->n + 1));
	for (i = 0; i < G->n; i++) {
		rvertices[i] = i + 1;
	}
	for (i = 0; i < G->n - 1; i++) {
		r = (rand() % G->n);
		temp = rvertices[r];
		rvertices[r] = rvertices[i];
		rvertices[i] = temp;

	}

	for (i = 0; i <= G->n; i++) {
		if (G->vwgt[i] > maxwgt) {
			maxwgt = G->vwgt[i];
		}
	}
	//for each u in graph find its corresponding bin
	for (u = 0; u < G->n; u++) {
		bi = 0;
		//for each v in adj(u) find its affinity to u
		kstart = G->vtx[rvertices[u]];
		kend = G->vtx[rvertices[u] + 1] - 1;
		for (k = kstart; k <= kend; k++) {
			v = G->adjncy[k];

			b = binha->Bin[v];
			if (b == 0) {
				continue;
			}
			else if (BinAff[b] == 0) {
				bi++;
				Blist[bi] = b;

			}
			BinAff[b] = BinAff[b] + G->adjwgt[k];

		}

		bmax = 0;
		float max = 0.0;
		//no neighbor of vertex u is assigned to any bin
		if (bi != 0) {
			bmax = 0;
			for (i = 1; i <= bi; i++) {
				b = Blist[i];
				if (((binha->binwgt[b] + G->vwgt[rvertices[u]]) <= capacity) && (float)(BinAff[b] * (1 - binha->binwgt[b] / capacity)) > max) {
					max = (float)BinAff[b] * (1 - binha->binwgt[b] / capacity);
					bmax = b;

				}
				BinAff[b] = 0;
			}
		}
		//when the vertex did not assign to the bin of its neighbor because of high capacity we assign it randomly
		int minimum = capacity;
		int minindez = -1;
		if (bi == 0 || bmax == 0) {

			//if whole vertices will fit in the bins with given capacity
			int inbins = (G->n - ((G->n / K) * K));
			int left = inbins / K;
			if (left < 1) {
				left = 1;
			}

			if (((left + (G->n / K)) * maxwgt) <= capacity) {

				while (binha->Bin[rvertices[u]] == 0) {
					bmax = (rand() % K) + 1;
					if (binha->binwgt[bmax] + G->vwgt[rvertices[u]] <= (capacity))
					{
						binha->Bin[rvertices[u]] = bmax;
						continue;
					}

				}
			}
			else
			{
				while (binha->Bin[rvertices[u]] == 0) {
					bmax = (rand() % K) + 1;
					if (binha->binwgt[bmax] + G->vwgt[rvertices[u]] <= (capacity))
					{
						binha->Bin[rvertices[u]] = bmax;
						continue;
					}
					else {
						for (i = 1; i <= K; i++) {
							if (binha->binwgt[i] < minimum) {
								minimum = binha->binwgt[i];
								minindez = i;
							}
						}
						bmax = minindez;
						binha->Bin[rvertices[u]] = bmax;

						continue;
					}

				}

			}

		}
		binha->Bin[rvertices[u]] = bmax;

		binha->binwgt[bmax] = binha->binwgt[bmax] + G->vwgt[rvertices[u]];

	}


	free(Blist);
	free(BinAff);
	free(rvertices);
	return binha;

}

/*
 * Function:  coarsening
 * --------------------
 *takes bin information of previous level as a parameter and create a new coarse graph using this bin information
 * such that previous bin becomes new graph with vertice weight equal to bin weight of the previous level and edges between vertices are edgecuts between bins
 *
 *
 *  K: number of partitions
 *  graph*G: previous level graph
 *	Bin: array of the bin information of each vertex
 *	binwgt: array for the bin weight
 *  returns: coarser graph
 *
 *
 */
struct graph* coarsening(int *Bin, int K, uint32_t *binwgt, struct graph *G) {
	uint32_t *bvix;
	bvix = (uint32_t*)calloc(K + 2, sizeof(uint32_t));
	uint32_t *BV;
	int v;
	int u;
	struct graph *Gprime = graph_coarse(K);
	int *Bmap;
	Bmap = (int*)calloc(K + 1, sizeof(int));
	BV = malloc(sizeof(uint32_t) * (G->n + 1));
	int b = 0;
	Gprime->vwgt[0] = 0;
	for (b = 1; b <= K; b++) {
		Gprime->vwgt[b] = binwgt[b];
	}

	for (v = 1; v <= G->n; v++) {
		int t;
		t = Bin[v];
		bvix[t] += 1;
	}
	//prefix sum for bvix
	for (b = 2; b <= K; b++) {
		bvix[b] += bvix[b - 1];
	}
	bvix[K + 1] = bvix[K];

	for (v = 1; v <= G->n; v++) {
		b = Bin[v];
		int t;
		t = bvix[b] - 1;
		BV[t] = v;
		bvix[b] -= 1;
	}

	for (b = 1; b <= K; b++) {
		int kstart;
		int kadjst;
		int kadjend;
		int kadj;
		int kend;
		kstart = bvix[b];
		kend = bvix[b + 1] - 1;
		int kix;
		for (kix = kstart; kix <= kend; kix++) {
			u = BV[kix];
			kadjst = G->vtx[u];
			kadjend = G->vtx[u + 1] - 1;

			for (kadj = kadjst; kadj <= kadjend; kadj++) {
				v = G->adjncy[kadj];
				int bt;
				bt = Bin[v];
				if (bt != b && (Bmap[bt] == 0 || Bmap[bt] != b)) {
					Gprime->vtx[b] += 1;
					Bmap[bt] = b;
				}
			}

		}
	}
	Bmap[1] = 0;
	for (b = 2; b <= K; b++) {
		Gprime->vtx[b] += Gprime->vtx[b - 1];
		Bmap[b] = 0;
	}

	Gprime->vtx[K + 1] = Gprime->vtx[K];

	int cadjsize = Gprime->vtx[K + 1];
	Gprime->m = cadjsize;
	Gprime->adjncy = malloc(sizeof(uint32_t) * cadjsize + 1);
	Gprime->adjwgt = malloc(sizeof(uint32_t) * cadjsize + 1);
	for (b = K; b >= 1; b--) {

		int kstart;
		int kend;
		int kadjst;
		int kadjend;
		int kadj;
		kstart = bvix[b + 1];
		kend = bvix[b] + 1;
		int kix;
		for (kix = kstart; kix >= kend; kix--) {
			u = BV[kix - 1];

			kadjst = G->vtx[u];
			kadjend = G->vtx[u + 1] - 1;
			for (kadj = kadjst; kadj <= kadjend; kadj++) {
				v = G->adjncy[kadj];

				int bt = Bin[v];

				if (bt != b && Bmap[bt] == 0) {
					int t = Gprime->vtx[b] - 1;
					Gprime->adjncy[t] = bt;
					Bmap[bt] = G->adjwgt[kadj];
					Gprime->vtx[b] -= 1;

				}
				else if (bt != b) {
					Bmap[bt] += G->adjwgt[kadj];

				}

			}
		}
		int bst = Gprime->vtx[b];
		int bend = Gprime->vtx[b + 1] - 1;
		int bx;
		for (bx = bst; bx <= bend; bx++) {

			Gprime->adjwgt[bx] = Bmap[Gprime->adjncy[bx]];
			;
			Bmap[Gprime->adjncy[bx]] = 0;
		}

	}

	free(BV);
	free(bvix);
	free(Bmap);
	return Gprime;
}


/*
 * Function:  refLDG
 * --------------------
 * makes another pass from partitioned graph and tries to find the best location for assigned vertices. 
 * if the previously assigned vertice finds a bin that gives higher score in terms of adjacent vertices,
*  it will move to that partition or bin if the bin capacity of the destination bin is low enough to allow this transition.
 *
 *
 *  K: number of partitions
 *  graph*G: graph that needs to be partitioned
 *	Bins*bin: information of the previous partition that needs to be changed
 *	capacity: capacity of the bins
 *  returns: bin information
 *
 *
 */
struct Bins* refLDG(int K, struct graph *G, struct Bins *bin, int capacity) {
	//printf("capacity : %d \n",capacity);
	int rpt;
	float edgecutsum = 0.0;
	float imbalancesum = 0.0;
	int i;
	int j;
	int u;
	int v;
	int b;
	//srand (time(NULL));
	int bmax;
	uint32_t r;
	int kstart;
	int kend;
	int k;
	int temp;
	int cutcount = 0;
	//initialize bin ,bin affinity and bin weight to zero

	uint32_t *BinAff;
	BinAff = (uint32_t*)calloc(K + 1, sizeof(uint32_t));

	int *Blist;
	Blist = (uint32_t*)calloc(K + 1, sizeof(uint32_t));

	//for each u in V  in some order do
	int bi;
	//for random number generation
	uint32_t* rvertices;
	rvertices = malloc(sizeof(uint32_t) * (G->n));
	for (i = 0; i < G->n; i++) {
		rvertices[i] = i + 1;
	}
	for (i = 0; i < G->n - 1; i++) {
		r = (rand() % G->n);
		temp = rvertices[r];
		rvertices[r] = rvertices[i];
		rvertices[i] = temp;

	}
	int flipped = 0;

	for (u = 0; u < G->n; u++) {

		bi = 0;
		kstart = G->vtx[rvertices[u]];
		kend = G->vtx[rvertices[u] + 1] - 1;
		for (k = kstart; k <= kend; k++) {

			v = G->adjncy[k];

			b = bin->Bin[v];

			if (BinAff[b] == 0) {
				bi++;
				Blist[bi] = b;

			}
			BinAff[b] = BinAff[b] + G->adjwgt[k];

		}
		//no neighbor of vertex u is assigned to any bin
		bmax = 0;
		float max = 0.0;
		for (i = 1; i <= bi; i++) {
			b = Blist[i];

			if (((bin->binwgt[b] + G->vwgt[rvertices[u]]) <= capacity) && (float)(BinAff[b] * (1 - bin->binwgt[b] / capacity)) > max) {
				max = (float)BinAff[b] * (1 - bin->binwgt[b] / capacity);
				bmax = b;


			}
			BinAff[b] = 0;
		}

		//when the vertex did not assign to the bin of its neighbor because of high capacity we assign it randomly

		if (bmax != 0 && bmax != bin->Bin[rvertices[u]]) {

			bin->binwgt[bin->Bin[rvertices[u]]] = bin->binwgt[bin->Bin[rvertices[u]]] - G->vwgt[rvertices[u]];
			bin->binwgt[bmax] = bin->binwgt[bmax] + G->vwgt[rvertices[u]];

		}
		if (bmax == 0) {
			bmax = bin->Bin[rvertices[u]];
		}
		bin->Bin[rvertices[u]] = bmax;



	}


	free(Blist);
	free(BinAff);
	free(rvertices);
	return bin;

}

/*
 * Function:  uncoarsen
 * --------------------
 * uncoarsens the graph into its original version such that it gets bin information of current and previous bin assignment and it 
 * changes bin information of current bin in terms of previous bin. 
 * initbins->Bin[1]=2 and initbins->Bin[2]=2
 * binprime->Bin[2]=3 
 * uncoarsen will return binprime->Bin[1]=3; and binprime->Bin[2]=3
 *
 *
 *  K: number of partitions
 *  graph*G: graph information
 *	Bins*binprime: information of the latest partition
 *	Bins*initbins: information of the pervious partition
 *  returns: uncoarsened bin value of the latest partition
 *
 *
 */
struct Bins* uncoarsen(struct Bins *binprime, struct Bins *initbins, struct graph* G, int K) {
	int i = 0;
	struct Bins *bin = Bins(G->n, K);
	for (i = 1; i <= G->n; i++) {
		bin->Bin[i] = binprime->Bin[initbins->Bin[i]];
	}
	bin->binwgt[0] = 0;
	for (i = 1; i <= K; i++) {
		bin->binwgt[i] = binprime->binwgt[i];

	}
	return bin;

}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int i;
	int j;
	int u;
	int rpt;
	uint64_t l;
	int v;
	int kstart;
	int kend;
	int k;
	float edgecutsum = 0.0;
	float imbalancesum = 0.0;
	char *fname = argv[1];
	int K = atoi(argv[2]);
	int lvl = 0;
	int beta = atoi(argv[3]);
	int setting = 1;

	struct graph* G = readgraphfromfile(fname);

	// running the method 5 times
	for (rpt = 0; rpt < 5; rpt++) {
		float eps = (float)0.0;
		uint32_t adjwgttotal = 0;
		for (l = 1; l < G->vtx[G->n + 1]; l++) {
			adjwgttotal += G->adjwgt[l];

		}

		int id = 0;
		G->id = id;
		struct graph* Gprime = G;

		int wtotal = 0;
		double wavg = 0.0;
		wtotal = G->n;
		eps = 0.300000;
		float temp = (float)eps - 0.100000;
		double refcapacity = 0.0;
		struct Bins* partitions;
		int capacity = 0;
		int i = 1;
		int ki = 0;
		//finding beta value 
		while (G->n / pow(beta, i) > K * beta) {

			i++;
		}
		lvl = i + 1;
		struct Bins* binarray[lvl + 2];
		struct graph* grapharray[lvl + 2];
		int karray[lvl + 2];
		karray[0] = 0;
		for (i = 1; i < lvl; i++) {
			karray[i] = G->n / pow(beta, i);
		}
		karray[lvl] = K;
		for (i = 1; i <= lvl; i++) {

			float range = (float)0.0;
			if (lvl > 1)
				range = temp / (lvl - 1);
			else
				range = 0;
			wavg = ((double)wtotal / karray[i]);
			capacity = ceil((wavg * (1 + eps)));
			grapharray[i] = Gprime;

			partitions = LDGpartitioning(karray[i], grapharray[i], capacity);

			if (setting == 1) {
				partitions = refLDG(karray[i], grapharray[i], partitions, capacity);
			}
			id++;
			binarray[i] = partitions;

			Gprime = coarsening(partitions->Bin, karray[i], partitions->binwgt, grapharray[i]);
			Gprime->m = Gprime->vtx[K + 1];
			Gprime->id = id;
			eps -= range;

		}



		uint32_t maximb = 0;
		double imbalance;


		wavg = ((double)wtotal / K);

		refcapacity = ceil(wavg*(1.1));
		int templvl = lvl;

		while (templvl >= 2) {
			binarray[templvl - 1] = uncoarsen(binarray[templvl], binarray[templvl - 1], grapharray[templvl - 1], K);

			binarray[templvl - 1] = refLDG(K, grapharray[templvl - 1], binarray[templvl - 1], refcapacity);
			templvl--;
		}

		int cutcount = 0;
		for (u = 1; u <= G->n; u++) {
			//for each v in adj(u)
			kstart = G->vtx[u];
			kend = G->vtx[u + 1] - 1;
			for (k = kstart; k <= kend; k++) {
				v = G->adjncy[k];
				if (binarray[1]->Bin[u] != binarray[1]->Bin[v]) {
					cutcount += G->adjwgt[k];
				}
			}
		}

		cutcount = cutcount / 2;

		float edgecut = 0.0;
		edgecut = (float)cutcount / (adjwgttotal / 2);
		edgecutsum += edgecut;
		imbalance = 0;
		maximb = 0;
		for (i = 0; i < K + 1; i++) {
			if (binarray[1]->binwgt[i] > maximb) {
				maximb = binarray[1]->binwgt[i];
			}

		}

		imbalance = (float)maximb / wavg;

		imbalancesum += imbalance;


	}//end of five run


	printf("number of levels : %d \n", lvl);
	printf("avg imbalance : %f\n", imbalancesum / 5);
	printf("avg edgecut of all : %f\n", edgecutsum / 5);
	return 0;
}

