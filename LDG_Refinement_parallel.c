#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <omp.h>
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
	struct Bins* LDGpartitioning(int K, struct graph *G, int capacity, int threadnum) {
		int i;
		int u;
		int v;
		int b;
		int bmax;
		uint32_t r;
		int kstart;
		int kend;
		int k;
		uint32_t temp;
		struct Bins *binha;
		uint32_t *BinAff;
		int *Blist;
		int bi;
		float max;
		int bwt = 0;
		int j;
		binha = Bins(G->n, K);
		#pragma omp parallel  num_threads(threadnum) private(u,j,k,i,bi,BinAff,Blist,kstart,kend,bmax,bwt,v,b,max) shared(capacity,G,K,binha)
		{
			BinAff = (uint32_t*)calloc(K + 1, sizeof(uint32_t));
			Blist = (int*)calloc(K + 1, sizeof(int));

			double start = omp_get_wtime();
			int last_bin;
			last_bin = 1;
			#pragma omp for schedule (dynamic,1024) nowait
			for (u = 1; u <= G->n; u++) {

				bi = 0;


				kstart = G->vtx[u];

				kend = G->vtx[u + 1] - 1;

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
				max = (float)0.0;
				if (bi != 0) {
					bmax = 0;
					for (i = 1; i <= bi; i++) {
						b = Blist[i];
						bwt = binha->binwgt[b];
						if (((bwt + G->vwgt[u]) <= capacity) && (float)(BinAff[b] * (1 - bwt / capacity)) > max) {
							max = (float)BinAff[b] * (1 - bwt / capacity);
							bmax = b;
						}
						BinAff[b] = 0;
					}
				}


				if (bi == 0 || bmax == 0) {

					while ((binha->binwgt[last_bin] + G->vwgt[u] > (capacity))) {


						last_bin = ((last_bin % K) + 1);


					}


					bmax = last_bin;
					last_bin = ((last_bin % K) + 1);

				}

				binha->Bin[u] = bmax;
				binha->binwgt[bmax] = binha->binwgt[bmax] + G->vwgt[u];


			}
			double end = omp_get_wtime();
			printf("Time taken by thread %d is %f\n", omp_get_thread_num(), end - start);
			//barrier is needed to synchronize threads
			#pragma omp barrier
			double startr = omp_get_wtime();

			BinAff = (uint32_t*)calloc(K + 1, sizeof(uint32_t));
			Blist = (int*)calloc(K + 1, sizeof(int));
			//refinement partitioning starts here
			#pragma omp for schedule (dynamic,1024) nowait
			for (u = G->n; u > 0; u--)
			{

				bi = 0;
				kstart = G->vtx[u];
				kend = G->vtx[u + 1] - 1;

				for (k = kend; k >= kstart; k--) {

					v = G->adjncy[k];
					b = binha->Bin[v];
					if (BinAff[b] == 0) {
						bi++;
						Blist[bi] = b;

					}
					BinAff[b] = BinAff[b] + G->adjwgt[k];
				}
				bmax = 0;
				max = 0.0;
				for (i = 1; i <= bi; i++) {
					b = Blist[i];
					//printf("b : %d \n",b);
					bwt = binha->binwgt[b];
					if (((bwt + G->vwgt[u]) <= capacity) && (float)(BinAff[b] * (1 - bwt / capacity)) > max) {
						max = (float)BinAff[b] * (1 - bwt / capacity);
						bmax = b;


					}
					BinAff[b] = 0;
				}
				if (bmax != 0 && bmax != binha->Bin[u]) {

					binha->binwgt[binha->Bin[u]] = binha->binwgt[binha->Bin[u]] - G->vwgt[u];
					binha->binwgt[bmax] = binha->binwgt[bmax] + G->vwgt[u];

				}
				if (bmax == 0) {
					bmax = binha->Bin[u];
				}

				binha->Bin[u] = bmax;


			}
			double endr = omp_get_wtime();
			printf("refinedTime taken by thread %d is %f\n", omp_get_thread_num(), endr - startr);

		}


		return binha;

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
		float edgecutsum = (float)0.0;
		float second_edgecutsum = (float)0.0;
		float second_refcutsum = (float)0.0;
		float imbalancesum = 0.0;
		char *fname = argv[1];
		int K = atoi(argv[2]);
		int threadnum = atoi(argv[3]);
		struct graph* G = readgraphfromfile(fname);

		float eps;
		eps = 0.1;
		float beforeedgecutsum = (float)0.0;
		//five time runs
		for (rpt = 0; rpt < 5; rpt++) {

			//struct Bins* bin=Bins(G->n,K);
			uint32_t cutcount = 0;
			int capacity = 0;
			int wtotal = 0;
			float wavg = 0.0;

			for (j = 1; j < G->n + 1; j++) {
				wtotal += G->vwgt[j];

			}

			wavg = (float)(wtotal / K);


			capacity = ceil(wavg * (1 + eps));


			clock_t start, end;
			double cpu_time_used;
			struct Bins* initbins = LDGpartitioning(K, G, capacity, threadnum);
			uint32_t adjwgttotal = 0;
			for (u = 1; u <= G->n; u++) {
				//for each v in adj(u)
				kstart = G->vtx[u];
				kend = G->vtx[u + 1] - 1;
				for (k = kstart; k <= kend; k++) {
					v = G->adjncy[k];
					if (initbins->Bin[u] != initbins->Bin[v]) {

						cutcount += G->adjwgt[k];

					}
				}
			}

			for (l = 1; l < G->vtx[G->n + 1]; l++) {
				adjwgttotal += G->adjwgt[l];

			}

			cutcount = cutcount / 2;

			float beforeedgecut = 0.0;
			beforeedgecut = (float)cutcount / (adjwgttotal / 2);

			beforeedgecutsum += beforeedgecut;



			float imbalance;
			uint32_t maximb = 0;
			for (i = 0; i < K + 1; i++) {
				if (initbins->binwgt[i] > maximb) {
					maximb = initbins->binwgt[i];
				}

			}

			imbalance = (float)maximb / wavg;
			imbalancesum += imbalance;
			printf("imbalance : %f\n", imbalance);

		}

		float avgimbalance = imbalancesum / 5;
		float avgedgecutbefore = beforeedgecutsum / 5;
		printf("beforere finement edgecut: %f \n", avgedgecutbefore);
		printf("avg imbalance : %f\n", avgimbalance);



		return 0;
	}
