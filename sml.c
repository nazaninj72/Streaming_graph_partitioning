#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <malloc.h>
#include <assert.h>
#include <inttypes.h>
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
    uint32_t bigm=2 *m + 2;
    printf("%" PRIu32"\n",bigm);
	G->vtx = (uint64_t*)calloc(n + 2, sizeof(uint64_t));
	G->adjncy = malloc(sizeof(uint32_t) * (bigm));
	G->adjwgt = malloc(sizeof(uint32_t) * (bigm));
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
	uint32_t m;
	int fmt;
	char ch = 'a';
	ifp = fopen(fname, "r");
	if (ifp == NULL)
	{
		printf("\nFailed to open file.\n");
		exit(1);
	}

	fscanf(ifp, "%d %d %d%c", &n, &m, &fmt, &ch);
    printf("here blah blah : n, m , fmt: %d %d %d \n",n,m, fmt);

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


void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
struct graph* randomization(struct graph* G,int perm_capacity){

    int r;
    int i;
    uint32_t j;
    uint32_t k;
    uint32_t temp;
    uint32_t* permutes;

    int newid;
    int neighbor;
    int neigh_weight;
    struct graph *new_G = graph_new(G->n,G->m,G->fmt);
    permutes=malloc(sizeof(uint32_t) * (G->n + 2));
    //temp_vtx[1] = 1;
    new_G->vtx[0] = 0;
    new_G->vtx[1] = 1;
    permutes[0]=-1;
   for (i=1;i<=G->n;i++){
        permutes[i]=i;
    }
   // srand(1);
    for (i = G->n; i>=perm_capacity; i--){
        r= (rand() % i) + 1;
        //printf("R : %d \n",r);
        swap(&permutes[i], &permutes[r]);

    }

 /* for (i=1;i<=G->n;i++){

        printf("permutes[%d]= : %d \n",i,permutes[i]);



    }*/
   /* for (i=1;i<=G->n;i++)
        printf("G->reversepermutes[%d] : %d \n",i,G->reversepermutes[i]);
        */
    new_G->vwgt=G->vwgt;
    for (i=1; i<=G->n; i++){
       new_G->vtx[permutes[i] + 1]=G->vtx[(i + 1)] - G->vtx[i];

     //  printf("temp_vtx[%d]= : %d \n",i,temp_vtx[i]);
      // new_G->vtx[i+1]=new_G->vtx[i]+temp_vtx[i+1];
     //  printf("new_G->vtx[%d]= : %d \n",i,new_G->vtx[i]);
     //  printf("new_G->vwgt[%d]= : %d \n",i,new_G->vwgt[i]);
    }
    for (i=2; i<=G->n + 1;i++){
            new_G->vtx[i]+= new_G->vtx[i-1];
          //  printf("new_G->vtx[%d]= : %d \n",i,new_G->vtx[i]);
    }
    new_G->adjncy[0]=0;
    new_G->adjwgt[0]=0;
    for (i = 1; i <= G->n; i++){
        newid = permutes[i];
        j=new_G->vtx[newid];
        assert(G->vtx[i + 1 ] - G->vtx[i] == new_G->vtx[newid + 1] - new_G->vtx[newid] );
        for (k=G->vtx[i]; k<G->vtx[i+1]; k++){

            neighbor = G->adjncy[k];

            neigh_weight=G->adjwgt[k];
            new_G->adjncy[j + (k - G->vtx[i])] = permutes[neighbor];
            new_G->adjwgt[j + (k - G->vtx[i])]=neigh_weight;


        }

    }


   /* printf(" ----------------adjncy------------------\n");
    for (i = 1; i <= 2 * G->m; i++){
        printf("new_G->adjncy[%d] : %d \n",i , new_G->adjncy[i]);
        printf("new_G->adjwgt[%d] : %d \n",i , new_G->adjwgt[i]);
    }*/


    return new_G;
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
struct Bins* LDGpartitioning(int K, struct graph *G, int capacity,int threadnum)
{
	int i;
	int j;
	int u;
	uint32_t v;
	int b;
	int bmax;
	uint32_t r;
	uint32_t kstart;
	int64_t kend;
	uint32_t k;
	int cutcount = 0;
	int temp;
	int maxwgt = 0;
	int bwt;
    float max = 0.0;
	//initialize bin ,bin affinity and bin weight to zero
	struct Bins *binha = Bins(G->n, K);
	uint32_t *BinAff;
	//BinAff = (uint32_t*)calloc(K + 1, sizeof(uint32_t));
	int *Blist;
	//Blist = (uint32_t*)calloc(K + 1, sizeof(uint32_t));

	int bi;

	for (i = 0; i <= G->n; i++) {
		if (G->vwgt[i] > maxwgt) {
			maxwgt = G->vwgt[i];
		}
	}
    int left = 0;
	int inbins = G->n % K;
    if (inbins != 0 )
       left = 1;
	
    int excess = ((left + (G->n / K)) * maxwgt);

    printf("excess : %d capacity : %d\n",excess,capacity);
    int chunk = 512;
	//for each u in graph find its corresponding bin
    #pragma omp parallel num_threads(threadnum) private(u,j,k,i,bi,BinAff,Blist,kstart,kend,bmax,v,b,max) shared(capacity,G,K,binha,excess) if (G->n >=1000)
    {
        if (18 * chunk > G->n)
            chunk = G->n/threadnum;

        BinAff=(uint32_t*)calloc(K+1,sizeof(uint32_t));
        Blist=(int*)calloc(K+1,sizeof(int));
        int last_bin;
        last_bin = (rand() % K) + 1;
        double start=omp_get_wtime();

        int bizeros=0;
        #pragma omp for schedule (dynamic,chunk) nowait
        for (u = 1; u <= G->n; u++) {

            bi = 0;

            //for each v in adj(u) find its affinity to u
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
            //no neighbor of vertex u is assigned to any bin
            if (bi != 0) {

                
               // norandom++;

                bmax = 0;
                for (i=1; i <= bi; i++){

                    b = Blist[i];
                    bwt = binha->binwgt[b];
                    if (((bwt + G->vwgt[u]) <= capacity )&&(float)(BinAff[b] * (1 - ((float)bwt / capacity)))>max){
                        max = (float)BinAff[b] * (1 - (float)bwt / capacity);
                        bmax = b;

                    }
                    BinAff[b] = 0;
                }
            }
            
            int minimum = capacity;
            int minindez = -1;

            //when the vertex did not assign to the bin of its neighbor because of high capacity we assign it randomly

            if (bi == 0 || bmax == 0) {

              //  printf("U : %d , bi : %d , randomness : %d \n",u,bi,randomness );
                //if whole vertices will fit in the bins with given capacity

                if (excess <= capacity) {

                        while ( ( binha->binwgt[last_bin]+ G->vwgt[u] > (capacity)) ){
                            bizeros++;

                           // printf("Bwt : %d vwgt : %d capacity : %d")
                          last_bin = ((last_bin % K) + 1);
                          bizeros++;


                        }


                       bmax = last_bin;
                       last_bin = ((last_bin % K) + 1);

                }
                else
                {
                    while (binha->Bin[u] == 0) {
                        
                       

                        if (binha->binwgt[last_bin] + G->vwgt[u] <= (capacity))
                        {
                            bmax = last_bin;
                            binha->Bin[u] = bmax;
                            last_bin = ((last_bin % K) +1);
                            bizeros++;
                            continue;
                        }else {
                            for (i = 1; i <= K; i++) {
                                bizeros++;
                                if (binha->binwgt[i] < minimum) {
                                    minimum = binha->binwgt[i];
                                    minindez = i;
                                }
                            }
                            bmax = minindez;
                            binha->Bin[u] = bmax;

                            continue;
                        }

                    }

                }

            }

            binha->Bin[u] = bmax;
          // printf("bin[%d] : %d by thread : %d \n",u,binha->Bin[u],omp_get_thread_num());
            #pragma omp atomic
            binha->binwgt[bmax] = binha->binwgt[bmax] + G->vwgt[u];


        }

        double end=omp_get_wtime();
        
	printf( "Time taken by thread %d is %f\n", omp_get_thread_num(), end-start );
  //  printf("bizeros : %d \n",bizeros );

    }

     
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
struct graph* coarsening(int *Bin, int K, uint32_t *binwgt, struct graph *G, int threadnum) {
	uint32_t *bvix;

	uint32_t *BV;
	int v;
	int u;
	struct graph *Gprime = graph_coarse(K);
	int *Bmap;
    int i;
    uint32_t kstart;
    uint32_t kadjst;
    int64_t kadjend;
    uint32_t kadj;
    int64_t kend;
    int64_t kix;
	int b = 0;
	Gprime->vwgt[0] = 0;
	int *suma;
	bvix = (uint32_t*)calloc(K + 2, sizeof(uint32_t));

    uint64_t *sumvtx;
    int bmpt;
    uint32_t bx;
    uint32_t bst;
    int64_t bend;
    uint32_t * priVtx = (uint32_t*)calloc(K + 3, sizeof(uint32_t));
	#pragma omp parallel num_threads(threadnum) private(u,v,b,kstart,kend,kadjst,kadjend,kadj,kix,i,bmpt,Bmap,bx,bst,bend) shared(binwgt,bvix,G,K,Bin, Gprime,BV,priVtx) if (K >=100)
    {
        double start=omp_get_wtime();

        Bmap = (int*)calloc(K + 1, sizeof(int));

        BV = malloc(sizeof(uint32_t) * (G->n + 1));
        int ithread = omp_get_thread_num();    //label of classmate
        int nthreads = omp_get_num_threads();
        suma = malloc(sizeof(uint64_t)* (nthreads+2));
        sumvtx = malloc(sizeof(uint64_t) * (nthreads+2));
        suma[0] = 0;
        sumvtx[0] = 0;

        #pragma omp for schedule(dynamic,512)
        for (b = 1; b <= K; b++) {
            Gprime->vwgt[b] = binwgt[b];

        }

        #pragma omp master
        {
            for (v = 1; v <= G->n; v++) {
            int t;
            t = Bin[v];
            bvix[t] += 1; //this can be problematic

           // printf( "bvix updated by thread %d is %d %d\n", omp_get_thread_num(),t,bvix[t] );
            }
        }
        
        #pragma omp barrier
        //prefix sum for bvix
        int s = 0;
        if (ithread == 0){
            s = bvix[1];
            //printf("s by thread 0 is : %d\n", s);
        }
        #pragma omp for schedule(static) nowait
        for (b = 2; b <= K; b++) {


           //  bvix[b] = bvix[b] + bvix[b - 1];
           s += bvix [b];
           bvix[b] = s;

        }
        suma[ithread + 1] = s;
       // printf("suma[%d] : %d\n", ithread + 1, suma[ithread + 1]);

        #pragma omp barrier
        int offset = 0;
        for(i=0; i<(ithread+1); i++) {
            offset += suma[i];
        }
       // printf("offset for thread %d : %d\n",omp_get_thread_num(),offset );

        #pragma omp for schedule(static)
        for (b = 2 ; b <= K ; b++) {
           bvix[b] += offset;
        //   printf( "bvix[%d] : FINAL cumulated by thread %d is %d\n",b, omp_get_thread_num(),bvix[b] );
        }

        //
   //     printf("after forth pragma: \n");
        bvix[K + 1] = bvix[K];
        #pragma omp barrier
        //printf("\nbvix[%d] = %d and bvix[%d] = %d in thread : %d\n", K, bvix[K], K+1, bvix[K+1], omp_get_thread_num());
        #pragma omp master
        {
        for (v = 1; v <= G->n; v++) {
            b = Bin[v];
            int f;
            f = bvix[b] - 1;
            BV[f] = v;
           // printf( "BV[%d] : by thread %d is %d\n",f, omp_get_thread_num(),BV[f] );
            bvix[b] -= 1;

            }
        }
      //  printf("\nAFTER FOR bvix[%d] = %d and bvix[%d] = %d in thread : %d\n", K, bvix[K], K+1, bvix[K+1], omp_get_thread_num());
        #pragma omp barrier

        #pragma omp for schedule (static)
        for (b = 1; b <= K; b++) {

            kstart = bvix[b];
           // printf("kstart for thread : %d is %d \n", omp_get_thread_num(),kstart);

            kend = bvix[b + 1] - 1;
          //  printf("kend for thread : %d is %d \n", omp_get_thread_num(),kend);

            for (kix = kstart; kix <= kend; kix++) {
                u = BV[kix];
                kadjst = G->vtx[u];
                kadjend = G->vtx[u + 1] - 1;

                for (kadj = kadjst; kadj <= kadjend; kadj++) {
                    v = G->adjncy[kadj];
                    int bt;
                    bt = Bin[v];
                    bmpt = Bmap[bt];
                    if (bt != b && (bmpt == 0 || bmpt != b)) {
                        Gprime->vtx[b] += 1;
                   //    printf( "Gprime->vtx[%d] : by thread %d is %d\n",b, omp_get_thread_num(), Gprime->vtx[b]);
                        Bmap[bt] = b;
                     //  printf( "Bmap[%d] : by thread %d is %d\n",bt, omp_get_thread_num(), Bmap[bt]);
                        }
                    }

            }
        }

        
        for (i = 1; i <= K + 1; i++){
            priVtx[i + 1] = Gprime->vtx[i] + priVtx [i];
          //  printf("priVtx[%d] : %d \n",i,priVtx[i]);
        }
      
       

      #pragma omp barrier
       Bmap = (int*)calloc(K + 1, sizeof(int));
        //cumulative sum

        uint64_t l = 0;
        if (ithread == 0)
            l = Gprime->vtx[1];

        #pragma omp for schedule(static) nowait
        for (b = 2; b <= K; b++) {


           //  bvix[b] = bvix[b] + bvix[b - 1];
           l += Gprime->vtx[b];
           Gprime->vtx[b] = l;

        }
        sumvtx[ithread + 1] = l;
       
        #pragma omp barrier
        uint64_t off = 0;
        for(i=0; i<(ithread+1); i++) {
            off += sumvtx[i];
        }


      //  printf( "Gprime->vtx[1] : FINAL cumulated by thread %d is %d\n", omp_get_thread_num(),Gprime->vtx[1] );
        #pragma omp for schedule(static)
        for (b = 2 ; b <= K ; b++) {
           Gprime->vtx[b] += off;
         //  printf( "Gprime->vtx[%d] : FINAL cumulated by thread %d is %d\n",b, omp_get_thread_num(),Gprime->vtx[b] );

        }

        Gprime->vtx[K + 1] = Gprime->vtx[K];

        int cadjsize = Gprime->vtx[K + 1];
       // printf("cadjsize : %d \n",cadjsize);
        Gprime->m = cadjsize;
        Gprime->adjncy = malloc(sizeof(uint32_t) * cadjsize + 1);
        Gprime->adjwgt = malloc(sizeof(uint32_t) * cadjsize + 1);

        #pragma omp barrier
        /*for (i = 0 ; i <= K+1; i++){
            printf("bvix[%d] : %d by thread :%d \n", i,bvix[i],omp_get_thread_num());
        }*/
        #pragma omp for schedule (static)
        for (b = K; b >= 1; b--) {

            kstart = bvix[b + 1];
            //printf("kstart for thread : %d is %d \n", omp_get_thread_num(),kstart);
            kend = bvix[b] + 1;
            //printf("kend for thread : %d is %d \n", omp_get_thread_num(),kend);
            for (kix = kstart; kix >= kend; kix--) {
                u = BV[kix - 1];
                kadjst = G->vtx[u];
                kadjend = G->vtx[u + 1] - 1;
                //printf("kadjst for thread : %d is %d \n", omp_get_thread_num(),kadjst);
                //printf("kadjend for thread : %d is %d \n", omp_get_thread_num(),kadjend);
                for (kadj = kadjst; kadj <= kadjend; kadj++) {
                    v = G->adjncy[kadj];
                    int bt;
                    bt = Bin[v];
                    //printf("bt : %d by thread %d\n",bt,omp_get_thread_num());
                    if (bt != b && Bmap[bt] == 0) {
                        int t = Gprime->vtx[b] - 1;
                        Gprime->adjncy[t] = bt;
                        //printf("t : %d by thread :%d\n", t,omp_get_thread_num());
        //                printf( "Gprime->adjncy[%d] : by thread %d is %d\n",t, omp_get_thread_num(), Gprime->adjncy[t]);
                        Bmap[bt] = G->adjwgt[kadj];
                        //printf( "Bmap[%d] : by thread %d is %d\n",bt, omp_get_thread_num(), Bmap[bt]);
                        Gprime->vtx[b] -= 1;
       //                 printf( "REDUCED VTX [%d] : by thread %d is %d\n",b, omp_get_thread_num(), Gprime->vtx[b]);

                    }
                    else if (bt != b) {
                        Bmap[bt] += G->adjwgt[kadj];
                //        printf( "ADDED Bmap[%d] : by thread %d is %d\n",bt, omp_get_thread_num(), Bmap[bt]);

                    }

                }
            }
            bst = priVtx[b];
        //    printf("bst for thread : %d is %d \n", omp_get_thread_num(),bst);
            bend = priVtx[b + 1] - 1;
         //   printf("bend for thread : %d is %d \n", omp_get_thread_num(),bend);

            for (bx = bst; bx <= bend; bx++) {
                int p;
                p = Gprime->adjncy[bx];
            //    printf( "Bmap[%d] : by thread %d is %d\n",p, omp_get_thread_num(), Bmap[p]);
                Gprime->adjwgt[bx] = Bmap[p];

             //   printf( "Gprime->adjwgt[%d] : by thread %d is %d\n",bx, omp_get_thread_num(), Gprime->adjwgt[bx]);
                Bmap[p] = 0;
            }
        }
        
        

        
        double end=omp_get_wtime();
        printf( "Coarsening time taken by thread %d is %f\n", omp_get_thread_num(), end-start );
    }

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
struct Bins* refLDG (int K,struct graph *G,struct Bins *binha,int capacity,int threadnum ){

    int i;
    int j;
    int u;
    int v;
    int b;
    int bmax;
    uint32_t r;
    uint32_t kstart;
    int64_t kend;
    uint32_t k;
    //initialize bin ,bin affinity and bin weight to zero
    uint32_t *BinAff;
    int *Blist;
    //for each u in V  in some order do
    int bi;
    //for random number generation
    float max=(float)0.0;
    int bwt=0;
    int currentbin;
    int chunk = 512;

#pragma omp parallel num_threads(threadnum) private(u,k,i,bi,b,kstart,kend,v,bmax,max,bwt,BinAff,Blist,currentbin) shared(capacity,G,binha,K,chunk) if (G->n >= 1000)
 {
            if (18 * chunk > G->n)
                chunk = G->n/threadnum;
            double startr=omp_get_wtime();

            BinAff=(uint32_t*)calloc(K+1,sizeof(uint32_t));
            Blist=(int*)calloc(K+1,sizeof(int));
            int zeros = 0;
            #pragma omp for schedule (dynamic,512) nowait
            for (u=1;u<=G->n;u++)
            {

                currentbin=binha->Bin[u];
                bi=0;
                kstart=G->vtx[u];
                kend=G->vtx[u+1]-1;

                for (k=kstart;k<=kend;k++){

                    v=G->adjncy[k];
                    b=binha->Bin[v];
                    if(BinAff[b]==0){
                        bi++;
                        Blist[bi]=b;

                    }
                    BinAff[b]=BinAff[b]+G->adjwgt[k];
                }
                bmax=0;
                max=0.0;
                for (i=1;i<=bi;i++){
                    b=Blist[i];
                    //printf("b : %d \n",b);
                    bwt=binha->binwgt[b];

                    if (currentbin==b)
                        bwt=bwt-G->vwgt[u];

                    if (((bwt+ G->vwgt[u])<=capacity )&&(float)(BinAff[b]*(1-((float)bwt/capacity)))>max){
                        max=(float)BinAff[b]*(1-(float)bwt/capacity);
                        bmax=b;


                    }
                    BinAff[b]=0;
                }
                if (bmax!=0 && bmax!=binha->Bin[u]){
                    #pragma omp atomic
                    binha->binwgt[binha->Bin[u]]=binha->binwgt[binha->Bin[u]]-G->vwgt[u];
                    #pragma omp atomic
                    binha->binwgt[bmax]=binha->binwgt[bmax]+G->vwgt[u];

                }
                if (bmax==0){
                    zeros++;
                    bmax=binha->Bin[u];
                }

                binha->Bin[u]=bmax;
                
               // printf("refine bin[%d] : %d  by thread : %d \n",u,binha->Bin[u],omp_get_thread_num());
                //printf("refine bmax : %d by thread : %d \n",bmax,omp_get_thread_num());
               // printf("refine binwgt[%d] : %d by thread : %d\n",bmax,binha->binwgt[bmax],omp_get_thread_num());


            }
         
            double endr=omp_get_wtime();
           printf( "refinedTime taken by thread %d is %f\n", omp_get_thread_num(), endr-startr );

    }

      //  printf("weight of all bins : %d \n", weights);

        /*for (i=1;i<=G->n;i++){
            binha->binwgt[binha->Bin[1]]++;
        }*/
    return binha;

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
struct Bins* uncoarsen(struct Bins *binprime, struct Bins *initbins, struct graph* G, int K,int threadnum) {

	int i = 0;
	struct Bins *bin = Bins(G->n, K);
	#pragma omp parallel num_threads(threadnum) private(i) shared(bin,binprime,initbins,K,G) if (G->n >= 1000)
    {
        double startr=omp_get_wtime();
        #pragma omp for schedule (static)
        for (i = 1; i <= G->n; i++) {
            bin->Bin[i] = binprime->Bin[initbins->Bin[i]];
        }
        bin->binwgt[0] = 0;
        #pragma omp for schedule (static)
        for (i = 1; i <= K; i++) {
            bin->binwgt[i] = binprime->binwgt[i];

        }
        double endr=omp_get_wtime();
        printf( "uncoarsening taken by thread %d is %f\n", omp_get_thread_num(), endr-startr );  
    }
	return bin;

}
float imbalance (int K, double wavg, struct Bins* bins){
	uint32_t maximb = 0;
    double imbalance = 0.0;
    int i;
    for (i = 0; i < K + 1; i++) {
			if (bins->binwgt[i] > maximb) {
				maximb = bins->binwgt[i];
			}

		}
	//	printf("maximb : %d\n", maximb);
	//	printf("wavg : %f\n", wavg);

		imbalance = (float)maximb / wavg;
		return imbalance;
}
float edgecuts(int K, struct graph* G,struct Bins* bins, uint32_t adjwgttotal){
        int cutcount = 0;
		int u;
        uint32_t kstart;
        int64_t kend;
        uint32_t k;
        int v;
		for (u = 1; u <= G->n; u++) {
			//for each v in adj(u)
			kstart = G->vtx[u];
			kend = G->vtx[u + 1] - 1;
			for (k = kstart; k <= kend; k++) {
				v = G->adjncy[k];
				if (bins->Bin[u] != bins->Bin[v]) {
					cutcount += G->adjwgt[k];
				}
			}
		}

		cutcount = cutcount / 2;

		float edgecut = 0.0;
		edgecut = (float)cutcount / (adjwgttotal / 2);
		return edgecut;
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
	float edgecutsum = 0.0;
	float imbalancesum = 0.0;
	char *fname = argv[1];
	int K = atoi(argv[2]);
	int lvl = 0;
	int beta = atoi(argv[3]);
	double alfa=atof(argv[4]);
   	int threadnum = atoi(argv[5]);
	int setting = 1;

	struct graph* G = readgraphfromfile(fname);
    //showgraph(G);
    //printf("here\n");
    int perm_capacity=ceil(alfa * G->n);
    struct graph* randG;
	printf("Vertices and Edges: %d %d\n",G->n,G->m);
	// running the method 5 times
    for (rpt = 0; rpt < 5; rpt++) {
        randG=randomization(G,perm_capacity);
        float eps = (float)0.0;
        uint32_t adjwgttotal = 0;
        for (l = 1; l < G->vtx[G->n + 1]; l++) {
            adjwgttotal += G->adjwgt[l];

        }

        int id = 0;
        G->id = id;
        struct graph* Gprime = randG;

        int wtotal = 0;
        double wavg = 0.00000;
        wtotal = G->n;
        eps = (float)0.300000;
        float temp = (float)eps - 0.100000;
        int refcapacity = 0;
        struct Bins* partitions;
        int capacity = 0;
        i = 1;
        int ki = 0;
        //finding beta value
        while (G->n / pow(beta, i) > K * beta) {

            i++;
        }
        lvl = i;
        struct Bins* binarray[lvl + 2];
        struct graph* grapharray[lvl + 2];
        int karray[lvl + 2];
        karray[0] = 0;
        for (i = 1; i < lvl; i++) {
            karray[i] = G->n / pow(beta, i);
        }
        karray[lvl] = K;
		time_t start = time(NULL);
        for (i = 1; i <= lvl; i++) {
          //  printf("LEVEL : %d===================================================================\n",i);

            float range = (float)0.0;
            if (lvl > 1)
                range = temp / (lvl - 1);
            else
                range = 0;
            wavg = ((double)wtotal / karray[i]);
            capacity = ceil((wavg * (1 + eps)));

            grapharray[i] = Gprime;
          //  if (i >= lvl)
               // showgraph(grapharray[i]);
          //  printf("before LDG\n");
            partitions = LDGpartitioning(karray[i], grapharray[i], capacity,threadnum);
          
         
           
            
                


            if (setting == 1) {
                partitions = refLDG(karray[i], grapharray[i], partitions, capacity,threadnum);
             

            }
            id++;
            binarray[i] = partitions;
           
            
      //      Gprime = coarsening(partitions->Bin, karray[i], partitions->binwgt, grapharray[i],1);
      //      printf("ONE THREAD COARSE : \n");
      //      showgraph(Gprime);
            Gprime = coarsening(partitions->Bin, karray[i], partitions->binwgt, grapharray[i],threadnum);
      //      printf("MORE THREAD COARSE : \n");
        //    showgraph(Gprime);
            Gprime->m = Gprime->vtx[K + 1];
            Gprime->id = id;
            eps -= range;

        }




            //printf("wtotal : %d\n", wtotal);

        wavg = ((double)wtotal / K);

        refcapacity = ceil(wavg*(1.1));
        //printf("refcapacity : %d\n", refcapacity);
        int templvl = lvl;

        while (templvl >= 2) {
            printf("UNLEVELING :%d ============================\n", templvl);
            binarray[templvl - 1] = uncoarsen(binarray[templvl], binarray[templvl - 1], grapharray[templvl - 1], K,threadnum);


            binarray[templvl - 1] = refLDG(K, grapharray[templvl - 1], binarray[templvl - 1], refcapacity,threadnum);
             //   printf("edgecuts of last ref: %f\n", edgecuts(K,grapharray[templvl - 1],binarray[templvl - 1],adjwgttotal));
            //    printf("imbalance of last ref: %f\n", imbalance(K,wavg,binarray[templvl - 1]));
            templvl--;
        }
        int m = 2;
        while (m > 0){
            printf("LASTREFINEMENTS :%d================\n", m);
            binarray[1] = refLDG(K, grapharray[1], binarray[1], refcapacity,threadnum);
            m--;
        }
		
		printf("WHOLE WALL TIME %.2f\n", (double)(time(NULL) - start));

        edgecutsum += edgecuts(K,grapharray[1],binarray[1],adjwgttotal);

        imbalancesum += imbalance(K,wavg,binarray[1]);
    		//printf("imbalance : %f\n", imbalance(K,wavg,binarray[1]));


	}//end of five run


	printf("number of levels : %d \n", lvl);
	printf("avg imbalance : %f\n", imbalancesum / 5);
	printf("avg edgecut of all : %f\n", edgecutsum / 5 );
	return 0;
}

