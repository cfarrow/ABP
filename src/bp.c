/* ------------------------------------------------------------------------------- 

 ---------------------------------------------------------------------------------- */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

//#define NSITES 100000
#define MAX_NEIGHBORS 50
#define MAX_BONDS (MAX_NEIGHBORS*NSITES/2)
#define RNG_N_OF_BITS 32
#define MAXGVALUES  510
#define MAXLAYERS (5000)
#define GO_STEP 1
#undef  RANDOM_BC 

#undef DEBUG_CUT_BOND
#undef CHECK_CULLING
#undef CHECK_QUEUEING
#undef MAKE_CHECKS

#define ENABLE_ASSERTS

#define RANDOM_BOND 0
#define LR_NETWORK 1

#define MAX_DIMENSION 6
#define MAX_NUMBER_OF_CORNERS 64 // = 2^MAX_DIMENSION
#define MAX_N_ALPHASTEPS 100

// Separation of bonds into short and long ones.
#define SHORT_BOND_LENGTH  ((double)2.3)
#define MAX_NUMBER_OF_SHORT_BONDS 100



FILE *out0_ptr;
FILE *out1_ptr;
FILE *out2_ptr;

int mcull=2;/* minimal n. of neighbors for culling */

// Basic arrays (network, k-core identity, etc.)
unsigned long neighbor[NSITES][MAX_NEIGHBORS]; /* index of neighbor */
unsigned long N_neighbors[NSITES];             /* number of connected neighbors */
unsigned long N_mcore_neighbors[NSITES];      /* n. of mcore neighbors */

short mcore[NSITES];        /* 1: mcore  0: inactive */
short present[NSITES];       /* 1: present 0: absent */

// a bond links site i to site j
unsigned long ivalue[MAX_BONDS]; /* site i of a present bond */
unsigned long kvalue[MAX_BONDS]; /* site j of a present bond */

// Variables for random link generation
double bond_weight[NSITES]; // One weight for each value of the difference of labels.

double site_weight[NSITES]; // This holds the accumulated normalized weight of
			    // all this site's yet unconnected bonds. Sites
			    // for connection are chosen with this
			    // probability.

unsigned long number_of_sites; /* total sites of the network, normally equal to NSITES */
unsigned long mcore_sites;  /* n. of mcore sites */
unsigned long present_sites; /* n. of present sites */
unsigned long corona_sites;  /* n. of sites with exactly mcull neigbors */
unsigned long mcore_bonds;  /* n. of mcore bonds */
unsigned long present_bonds; /* n. of present bonds */


/* for creation of d-dim. LR-Networks: */
int     space_dimension;
long    number_of_corners;              // Number of corners of a hypercube (2^d)
unsigned long jump[MAX_DIMENSION];	// contains periodicity-information
                                        //     for labelling regular 
                                        //     d-dim. lattice
long direction;				// loop-index for jump[]
long coordinate[MAX_DIMENSION];         // Cartesian coordinates of a site.

// Used to build coordinates of plaquette corners
long plaquette_limit[MAX_DIMENSION][2];

unsigned long lattice_length, half_lattice_length;  /* for regular lattice of
						     * LR-Network */
double lambda, dexponent;		/* constants for calculation of the
					 * distribution (inversion method) */
double alpha;				/* probablity distribution decay exponent */
double alpha_intpart, alpha_floatpart;
double alpha_vector[MAX_N_ALPHASTEPS];
double alpha_min, alpha_max, alpha_step;  /* array and frontiers of link
					   * probability decay exponents */

unsigned long n_alphaintervals;		  /* # of alpha values */
unsigned long n;			  /* loop-index for alphas */


int current_gamma_index;  // routine label_sites works with variable k
			  // (gamma-value-index) from main, passed via this
			  // variable

// Separation of bonds into short and long ones.
unsigned long number_of_short_bonds;
unsigned long short_bond_label[MAX_NUMBER_OF_SHORT_BONDS]; // label of "short bonds"
double        p_s;  // Accumulated weight of all short bonds



/* -------------------------------------------------------- *\
 *	   Definitions for shift-register RNG	
\* -------------------------------------------------------- */
#define 	JJQ		32
#define 	JJP		521
#define		MODUL		JJP
#define 	JJD		(JJP-JJQ)

unsigned long 		ibits[MODUL];


#define         UL_MAX          4294967295U
#define         L_MAX           4294967296
#define         DRNG_FACTOR     2.328306437e-10



/* -------------------------------------------------------- *\
		 Some macros 
\* -------------------------------------------------------- */

#define DEFAULT_LRANDOM_NUMBER  lrand48()
#define DEFAULT_DRANDOM_NUMBER  drand48()
#define DEFAULT_SET_SEED_RNG(a) srand48(a)

#define LRANDOM_NUMBER          my_lrand()
#define DRANDOM_NUMBER          (LRANDOM_NUMBER*DRNG_FACTOR)
#define SET_SEED_RNG(a)         my_RNG_seed_init(a)

#define MSG(s1)			{printf(s1);fflush(stdout);}
#define DIE(s1,s2)		{printf(s1);printf(s2);exit(1);}
#define MEMERR			" memory allocation error\n"
#define FOPERR			" error opening file\n"
#define MIN(a,b)		((a<b)? a : b )
#define MAX(a,b)		((a>b)? a : b )
#define ABS(a)		        ((a>0)? a :-a )
#define SGN(x)          (x>0.0?1.0:-1.0)
#define HERE  "undefined"
#define TELL(s1)   {printf(HERE);MSG(s1)}

#define SWAP(a,b)  {unsigned long temp;temp=a;a=b;b=temp;}

#define         FOPEN(fptr,file_name)                              \
    if( !(fptr = fopen(file_name,"w"))){                           \
  printf(":( cannot open file %s for output.\n",file_name);exit(1);}


/* Seed for random numbers */
unsigned long int seed = 1234567;

/* Queue definitions */
#define QLEN    (NSITES+1)
unsigned long front,rear;
unsigned long queue[QLEN];



/* Prototypes */
void add_to_queue(unsigned long);
unsigned long  extract_from_queue(void);
void my_RNG_seed_init(unsigned long);
unsigned long my_lrand (void);

void cut_one_bond(int);
void cull_sites(int);
void init_random_sites(float);
void check_connectivity(void);
int check_activity(int);
unsigned long count_mcore_neighbors(unsigned long);
void label_corona(void);
void collect_to_queue(int);
int  are_neighbors(long, long);
void d_dim_lr_network(double, float);
void calculate_bond_weights(void);

void average_neighbors(void);
void plot_links(float);


int main(void){
    
    unsigned long sample,nsamples;
    unsigned long k,kmax;
    
    float gamma0,psite=1.;
    float gamma,gammaf;
    float ming,maxg;
    int   intervals;                  // the number of gamma-steps..
    float gamma_vector[MAXGVALUES];   // ..which are stored in this array

    
    
    printf("\nNSITES = %u\n", NSITES);    
    printf("Enter base dimension d: ");
    scanf("%u", &space_dimension);
    assert(space_dimension<=MAX_DIMENSION);
    printf("%u\n",space_dimension);
    number_of_corners=(long)pow(2,space_dimension);
    assert(number_of_corners<=MAX_NUMBER_OF_CORNERS);    
    printf("LR-Network embedded in %u dimensions\n", space_dimension);
    // number_of_sites stores the real number of sites, which may be smaller
    // than NSITES.
    lattice_length = (unsigned long) pow(NSITES, 1.0/space_dimension);
    half_lattice_length = (unsigned long) lattice_length/2;
    number_of_sites=(int)pow((double)lattice_length,(double)space_dimension);
    printf("with lattice_length=%lu => half_lattice_length=%lu\n", lattice_length, half_lattice_length);
    printf("and consist of %lu sites.\n\n", number_of_sites);
    
    /* Calculate "coordinate jumps" for labelling of regular d-dim. lattice with helical BC */
    jump[0]=(unsigned long) 1;
    for(direction=1;direction<(unsigned long)space_dimension;direction++){
	jump[direction]=lattice_length*jump[direction-1];
    }
    
    printf("Enter min. and max. decay exponent, and # of alpha-values (alphamin, alphamax, nintervals): ");
    scanf("%lf %lf %lu",&alpha_min, &alpha_max, &n_alphaintervals); 
    printf("%lf %lf %lu\n",alpha_min,alpha_max,n_alphaintervals); 
    assert(n_alphaintervals<MAX_N_ALPHASTEPS && n_alphaintervals>=0);
    

    // Generate alpha values
    if(n_alphaintervals==0){ // do a single alpha value
	alpha_step = (double)0;
    }else{  // do many alpha values
	alpha_step = (alpha_max-alpha_min)/n_alphaintervals;
    }
    for(n=0;n<=n_alphaintervals;n++)alpha_vector[n]=(alpha_min + n*alpha_step);
    
    
    printf("Enter value of culling parameter k: ");
    scanf("%d",&mcull);
    printf("%d\n",mcull);
    printf("Initial value of coordination, final value and # of intervals (gammai,gammaf,nint): ");
    scanf("%f %f %d",&gamma0,&gammaf,&intervals);
    ming=MIN(gamma0,gammaf);
    maxg=MAX(gamma0,gammaf);
    printf("%f %f %d\n",gamma0,gammaf,intervals);
    assert((intervals<MAXGVALUES)&&(intervals>=0));
    printf("Enter integer seed for RNG: ");
    scanf("%lu",&seed);
    printf("%lu\n",seed);
    printf("Enter # of samples: ");
    scanf("%lu",&nsamples);
    printf("%lu\n",nsamples);
    // Initialize RNG
    seed|=1L;(void) SET_SEED_RNG(seed);

    
    
    
    /* BEGIN ALPHA-LOOP */
    for(n=0;n<=n_alphaintervals;n++){
	
	alpha=alpha_vector[n];
	printf("Alpha = %2.3lf (value %lu of %lu)\n",alpha,n+1,n_alphaintervals+1);

	// Calculate weights needed by routine d_dim_lr_network(). This must
	// be called once for each value of alpha, BEFORE generating any
	// networks.
	(void)calculate_bond_weights();
	
	for(k=0;k<intervals+1;k++){
	    if(intervals){
		gamma_vector[k]=maxg-k*(maxg-ming)/intervals;
	    }else{ /* 0 intervals, just one point */
		gamma_vector[0]=maxg;
	    }
	}
	
	kmax=intervals+1;
	
	/* ----------------- */
	/* BEGIN SAMPLE-LOOP */
	/* ----------------- */
	for(sample=0;sample<nsamples;sample++){

	    if(( !((sample+1)%(MAX(1,nsamples/10)))))printf("Sample %lu\n",sample+1);
	    
	    /* init queue */
	    front=rear=0; 
	    
	    /* Build the network */
	    (void) d_dim_lr_network((double) alpha,(float) maxg);



	    // NOTE: Comment this line out if you don't need to visualize the
	    // network
	    (void) plot_links(maxg);

	    
	    /* Initiallize site dilution with probability psite, and check for
	     * mcull neigbours */
	    (void) init_random_sites(psite);
	    (void) collect_to_queue(mcull); 
	    (void) average_neighbors();
	    
	    
	    /* ---------------- */
	    /* BEGIN GAMMA-LOOP */
	    /* ---------------- */
	    for(k=0;k<kmax;k++){
		
		current_gamma_index=k;
		gamma=gamma_vector[k];
		
		/* ------------ 
		   BOND-CUTTING 
		   ------------ */
		while(present_bonds>(unsigned long)(gamma*number_of_sites/2)){
		    (void)cut_one_bond(mcull);
		}
		
		/* ------- 
		   CULLING    
		   ------- */
		(void) cull_sites(mcull);
		
		
		/*--------------------------------------------------------------*\
		  CLUSTER ANALYSIS
		  \*--------------------------------------------------------------*/
		
		// Code has been removed 
		
	    } /* END GAMMA LOOP */
	    
	}/* END SAMPLE LOOP */
	
    }  /* end of alpha-loop */
    
    
    return(0);
    
    
} /* END OF MAIN */


/*----------------------------------------------------------------*\
  Queue management routines
  
  front=rear means queue is empty
  front=rear+1 means array capacity is full
  
  to add an element do:   
  
  if(front+QLEN> rear){
  data -> queue[rear%QLEN];
  rear++;
  }else{queue full}
  

   to remove an element do:
   
   if(rear > front){
   queue[front%QLEN] -> data;
   front++;
   }else{ queue empty}

\*----------------------------------------------------------------*/


void add_to_queue(unsigned long data){
  int i;
  if(front+QLEN > rear){
    queue[rear%QLEN]=data;
    rear++;
  }else{
    printf("Queue full. Quitting.../n");
    printf("front=%lu, rear=%lu\n",front,rear);
    for(i=0;i<QLEN;i++) printf("%lu\n",queue[i]);
    fflush(0);
    exit(1);
  }
}

unsigned long extract_from_queue(void){
  unsigned long data;
  
  if(rear>front){
    data=queue[front%QLEN];
    front++;
    return(data);
  }else{
    return(-1);
  }
}

/* ---------------------------------------------------------------*\
   
   Set sites present and mcore with probability psite.
   
\* ---------------------------------------------------------------*/
void init_random_sites(float psite){
#undef HERE
#define HERE "init_random_sites:"
    
    unsigned long i;
    
    memset(N_mcore_neighbors,0,sizeof(N_mcore_neighbors));
    memset(present,0,sizeof(present));
    memset(mcore,0,sizeof(mcore));
    
    
    mcore_sites=present_sites=mcore_bonds=0;
    
    /* set sites present and mcore with probability p*/
    for(i=0;i<number_of_sites;i++){
	if( DRANDOM_NUMBER < psite ){
	    mcore[i]=present[i]=1;
	    mcore_sites++;
	    present_sites++;
	}
    }
    
    /* count mcore neighbors */
    for(i=0;i<number_of_sites;i++){
	if(mcore[i]==1){
	    N_mcore_neighbors[i]=count_mcore_neighbors(i);
	    mcore_bonds+=N_mcore_neighbors[i];      
	}
    }
    mcore_bonds/=2;   /* the last routine counts every bond twice */
    return;
}

/* ---------------------------------------------------------------*\
   
Send mcore sites with less than cull neighbors to culling queue.

\* ---------------------------------------------------------------*/
void collect_to_queue(int cull){
#undef HERE
#define HERE "collect_to_queue:"
    
    unsigned long i;
    
    /* count mcore neighbors and queue */
    for(i=0;i<number_of_sites;i++){
	if((mcore[i]==1)&&(count_mcore_neighbors(i)<cull)){
	    mcore[i]=2;
	    queue[rear%QLEN]=i;
	    rear++;
#ifdef  CHECK_QUEUEING
	    TELL(" "); 
	    printf("%lu -> QUEUE\n",i);
#endif  /* CHECK_QUEUEING*/
	}
    }
    return;
}

/* ---------------------------------------------------------------*\
   
Cut a present bond.

This routine now works by keeping a list of all present bonds, and choosing
one for deletion. The list is stored as two arrays ivalue[] and kvalue[],
containing the sites linked by this bond.

\* ---------------------------------------------------------------*/
void cut_one_bond(int cull){
#undef HERE
#define HERE "cut_one_bond:"
    
    unsigned long nnbr1,nnbr2;
    unsigned long i,j,k,l;
    unsigned long bond;
    
#ifdef ENABLE_ASSERTS
    assert(present_bonds>0);
#endif /* ENABLE_ASSERTS */
    
    /* pick a bond */
#ifdef RANDOM_BC
    bond=DRANDOM_NUMBER*present_bonds; // pick at random
#else
    bond=(unsigned long)(present_bonds-1); // pick the "newest"
#endif
    
    /* read list at position bond */
    i=ivalue[bond];
    k=kvalue[bond];
    /* get number of neighbors of the linked sites */
    nnbr1=N_neighbors[i];
    nnbr2=N_neighbors[k];
    
#ifdef ENABLE_ASSERTS
    assert(nnbr1>0); 
    assert(nnbr2>0); 
#endif /* ENABLE_ASSERTS */
    
    /*find index j at i that points to k */
    for(j=0;j<nnbr1;j++){
	if(neighbor[i][j]==k)break;   /* thus neigbor j of i is k */
    }
#ifdef ENABLE_ASSERTS
    assert(neighbor[i][j]==k);
#endif /* ENABLE_ASSERTS */
    
    
    /*find index l at k that points to i */
    for(l=0;l<nnbr2;l++){
	if(neighbor[k][l]==i)break;   /* thus neigbor l of k is i */
    }
    
    /* delete bond from on-site list*/
    neighbor[i][j]=neighbor[i][N_neighbors[i]-1];
    neighbor[k][l]=neighbor[k][N_neighbors[k]-1];
    N_neighbors[i]--;
    N_neighbors[k]--;
    
    /* delete bond from global list */
    ivalue[bond]=ivalue[present_bonds-1];
    kvalue[bond]=kvalue[present_bonds-1];
    present_bonds--;
    
    
    /* test and send to queue if numbers of mcore neighbors < mcull*/
    if((mcore[i]>0)&&(mcore[k]>0)){
	N_mcore_neighbors[i]--;
	N_mcore_neighbors[k]--;
	mcore_bonds--;
    }
    if((mcore[i]==1)&&(N_mcore_neighbors[i]<cull)){   /* cull=mcull */
	mcore[i]=2;
	queue[rear%QLEN]=i;
	rear++;
#ifdef  CHECK_QUEUEING
	TELL(" "); 
	printf("%lu -> QUEUE\n",i);
#endif  /* CHECK_QUEUEING*/
    }
    if((mcore[k]==1)&&(N_mcore_neighbors[k]<cull)){
	mcore[k]=2;
	queue[rear%QLEN]=k;
	rear++;
#ifdef  CHECK_QUEUEING
	TELL(" "); 
	printf("%lu -> QUEUE\n",k);
#endif  /* CHECK_QUEUEING*/
    }
    
    return;
}

/* ---------------------------------------------------------------*\
   
   Recursively cull sites, starting from those in Queue.
   
\* ---------------------------------------------------------------*/

void cull_sites(int cull){
#undef HERE
#define HERE "cull_sites:"

  
  unsigned long i,j,k,culled;

  culled=0;
  while (rear>front){ /* while queue is nonempty */
    i=queue[front%QLEN];  /* pick and deactivate first site in culling-queue */
    front++;
    mcore_sites--;
    mcore[i]=0;
    culled++;

    for(j=0;j<N_neighbors[i];j++){  /* loop thru all neigbors j of i */
      k=neighbor[i][j];
      if(mcore[k]==1){
	N_mcore_neighbors[k]--;
	mcore_bonds--;   /* reduce mcore-values for k due to i being deactivated */
	if(N_mcore_neighbors[k]<cull){
	  mcore[k]=2;
	  queue[rear%QLEN]=k;
	  rear++;   /* add also k to queue, if too few mcore neigbors */
#ifdef  CHECK_QUEUEING
	  TELL(" "); 
	  printf("%lu -> QUEUE\n",k);
#endif  /* CHECK_QUEUEING*/
	}
      }
    }
  }
#ifdef CHECK_CULLING
  TELL(" "); 
  printf("Culled %lu sites.\n",culled); 
#endif /* CHECK_CULLING */
  return;
}

/* ---------------------------------------------------------------*\
   
   Check lists of neighbors for consistency.
   
\* ---------------------------------------------------------------*/

void check_connectivity(void){
#undef HERE
#define HERE "check_connectivity:"

  unsigned long nnbr1,nnbr2;
  unsigned long i,j,k,l;
  
  if(present_bonds==0)return;

  for(i=0;i<number_of_sites;i++){

    nnbr1=N_neighbors[i];
    if(nnbr1>0){
      for(j=0;j<nnbr1;j++){
	k=neighbor[i][j];
	assert(k!=i); /* no self-connections */
	nnbr2=N_neighbors[k];
	assert(nnbr2>0);
	/*find index l at k that points to i */
	for(l=0;l<nnbr2;l++){
	  if(neighbor[k][l]==i)break;
	}	
	assert(neighbor[k][l]==i);

      }
    }
  }
  return;

}


/* ---------------------------------------------------------------*\
   
   Check N_mcore_neighbors for consistency.
   
\* ---------------------------------------------------------------*/
int check_activity(int cull){
#undef HERE
#define HERE "check_activity:"

  unsigned long i,anbrs;
  int error;
  
  if(present_bonds==0)return(0);
  error=0;

  for(i=0;i<number_of_sites;i++){
    if(mcore[i]==1){
      if(!((anbrs=count_mcore_neighbors(i))==N_mcore_neighbors[i])){
	TELL(" ");
	printf("Mcore neighbors mismatch: Site %lu has %lu, claims %lu.\n",
	       i,anbrs,N_mcore_neighbors[i]);
	error=1;
      }
      if(anbrs<cull){
	TELL(" ");
	printf("Mcore site with %lu neighbors: %lu\n",anbrs,i);
	error=1;
      }
    }else if ((present[i]==1)&&(mcore[i]==0)){
      anbrs=count_mcore_neighbors(i);
      if(anbrs>=cull){
	TELL(" ");
	printf("Dead site with %lu neighbors: %lu\n",anbrs,i);
	error=1;
      }
    }
  }
  return(error);

}

/* ---------------------------------------------------------------*\
   
   Return number of mcore neighbors of a site. 
   
\* ---------------------------------------------------------------*/

unsigned long count_mcore_neighbors(unsigned long i){
#undef HERE
#define HERE "count_mcore_neighbors:"

  unsigned long j,anbrs;
  anbrs=0;
  for(j=0;j<N_neighbors[i];j++)if(mcore[neighbor[i][j]]==1)anbrs++;
  return(anbrs);

}


/*---------------------------------------------------------------------------*\
 
  lr_network: 

  Sites are those of a d-dimensional hypercubic network, linked by LR bonds
  whose lengths in d-space are power-law distributed with an exponent alpha.

  P(l) ~ 1/l^{alpha}

  For each link to be connected first a site (isite) is chosen at random, next
  another site (ksite) is chosen for linking.

  Multiple neighbor-linking is not accepted, i.e. if two sites are already
  linked, they are never linked again. Beware of too large values of alpha and
  gamma, the bond-linking routine might never end if it fails to find a pair
  of linkable sites.


  In order to choose the far site (ksite), the following is done:

  Bonds are separated into SHORT and LONG ones, and a different method is used
  for each class.
  
  First the normalized probabilities of all possible bonds stemming from a
  site is calculated. Then the probability p_s that a bond is SHORT (i.e. its
  length is smaller than a given prespecified limit SHORT_BOND_LENGTH) is
  calculated from these numbers. Then, for each site (isite) to which a link
  is to be connected, first it is decided whether the ensuing link is to be
  short or long. This is simply done by drawing a random number and comparing
  it to p_s.


  a) If the resulting link is LONG, the following is done: First a value of
  the radius is generated with power-law distribution (with an exponent
  alpha-(dimension-1)).  Next a random direction is drawn, and a vector is
  built that has this direction and teh previously obtained length. Finally a
  neighbor (ksite) is found by taking the lattice point closest to the tip of
  this vector.

  b) If the resulting link is SHORT (this happens more often for large alpha),
  one instead chooses from a list of short bonds according to their relative
  probabilities. If the number of short bonds is kept small (this is
  controlled by parameter SHORT_BOND_LENGTH) this second method is not too
  inefficient.


\*----------------------------------------------------------------------------*/
void d_dim_lr_network(double alpha, float gamma){
#undef HERE
#define HERE "d_dim_lr_network1:"


    long inbr;
    unsigned long isite,ksite,b_label;
    unsigned long direction,kounter;
    unsigned long total_bonds;
    unsigned long bond;

    // Separation of bonds into short and long ones.
    double        short_bond_weight[MAX_NUMBER_OF_SHORT_BONDS]; // label of "short bonds"


    double bond_len, bond_mod;
    double vctr_coord[space_dimension+1];

    double  total_weight,psum,rweight;

    // Parameters for creation of the bond length probability distribution obtained via inversion method.
    // half_lattice_length appears because of the distribution being finite. Else, the integrals would 
    // not converge in case of exponent being >= (-1).
    double dexponent = 1.0/( (double)space_dimension - (double)alpha);
    double lambda = pow( (1.0*half_lattice_length + 1.0), (double)space_dimension - (double)alpha )  - 1.0;
    
    memset(N_neighbors,0,sizeof(N_neighbors));

    //
    // The AVERAGE coordination is 2*bonds per site. Individual sites have a
    // Poisson-distributed number of links (for large systems).

    total_bonds=(unsigned long) (gamma * number_of_sites / 2.0);
    present_bonds=0;


    // Now connect a total of total_bonds bonds
    for(bond=0;bond<total_bonds;bond++){ 

	kounter=(unsigned long)0;
	do{ // this loop is repeated until a pair of neighbors is found that
	    // can be connected.
	    kounter++;

	    // NOTE: A new starting site is chosen at each trial 
	    isite=DRANDOM_NUMBER*number_of_sites;

	    // Now choose the OTHER site.
	    if(DRANDOM_NUMBER<p_s){ // use method for short bonds
		// Loop over short bonds and pick one with probability
		// short_bond_weight[]
		total_weight=(double)0;
		for(b_label=0;b_label<number_of_short_bonds;b_label++){
//		    printf("Short bond %ld:",short_bond_label[b_label]);
		    ksite=(isite+short_bond_label[b_label]+number_of_sites)%number_of_sites;
		    // Ignore if this site is already connected to isite. 
		    for(inbr=0;inbr<N_neighbors[isite];inbr++)
			if(neighbor[isite][inbr]==ksite)break;
		    if(inbr<(long)N_neighbors[isite]){
			// sites are connected. Ignore
//			printf("\tis connected. (discarded) \tW=0\n");
			short_bond_weight[b_label]=(double)0;
		    }else{ // site is elligible. 
			short_bond_weight[b_label]=bond_weight[short_bond_label[b_label]];
//			printf("\tis elligible for linking. \tW=%e\n",short_bond_weight[b_label]);
		    }
		    total_weight+=short_bond_weight[b_label];
		} // end loop over short bonds
//		printf("total weight of linkable bonds: %12.4e\n",total_weight);

		if(total_weight>(double)0){ 
		    // Normalize and choose one
		    rweight=DRANDOM_NUMBER;
		    psum=(double)0;
		    for(b_label=0;b_label<number_of_short_bonds;b_label++){
			short_bond_weight[b_label]/=total_weight;
			psum+=short_bond_weight[b_label];
			if(rweight<psum){
			    ksite=(isite+short_bond_label[b_label]+number_of_sites)%number_of_sites;
			    break;
			}
		    }
		    assert(b_label<number_of_short_bonds); 
		    break;
		}else{// no linkable site found
//		    printf("no linkable site found - continuing... \n");
		    continue; // jump to end of loop - a new isite will be chosen
		}	


	    }else{ // use method for long bonds
		// The LONG-BOND method consists of the folowing steps: 1) picking
		// a random distance, 2) picking a random direction, and then
		// choosing the site closest to the point so generated.
		
		// 1) Generate bond_len, a real random variable between 1 and
		// L/2, with PowerLaw distribution.
		do {
		    if((float)alpha==(float)space_dimension){ // special case: P ~ 1/r => Integral ~ ln(r)
			bond_len=(double) pow((double)half_lattice_length+1.0,DRANDOM_NUMBER);
		    }else {                                 // normal case: P ~ 1/r^kappa with kappa!=1
			bond_len=(double) pow(lambda*DRANDOM_NUMBER+1.0,dexponent);
		    }
		    
		    // 2) Generate random d-dimensional versor
		    do{
			bond_mod=(double)0;
			for(direction=0;direction<(unsigned long)space_dimension;direction++){
			    vctr_coord[direction]=2*DRANDOM_NUMBER-(double)1;
			    bond_mod+=vctr_coord[direction]*vctr_coord[direction];
			}
		    }while (bond_mod>1.0);    // Within unit sphere? - If not, discard!
		    bond_mod=sqrt(bond_mod);
		    
		    // Build vector with random direction and length equal to bond_mod
		    for(direction=0;direction<(unsigned long)space_dimension;direction++){
			vctr_coord[direction]=bond_len*vctr_coord[direction]/bond_mod;
		    } 
		    // Find site closest to this point. Recalculate length.
		    bond_mod=(double)0;
		    ksite=(long)0;
		    for(direction=0;direction<space_dimension;direction++){
			ksite+=(long)nearbyintl(vctr_coord[direction])*jump[direction];
			bond_mod+=nearbyintl(vctr_coord[direction])*nearbyintl(vctr_coord[direction]);
		    } 
		}while (bond_len<SHORT_BOND_LENGTH*SHORT_BOND_LENGTH); // make sure that bond is LONG
		ksite=(isite+ksite+number_of_sites)%number_of_sites;
		break;
	    } // End of long-bond method
//	    break;
	}while(1); // repeat until a break is executed

	
/* 	printf("For isite=%ld with %ld neighbors, chosen ksite=%ld with %ld neighbors\n", */
/* 	       isite,N_neighbors[isite],ksite,N_neighbors[ksite]); */

#ifdef ENABLE_ASSERTS
	assert((isite<number_of_sites)&&(isite>=0));
	assert((ksite<number_of_sites)&&(isite>=0));
#endif

	/* Final step: Modify neighbor list to include bond ij and count new bond */
	neighbor[isite][N_neighbors[isite]]=ksite;
	neighbor[ksite][N_neighbors[ksite]]=isite;
	N_neighbors[isite]++;
	N_neighbors[ksite]++;
	
#ifdef ENABLE_ASSERTS	
	assert(N_neighbors[isite]<MAX_NEIGHBORS);
	assert(N_neighbors[ksite]<MAX_NEIGHBORS);
#endif
	
	ivalue[bond]=isite;
	kvalue[bond]=ksite;
	present_bonds++;

    } /* END BOND LOOP */
    
    return;
}

  


/* ---------------------------------------------------------------*\
   are_neighbors: checks, whether j is neighbor of i
\* ---------------------------------------------------------------*/
int are_neighbors(long i,long j){
#undef HERE
#define HERE "are_neighbors:"

   unsigned long n;
   int answer = 0;
   
   for(n=0;n<=N_neighbors[i];n++) { if(neighbor[i][n]==j) answer=1; };
   return answer;
   
}



/* ---------------------------------------------------------------*\
   average_neighbors:   Calculate the average neighbors per site
                        and output it
\* ---------------------------------------------------------------*/
void average_neighbors(){
#undef HERE
#define HERE "average_neighbors:"

	int i, imax;
	unsigned long tot_N_neighbors=0;
	double avg_N_neighbors;
	imax=number_of_sites;
	
	for(i=0;i<imax;i++) tot_N_neighbors+=N_neighbors[i];
	avg_N_neighbors=1.0*tot_N_neighbors/imax;
	
//	TELL(" ");
//	printf("Average neighbors per site = %f.\n", avg_N_neighbors);
	
	return;
}



/* ---------------------------------------------------------------*\
   
   Plot links
   
\* ---------------------------------------------------------------*/

void plot_links(float maxg){
#undef HERE
#define HERE "plot_links:"

    FILE  *out0_ptr;
    FOPEN(out0_ptr,"lattice.links");
    
    unsigned long bond,sa,sb,ia,ib,ja,jb;
    unsigned long i;

    for (bond=0;bond<present_bonds;bond++){

	sa=ivalue[bond]; 
	sb=kvalue[bond]; 
	
	ia=(sa%lattice_length);
	ja=(sa/lattice_length);
	
	ib=(sb%lattice_length);
	jb=(sb/lattice_length);

	if((unsigned long) sqrt((ia-ib)*(ia-ib)+(ja-jb)*(ja-jb))<half_lattice_length)
	    fprintf(out0_ptr,"\n%ld %ld\n%ld %ld\n",ia,ja,ib,jb);
	
    } 
    fclose(out0_ptr);



    FOPEN(out0_ptr,"lattice.strong.sites");
    for(i=0;i<number_of_sites;i++){
	if (N_neighbors[i]>=mcull){
	    ia=(i%lattice_length);
	    ja=(i/lattice_length);
	    fprintf(out0_ptr,"%ld %ld %ld\n",ia,ja,N_neighbors[i]);
	}
    }
    fclose(out0_ptr);




    return;
}

/*---------------------------------------------------------------------------*\
 
  Calculates all bond weights, and in particular the weights and identities of
  all short bonds, i.e. those with length < SHORT_BOND_LENGTH.

\*----------------------------------------------------------------------------*/
void calculate_bond_weights(void){
#undef HERE
#define HERE "calculate_bond_weights:"

    unsigned long isite,direction;
    unsigned long site_label,center_label,new_label;


    double bond_len;
    double  total_sweight;

    //
    // Calculate bond weights:
    // ----------------------
    //
    // Bond weights are defined as 1/(len)^{alpha}, where len is the length of
    // the bond, i.e. its distance to the origin. We are using helical BC's
    // here so there is a slight complication in order to define bond
    // lenghts. For a given site we will define its distance to the origin as
    // the minimum distance among all possible images of this site. This in
    // fact makes all distances smaller than approximately sqrt(d)*L/2.
    // 
    // In order to enforce this condition, the origin is defined to be at
    // {L/2,L/2,...,L/2} and bond lengths are calculated as distances to this
    // point (whose label is center_label).  

    // Calculate center_label
    center_label=(unsigned long)0;
    for(direction=0;direction<(unsigned long)space_dimension;direction++)
	center_label+=half_lattice_length*jump[direction];
//    printf("Center label=%ld\n",center_label);

    // Calculate bond weights bond_weight[], and short-bond weight p_s
    total_sweight=(double)0;
    bond_weight[0]=(double)0;
    p_s=(double)0;
    number_of_short_bonds=(unsigned long)0;
    for (isite=0;isite<number_of_sites;isite++){// loop over all possible links
	site_label=isite;
	// Label with respect to central point at {L/2,L/2,...,L/2}
	new_label=(isite-center_label+number_of_sites)%(number_of_sites);
	bond_len=(double)0;
//	printf("Site %4ld : (new label %4ld)  (",isite,new_label);
	// calculate d cartesian coordinates
	for(direction=0;direction<(unsigned long)space_dimension;direction++){
	    coordinate[direction]=(site_label)%(lattice_length);
	    site_label/=lattice_length;
//	    printf("%4ld ",coordinate[direction]);
	    bond_len+=(coordinate[direction]-half_lattice_length)*(coordinate[direction]-half_lattice_length);
	}
	if(new_label!=0){// ignore central site
	    bond_weight[new_label]=pow(bond_len,-alpha/2);
	    // If bond is "short", add its weight to that of short bonds, and store its identity.
	    if(bond_len < SHORT_BOND_LENGTH*SHORT_BOND_LENGTH ){
		p_s+=bond_weight[new_label];
		short_bond_label[number_of_short_bonds]=new_label;
		number_of_short_bonds++;
	    }
	    total_sweight+=bond_weight[new_label];
	}
//	printf(") L=%lf   W=%16.6e\n",sqrt(bond_len),bond_weight[new_label]);
    } // end of loop over isites (i.e. over all possible bonds stemming from site at center_label)


    // Normalize bond weights
    for (isite=0;isite<number_of_sites;isite++){
	bond_weight[isite]/=total_sweight;
    }
    p_s/=total_sweight;


/*
    for(isite=0;isite<number_of_short_bonds;isite++){ 
	printf("%ld: short bond %ld, W=%12.4e\n",
	       isite,short_bond_label[isite],bond_weight[short_bond_label[isite]]);
    }
*/

    printf("There are %ld SHORT bonds (len < %lf). Their relative weight is %16.6e\n\n",
	   number_of_short_bonds,SHORT_BOND_LENGTH,p_s);
    return;
}

  


/*--------------------------------------------------------------------*\
 +
 + unsigned long my_lrand( void ): A Kirkpatrick-Stoll shift-register RNG.
 +
\*--------------------------------------------------------------------*/
unsigned long my_lrand (void){
  
#undef  HERE
#define HERE "\nmy_lrand: "
  
  static unsigned long	cnt=0;
  unsigned long               rnd_word;
  int                 i0,i1;
  
  i0 = cnt%MODUL;
  i1 = (cnt+JJD)%MODUL;
  rnd_word = ibits[i0] = ibits[i0]^ibits[i1];
  cnt++;
  
  return(rnd_word);
}



/* ------------------------------------------------------------------------- *\
 |
 |	my_RNG_seed:
 |	
 |      fills in an array of random bits for the use of my_RNG.
 |
\* ------------------------------------------------------------------------- */

void my_RNG_seed_init(unsigned long seed){
  
  int	i0,bit;
  
  printf("Calling default seed initializer with seed=%lu.\n",
	 seed);
  /* Initialize default RNG */
  DEFAULT_SET_SEED_RNG(seed);
  
  printf("Initializing arrays for shift-register RNG.");
  
  for(i0=0;i0<MODUL;i0++){
    
    ibits[i0] = 0L;
    
    for(bit=0;bit<RNG_N_OF_BITS;bit++){            
      if( DEFAULT_DRANDOM_NUMBER  > 0.5 ) ibits[i0] |= (1L<<bit);      
    }
  }
  
  printf("\tShift-register arrays initialized.\n");
  
  return;
  
}


