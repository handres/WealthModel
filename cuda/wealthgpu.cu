/*----------------------------------------------------------------------------*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
#include <cuda.h>
#include <curand_kernel.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;
using namespace thrust::placeholders;

/*----------------------------------------------------------------------------*/
#define EPS 1e-5
#define ISSMALL(_e) (-EPS<=(_e) && (_e)<= EPS)
#define EQUAL(_a,_b) (((_a)+(_b)==0)?true:ISSMALL(((_a)-(_b))/((_a)+(_b))))
/*----------------------------------------------------------------------------*/
typedef float WTYPE;
#define WTSZ sizeof(WTYPE)

/*----------------------------------------------------------------------------*/
struct MyTimer { cudaEvent_t start, stop; float dt; };
#define TIMER_TYPE MyTimer
#define TIMER_START(_x) do{cudaEventCreate(&(_x).start);cudaEventCreate(&(_x).stop);cudaEventRecord((_x).start,0);}while(0)
#define TIMER_END(_x) do{cudaEventRecord((_x).stop,0);cudaEventSynchronize((_x).stop); cudaEventElapsedTime(&((_x).dt), (_x).start, (_x).stop);}while(0)
#define TIMER_SECS(_x) (((_x).dt)*1e-3)

/*----------------------------------------------------------------------------*/
static int dbglvl = 1;
#define DBG(_l,_s) if((_l)<=dbglvl){cout<<_s;}else
#define DBGNL(_l,_s) DBG(_l,_s<<endl)

/*----------------------------------------------------------------------------*/
bool normalize = false, verify = true, writeresults = true;

/*----------------------------------------------------------------------------*/
WTYPE initw( WTYPE *w, long n, WTYPE iw )
{
    WTYPE wsum = 0;
    srand48(345);
    for( long i = 0; i < n; i++ )
    {
        w[i] = (iw >= 0 ? iw : (-iw)*drand48());
        wsum += w[i];
    }
    return wsum;
}

/*----------------------------------------------------------------------------*/
void printw( int d, const char *s, WTYPE *w, long n )
{
    DBG(d,"** "<<s<<": ");
    for( long i = 0; i < n; i++ )
    {
        DBG(d,(i>0?", ":"")<<w[i]);
        if( i >= 4 && n > 10 && i < n-6 )
        {
            DBG(d,", ...");
            i = n-6;
        }
    }
    DBGNL(d,"");
}

/*----------------------------------------------------------------------------*/
void writew( WTYPE *w, long n, long ts )
{
    char fname[100];
    sprintf( fname, "wealth-npersons-%ld-at-time-%ld.bin", n, ts );
    DBGNL(0,"Writing wealth distribution to file \""<<fname<<"\"");
    ofstream fout;
    fout.open( fname, ios::binary | ios::out );
    fout.write( (char*) w, n*sizeof(WTYPE) );
    fout.close();
    DBGNL(2,"Writing done");
}

/*----------------------------------------------------------------------------*/
long nthr = 16*32, nblk = 32;
long totthr = 0;
curandState *randstate = 0;

/*----------------------------------------------------------------------------*/
#define randnext( _rs, _i ) curand_uniform( &(_rs)[(_i)] )

/*----------------------------------------------------------------------------*/
__global__
void k_enrich( WTYPE *w, WTYPE *wdf,
               long n, long npt,
               WTYPE sumbeta, WTYPE dw,
               long ndw, curandState *rs )
{
    long nthr = (gridDim.x * blockDim.x);
    long tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    long nr = ndw/nthr;
    if( tid < ndw%nthr ) nr++;

    for( long r = 0; r < nr; r++ )
    {
      float rng = randnext( rs, tid );
      WTYPE val = rng*sumbeta;
      long si = 0, ei = n-1;
      //long si = tid*npt, ei = si+npt-1;
      if( ei >= n )
      {
        ei = n-1;
      }
      if( si < n )
      {
        long mi = -1;
        if( val < wdf[si] )
        {
            WTYPE prev = (si > 0 ? wdf[si-1] : 0.0);
            if( prev < val )
            {
                mi = si;
            }
        }
        else if( val <= wdf[ei] )
        {
            mi = si+(ei-si)/2;
            while( si < ei )
            {
                if( wdf[mi] == val )
                {
                    break;
                }
                else if( val < wdf[mi] )
                {
                    ei = mi-1;
                }
                else
                {
                    si = mi+1;
                }
                mi = si+(ei-si)/2;
            }
        }
        if( mi >= 0 )
        {
            atomicAdd( &w[mi], dw );
        }
      }
    }
}

/*----------------------------------------------------------------------------*/
void enrich( WTYPE *w, WTYPE *wdf, long n, WTYPE sumw2, WTYPE dw, long ndw )
{
    long npt = (n <= totthr) ? 1 : (1+(n-1)/totthr);
    DBGNL(2,"npt="<<npt);
    DBGNL(1,"val="<<sumw2);
    k_enrich<<<nblk,nthr>>>( w, wdf, n, npt, sumw2, dw, ndw, randstate );
}

/*----------------------------------------------------------------------------*/
__global__
void k_initrand( curandState *randstate )
{
    long tid = (blockIdx.x * blockDim.x) + threadIdx.x;
    long rseed = 1234;
    curand_init( rseed, tid, 0, &randstate[tid] );
}

/*----------------------------------------------------------------------------*/
void initrand( void )
{
    totthr = nblk*nthr;
    cudaMalloc( &randstate, totthr*sizeof(randstate[0]) );
    k_initrand<<<nblk,nthr>>>( randstate );
}

/*----------------------------------------------------------------------------*/
class PowerFunctor
{
    WTYPE e;
    public: PowerFunctor(WTYPE _e){e=_e;}
    public: __host__ __device__ WTYPE operator()(WTYPE x)const{return pow(x,e);}
};

/*----------------------------------------------------------------------------*/
int main( int ac, char *av[] )
{
    long np = 10, niter = 10, ndw = 1024;
    float beta = 1.36;
    WTYPE iw = 1.0, dw = 0.1;
    if(getenv("VERIFY"))verify=!strcmp(getenv("VERIFY"),"TRUE");
    if(getenv("NORMALIZE"))normalize=!strcmp(getenv("NORMALIZE"),"TRUE");
    if(getenv("WRITERESULTS"))writeresults=!strcmp(getenv("WRITERESULTS"),"TRUE");
    if(getenv("DEBUG"))dbglvl=atoi(getenv("DEBUG"));
    if(getenv("BETA"))beta=atof(getenv("BETA"));
    if( ac <= 1 )
    {
        DBGNL(0,"Usage: "<<av[0]<<" npersons niterations initwealth deltawealth ndw beta");
        DBGNL(0,"All persons have same given init wealth if initwealth > 0 ");
        DBGNL(0,"else (if initwealth < 0), have random wealth with mean=-initwealth");
        exit(1);
    }
    if( ac >= 2 ) np = atoi( av[1] );
    if( ac >= 3 ) niter = atoi( av[2] );
    if( ac >= 4 ) iw = atof( av[3] );
    if( ac >= 5 ) dw = atof( av[4] );
    if( ac >= 6 ) ndw = atoi( av[5] );
    if( ac >= 7 ) beta = atof( av[6] );

    DBGNL(0,"dbglvl="<<dbglvl);
    DBGNL(0,"verify="<<verify);
    DBGNL(0,"normalize="<<normalize);
    DBGNL(0,"writeresults="<<writeresults);
    DBGNL(0,"npersons="<<np<<" niterations="<<niter<<" iw="<<iw<<" dw="<<dw<<" beta="<<beta);

    long nbytes = np*WTSZ;
    WTYPE *h_w = (WTYPE*)malloc( nbytes ); DBGNL(1,"mallocated w"<<nbytes);
    WTYPE *h_wdf = (WTYPE*)malloc( nbytes ); DBGNL(2,"mallocated wdf"<<nbytes);
    WTYPE h_totw = initw( h_w, np, iw );
    printw( 0, "Initialized", h_w, np );

    initrand(); DBGNL(0,"RNG initialized");
    WTYPE *d_w = 0, *d_wdf = 0;
    cudaMalloc( &d_w, nbytes ); DBGNL(2,"w cudaMallocated "<<nbytes);
    cudaMalloc( &d_wdf, nbytes ); DBGNL(3,"wdf cudaMallocated "<<nbytes);
    cudaMemcpy( d_w, h_w, nbytes, cudaMemcpyHostToDevice ); DBGNL(3,"H2D "<<nbytes);
    thrust::device_ptr<WTYPE> td_w(d_w), td_wdf(d_wdf);

    #define GPW(_d,_s) if((_d)<=dbglvl){ \
        cudaMemcpy( h_w, d_w, nbytes, cudaMemcpyDeviceToHost ); DBGNL(10,"D2H W"<<nbytes); \
        printw( _d, _s, h_w, np ); \
      }else
    #define GPWDF(_d,_s) if((_d)<=dbglvl){ \
        cudaMemcpy( h_wdf, d_wdf, nbytes, cudaMemcpyDeviceToHost ); DBGNL(10,"D2H WDF"<<nbytes); \
        printw( _d, _s, h_wdf, np ); \
      }else

    long nprints = getenv("NPRINTS") ? atoi(getenv("NPRINTS")) : 100;
    long printevery = niter<nprints?1:(niter/nprints);
    DBGNL(0,"nprints="<<nprints);
    long nwrites = getenv("NWRITES") ? atoi(getenv("NWRITES")) : 10;
    long writeevery = niter<nwrites?1:(niter/nwrites);
    DBGNL(0,"nwrites="<<nwrites);

    DBGNL(0,"Starting simulation...");
    TIMER_TYPE timer;
    TIMER_START(timer);
    thrust::sort( td_w, td_w+np ); GPW(3,"Sorted");
    for( long iter = 0; iter < niter; iter++ )
    {
        if((iter+1)%printevery==0)DBGNL(2,"------ Iteration "<<iter<<" Start ------");
        thrust::transform( td_w, td_w+np, td_wdf, PowerFunctor(beta) ); GPWDF(4,"Exponentiated");
        WTYPE h_totwdf = thrust::reduce( td_wdf, td_wdf+np ); DBGNL(1,"Totwdf="<<h_totwdf);
        if( normalize )
        {
            WTYPE h_totwdf_inv = ( (h_totwdf == 0) ? 0.0 : (1.0 / h_totwdf));
            thrust::transform( td_wdf, td_wdf+np, td_wdf, _1*h_totwdf_inv ); GPWDF(5,"Normalized");
            if(verify)
            {
                WTYPE hsum = thrust::reduce( td_wdf, td_wdf+np ); DBGNL(1,"Totwdfnorm="<<hsum);
                if( !EQUAL(hsum,1.0) ) { cerr<<"SUM MISMATCH "<<(hsum-1.0)<<endl; exit(1); }
            }
        }
        thrust::inclusive_scan( td_wdf, td_wdf+np, td_wdf ); GPWDF(6,"Scanned");
        enrich( d_w, d_wdf, np, h_totwdf, dw, ndw ); GPW(2,"Enriched");
        thrust::sort( td_w, td_w+np ); GPW(3,"Sorted");
        if( verify || iter==niter-1 )
        {
            WTYPE hsum = thrust::reduce( td_w, td_w+np ); DBGNL(1,"Totw="<<hsum);
            WTYPE isum = h_totw + (iter+1)*ndw*dw;
            if( !EQUAL(hsum,isum) ) { cerr<<"SUM MISMATCH "<<hsum<<" != "<<isum<<" ("<<(hsum-isum)<<")"<<endl; if(verify)exit(1); }
            DBGNL(0,"Total wealth="<<isum);
        }
        if( (iter+1)%writeevery==0 || iter==niter-1 )
        {
            if( writeresults )
            {
                GPW(0,"Writing");
                writew( h_w, np, iter+1 );
            }
        }
        if((iter+1)%printevery==0)DBGNL(0,"------ Iteration "<<iter<<" ("<<((iter+1)*100.0/niter)<<"%) End ------");
    }
    TIMER_END(timer);

    float dt = TIMER_SECS(timer);
    float million = 1e6;
    float thousand = 1e3;
    DBGNL(0,"           Iterations = "<<niter/million<<" million");
    DBGNL(0,"              Persons = "<<np/million<<" million");
    DBGNL(0," #Wealth packets/Iter = "<<ndw/thousand<<" thousand");
    DBGNL(0,"Total #wealth packets = "<<ndw*niter/million<<" million");
    DBGNL(0,"                 Time = "<<dt<<" secs");
    DBGNL(0,"       Time/Iteration = "<<(dt/niter)*1e3<<" millisecs");
    DBGNL(0,"Time/Iteration/Person = "<<(dt/niter/np)*1e6<<" microsecs");

    GPW(0,"Final enrichment");

    DBGNL(0,"End");
}

/*----------------------------------------------------------------------------*/
