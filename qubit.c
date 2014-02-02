#include "cpuminer-config.h"
#include "miner.h"

#include <string.h>
#include <stdint.h>

#include "x5/luffa_for_sse2.h" //sse2 opt
//--ch h----
#include "x5/cubehash_sse2.h" //sse2
//----------
#include "x5/sph_shavite.h"
//#include "x5/low-mem/SHA3api_ref.h"
//-----simd vect128---------
#include "x5/vect128/nist.h"
//---echo ----------------
#include "x5/sph_echo.h"


/* Move init out of loop, so init once externally, and then use one single memcpy with that bigger memory block */
typedef struct {
	sph_shavite512_context  shavite1;
	sph_echo512_context		echo1;
} qubithash_context_holder;

qubithash_context_holder base_contexts;
//--luffa var --
hashState base_context_luffa;
//--cubehash var --
cubehashParam base_context_cubehash;

void init_qubithash_contexts()
{
  //-- luffa init --
  init_luffa(&base_context_luffa,512);
  //--cubehash init--
  cubehashInit(&base_context_cubehash,512,16,32);
  //---------------
  sph_shavite512_init(&base_contexts.shavite1);
   //-------------------------------
  sph_echo512_init(&base_contexts.echo1);
}

static void qubithash(void *state, const void *input)
{
	qubithash_context_holder ctx;
	//-- local luffa var --
	hashState			 ctx_luffa;
	//---local cubehash var ---
	cubehashParam		 ctx_cubehash;
	//---local simd var ---
	hashState_sd *	     ctx_simd1;
	
    uint32_t hashA[16], hashB[16];	
	
	memcpy(&ctx, &base_contexts, sizeof(base_contexts));
	//--luffa copy date --
	memcpy(&ctx_luffa,&base_context_luffa,sizeof(hashState));
	//--cubehash copy date----
	memcpy(&ctx_cubehash,&base_context_cubehash,sizeof(cubehashParam));
		
	//-------luffa sse2--------
   update_luffa(&ctx_luffa,(const BitSequence *)input,640);
   final_luffa(&ctx_luffa,(BitSequence *)hashA);	
    //---cubehash sse2---    
	cubehashUpdate(&ctx_cubehash,(const byte *)hashA,64);
	cubehashDigest(&ctx_cubehash,(byte *)hashB);
  //------shavite  ------	
    sph_shavite512 (&ctx.shavite1, hashB, 64);   
    sph_shavite512_close(&ctx.shavite1, hashA);  
 //Hash_sh(512,(const BitSequence *)hashB,512,(BitSequence *)hashA);
//-------simd512 vect128 --------------	
	ctx_simd1=malloc(sizeof(hashState_sd));
	Init(ctx_simd1,512);
	Update(ctx_simd1,(const BitSequence *)hashA,512);
	Final(ctx_simd1,(BitSequence *)hashB);
	free(ctx_simd1->buffer);
    free(ctx_simd1->A);
	free(ctx_simd1);
//-----------------	
	sph_echo512 (&ctx.echo1, hashB, 64);   
    sph_echo512_close(&ctx.echo1, hashA); 


	memcpy(state, hashA, 32);
	
}

int scanhash_qubit(int thr_id, uint32_t *pdata, const uint32_t *ptarget,
	uint32_t max_nonce, unsigned long *hashes_done)
{
	uint32_t n = pdata[19] - 1;
	const uint32_t first_nonce = pdata[19];
	const uint32_t Htarg = ptarget[7];

	uint32_t hash64[8] __attribute__((aligned(32)));
	uint32_t endiandata[32];
	
	
	int kk=0;
	for (; kk < 32; kk++)
	{
		be32enc(&endiandata[kk], ((uint32_t*)pdata)[kk]);
	};

	if (ptarget[7]==0) {
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (((hash64[7]&0xFFFFFFFF)==0) && 
					fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	
	} 
	else if (ptarget[7]<=0xF) 
	{
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (((hash64[7]&0xFFFFFFF0)==0) && 
					fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	
	} 
	else if (ptarget[7]<=0xFF) 
	{
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (((hash64[7]&0xFFFFFF00)==0) && 
					fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	
	} 
	else if (ptarget[7]<=0xFFF) 
	{
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (((hash64[7]&0xFFFFF000)==0) && 
					fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	

	} 
	else if (ptarget[7]<=0xFFFF) 
	{
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (((hash64[7]&0xFFFF0000)==0) && 
					fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	

	} 
	else 
	{
		do {
			pdata[19] = ++n;
			be32enc(&endiandata[19], n); 
			qubithash(hash64, &endiandata);
			if (fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		} while (n < max_nonce && !work_restart[thr_id].restart);	
	}
	
	
	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;
	return 0;
}