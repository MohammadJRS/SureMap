#include "CompressedString.h"
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
using namespace std;

#define size_uchar 200
#define MAXALPH 10

class BWT {
public:
	char dna[4];
	
	int mode;
	uint32_t zeroPos;
	uint8_t remap[size_uchar];
	uint32_t sigma;
	uint8_t* remap_reverse;
	uint32_t n;
	uint32_t C[size_uchar+1];
	uint32_t* suffixes;
	BitArray bit[MAXALPH];
	CompressedString bwt;
	
	BWT(){ zeroPos = sigma = n = 0; dna[0] = 'a', dna[1] = 'c', dna[2] = 'g', dna[3] = 't';
		remap_reverse = NULL;
		suffixes = NULL;
	}
	
    BWT(uint8_t* T,uint32_t _n, int _m ){
		mode = _m;
		n = _n;	
		build(T,n);
		dna[0] = 'a', dna[1] = 'c', dna[2] = 'g', dna[3] = 't';
	}
	
	uint8_t* remap0(uint8_t* T,uint32_t n) {
		uint8_t* X;
		uint32_t i,j,size=0;
		uint32_t freqs[size_uchar];
		
		for(i=0;i<size_uchar;i++) freqs[i]=0;
		for(i=0;i<n;i++) if(freqs[T[i]]++==0) size++;
		
		sigma=size;
		
		// remap alphabet
		if (freqs[0]>1) {i=1;sigma++;} //test if some character of T is zero, we already know that text[n-1]='\0'
		else i=0;

		remap_reverse = (uint8_t*) malloc(size*sizeof(uint8_t));
		for(j=0;j<size_uchar;j++) {
		  if(freqs[j]!=0) {
			remap[j]=i;
			remap_reverse[i++]=j;
		  }
		}
		// remap text
		X = new uint8_t[n];
		for(i=0;i<n-1;i++) // the last character must be zero
		  X[i]=remap[T[i]];
		X[n-1] = 0;    
		return X;
	}
	
    void build(uint8_t* T,uint32_t n){
		uint8_t* X;
		unsigned int* SA;
		uint32_t prev,tmp,start,stop;
		float elapsed;
		X = remap0(T,n);
		for( int i = 0; i < sigma; i++ ){
			uint8_t ch = remap_reverse[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' ){
				bit[i].reset( n );
			}
		}
		unsigned int* xx = new unsigned int[n + 100];
		for( unsigned int i = 0; i < n; i++ )
			xx[i] = X[i];
		xx[n] = xx[n + 1] = xx[n + 2] = 0;
		for (unsigned int i=0;i<size_uchar+1;i++) C[i]=0;
		for (unsigned int i=0;i<n;++i) C[X[i]]++;
		prev=C[0];C[0]=0;
		for (unsigned int i=1;i<size_uchar+1;i++) {
			tmp = C[i];
			C[i]=C[i-1]+prev;
			prev = tmp;
		}
		SA = new unsigned int[n + 100];
		string saAdr = ( mode == 0 ) ? "/home/javad/HgSuffixArray/hgSuff.bin" : "/home/javad/HgSuffixArray/hgInvSuff.bin";
		FILE* fin = fopen( saAdr.c_str(), "rb" );
		const int maxR = 100000;
		unsigned int buffer[maxR];
		unsigned int idx = 0;
		while( true ){
			int rd = fread( buffer, sizeof( unsigned int ), maxR, fin );
			for( int i = 0; i < rd; i++ ){
				SA[idx++] = buffer[i];
			}
			if( rd == 0 )
				break;
		}
		for( unsigned int i = 0; i < n; i++ ){
			if( SA[i] == 0 )
				T[i] = 0;
			else T[i] = remap_reverse[X[SA[i] - 1]];
			if( T[i] == 0 )
				zeroPos = i;
		}
		for( uint32_t i = 0; i < n; i++ ){
			uint8_t ch = T[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' ){
				bit[remap[T[i]]].setBit( i );
			}
		}
		for( int i = 0; i < sigma; i++ ){
			uint8_t ch = remap_reverse[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' ){
				bit[i].setSum();
			}
		}
		bwt.set( T, n );
		delete[] T;
		suffixes = new uint32_t[((n + SA_SAMPLERATE - 1)/SA_SAMPLERATE)];
		for(unsigned int i=0;i<n;i++) {
			if( i % SA_SAMPLERATE == 0)
				suffixes[i / SA_SAMPLERATE] = SA[i];
		}
	}
	
	
    void load(string filename){
		FILE* fin = fopen( filename.c_str(), "rb" );
		fread( remap, sizeof( uint8_t ), size_uchar, fin );
		fread( &sigma, sizeof( uint32_t ), 1, fin );
		remap_reverse = new uint8_t[sigma];
		fread( remap_reverse, sizeof( uint8_t ), sigma, fin );
		fread( &n, sizeof( uint32_t ), 1, fin );
		fread( C, sizeof( uint32_t ), size_uchar, fin );
		suffixes = new uint32_t[( n + SA_SAMPLERATE - 1 ) / SA_SAMPLERATE];
		fread( suffixes, sizeof( uint32_t ), ( n + SA_SAMPLERATE - 1 ) / SA_SAMPLERATE, fin );
		fread( &zeroPos, sizeof( uint32_t ), 1, fin );
		for( int i = 0; i < MAXALPH; i++ )
			bit[i].load( fin );
		bwt.load( fin );
		fclose( fin );
	}
	
    int32_t save(string filename){
		FILE* fout = fopen( filename.c_str(), "wb" );
		for( int i = 0; i < size_uchar; i++ )
			fwrite((const void*)& remap[i],sizeof(uint8_t),1,fout);
		fwrite((const void*)& sigma,sizeof(uint32_t),1,fout);
		for( int i = 0; i < sigma; i++ )
			fwrite((const void*)& remap_reverse[i],sizeof(uint8_t),1,fout);
		fwrite((const void*)& n,sizeof(uint32_t),1,fout);
		for( int i = 0; i < size_uchar; i++ )
			fwrite((const void*)& C[i],sizeof(uint32_t),1,fout);
		for( int i = 0; i < ( n + SA_SAMPLERATE - 1 ) / SA_SAMPLERATE; i++ )
			fwrite((const void*)& suffixes[i],sizeof(uint32_t),1,fout);
		fwrite((const void*)& zeroPos,sizeof(uint32_t),1,fout);
		for( int i = 0; i < MAXALPH; i++ )
			bit[i].save( fout );
		bwt.save( fout );
		fclose( fout );
	}
	
	inline uint8_t charAt( uint32_t pos ){
		if( pos == zeroPos )
			return 0;
		if( bwt.sz )
			return bwt.charAt( pos );
		for( int i = 0; i < 4; i++ ){
			if( bit[remap[dna[i]]].getPos( pos ) )
				return dna[i];
		}
		return 'n';
	}
	
	uint32_t locateRow( uint32_t idx ){
		int dist = 0;
		//while( idx % SA_SAMPLERATE != 0 ){
		while( idx & ( SA_SAMPLERATE - 1 ) ){
			char ch = charAt( idx );
			if( ch == 0 ){
				return dist;
			}
			int id = remap[ch];
			idx = C[id] + bit[id].getRank( idx ) - 1;
			dist++;
		}
		uint32_t ret = ( suffixes[idx / SA_SAMPLERATE] + dist );
		if( ret >= n )
			ret %= n;
		return ret;
	}
	
	void updateInterval( uint32_t& l, uint32_t& r, uint8_t c ){
		if( l == 0 )
			l = C[c];
		else{
			l = C[c] + bit[c].getRank( l - 1 );
		}
		r = C[c] + bit[c].getRank( r ) - 1;
	}
};

