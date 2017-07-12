/**
 * @file    BWT.h
 * @author  MohammadJavad Rezaei Seraji <mjrezaei (at) ce.sharif.edu>
 *
 * @section LICENCE
 *
 *
 * Copyright (C) 2017-2020
 *   MohammadJavad Rezaei Seraji <mjrezaei (at) ce.sharif.edu>
 *	 Seyed Abolfazl Motahari <motahari (at) sharif.edu
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>
#include "psascan_src/psascan.h"
#include "CompressedString.h"

//using namespace seqan;
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
	unsigned char* X;
	
	BWT(){ zeroPos = sigma = n = 0; dna[0] = 'a', dna[1] = 'c', dna[2] = 'g', dna[3] = 't';
		remap_reverse = NULL;
	}
	
    BWT(uint8_t* T,uint32_t _n, int _m ){
		mode = _m;
		n = _n;	
		build(T,n);
		dna[0] = 'a', dna[1] = 'c', dna[2] = 'g', dna[3] = 't';
	}
	
	void remap0(uint8_t* T,uint32_t n) {
		//uint8_t* X;
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
		//resize( X, n );
		for(i=0;i<n-1;i++) // the last character must be zero
		  X[i]=remap[T[i]];
		X[n-1] = 0;
		//cout << "HERE " << (int)remap['n'] << endl;
	}
	
    void build(uint8_t* T,uint32_t n){
		uint32_t prev,tmp,start,stop;
		float elapsed;
		remap0(T,n);
		for( int i = 0; i < sigma; i++ ){
			uint8_t ch = remap_reverse[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' || ch == 'n' ){
				bit[i].reset( n );
			}
		}
		/*unsigned int* xx = new int[n + 100];
		for( unsigned int i = 0; i < n; i++ ){
			xx[i] = X[i];
		}
		xx[n] = xx[n + 1] = xx[n + 2] = 0;
		*/
		for (unsigned int i=0;i<size_uchar+1;i++) C[i]=0;
		for (unsigned int i=0;i<n;++i) C[X[i]]++;
		prev=C[0];C[0]=0;
		for (unsigned int i=1;i<size_uchar+1;i++) {
			tmp = C[i];
			C[i]=C[i-1]+prev;
			prev = tmp;
		}
		string adr = "temp.ref";
		string pAdr = "sa_" + adr + ".sdsl";
		ofstream fout( adr );
		for( unsigned int i = 0; i < n; i++ )
			fout << (int)X[i];
		//fout << '0';
		fout.flush();
		fout.close();
		
		/*csa_bitcompressed<> SA;
		construct(SA, adr, 1);
		*/
		
		long ram_use = 3072L << 20;
		std::string text_fname = adr;
		std::string out_fname("");
		std::string gap_fname("");

		// Parse the text filename.
		
		out_fname = text_fname + ".sa5";
		gap_fname = out_fname;
		// Find the number of (logical) cores on the machine.
		long max_threads = (long)omp_get_max_threads();
		cerr << "Building BWT and fmidx ...\n";
		// Run pSAscan.
		pSAscan(text_fname, out_fname, gap_fname,
		  ram_use, max_threads, true);
		FILE* fin = fopen( out_fname.c_str(), "rb" );
		suffixes = new uint32_t[((n + SA_SAMPLERATE - 1)/SA_SAMPLERATE)];
		int bufSize = (1LL << 20) - 1;
		unsigned char* buff = new unsigned char[bufSize];
		long long curIdx = 0;
		while( true ){
			int readBytes = fread( buff, sizeof( uint8_t ), bufSize, fin );
			if( readBytes == 0 )
				break;
			for( int i = 0; i < readBytes; i += 5 ){
				long long val = 0;
				for( int j = i + 4; j >= i; j-- ){
					val = ( val << 8 ) + buff[j];
					//cerr << (int)buff[j] << ' ';
				}
				//if( curIdx < 1000 )
				//	cerr << curIdx << ' ' << val << endl;
				//cerr << endl << curIdx << ' ' << val << endl;
				if( val == 0 )
					T[curIdx] = 0;
				else T[curIdx] = remap_reverse[X[val - 1]];
				if( T[curIdx] == 0 )
					zeroPos = curIdx;
				if( curIdx % SA_SAMPLERATE == 0)
					suffixes[curIdx / SA_SAMPLERATE] = val;
				curIdx++;
			}
		}
		cerr << "done\n";
		
		for( uint32_t i = 0; i < n; i++ ){
			uint8_t ch = T[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' || ch == 'n' ){
				bit[remap[T[i]]].setBit( i );
			}
		}
		
		for( int i = 0; i < sigma; i++ ){
			uint8_t ch = remap_reverse[i];
			if( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' || ch == 'n' ){
				bit[i].setSum();
			}
		}
		bwt.set( T, n );
		delete[] T;
		
		
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
		//cout << "TTTTest " << suffixes[0] << endl;
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
		//cout << charAt( 0 ) << endl;
		//cout << "start = " << idx << ' ' << charAt( idx ) << ' ' << bit[remap[charAt( idx )]].getRank( idx ) << endl;
		while( idx % SA_SAMPLERATE != 0 ){
			char ch = charAt( idx );
			
			if( ch == 0 ){
				return dist;
			}
			int id = remap[ch];
			//updateInterval( idx, idx, id );
			//cout << "fff " << idx << ' ' << ch << ' ' << id << ' ' << C[id] << ' ' << ( idx - 1 ) << endl;
			idx = C[id] + bit[id].getRank( idx ) - 1;
			//cout << "sss " << idx << ' ' << ch << endl;
			//idx = C[id] + bit[id].getRank( idx );
			//cout << "HERE " << id << ' ' << idx << ' ' << dist << ' ' << SA_SAMPLERATE << ' ' << ( idx % SA_SAMPLERATE ) << endl;
			dist++;
		}
		return ( suffixes[idx / SA_SAMPLERATE] + dist ); //% n;
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

