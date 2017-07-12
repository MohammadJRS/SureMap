#include <string>
#include <cstdlib>
#include <cstdio>
#include "BitArray.h"

using namespace std;

struct CompressedString{
	char dna[4];
	BitArray invalidPos;
	uint32_t sz;
	uint8_t* encodedStr;
	uint32_t zeroPos;
	
	CompressedString(){
		sz = 0;
		dna[0] = 'a', dna[1] = 'c', dna[2] = 'g', dna[3] = 't';
		encodedStr = NULL;
	}
	
	void save( FILE* fout ){
		fwrite((const void*)& sz,sizeof(uint32_t),1,fout);
		if( sz > 0 ){
			invalidPos.save( fout );
			for( uint32_t i = 0; i < ( sz + 3 ) / 4; i++ )
				fwrite((const void*)& encodedStr[i],sizeof(uint8_t),1,fout);
			fwrite((const void*)& zeroPos,sizeof(uint32_t),1,fout);
		}
	}
	
	void load( FILE* fin ){
		fread( &sz, sizeof( uint32_t ), 1, fin );
		if( sz > 0 ){
			invalidPos.load( fin );
			encodedStr = new uint8_t[( sz + 3 ) / 4];
			fread( encodedStr, sizeof( uint8_t ), ( sz + 3 ) / 4, fin );
			fread( &zeroPos, sizeof( uint32_t ), 1, fin );
		}
	}
	
	void set( const uint8_t* T, uint32_t n ){
		sz = n;
		invalidPos.reset( n );
		encodedStr = new uint8_t[( n + 3 ) / 4];
		for( uint32_t i = 0; i < ( n + 3 ) / 4; i++ )
			encodedStr[i] = 0;
		for( unsigned int i = 0; i < n; i++ )
			if( T[i] == 'n' )
				invalidPos.setBit( i );
			else if( T[i] == 0 )
				zeroPos = i;
			else{
				int idx = 0;
				while( dna[idx] != T[i] )
					idx++;
				encodedStr[i >> 2] |= ( idx << ( ( i & 3 ) << 1 ) );
			}
	}
	
	uint8_t charAt( uint32_t pos ){
		if( pos == zeroPos )
			return 0;
		if( invalidPos.getPos( pos ) )
			return 'n';
		return dna[( encodedStr[pos>>2] >> ( ( pos & 3 ) << 1 ) ) & 3];
	}
	
	uint32_t getBytes(){
		long long ret = 0;
		ret += sz / 4;
		ret += invalidPos.getBytes();
		return ret;
	}
};