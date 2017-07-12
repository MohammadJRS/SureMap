#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#define SA_SAMPLERATE	128
#define BT_SAMPLERATE	8

using namespace std;

struct BitArray{
	static const int blockLog2 = 6;
	static const int block = ( 1 << blockLog2 ) ;
	uint32_t smRate;
    uint64_t* arr;
	uint8_t* pc;
    uint32_t* sum;
    uint32_t sz;
	uint64_t needMask[64];
	BitArray(){
		smRate = 0;
		needMask[0] = 1;
		arr = NULL, sum = NULL, pc = NULL;
		for( int i = 1; i < 64; i++ )
			needMask[i] = needMask[i - 1] + ( 1ULL << i );
		sz = 0;
	}
	
	int myPopcount(unsigned long long x)
	{
		x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
		x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
		x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
		return (x * 0x0101010101010101ULL) >> 56;
	}
	
    void reset( uint64_t aS ){
        sz = ( aS + block - 1 ) / block;
        arr = new uint64_t[sz];
        for( uint64_t i = 0; i < sz; i++ )
			arr[i] = 0;
    }
	
    void setBit( uint32_t pos ){
        uint32_t idx = pos >> blockLog2;
        arr[idx] |= ( 1ULL << ( pos & ( block - 1 ) ) );
    }
	
    void setSum(){
		smRate = BT_SAMPLERATE;
		unsigned int smSize = ( sz + smRate - 1 ) / smRate;
        sum = new unsigned int[smSize];
		pc = new uint8_t[sz];
		for( uint32_t i = 0; i < sz; i++ )
			pc[i] = myPopcount( arr[i] );
        for( uint64_t i = 0; i < smSize; i++ )
			sum[i] = 0;
        sum[0] = myPopcount(arr[0]);
        for( long long i = smRate; i < sz; i += smRate ){
			uint32_t idx = i / smRate;
			sum[idx] = sum[idx - 1];
			for( unsigned int j = ( idx - 1 ) * smRate + 1; j <= i; j++ )
				sum[idx] += myPopcount(arr[j]);
		}
    }
	
    int getPos( uint32_t pos ){
        uint32_t idx = pos >> blockLog2;
        return ( arr[idx] & ( 1ULL << ( pos & ( block - 1 ) ) ) ) ? 1 : 0;
    }
	
    inline uint32_t getRank( uint32_t pos ){
        int idx = ( pos >> blockLog2 );
		if( sz == 0 )
			return 0;
		uint32_t ret;
		try{
			ret = myPopcount( arr[idx] & needMask[( pos & ( block - 1 ) )] );
		}catch( const std::exception& e ){
			 std::cout << e.what() << endl;
		}
		if( idx == 0 )
			return ret;
		int pre = idx - 1;
		while( pre >= 0 && ( pre & ( smRate - 1 ) ) ){
			//ret += myPopcount( arr[pre] );
			//cout << pre << ' ' << ret << endl;
			ret += pc[pre];
			pre--;
		}
		//cout << "fuck " << pre << ' ' << ret << endl;
		if( pre >= 0 && ( pre & ( smRate - 1 ) ) == 0 ){
			ret += sum[pre / smRate];
		}
        return ret;
    }
	
	void save( FILE* fout ){
		fwrite((const void*)& sz,sizeof(uint32_t),1,fout);
		if( sz > 0 ){
			for( uint32_t i = 0; i < sz; i++ )
				fwrite((const void*)& arr[i],sizeof(uint64_t),1,fout);
			fwrite((const void*)& smRate,sizeof(uint32_t),1,fout);
			if( smRate > 0 ){
				unsigned int smSize = ( sz + smRate - 1 ) / smRate;
				for( uint32_t i = 0; i < smSize; i++ )
					fwrite((const void*)& sum[i],sizeof(uint32_t),1,fout);
			}
		}
	}
	
	void load( FILE* fin ){
		pc = NULL;
		fread( &sz, sizeof( uint32_t ), 1, fin );
		if( sz > 0 ){
			arr = new uint64_t[sz];
			fread( arr, sizeof( uint64_t ), sz, fin );
			fread( &smRate, sizeof( uint32_t ), 1, fin );
			if( smRate > 0 ){
				unsigned int smSize = ( sz + smRate - 1 ) / smRate;
				sum = new uint32_t[smSize];
				fread( sum, sizeof( uint32_t ), smSize, fin );
				pc = new uint8_t[sz];
				for( uint32_t i = 0; i < sz; i++ )
					pc[i] = myPopcount( arr[i] );
			}
		}
	}
	
	long long getBytes(){
		long long ret = sz / 8;
		if( sum != NULL )
			ret += ( sz / smRate ) * 4;
		return ret;
	}
};