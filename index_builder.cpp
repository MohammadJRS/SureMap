/*  FM-Index - Text Index
 *  Copyright (C) 2011  Matthias Petri
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string>
#include <vector>
#include <fstream>
#define puu pair< unsigned int, unsigned int >
#include "BWT.h"

static void
print_usage(const char *program)
{
    fprintf(stderr, "USAGE: %s <input> <outputFolder>\n", program);
	fprintf(stderr, "  input : reference file to be indexed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EXAMPLE: %s alice29.fa <file path>\n",program);
    fprintf(stderr, "\n");
    return;
}

char toLower( char ch ){
	if( ch >= 'A' && ch <= 'Z' )
		return 'a' + ch - 'A';
	return ch;
}

int main(int argc, char** argv) {
    int32_t opt, samplerate;
    string inname;
    uint8_t* T;
    uint32_t n;
    //cout << argc << endl;
    /* parse command line parameter */
    if (argc != 3) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    string pAdr = argv[2];
	if( pAdr[pAdr.length() - 1] != '/' )	
		pAdr += '/';
	opt = -1;
    inname = "";
	
	/* read filenames */
	inname = argv[1];
	int pt = inname.length() - 1;
	while( inname[pt] != '/' )
		pt--;
	//inname = inname.substr( pt + 1 );
	string refPrefix = pAdr + inname.substr( pt + 1 );
	cout << refPrefix << endl;
	string idxname;
	
	/* load input file */
	
	unsigned int offSet = 0;
	vector< string > refNames;
	vector< unsigned int > refOffset;
	ifstream fin( inname );
	string line, Ref;
	while( fin >> line ){
		if( line[0] == '>' ){
			string refName = line.substr( 1 );
			refNames.push_back( refName );
			refOffset.push_back( offSet );
		}
		else{
			bool valid = true;
			for( unsigned int i = 0; i < line.length(); i++ ){
				line[i] = toLower( line[i] );
				if( line[i] != 'a' && line[i] != 'c' && line[i] != 'g' && line[i] != 't' && line[i] != 'n' ){
					valid = false;
				}
			}
			if( !valid )
				continue;
			Ref += line;
			offSet += line.length();
		}
	}
	{
		uint8_t* T = new uint8_t[Ref.length()];
		for( uint32_t i = 0; i < Ref.length(); i++ )
			T[i] = Ref[i];
		FILE* fout = fopen( ( refPrefix + ".rinfo" ).c_str(), "wb" );
		ofstream fout2( ( refPrefix + ".cinfo" ) );
		CompressedString ref;
		ref.set( T, Ref.length() );
		ref.zeroPos = Ref.length();
		ref.save( fout );
		fclose( fout );
		delete[] T;
		for( int i = 0; i < refNames.size(); i++ )
			fout2 << refNames[i] << ' ' << refOffset[i] << endl;
	}

	n = Ref.length();
	cerr << "ref len = " << n << endl;
	T = new uint8_t[n + 100];
	// buiding forward bwt
	{
		idxname = refPrefix + ".fm";
		for( unsigned int i = 0; i < n; i++ )
			T[i] = Ref[i];
		T[n] = 0;
		cerr << "Building Forward BWT ...\n";
		BWT fm( T, n + 1, 0 );
		fm.save( idxname );
		//cout << "done\n";
    }
	cerr << "Forward bwt building is done!\n\n";
	// buiding reverse bwt
	{
		//delete[] T;
		T = new uint8_t[n + 100];
		idxname = refPrefix + ".rev.fm";
		unsigned int ptr = 0;
		for( long long i = n - 1; i >= 0; i--, ptr++ ){
			//cerr << ptr << endl;
			T[ptr] = Ref[i];
		}
		T[n] = 0;
		cerr << "Building Reverse BWT ...\n";
		BWT rev( T, n + 1, 1 );
		rev.save( idxname );
    }
	cerr << "Reverse bwt building is done!\n";
    return (EXIT_SUCCESS);
}

