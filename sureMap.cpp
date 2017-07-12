
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <memory.h>
#include <chrono>
#include <set> 
#include <fstream>
#include <unordered_map>
#include <queue>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <boost/atomic.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <sstream>
#include "BWT2.h"
#define pci pair< char, int >
#define NV N[v]
#define MAXTHREADS 40
#define MAXSTATE 100000
#define MAXTOT MAXSTATE
#define MAXINROW MAXSTATE
#define MAXINPLEN 100000
#define MAXINDXES 1000000
#define pll pair< long long, long long >
#define puu std::pair< unsigned int, unsigned int >
#define shift 21
#define seed1 1000000007
#define seed2  1000019
#define seed3 100000007 
#define P pair< int, int >
#define PLI pair< long long, int >
#define MAXREADSIZE 100000
#define MAXREPORTSIZE 409600
#define MAXSHORTREADLEN 300 
#define MAXMASK 500
#define psu pair< string, uint32_t >
#define MAX_WRITER_THREAD 1
#define INF ( 1 << 22 )
#define pInterval pair< intervalNode,  intervalNode >
#define retVal pair< bool, bwtNode >
#define pss pair< string, string >

#pragma comment(linker, "/STACK:500000000")

 using namespace std;
 
struct compressedArray{
    unsigned long long* row = NULL; 
    unsigned long long HASHCode;
    unsigned short rowSize;
	compressedArray(){}
    compressedArray( vector<unsigned short> _row, unsigned short _base ){
        unsigned short base = _base;
        unsigned char cellWidth = 0;
        unsigned char cntInRow = 0;
        while( ( 1 << cellWidth ) - 1 < base + 1 )
            cellWidth++;
        HASHCode = 0;
        cntInRow = 64 / cellWidth;
        rowSize = ( _row.size() + cntInRow - 1 ) / ( cntInRow );
        row = new unsigned long long[rowSize];
        for( int i = 0; i < rowSize; i++ )
            row[i] = 0;
        for( int i = 0; i < _row.size(); i += cntInRow ){
            for( int j = min( (int)_row.size() - 1, (int)( i + cntInRow - 1 ) ); j >= i; j-- ){
                int idx = i / ( cntInRow );
                row[idx] <<= ( cellWidth );
                row[idx] += _row[j];
            }
        }
		for( int i = 0; i < rowSize; i++ )
			HASHCode = HASHCode * seed3 + row[i];
    }
    compressedArray( const compressedArray& m ){
    	HASHCode = m.HASHCode;
    	rowSize = m.rowSize;
    	row = new unsigned long long[rowSize];
    	for( int i = 0; i < rowSize; i++ )
    		row[i] = m.row[i];
    }
    unsigned short getCell( int pos, unsigned short base ){
        unsigned char cellWidth = 0;
        unsigned char cntInRow = 0;
        while( ( 1 << cellWidth ) - 1 < base + 1 )
            cellWidth++;
        cntInRow = 64 / cellWidth;
        int idx = pos / ( cntInRow );
        return ( row[idx] >> (cellWidth * ( pos % cntInRow ) ) ) & ( ( 1 << cellWidth ) - 1 ) ;
    }
    vector<unsigned short> getRow( unsigned short base ){
        vector<unsigned short> ret( 2 * base + 1 );
        for( int i = 0; i < ret.size(); i++ )
            ret[i] = getCell( i, base );
        return ret;
    }
    bool operator == ( const compressedArray& m ) const{
        return HASHCode == m.HASHCode;
    }

    compressedArray& operator=(compressedArray arg) // copy/move constructor is called to construct arg
	{
		HASHCode = arg.HASHCode;
		rowSize = arg.rowSize;
		row = new unsigned long long[rowSize];
		for( int i = 0; i < rowSize; i++ )
			row[i] = arg.row[i];
	    return *this;
	} 

    int getHachCode(){
        return HASHCode;

    }
    ~compressedArray(){
    	delete[] row;
    }
};

struct HASHNode{
    compressedArray to;
    unsigned long long mask;
	unsigned long long m2;
    bool isValid;
    int base;
    HASHNode( ){}
    HASHNode( compressedArray _to, unsigned long long _mask, int _base, bool _isValid ){
        to = _to, mask = _mask, base = _base, isValid = _isValid;
    }
};
 
 struct bwtNode{
    unsigned int first, second, len;
	int acceptedValue;
    uint32_t refRow;
    char flag;
    bwtNode( unsigned int _f, unsigned int _s, unsigned int _l, int _ac ){
        first = _f, second = _s, len = _l, acceptedValue = _ac;
    }
    bwtNode(){}
    bool operator < ( const bwtNode& m ) const{
        if( first != m.first )
                return first < m.first;
        if( second != m.second )
                return second < m.second;
        return len < m.len;
    }
    bool operator == ( const bwtNode& m ) const{
            return first == m.first && second == m.second && len == m.len;
    }
    void print(){
        cerr << first << ' ' << second << ' ' << len << ' ' << acceptedValue << ' ' << endl; 
    }
};

struct mappingInfo{
    string qName, read, QC;
    uint32_t refPos;
    psu cigar;
    int len, acceptedValue, flag;
    int startPos;
    int order;
    mappingInfo(){}
    mappingInfo( string _q, string _r, uint32_t _rf, int _len, int _ac, int _flag, string _QC, psu _c = psu( "", 0 ) ){
        qName = _q, read = _r, refPos = _rf, len = _len, acceptedValue = _ac, flag = _flag, QC = _QC, cigar = _c;
    }
	void print(){
		cerr << qName << endl << read << ' ' << QC << endl;
		cerr << refPos << ' ' << len << ' ' << acceptedValue << ' ' << flag << endl;
	}
	bool operator < ( const mappingInfo& m ) const{
		return order < m.order;
	}
};

struct intervalNode{
    unsigned int left, right;
    int editDistance, leftToRight;
	int mask;
    intervalNode(){}
    intervalNode( unsigned int _l, unsigned int _r,  int _ed, int _lr, int _m ){
        left = _l, right = _r, editDistance = _ed, leftToRight = _lr, mask = _m;
    }
};

struct nextState{
    int acceptedValue;
    vector<unsigned short> nextRowDp;
    bool valid;
    puu nxInterval;
    nextState( int edit ){
        nextRowDp.resize( 2 * edit + 1, edit + 1 );
        valid = false;
        acceptedValue = -1;
    }
    nextState(){}
    bool operator < ( const nextState& m ) const{
		if( acceptedValue == -1 )
			return false;
		if( m.acceptedValue == -1 )
			return true;
        return acceptedValue < m.acceptedValue;
    }
};

struct nextState3{
    int acceptedValue;
    compressedArray nextRowDp;
    bool valid;
    puu nxInterval;
    nextState3( int edit ){
        valid = false;
        acceptedValue = -1;
    }
    nextState3(){}
    bool operator < ( const nextState3& m ) const{
		if( acceptedValue == -1 )
			return false;
		if( m.acceptedValue == -1 )
			return true;
        return acceptedValue < m.acceptedValue;
    }
};

struct nextState2{
    int acceptedValue, error;
    bool valid;
    puu nxInterval;
    nextState2( int edit ){
        valid = false;
        acceptedValue = -1;
    }
    nextState2(){}
    bool operator < ( const nextState2& m ) const{
        return error < m.error;
    }
};

struct readMappingInfo{
	string read, qual, readName;
	int maxReport, bestOne, maxDiffMismatch, maxDiffEdit, gap, uniqeOption;
	double noisePercent = -1.0;
	readMappingInfo( string _read, string _qual, string _readName, int _maxReport, int _bestOne, int _maxDiffMismatch, int _maxDiffEdit,
	int _gap, int _uniqeOption, double _noicePercent ){
		read = _read, qual = _qual, readName = _readName, maxReport = _maxReport, bestOne = _bestOne, maxDiffMismatch = _maxDiffMismatch,
		maxDiffEdit = _maxDiffEdit, gap = _gap, uniqeOption = _uniqeOption, noisePercent = _noicePercent;
	}
};
 
BWT fmIdx, revIdx;
CompressedString Ref;
mutex mu;

/* runing options */
int mc = 10;
int core = 1;
int globalMaxReport = 1;
int maxReport[MAXTHREADS];
int globalBestOne = 0;
string globalMode = "normal";
int bestOne[MAXTHREADS];
double globalNoisePercent = -0.1;
double noisePercent[MAXTHREADS];
int globalMaxDiffMismatch = 1;
int maxDiffMismatch[MAXTHREADS];
int globalMaxDiffEdit = 1;
int maxDiffEdit[MAXTHREADS];
int globalGap = 0;
int gap[MAXTHREADS];
int globalUniqeOption = 0;
int globalLongReadNoice = 30;
int uniqeOption[MAXTHREADS];
string outputAdr = "report.sam";
bool longRead = false;
int mxFailed[MAXTHREADS];
int minVal[MAXTHREADS];
int hardToMap = 0;
int Mode = 20;
int fragLen = 500;

/********** GLOBAL VARIABLES ********************/
bwtNode indexes[MAXTHREADS][MAXINDXES];
vector< bwtNode > threadResults[MAXTHREADS];

boost::atomic<bool> readVector[MAXTHREADS], writeVector[MAXTHREADS];
boost::atomic<int> readAccess[MAXTHREADS];
boost::atomic<bool> finished;
boost::atomic<bool> reading[seed2], writing[seed2];

HASHNode HASH[seed2];
bool isSetHASH[seed2];

string reads[MAXREADSIZE + 1000];
string qc[MAXREADSIZE + 1000];
string readsName[MAXREADSIZE + 1000];
string readsPerThread[MAXTHREADS];
unsigned char nonComplete[MAXTHREADS][MAXINPLEN * 4];
int leftIndex[MAXTHREADS], rightIndex[MAXTHREADS];
char dna[] = {'a', 'c', 'g', 't'};
char tst[100];
map< P, int > failedTry[MAXTHREADS][MAXMASK];
vector< mappingInfo > toWrite;

vector< string > refNames;
vector< string > findAllStrings[MAXTHREADS];

vector< uint32_t > refOffSets;
vector< long long > ms;

vector<mappingInfo> results[MAXTHREADS];


long long cntReads;
long long sptr = 0;
long long vis = 0;
long long mch[MAXTHREADS] = {0};
long long lp[MAXTHREADS];
long long totReads = 0;
unsigned long long needMask[64];
long long zeroGen[64][64];
long long windowMask[MAXTHREADS][6][MAXMASK][MAXMASK];
int Flag[MAXTHREADS];
int dbound[MAXTHREADS][MAXMASK][MAXMASK];
int gg = 1;
int sampleRate = 11;
int MAXD = 0;
int debg = 0;
int cntDFS = 0;
int cntBut = 0;
int hit = 0;
int unhit = 0;
long long rc = 0;
int drX[3] = { 0, -1, -1 };
int drY[3] = { -1, -1, 0 };
int totChar;
int indexesPtr[MAXTHREADS];
int dbg = 0;
int loopTh = 100;
int resPtr[MAXTHREADS];
int direction[MAXTHREADS];
int readPtr = 0;
string running_mode;

unsigned int refLen;
unsigned int preSearched[MAXTHREADS][MAXMASK * 4 + 1];

bool forceStop[MAXTHREADS];
bool onlyMis = true;
bool threadJobDone[100];
/*****************************************************************/

bool dfsButtumUpWithGap( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, vector<unsigned short>& curRowDp, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked );
	
bool CompressedDfsButtumUpWithGap( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, compressedArray& dpRow, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked );

bool dfsButtumUp( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, int error, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked  );


template<typename T>
void print( vector<T> v){
	for( int i = 0; i < v.size(); i++ )
		cerr << v[i] << ' ';
	cerr << endl;
}


inline bool isValidDnaCharacter( char ch ){
    return ( ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' );
}

void toLower( char& ch ){
    if( ch >= 'A' && ch <= 'Z' )
            ch = ( 'a' + ch - 'A' );
}

inline unsigned int locateBwt( unsigned int idx ){
    return fmIdx.locateRow( idx );
}

inline unsigned int locateInvBwt( unsigned int idx ){
    return revIdx.locateRow( idx );
}

bool newCmp( const bwtNode a, const bwtNode b ){
    return a.acceptedValue < b.acceptedValue;
}

psu getPosition( uint32_t offSet ){
	int idx = 0;
	for( idx = 0; idx < refNames.size() - 1; idx++ )
		if( offSet >= refOffSets[idx] && offSet < refOffSets[idx + 1] )
			break;
	return psu( refNames[idx], offSet - refOffSets[idx] + 1 );
}

void resetControllers(){
    for( int i = 0; i < core; i++ ){
        readVector[i] = writeVector[i] = false;
        readAccess[i] = -1;
    }
	for( int i = 0; i < seed2; i++ )
		reading[i] = writing[i] = false;
    finished = false;
}

string reverseComplement( string s ){
    reverse( s.begin(), s.end() );
    for( int i = 0; i < s.length(); i++ ){
        if( s[i] == 'n' )
			continue;
        if( s[i] == 'a' )	s[i] = 't';
		else 
			if( s[i] == 't' )	s[i] = 'a';
		else
			if( s[i] == 'c' )	s[i] = 'g';
        else
			if( s[i] == 'g' )	s[i] = 'c';
    }
    return s;
}

int getMinEdit( int leftIdx, int rightIdx, int tIdx ){
    unsigned int lb = 0, rb = fmIdx.n - 1;
    int z = 0, j = 0;
    for( int i = leftIdx; i <= rightIdx; i++ ){
        int id = fmIdx . remap[readsPerThread[tIdx][i]];
		revIdx . updateInterval( lb, rb, id );
        if( lb > rb ){
            z++;
            lb = 0, rb = fmIdx.n - 1; 
        }
    }
    return z;
} 

puu getInvToBwt( int tIdx ){
    unsigned int tl = 0, tr = fmIdx.n - 1, iter = 0;
	for( int i = rightIndex[tIdx] - 1; i >= leftIndex[tIdx]; i-- ){
		int id = fmIdx.remap[nonComplete[tIdx][i]];
		fmIdx.updateInterval( tl, tr, id );
		//cerr << iter++ << ' ' << tl << ' ' << tr << ' ' << (char)nonComplete[tIdx][i] << ' ' << id << endl;
	}
    return puu( tl, tr );
}

puu getBwtToInv( int tIdx ){
	unsigned int tl = 0, tr = fmIdx.n - 1;
	for( int i = leftIndex[tIdx]; i < rightIndex[tIdx]; i++ ){
		int id = fmIdx.remap[nonComplete[tIdx][i]];
		revIdx.updateInterval( tl, tr, id );
	}
    return puu( tl, tr );
}

puu getInterVal( string& q, BWT& fm ){
    unsigned int tl = 0, tr = fm.n - 1;
    for( int i = 0; i < q.length(); i++ ){
        int id = fmIdx . remap[q[i]];
		fm . updateInterval( tl, tr, id );
    }
    return puu( tl, tr );
}


puu getNextInterval( bwtNode& cur, BWT& fm, int id ){
    unsigned int tl = cur.first;
    unsigned int tr = cur.second;
	fm . updateInterval( tl, tr, id );
    return puu( tl, tr );
}

void buildLevstein( string q, bool REV, int tId, BWT& fm ){
	if( REV )
		reverse( q.begin(), q.end() );
	unsigned int lb = 0, rb = fm.n - 1;
	for( int i = 0; i < q.length(); i++ ){
		fm.updateInterval( lb, rb, fm.remap[q[i]] );
		if( lb > rb )
			break;
	}
	if( lb <= rb )
		indexes[tId][indexesPtr[tId]++] = (bwtNode( lb, rb, q.length(), 0 ) );
}

void newBuildDArray( int leftIdx, int rightIdx, int st, int dir, int tIdx, int DArr[], BWT& fm,  int D ){
    unsigned int lb = 0, rb = fm.n - 1;
    int len = rightIdx - leftIdx + 1;
    for( int i = 0; i < len; i++ )
        DArr[i] = D + 1;
    int z = 0, j = 0, idx = 0;
    for( int I = st; ( I >= leftIdx && I <= rightIdx ); I += dir, idx++ ){
        int i = fmIdx . remap[readsPerThread[tIdx][I]];
		fm.updateInterval( lb, rb, i );
        if( lb > rb ){
            z++;
            j = i + 1;
            lb = 0, rb = fm.n - 1; 
        }
        DArr[idx] = z;
    }
}



vector< int > getLastRow( string s1, string s2, int maxDiff, int& bestEnd, int mode = 0 ){
	vector< int > lastRows[2];
	lastRows[0] = vector< int >( s2.length() + 1, 10 * maxDiff + 1 );
	lastRows[1] = vector< int >( s2.length() + 1, 10 * maxDiff + 1 );
	int turn = 1;
	int curBest = 10000000;
	for( int i = 0; i < lastRows[turn].size(); i++ )
		lastRows[1 - turn][i] = i;
	for( int i = 1; i <= s1.length(); i++ ){
		int preRow = 1 - turn;
		for( int j = 0; j <= s2.length(); j++ ){
			if( j < 0 || j > s2.length() )
				continue;
			if( j == 0 )
				lastRows[turn][j] = ( mode == 0 ) ? 0 : i;
			else{
				lastRows[turn][j] = 10 * maxDiff + 1;
				lastRows[turn][j] = min( lastRows[preRow][j] + 1, lastRows[turn][j - 1] + 1 );
				lastRows[turn][j] = min( lastRows[turn][j], lastRows[preRow][j - 1] + 1 );
				if( s1[i - 1] == s2[j - 1] )
					lastRows[turn][j] = min( lastRows[turn][j], lastRows[preRow][j - 1] );
			}
		}
		if( lastRows[turn][s2.length()] < curBest ){
			bestEnd = i;
			curBest = lastRows[turn][s2.length()];
		}
		turn = 1 - turn;
	}
	return lastRows[1 - turn];
}

int getBestStartingPos( string& s1, string& s2, int maxDiff, bool canCutFirst = true, bool canCutLast = true ){
	int bestEnd;
	if( canCutLast == true ){
		getLastRow( s1, s2, maxDiff, bestEnd, 0 );
		s1 = s1.substr( 0, bestEnd );
	}
	if( canCutFirst == false ){
		return 0;
	}
	reverse( s1.begin(), s1.end() );
	reverse( s2.begin(), s2.end() );
	getLastRow( s1, s2, maxDiff, bestEnd, 0 );

	int ret = s1.length() - bestEnd;
	s1 = s1.substr( 0, bestEnd );
	reverse( s1.begin(), s1.end() );
	reverse( s2.begin(), s2.end() );
	return ret;
}

string Rev( string s1 ){
	reverse( s1.begin(), s1.end() );
	return s1;
}

vector< int > Rev( vector< int > v ){
	reverse( v.begin(), v.end() );
	return v;
}

// return pair( Z, W )
pss Hirschberg( string s1, string s2, int maxDiff, int depth = 0 ){
	pss ret;
	string Z, W;
	int tmp;
	if( s1.length() == 0 ){
		for( int i = 0; i < s2.length(); i++ ){
			Z += s2[i];
			W += '-';
		}
		return pss( Z, W );
	}
	else if( s2.length() == 0 ){
		for( int i = 0; i < s1.length(); i++ ){
			Z += '-';
			W += s1[i];
		}
		return pss( Z, W );
	}
	else if( s1.length() == 1 ){
		for( int i = 0; i < s2.length(); i++ ){
			Z += s2[i];
			W += '-';
		}
		if( s2.find( s1[0] ) != string::npos )
			W[s2.find( s1[0] )] = s1[0];
		else W[0] = s1[0];
		return pss( Z, W );
	}
	else if( s2.length() == 1 ){
		for( int i = 0; i < s1.length(); i++ ){
			Z += '-';
			W += s1[i];
		}
		if( s1.find( s2[0] ) != string::npos )
			Z[s1.find( s2[0] )] = s2[0];
		else Z[0] = s2[0];
		return pss( Z, W );
	}
	int xlen = s1.length();
	int xmid = xlen / 2;
	int ylen = s2.length();
	vector< int > ScoreL = getLastRow( s1.substr( 0, xmid ), s2, maxDiff, tmp, 1 );
	vector< int > ScoreR = Rev( getLastRow( Rev( s1.substr( xmid ) ), Rev(s2), maxDiff, tmp, 1 ) );
	int ymid = 0, best = 10000000;
	for( int i = 0; i <= s2.length(); i++ )
		if( ScoreL[i] + ScoreR[i] < best ){
			ymid = i;
			best = ScoreL[i] + ScoreR[i];
		}
	pss h1 = Hirschberg(s1.substr( 0, xmid ), s2.substr( 0, ymid ), maxDiff, depth + 1 );
	pss h2 = Hirschberg(s1.substr( xmid ), s2.substr( ymid ), maxDiff, depth + 1 );
	return pss( h1.first + h2.first, h1.second + h2.second );
}

void testHeirch( string s1, string s2, int maxDiff ){
	cerr << s1 << ' ' << s2 << endl;
	getBestStartingPos( s1, s2, maxDiff );
	cerr << s1 << ' ' << s2 << endl;
	cerr << Hirschberg( s1, s2, maxDiff ).first << endl << Hirschberg( s1, s2, maxDiff ).second << endl;
}

bool processResult( bwtNode& curMap, int tIdx, int preLeftToRight ){
	if( preLeftToRight == 1 ){
        puu nx = getInvToBwt( tIdx );
		curMap.first = nx.first, curMap.second = nx.second;
   	}
 	for( long long BwtIdx = curMap.first; BwtIdx <= curMap.second; BwtIdx++ ){
        bwtNode newCurMap = curMap;
		newCurMap.first = newCurMap.second = BwtIdx;
        newCurMap.refRow = fmIdx.locateRow( newCurMap.first );
		bool valid = true;
        for( int i = 0; i < resPtr[tIdx]; i++ ){
            if( abs( (long long)threadResults[tIdx][i].refRow - (long long)newCurMap.refRow ) <= readsPerThread[tIdx].length() ){ //max( threadResults[tIdx][i].acceptedValue, newCurMap.acceptedValue ) ){
                if( newCurMap.acceptedValue >= threadResults[tIdx][i].acceptedValue ){
                    valid = false;
                    break;
                }
                else{
                    swap( threadResults[tIdx][resPtr[tIdx] - 1], threadResults[tIdx][i] );
					threadResults[tIdx].pop_back();
                    resPtr[tIdx]--;
                    i--;
                }
            }
        }
        if( valid == false )
            continue;
        newCurMap.flag = ( direction[tIdx] == 0 ) ? 0 : 16;
        threadResults[tIdx].push_back( newCurMap );
        resPtr[tIdx]++;
        if( uniqeOption[tIdx] && resPtr[tIdx] > 1 ){
            forceStop[tIdx] = true;
            resPtr[tIdx] = 0;
            return true;
        }
        if( resPtr[tIdx] >= maxReport[tIdx] )
            return true;
    }
    return  false;
}



bool buttumUpMapping( bwtNode curMap, vector<pInterval>& path, int idxPath, int preLeftToRight, int tIdx ){
    if( lp[tIdx] > loopTh )
            return false;
	if( resPtr[tIdx] >= maxReport[tIdx] )
		return false;
    cntBut++;
    if( idxPath == -1 )
        return processResult( curMap, tIdx, preLeftToRight );
	//cerr << idxPath << endl;
    int curLeftToRight = ( path[idxPath].first.left < path[idxPath].second.left ) ? 1 : 0;
    if( curLeftToRight == 1 && preLeftToRight == 0 ){
        puu nx = getBwtToInv( tIdx );
        curMap.first = nx.first, curMap.second = nx.second;
    }
    if( curLeftToRight == 0 && preLeftToRight == 1 ){
        puu nx = getInvToBwt( tIdx );
        curMap.first = nx.first, curMap.second = nx.second;
    }
    int curMask = path[idxPath].first.mask;
    int edit = path[idxPath].second.editDistance;
    int len = path[idxPath].second.right - path[idxPath].second.left + 1;
    int* DArray = new int[len];
    if( curLeftToRight == 1 ){
        newBuildDArray( path[idxPath].second.left, path[idxPath].second.right, path[idxPath].second.left, 1, tIdx, DArray, revIdx, path[idxPath].second.editDistance );
    }
    else{
        newBuildDArray( path[idxPath].second.left, path[idxPath].second.right, path[idxPath].second.right, -1, tIdx, DArray, fmIdx, path[idxPath].second.editDistance );
    }
    
    if( DArray[path[idxPath].second.right - path[idxPath].second.left] > edit )
		return false;
	bool nx = false;
	int prePtVal = resPtr[tIdx];
    if( gap[tIdx] ){
            vector< unsigned short > root( 2 * edit + 1, edit + 1 );
            for( int i = edit; i < root.size(); i++ )
				root[i] = min( edit + 1, abs( i - edit ) + curMap.acceptedValue );
			int lower = (curLeftToRight == 0 ) ? ( curMap.acceptedValue + path[idxPath].second.editDistance / 2 ) : -1;
            curMap.acceptedValue = -1;
			
			if( edit <= 30 ){//&& readsPerThread[tIdx].length() < MAXMASK / 3 - 50 ){
				compressedArray tmp( root, edit );
				nx = CompressedDfsButtumUpWithGap(curMap, path[idxPath].second.left, path[idxPath].second.right, 0, tIdx, tmp, edit, curLeftToRight, ( curLeftToRight == 1 ) ? revIdx : fmIdx, 
				path, idxPath, DArray, edit + 10, lower, false );
			}
			else{
			
				nx = dfsButtumUpWithGap(curMap, path[idxPath].second.left, path[idxPath].second.right, 0, tIdx, root, edit, curLeftToRight, ( curLeftToRight == 1 ) ? revIdx : fmIdx, 
				path, idxPath, DArray, edit + 10, lower, false );
			
			}
    }
    else{
        int edit = path[idxPath].second.editDistance;
        nx = dfsButtumUp(curMap, path[idxPath].second.left, path[idxPath].second.right, 0, tIdx, curMap.acceptedValue, edit, curLeftToRight, ( curLeftToRight == 1 ) ? revIdx : fmIdx, 
        path, idxPath, DArray, edit + 10, (curLeftToRight == 0 ) ? ( curMap.acceptedValue + path[idxPath].second.editDistance / 2 ) : -1, false );
    }
    delete[] DArray;
    if( resPtr[tIdx] == prePtVal )
        preSearched[tIdx][curMask]++;
    if( nx )
        return true;
    return false;
}


bool dfsButtumUp( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, int curError, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked ){
	//cerr << leftIdx << ' ' << rightIdx << ' ' << curError << ' ' << curRow << ' ' << maxEditDistance << endl;
	cntDFS++;
	int mask = path[idxPath].first.mask;
	/*int tmask = mask;
    vector< int > ups;
	while( mask ){
		ups.push_back( mask );
		mask /= 2;
	}
	reverse( ups.begin(), ups.end() );
	bool inv = false;
	for( int i = 1; i < ups.size(); i++ ){
		if(  preSearched[tIdx][ups[i]] > ( 1 << ( i + Mode ) ) + i * mc )
			inv = true;
	}
	if( inv )
		return false;
	*/
	int cm = curError;//dpRow.getCell( maxEditDistance, maxEditDistance );
	if( failedTry[tIdx][mask][P( curRow,cm )] >= mxFailed[tIdx] && ( dbound[tIdx][mask][curRow] < cm ) ){
		return false;
	}
	if( resPtr[tIdx] >= maxReport[tIdx] )
		return false;
	
	//cerr << curMap.first << ' ' << curMap.second << ' ' << fm.bwt.charAt( curMap.first ) << endl;
    int len = ( rightIdx - leftIdx + 1 );
    int nextRowNum = curRow + 1;
    nextState2 nxState[10];
    for( int i = 0; i < totChar; i++ )
        nxState[i] = nextState2( maxEditDistance );
	if( curMap.first > curMap.second || fm.bwt.charAt( curMap.first ) == 0 )
		return false;
	char on = ( curMap.first == curMap.second ) ? fm.remap[fm.bwt.charAt( curMap.first )] : -1;
	int lch = 1, rch = totChar - 1;
	if( on != -1 )
		lch = rch = on;
	int preVal = resPtr[tIdx];
	bool useLess = true;
    for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( fm.remap_reverse[nextChar] == 'n' )
                continue;
        puu nxInterval = getNextInterval(curMap, fm, nextChar);
        if( nxInterval.first > nxInterval.second ){
            continue;
        } 
		int col = nextRowNum;
        nxState[nextChar].nxInterval = nxInterval;
        int minAcVal = maxEditDistance + 1, chk = maxEditDistance + 1;
        int ref = curError + 1;
        if( leftToRight == 1 )
            ref = min( ((( nextChar == fm.remap[readsPerThread[tIdx][leftIdx + col - 1]] ) ? 0 : 1 ) + curError), ref );
        else 
            ref = min( ((( nextChar == fm.remap[readsPerThread[tIdx][rightIdx - col + 1]] ) ? 0 : 1 ) + curError), ref );
        if( ref <= maxEditDistance ){
            nxState[nextChar].valid = true;
            nxState[nextChar].error = ref;
            if( col == len )
                    nxState[nextChar].acceptedValue = ref;
        }
    }
	if( on == -1 ){
        lch = 1, rch = totChar - 1;
	}
	for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( nxState[nextChar].valid == false )
            continue;
		if( nxState[nextChar].acceptedValue != -1
                && nxState[nextChar].acceptedValue > lower ){
			if( leftToRight )
				nonComplete[tIdx][rightIndex[tIdx]++] = fm.remap_reverse[nextChar];
			else nonComplete[tIdx][--leftIndex[tIdx]] = fm.remap_reverse[nextChar];
			
            bwtNode nextMap = curMap;
            nextMap.first = nxState[nextChar].nxInterval.first;
            nextMap.second = nxState[nextChar].nxInterval.second;
            nextMap.len++;
            nextMap.acceptedValue = nxState[nextChar].acceptedValue;
            bool nx = buttumUpMapping( nextMap, path, idxPath - 1, leftToRight, tIdx );
			if( leftToRight )
				rightIndex[tIdx]--;
			else leftIndex[tIdx]++;
            if( nx )
                return true;
        }
        {
			if( nxState[nextChar].acceptedValue != -1 )
				continue;
			if( leftToRight )
				nonComplete[tIdx][rightIndex[tIdx]++] = fm.remap_reverse[nextChar];
			else nonComplete[tIdx][--leftIndex[tIdx]] = fm.remap_reverse[nextChar];
            bwtNode nextMap = curMap;
            nextMap.first = nxState[nextChar].nxInterval.first;
            nextMap.second = nxState[nextChar].nxInterval.second;
			nextMap.len++;
			bool nx = dfsButtumUp( nextMap, leftIdx, rightIdx, nextRowNum, tIdx, nxState[nextChar].error, maxEditDistance, leftToRight, fm, path, idxPath,
			DArray, ( nxState[nextChar].acceptedValue == -1 ) ? lastAc : nxState[nextChar].acceptedValue, lower, isChecked );
			if( leftToRight )
				rightIndex[tIdx]--;
			else leftIndex[tIdx]++;
			if( nx )
				return true;
        }
        
    }
	if( resPtr[tIdx] > preVal )
		useLess = false;
	/*if( useLess ){
		dbound[tIdx][mask][curRow] = min( dbound[tIdx][mask][curRow], cm );
		failedTry[tIdx][mask][P( curRow,cm )]++;
	}*/
    return false;
}


bool dfsButtumUpWithGap( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, vector<unsigned short>& curRowDp, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked ){
	rc++;
	cntDFS++;
	int mask = path[idxPath].first.mask;
	int cm = curRowDp[maxEditDistance];
	/*if( failedTry[tIdx][mask][P( curRow,cm )] >= mxFailed[tIdx] && ( dbound[tIdx][mask][curRow] < cm ) ){
		return false;
	}*/
	vector< int > ups;
    if( resPtr[tIdx] >= maxReport[tIdx] )
        return true;
    int nextRowNum = curRow + 1;
    int len = rightIdx - leftIdx + 1;
    nextState nxState[6];
    for( int i = 0; i < totChar; i++ )
        nxState[i] = nextState( maxEditDistance );
    char on = ( curMap.first == curMap.second ) ? fm.remap[fm.bwt.charAt( curMap.first )] : -1;
	if( on == 0 )
            return false;
	int lch = 1, rch = totChar - 1;
	if( on != -1 )
		lch = rch = on;
    for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( fm.remap_reverse[nextChar] == 'n' )
                continue;
        puu nxInterval = getNextInterval(curMap, fm, nextChar);
        if( nxInterval.first > nxInterval.second ){
            continue;
        } 
        nxState[nextChar].nxInterval = nxInterval;
            int minAcVal = maxEditDistance + 1, chk = maxEditDistance + 1 ;
            for( int k = 0; k < 2 * maxEditDistance + 1; k++ ){
                int col = nextRowNum + ( k - maxEditDistance );
                unsigned short& ref = nxState[nextChar].nextRowDp[k];
                ref = maxEditDistance + 1;
                if( col >= 0 && col <= len ){
                    if( k > 0 )
                        ref = min( (unsigned short)( 1 + nxState[nextChar].nextRowDp[k - 1] ), ref );
                    if( k < curRowDp.size() && col > 0 ){
                        if( leftToRight == 1 )
                            ref = min( (unsigned short)((( nextChar == fm . remap[readsPerThread[tIdx][leftIdx + col - 1]] ) ? 0 : 1 ) + curRowDp[k]), ref );
                        else 
                            ref = min( (unsigned short)((( nextChar == fm . remap[readsPerThread[tIdx][rightIdx - col + 1]] ) ? 0 : 1 ) + curRowDp[k]), ref );
                    }
                    if( k + 1 < curRowDp.size() )
                        ref = min( (unsigned short)( 1 + curRowDp[k + 1] ), ref );
                    if( col < len )
                        minAcVal = min(  minAcVal, ref + DArray[len - col - 1] );
                    else
                        minAcVal = min( minAcVal, (int)ref );
                    if( col == len && ref <= maxEditDistance )
                        nxState[nextChar].acceptedValue = ref;
                }
            }
        nxState[nextChar].valid = minAcVal <= maxEditDistance;
    }
    if( on == -1 ){
        lch = 0, rch = totChar - 1;
    }
	int preVal = resPtr[tIdx];
	bool useLess = true;
    bool cn = false;
    bool shTry = false;
    for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( nxState[nextChar].valid == false )
            continue;
        {
            cn = true;
            bwtNode nextMap = curMap;
            nextMap.first = nxState[nextChar].nxInterval.first;
            nextMap.second = nxState[nextChar].nxInterval.second;
            nextMap.acceptedValue = nxState[nextChar].acceptedValue;
            nextMap.len++;
            bool cc = true;
            if( nxState[nextChar].acceptedValue != -1 && nxState[nextChar].acceptedValue >= lastAc ){
                    shTry = true;
            }
            if( nxState[nextChar].acceptedValue != -1 && nxState[nextChar].acceptedValue > lastAc ){
                    cc = false;
            }
            if( leftToRight )
				nonComplete[tIdx][rightIndex[tIdx]++] = fm.remap_reverse[nextChar];
			else nonComplete[tIdx][--leftIndex[tIdx]] = fm.remap_reverse[nextChar];
			
			bool nx = dfsButtumUpWithGap( nextMap, leftIdx, rightIdx, nextRowNum, tIdx, nxState[nextChar].nextRowDp, maxEditDistance, leftToRight, fm, path, idxPath,
			DArray, ( nxState[nextChar].acceptedValue == -1 ) ? lastAc : nxState[nextChar].acceptedValue, lower, cc );
			if( nx )
				return true;
			if( leftToRight )
				rightIndex[tIdx]--;
			else leftIndex[tIdx]++;
        }
    }
    if( isChecked && ( shTry || cn == false ) &&  curMap.acceptedValue != -1 && curMap.acceptedValue <= maxEditDistance && curMap.acceptedValue > lower ){
		bwtNode nextMap = curMap;
		bool nx = buttumUpMapping( nextMap, path, idxPath - 1, leftToRight, tIdx );
		if( nx )
			return true;
    }
	if( resPtr[tIdx] > preVal )
		useLess = false;
	if( useLess ){
		dbound[tIdx][mask][curRow] = min( dbound[tIdx][mask][curRow], cm );
		failedTry[tIdx][mask][P( curRow,cm )]++;
	}
    return false;
}

bool CompressedDfsButtumUpWithGap( bwtNode curMap, int leftIdx, int rightIdx, int curRow, int tIdx, compressedArray& dpRow, 
    int maxEditDistance, int leftToRight, BWT& fm, vector<pInterval>& path, int idxPath, int DArray[], int lastAc, int lower, bool isChecked ){
	rc++;
	cntDFS++;
	minVal[tIdx] = min( minVal[tIdx], idxPath );
	int mask = path[idxPath].first.mask;
	/*int cm = dpRow.getCell( maxEditDistance, maxEditDistance );
	if( failedTry[tIdx][mask][P( curRow,cm )] >= mxFailed[tIdx] && ( dbound[tIdx][mask][curRow] < cm ) ){
		return false;
	}*/
	int tmask = mask;
    vector< int > ups;
	while( tmask ){
		ups.push_back( tmask );
		tmask /= 2;
	}
	reverse( ups.begin(), ups.end() );
	bool inv = false;
	for( int i = 1; i < ups.size(); i++ ){
		if(  preSearched[tIdx][ups[i]] > ( i / 2 + Mode ) * mc )
			inv = true;
	}
	if( inv )
		return false;
	
	if( resPtr[tIdx] >= maxReport[tIdx] )
        return true;
    int nextRowNum = curRow + 1;
    int len = rightIdx - leftIdx + 1;
    nextState3 nxState[6];
    for( int i = 0; i < totChar; i++ ){
        nxState[i].valid = false;
        nxState[i].acceptedValue = -1;
    }
    char on = ( curMap.first == curMap.second ) ? fm.remap[fm.bwt.charAt( curMap.first )] : -1;
	if( on == 0 )
		return false;
	int lch = 1, rch = totChar - 1;
	if( on != -1 )
		lch = rch = on;
	int mid = ( leftToRight ) ? ( leftIdx + nextRowNum - 1 ) : ( rightIdx - nextRowNum + 1 );
    for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( fm.remap_reverse[nextChar] == 'n' || nextChar == 0 ){
        	nxState[nextChar].valid = false;
            continue;
        }
        puu nxInterval = getNextInterval(curMap, fm, nextChar);
        if( nxInterval.first > nxInterval.second ){
            continue;
        } 
        nxState[nextChar].nxInterval = nxInterval;
        long long realMask = windowMask[tIdx][nextChar][mask][mid + maxEditDistance];
        unsigned int index = (( ( dpRow.HASHCode * seed1 ) + realMask )) % seed2;
		reading[index] = true;
		while( writing[index] );
		
        if( !isSetHASH[index] || HASH[index].mask != dpRow.HASHCode || maxEditDistance != HASH[index].base || realMask != HASH[index].m2){
            unhit++;
            vector<unsigned short> curRowDp = dpRow.getRow(maxEditDistance);
            vector<unsigned short> nextRowDp( 2 * maxEditDistance + 1, maxEditDistance + 1 );
            int minAcVal = maxEditDistance + 1, chk = maxEditDistance + 1 ;
            for( int k = 0; k < 2 * maxEditDistance + 1; k++ ){
                int col = nextRowNum + ( k - maxEditDistance );
                unsigned short& ref = nextRowDp[k];
                ref = 1 + maxEditDistance;
                ref = min( (unsigned short)(maxEditDistance + 1), (unsigned short)(1 + curRowDp[k]) );
                if( k > 0 )
                    ref = min( (unsigned short)( 1 + nextRowDp[k - 1] ), ref );
                int rIdx = ( leftToRight ) ? ( leftIdx + col - 1 ) : ( rightIdx - col + 1 );
                if( rIdx >= leftIdx && rIdx <= rightIdx ){
                    ref = min( (unsigned short)((( nextChar == fm . remap[readsPerThread[tIdx][rIdx]] ) ? 0 : 1 ) + curRowDp[k]), ref );
                }
                if( k + 1 < curRowDp.size() )
                    ref = min( (unsigned short)( 1 + curRowDp[k + 1] ), ref );
                minAcVal = min( minAcVal, (int)ref );
                if( col == len && ref <= maxEditDistance )
                    nxState[nextChar].acceptedValue = ref;
            }
			
			while( writing[index] );
			writing[index] = true;
			
            isSetHASH[index] = true;
            HASH[index].mask = dpRow.HASHCode;
            HASH[index].to = compressedArray( nextRowDp, maxEditDistance );
            HASH[index].base = maxEditDistance;
            HASH[index].isValid = minAcVal <= maxEditDistance;
            HASH[index].m2 = realMask;
			writing[index] = false;
        }
        if( HASH[index].isValid == false )
            continue;
        nxState[nextChar].valid = true;
		nxState[nextChar].nextRowDp = compressedArray(HASH[index].to);
        if( abs( len - nextRowNum ) <= maxEditDistance ){
            int val = HASH[index].to.getCell( len - nextRowNum + maxEditDistance, maxEditDistance  );
            if( val <= maxEditDistance )
                nxState[nextChar].acceptedValue = val;
        }
		reading[index] = false;
    }
    if( on == -1 ){
        lch = 1, rch = totChar - 1;
    }
    bool cn = false;
    bool shTry = false;
	bool lq = false;
	int preVal = resPtr[tIdx];
	bool useLess = true;
    for( int nextChar = lch; nextChar <= rch; nextChar++ ){
        if( nxState[nextChar].valid == false )
            continue;
		
        {
			cn = true;
            bwtNode nextMap = curMap;
            nextMap.first = nxState[nextChar].nxInterval.first;
            nextMap.second = nxState[nextChar].nxInterval.second;
            nextMap.acceptedValue = nxState[nextChar].acceptedValue;
            nextMap.len++;
            bool cc = true;
            if( nxState[nextChar].acceptedValue != -1 && nxState[nextChar].acceptedValue >= lastAc ){
                    shTry = true;
            }
            if( nxState[nextChar].acceptedValue != -1 && nxState[nextChar].acceptedValue > lastAc ){
                    cc = false;
            }
			if( nxState[nextChar].acceptedValue != -1 && nxState[nextChar].acceptedValue < lastAc ){
                    lq = true;
            }
			if( leftToRight )
				nonComplete[tIdx][rightIndex[tIdx]++] = fm.remap_reverse[nextChar];
			else nonComplete[tIdx][--leftIndex[tIdx]] = fm.remap_reverse[nextChar];
			bool nx = CompressedDfsButtumUpWithGap( nextMap, leftIdx, rightIdx, nextRowNum, tIdx, nxState[nextChar].nextRowDp, maxEditDistance, leftToRight, fm, path, idxPath,
			DArray, ( nxState[nextChar].acceptedValue == -1 ) ? lastAc : nxState[nextChar].acceptedValue, lower, cc );
			if( leftToRight )
				rightIndex[tIdx]--;
			else leftIndex[tIdx]++;
			if( nx ){
				return true;
			}
			
			//}
        }
    }
	if( isChecked && ( shTry || cn == false ) &&  curMap.acceptedValue != -1 && curMap.acceptedValue <= maxEditDistance && curMap.acceptedValue > lower ){
		bwtNode nextMap = curMap;
		bool nx = buttumUpMapping( nextMap, path, idxPath - 1, leftToRight, tIdx );
		if( nx ){
			return true;
		}
	}
	if( resPtr[tIdx] > preVal )
		useLess = false;
	/*if( useLess ){
		dbound[tIdx][mask][curRow] = min( dbound[tIdx][mask][curRow], cm );
		failedTry[tIdx][mask][P( curRow,cm )]++;
	}
	else
		failedTry[tIdx][mask][P( curRow,cm )] = 0;
	*/
	return false;
}

bool topDownMapping( int leftIdx, int rightIdx, int editDistance, int leftToRight, int tIdx, vector<pInterval> path, int curMask ){
	if( lp[tIdx] > loopTh )
		return false;
    int mid = ( leftIdx + rightIdx ) / 2;
    int minEditLeft = getMinEdit( leftIdx, mid, tIdx );
    int minEditRight = getMinEdit( mid + 1, rightIdx, tIdx );
    int minEditAll = getMinEdit( leftIdx, rightIdx, tIdx );
    if( editDistance < minEditAll )
        return false;
    if( editDistance > MAXD ){
        {
            intervalNode left( leftIdx, mid, editDistance / 2, 1, curMask * 2 );
            intervalNode right( mid + 1, rightIdx, editDistance, 0, 0 );
            path.push_back( pInterval( left, right ) );
            if( editDistance / 2 >= minEditLeft ){
                bool leftRes = topDownMapping( left.left, left.right, left.editDistance, left.leftToRight, tIdx, path, curMask * 2 );
                if( leftRes )
                    return true;
            }
        }
		path.pop_back();
		{
            intervalNode left( leftIdx, mid, editDistance, 1, 0 );
            intervalNode right( mid + 1, rightIdx, editDistance / 2, 0, curMask * 2 + 1 );
            path.push_back( pInterval( right, left ) );
            if( editDistance / 2 >= minEditRight ){
                bool rightRes = topDownMapping( right.left, right.right, right.editDistance, right.leftToRight, tIdx, path, curMask * 2 + 1 );
                if( rightRes )
                    return true;
            }
        }
        return false;
    }
	indexesPtr[tIdx] = 0;
    int len = rightIdx - leftIdx + 1;
    string q = readsPerThread[tIdx].substr(leftIdx, rightIdx - leftIdx + 1 );
    if( leftToRight == 1 ){
		buildLevstein(q, false, tIdx, revIdx );
    }
    else{
		buildLevstein(q, true, tIdx, fmIdx );
    }
    for( int i = 0; i < indexesPtr[tIdx]; i++ ){
		for( int j = leftIdx; j <= rightIdx; j++ ){
			nonComplete[tIdx][rightIndex[tIdx]++] = readsPerThread[tIdx][j];
		}
        bool nx = buttumUpMapping(indexes[tIdx][i], path, path.size() - 1, leftToRight, tIdx);
        if( nx ){
            return true;
        }
		leftIndex[tIdx] = rightIndex[tIdx] = 2 * MAXINPLEN;
    }
    return false;
}

int errorThread[MAXTHREADS];

int getLenCigar( string cigar ){
	int ret = 0;
	for( int i = 0; i < cigar.length(); i++ ){
        int lenCheck = 0;
        int j = i;
        if( cigar[i] <= '9' && cigar[i] >= '0' ){
            while( j < cigar.length() && cigar[j] <= '9' && cigar[j] >= '0' ){
                lenCheck = 10 * lenCheck + ( cigar[j] - '0' );
                j++;
            }
        }
		i = j;
		if( cigar[i] == 'M' || cigar[i] == 'D' )
			ret += lenCheck;
	}
	return ret;
}

psu getCigar( long long refPos, string read, int diff, int tIdx, long long L = -1, long long R = -1, bool canCutFirst = true, bool canCutEnd = true ){

    /*if( gap[tIdx] == 0 ){
        ostringstream sout;
        sout << read.length() << 'M';
        return psu( sout.str(), 0 );
    }*/
	errorThread[tIdx] = 0;
    string refSub = "";
	long long leftPos = max( 0LL, refPos - diff );
	long long rightPos = min( (unsigned long long)Ref.sz, refPos + read.length() + diff );
	if( L != -1 )
		leftPos = L, rightPos = R;
	leftPos = max( leftPos, 1LL );
	rightPos = min( rightPos, (long long)Ref.sz );
    for( long long i = leftPos; i < rightPos; i++ )
	    refSub += Ref.charAt( i );
    int len = read.length();
    string ret;
	//cerr << refSub.length() << endl << read.length() << endl;
    int curRow = getBestStartingPos( refSub, read, diff, canCutFirst, canCutEnd );
	//cerr << refPos << ' ' << diff << ' ' << curRow << endl;
	//cerr << "HERE " << curRow << endl;
    pss hs = Hirschberg( refSub, read, diff );
    for( int i = 0; i < hs.first.length(); i++ ){
    	char c1 = hs.first[i];
    	char c2 = hs.second[i];
    	if( c1 == '-' )
    		ret += 'D';
    	else if( c2 == '-' )
    		ret += 'I';
    	else ret += 'M';
		//cerr << c1 << ' ' << c2 << endl;
		if( c1 != c2 )
			errorThread[tIdx]++;
    }
	/*int ps = 0;
	while( ps < ret.length() && ret[ps] != 'M' )
		ret[ps++] = 'M';
	ps = ret.length() - 1;
	while( ps >= 0 && ret[ps] != 'M' )
		ret[ps--] = 'M';
    */
	string fret;
    for( int i = 0; i < ret.length(); i++ ){
        int j = i;
        while( j < ret.length() && ret[j] == ret[i] )
            j++;
        fret += std::to_string( j - i );
        fret += ret[i];
        i = j - 1;
    }
    return psu( fret, curRow ); 
}

string samFormatString( string& qName, string& read, uint32_t refPos, int len, int acceptedValue, int flag, bool Cigar, string& QC, psu cigar = psu( "", 0 ) ){	
    string ret;
    ret += ( ( qName.length() > 0 ) ? qName : "*" ) + '\t';   // QNAME field
    ret += std::to_string( flag ) + '\t';	// FLAG field
    if( flag == 4 ){
        ret += ( "*\t0\t255\t*\t*\t0\t0\t" + read + "\t*" );
        return ret;
    }
    int dif = ( acceptedValue == -1 ) ? ((( globalLongReadNoice * read.length() ) / 100 ) ) : ( acceptedValue );
    
    string tp;
    //psu ps = getPosition( max( 0LL, (long long)refPos - (long long)dif ) + cigar.second );
	psu ps = getPosition( refPos );
    ret += ps.first + '\t';	// RNAME field
    ret += std::to_string( ps.second ) + '\t';	// POS field
    ret += "255\t";	// MAPQ field
    ret +=  cigar.first + '\t';	// CIGAR field
    
    ret += "*\t";	// RNEXT field
    ret += "0\t";	// PNEXT field
    ret += "0\t";	// TLEN field
    ret += read + '\t';		// SEQ field
    if( QC == "" )
            ret += "*";	// QUAL field
    else ret += QC;
    return ret;
}

void threadWriterJob(){
	std::ofstream fout;
  	fout.open (outputAdr, std::ofstream::out | std::ofstream::app);
  	for( int i = 0; i < toWrite.size(); i++ ){
  		mappingInfo mi = toWrite[i];
		//if( mi.flag & 4 )
		//	continue;
        string smf = samFormatString( mi.qName, mi.read, mi.refPos, mi.len, mi.acceptedValue, mi.flag, true, mi.QC, mi.cigar );
		fout << smf << '\n';
  	}
  	toWrite.clear();
  	fout.close();
}

void prepareForWriting(){
	for( int i = 0; i < MAXTHREADS; i++ ){
		for( int j = 0; j < results[i].size(); j++ )
			toWrite.push_back( results[i][j] );
		results[i].clear();
	}
	sort( toWrite.begin(), toWrite.end() );
	threadWriterJob();
}

void buildMaskArray( int idx, int allowedEdit, string curRead ){
	if( gap[idx] == 0 )
		return;
    vector< P > qu[2];
    vector< int > mm[2];
    qu[0].push_back( P( 0, curRead.length() - 1 ) );
    mm[0].push_back( 1 );
    int turn = 0;
    int mxEr = 2 * allowedEdit;
    while( mxEr ){
        for( int ci = 0; ci < qu[turn].size(); ci++ ){
            int cL = qu[turn][ci].first;
            int cR = qu[turn][ci].second;
            int cM = mm[turn][ci];
            int mid = ( cL + cR ) / 2;
            qu[1 - turn].push_back( P( cL, mid ) );
            mm[1 - turn].push_back( 2 * cM );
            qu[1 - turn].push_back( P( mid + 1, cR ) );
            mm[1 - turn].push_back( 2 * cM + 1 );
            if( mxEr == 2 * allowedEdit )
                continue;
            cL = ( ci % 2 == 0 ) ? qu[turn][ci + 1].first : qu[turn][ci - 1].first;
            cR = ( ci % 2 == 0 ) ? qu[turn][ci + 1].second : qu[turn][ci - 1].second;
            int len = curRead.length();
            if( ci % 2 == 0 ){
                for( int ii = cL - mxEr; ii <= cR + mxEr; ii++ ){
                    for( int j = 0; j < 6; j++ ){
                        long long& ref = windowMask[idx][j][cM][ii + mxEr];
                        ref = 0;
                        for( int k = ii - mxEr; k <= ii + mxEr; k++ ){
                            ref <<= 1;
                            if( k >= cL && k <= cR && fmIdx.remap[curRead[k]] == j )
                                ref++;
                        }
                    }
                }
            }
            else{
                for( int ii = cL - mxEr; ii <= cR + mxEr; ii++ ){
                    for( int j = 0; j < 6; j++ ){
                        long long& ref = windowMask[idx][j][cM][ii + mxEr];
                        ref = 0;
                        for( int k = ii + mxEr; k >= ii - mxEr; k-- ){
                            ref <<= 1;
                            if( k >= cL && k <= cR && fmIdx.remap[curRead[k]] == j )
                                ref++;
                        }
                    }
                }
            }
        }
        qu[turn].clear();
        mm[turn].clear();
        turn = 1 - turn;
        mxEr /= 2;
    }
    qu[turn].clear();
    mm[turn].clear();
}

inline void initDbound( int idx ){
	/*if( gap[idx] ){
		for( int i = 0; i < MAXMASK; i++ )
			for( int j = 0; j < MAXMASK; j++ )
				dbound[idx][i][j] = 1000000;//.clear();
		for( int i = 0; i < MAXMASK; i++ )
			failedTry[idx][i].clear();
		
	}
	*/
	leftIndex[idx] = 2 * MAXINPLEN;
	rightIndex[idx] = 2 * MAXINPLEN;
}

bool cmpBwt( bwtNode& a, bwtNode& b ){
	return a.acceptedValue < b.acceptedValue;
}

vector< mappingInfo > getAllMapping( int idx, readMappingInfo& rm ){
	/*readMappingInfo rm( reads[ii], qc[ii], readsName[ii], globalMaxReport, globalBestOne, globalMaxDiffMismatch, globalMaxDiffEdit, globalGap, 
		globalUniqeOption, globalNoisePercent );
	*/	
	maxReport[idx] = rm.maxReport, bestOne[idx] = rm.bestOne, maxDiffMismatch[idx] = rm.maxDiffMismatch, maxDiffEdit[idx] = rm.maxDiffEdit;
	gap[idx] = rm.gap, uniqeOption[idx] = rm.uniqeOption, noisePercent[idx] = rm.noisePercent;
	//cerr << maxReport[idx] << ' ' << uniqeOption[idx] << endl;
    vector< mappingInfo > ret;
    cntDFS = 0;
    hit = 0;
    unhit = 0;
    forceStop[idx] = false;
    string qual = rm.qual;
    for( int j = 0; j < rm.read.length(); j++ ){
        bool v = false;
        for( int k = 0; k < 4; k++ )
            if( rm.read[j] == dna[k] )
                v = true;
        if( v == false )
            rm.read[j] = 'a';
    }
    string curRead = rm.read;
    string curReadRev = reverseComplement( rm.read );
    string curName = rm.readName;
    readsPerThread[idx] = rm.read;
    vector<pInterval> path;
    lp[idx] = 0;
    resPtr[idx] = 0;
    threadResults[idx].clear();
    int allowedEdit = 0;
    for( int j = 0; j < rm.read.length() * 4 + 1; j++ )
        preSearched[idx][j] = 0;
    if( noisePercent[idx] >= 0 ){
        noisePercent[idx] = min( noisePercent[idx], 0.2 );
        allowedEdit = (int)( (double)rm.read.length() * noisePercent[idx] );
    }
    else if( maxDiffEdit[idx] >= 0 )
        allowedEdit = maxDiffEdit[idx];
    else allowedEdit = 1;
	
	buildMaskArray( idx, allowedEdit,  readsPerThread[idx] );
	initDbound( idx );
    direction[idx] = 0;
    topDownMapping( 0, rm.read.length() - 1, allowedEdit, 1, idx, path, 1 );
    if( resPtr[idx] < maxReport[idx] && forceStop[idx] == false ){
        debg = 0;
        direction[idx] = 1;
        for( int j = 0; j < rm.read.length() * 4 + 1; j++ )
            preSearched[idx][j] = 0;
        path.clear();
        readsPerThread[idx] = reverseComplement( readsPerThread[idx] );
		buildMaskArray( idx, allowedEdit,  readsPerThread[idx] );
		initDbound( idx );
        topDownMapping( 0, rm.read.length() - 1, allowedEdit, 1, idx, path, 1 );
        readsPerThread[idx] = reverseComplement( readsPerThread[idx] );
    }
	//cerr << rm.read << ' ' << allowedEdit << ' ' << readsPerThread[idx] << ' ' << gap[idx] << ' ' << resPtr[idx] << endl;
    if( forceStop[idx] )
        resPtr[idx] = 0;
   
    if( resPtr[idx] == 0 ){
        ret.push_back( mappingInfo( curName, curRead, 0, 0, 0, 4, qual ) );
    }
    else if( bestOne[idx] ){
        sort( threadResults[idx].begin(), threadResults[idx].end(), cmpBwt );
		//int ls = 0;
		for( int i = 0; i < threadResults[idx].size(); i++ )
			if( threadResults[idx][i].acceptedValue <= 7 && threadResults[idx][i].acceptedValue >= 5 )
				ret.push_back( mappingInfo( curName, (threadResults[idx][i].flag == 0 ) ? curRead : curReadRev, fmIdx.locateRow( threadResults[idx][i].first ), 
				threadResults[idx][i].len, threadResults[idx][i].acceptedValue, threadResults[idx][i].flag, qual ) );
		/*while( ls < threadResults[idx].size() && abs( threadResults[idx][ls].acceptedValue - threadResults[idx][0].acceptedValue ) <= 0 ){
			ret.push_back( mappingInfo( curName, (threadResults[idx][ls].flag == 0 ) ? curRead : curReadRev, fmIdx.locateRow( threadResults[idx][ls].first ), 
			threadResults[idx][ls].len, threadResults[idx][ls].acceptedValue, threadResults[idx][ls].flag, qual ) );
			ls++;
		}*/
    }
    else{
        for( int i = 0; i < resPtr[idx]; i++ )
            for( uint32_t j = threadResults[idx][i].first; j <= min( threadResults[idx][i].first + maxReport[idx] - 1, threadResults[idx][i].second ); j++ )
                ret.push_back( mappingInfo( curName, ( threadResults[idx][i].flag == 0 ) ? curRead : curReadRev, fmIdx.locateRow( j ), threadResults[idx][i].len, threadResults[idx][i].acceptedValue, threadResults[idx][i].flag, qual ) );
    }
    return ret;
}

void shortReadsMapping( int idx, readMappingInfo& rm, int ii ){
	string mos = rm.read;
	string rev = reverseComplement( rm.read );
	minVal[idx] = 1000;
	if( globalMode == "fast" )
		mxFailed[idx] = 1000, Mode = 4, mc =  4;
	else if( globalMode == "normal" )
		mxFailed[idx] = 1000, Mode = 3, mc = 10;
	else if( globalMode == "sensitive" )
		mxFailed[idx] = 1000, Mode = 3, mc = 10;
	else mxFailed[idx] = 1000, Mode = 1000, mc = 1000;
	/*
	string read, qual, readName;
	int maxReport, bestOne, maxDiffMismatch, maxDiffEdit, gap, uniqeOption;
	double noisePercent = -1.0;
	*/
	//cerr << globalUniqeOption << endl;
	vector< mappingInfo > mps = getAllMapping( idx, rm );
	if( resPtr[idx] ){
		mch[idx]++;
	}
	for( int j = 0; j < mps.size(); j++ ){
		if( mps[j].flag != 4 ){
			mps[j].cigar = getCigar( mps[j].refPos, ( mps[j].flag == 0 ) ? mos : rev, mps[j].acceptedValue, idx  );
			mps[j].refPos = mps[j].refPos - mps[j].acceptedValue + mps[j].cigar.second;
		}
		mps[j].order = ii;
		results[idx].push_back( mps[j] );
	}
	//clock_t end = clock();  
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cerr << rc << ' ' << mps[0].flag << ' ' << elapsed_secs << ' ' << uptT << ' ' << cntBut << endl;
	//for( int i = 0; i < 2 * mos.length(); i++ )
	//	cerr << i << ' ' << preSearched[idx][i] << endl;
}

string tokenizeCigar( string& cigar ){
	string ret;
	for( int i = 0; i < cigar.length(); i++ ){
        int lenCheck = 0;
        int j = i;
        if( cigar[i] <= '9' && cigar[i] >= '0' ){
            while( j < cigar.length() && cigar[j] <= '9' && cigar[j] >= '0' ){
                lenCheck = 10 * lenCheck + ( cigar[j] - '0' );
                j++;
            }
        }
        i = j;
		for( int k = 0; k < lenCheck; k++ )
			ret += cigar[i];
	}
	return ret;
}

string packCigar( string& cigar ){
	string ret;
	for( int i = 0; i < cigar.length(); i++ ){
        int lenCheck = 0;
        int j = i;
        while( j < cigar.length() && cigar[j] == cigar[i] )
			j++;
		lenCheck = j - i;
        ret += ( std::to_string( lenCheck ) + cigar[i] );
		i = j - 1;
	}
	return ret;
}

pll getPositions( string cigar, long long refPos, long long readPos ){
	for( int i = 0; i < cigar.length(); i++ ){
		char ch = cigar[i];
		if( ch == 'M' )
			readPos++, refPos++;
		else if( ch == 'D' )
			refPos++;
		else if( ch == 'I' )
			readPos++;
	}
	return pll( refPos, readPos );
}

pll longReadAlignment( long long appPos, string& mos, int FLAG, string& rev, int segLen, int repairLen, int idx, string& repairCigar ){
	if( appPos - ( globalLongReadNoice * mos.length() ) / 100 < 0 )
		appPos = ( globalLongReadNoice * mos.length() ) / 100;
	long long appSt = appPos - ( globalLongReadNoice * mos.length() ) / 100;
	int tol = ( globalLongReadNoice * mos.length() ) / 100;
	
	long long preR = -1;
	long long realPos = -1;
	string finCig = "";
	string realRead = ( FLAG == 0 ) ? mos : rev;
	
	vector< pll > segs;
	for( int startPos = 0; startPos < realRead.length(); startPos += segLen ){
		int endPos = min( (int)realRead.length(), startPos + segLen );
		if( startPos )
			segs.push_back( pll( max( 0, startPos - repairLen ), min( startPos + repairLen, (int)realRead.length() - 1 ) ) );
		int rem = realRead.length() - endPos;
		string curSubRead = realRead.substr( startPos, endPos - startPos );
		long long L = ( startPos == 0 ) ? appSt : preR;
		long long R = L + 4 * ( endPos - startPos + 2 ) / 3;
		if( startPos == 0 )
			R += ( globalLongReadNoice * mos.length() ) / 100;
		psu curCig = getCigar( appPos, curSubRead, ( globalLongReadNoice * mos.length() ) / 100, idx, L, R, ( startPos == 0 ) ? true : false );
		if( startPos == 0 )
			realPos = L + curCig.second;
		finCig += curCig.first;
		preR = L + getLenCigar( curCig.first );
		if( startPos == 0 )
			preR += curCig.second;
	}
	if( realPos == -1 || realPos >= Ref.sz  )
		return pll( realRead.length(), 0 );
	string tokCig = tokenizeCigar( finCig );
	int segPtr = 0;
	int cigPtr = 0;
	long long readPos = 0;
	long long refPos = realPos;
	//cerr << "real = " << realPos << endl;
	while( readPos < realRead.length() ){
		
		if( segPtr < segs.size() && readPos == segs[segPtr].first ){
			long long lRef = refPos;
			while( readPos <= segs[segPtr].second ){
				char ch = tokCig[cigPtr++];
				if( ch == 'M' )
					readPos++, refPos++;
				else if( ch == 'D' )
					refPos++;
				else if( ch == 'I' )
					readPos++;
			}
			long long rRef = refPos;
			string sread = realRead.substr( segs[segPtr].first, segs[segPtr].second - segs[segPtr].first + 1 );
			string curCig = getCigar( 0, sread, 0, idx, lRef, rRef, false, false ).first;
			repairCigar += tokenizeCigar( curCig );
			segPtr++;
			continue;
		}
		char ch = tokCig[cigPtr++];
		repairCigar += ch;
		if( ch == 'M' )
			readPos++, refPos++;
		else if( ch == 'D' )
			refPos++;
		else if( ch == 'I' )
			readPos++;
	}
	
	readPos = 0;
	refPos = realPos;
	/*cerr << getPosition( realPos ).first << ' ' << getPosition( realPos ).second << endl;
	for( long long i = 0; i < 200; i++ )
		cerr << Ref.charAt( refPos + i );
	cerr << endl;
	*/
	int mpError = 0;
	for( int i = 0; i < repairCigar.length(); i++ ){
		char ch = repairCigar[i];
		if( ch == 'M' ){
			if( realRead[readPos] != Ref.charAt( refPos ) )
				mpError++;
			readPos++, refPos++;
		}
		else if( ch == 'D' )
			refPos++, mpError++;
		else if( ch == 'I' )
			readPos++, mpError++;
	}	
	return pll( mpError, realPos );
}

void longReadsMapping( int idx, readMappingInfo& rm, int ii ){
	//if( totReads <= 700000 )
	//	return;
	//cerr << ii << ' ' << rm.read << endl;
	string tmp = rm.read;
	string tmpQual = rm.qual;
	string tmpName = rm.readName;
	string mos = rm.read;
	string rev = reverseComplement( rm.read );
	
	if( tmp.length() + ( 4 * tmp.length() + 9 ) / 10 >= MAXINPLEN ){
		results[idx].push_back( mappingInfo( readsName[ii], reads[ii], 0, 0, 0, 4, qc[ii] ) );
		return ;
	}

	
		/*reads[ii] = tmp.substr( T, MAXINPLEN );
		qc[ii] = tmpQual.substr( T, MAXINPLEN );
		readsName[ii] = ( T == 0 ) ? tmpName : ( tmpName + std::to_string( ( T * 2 ) / MAXINPLEN ) );
		*/
	int partSize = 32;
	int step = partSize / 4;
	int minMatched = 1;
	int segLen = fragLen;
	noisePercent[idx] = -1;
	maxReport[idx] = 3;
	uniqeOption[idx] = 0;
	gap[idx] = 1;
	Mode = 1000, mc = 1000;
	partSize = 20;
	maxDiffEdit[idx] = 1;
	mxFailed[idx] = 1000;
	minMatched = 2;
	step = 1;
	/*if( globalMode == "fast" ){
		segLen = 200;
	}
	else if( globalMode == "normal" ){
		segLen = 650;
	}
	else if( globalMode == "sensitive" ){
		segLen = 1100;
	}
	else{
		segLen = 3000;
	}*/
	int repairLen = segLen / 2;
	vector< int > startPos, tpos;
	for( int i = 0; i + partSize <= reads[ii].length(); i += step )
		startPos.push_back( i );
	/*int ll = 0, rr = startPos.size() - 1;
	while( ll <= rr ){
		tpos.push_back( startPos[ll] );
		if( rr > ll )
			tpos.push_back( startPos[rr] );
		ll++, rr--;
	}
	startPos = tpos;
	*/
	random_shuffle( startPos.begin(), startPos.end() );	
	vector< mappingInfo > curMapping;
	long long appPos = -1;
	set< long long > checkedPos;
	vector< long long > allP;
	int FLAG = 0;
	int maxCheck = 3000000;
	int curCheck = 0;
	int Lidx, Ridx;
	long long tapos = 0;
	for( int i = 0; i < startPos.size() && appPos == -1 && curCheck < maxCheck; i++ ){
		string subRead, subQc;
		curCheck++;
		for( int j = startPos[i]; j < startPos[i] + partSize; j++ )
			subRead += reads[ii][j], subQc += qc[ii][j];
		readMappingInfo rm( subRead, subQc, "11", maxReport[idx], 0, 0, maxDiffEdit[idx], gap[idx], 
		uniqeOption[idx], -1 );
		//readMappingInfo rm( subRead, subQc, "11", 1, 0, 0, maxDiffEdit[idx], 1, 
		//0, -1 );
		vector<mappingInfo> A = getAllMapping( idx, rm );
		for( int j = 0; j < A.size(); j++ ){
			if( A[j].flag & 4 )
				continue;
			A[j].startPos = startPos[i];
			/*set< long long >::iterator ps = checkedPos.lower_bound( (long long)A[j].refPos - 2LL * mos.length() );
			for( int k = 0; k < 3 && ps != checkedPos.end(); k++, ps++ ); 
			if( ps != checkedPos.end() && abs( (*ps) - (long long)A[j].refPos ) <= 2 * mos.length() )
				continue;
			*/
			int cnt = 0;
			int lidx = 1e9, ridx = 0;
			for( int k = 0; k < curMapping.size(); k++ ){
				if( curMapping[k].flag == A[j].flag && 
					abs( (long long)A[j].refPos - (long long)curMapping[k].refPos ) <=  ( 13 * reads[ii].length() ) / 10
					&& abs( A[j].startPos - curMapping[k].startPos ) >= partSize ){
						long long m1 = abs( (long long)A[j].refPos - (long long)curMapping[k].refPos );
						long long m2 = abs( A[j].startPos - curMapping[k].startPos );
					if(  m1 >= 2 * partSize && m2 >= 2 * partSize && (double)max( m1, m2 )/(double)min( m1, m2 ) <= 1.4 ){
						cnt++;
						lidx = min( lidx, min( A[j].startPos, curMapping[k].startPos ) );
						ridx = max( ridx, max( A[j].startPos, curMapping[k].startPos ) + partSize );
						if( cnt >= minMatched )
							break;
					}
				}
			}
			if( cnt >= minMatched ){
				Lidx = lidx;
				Ridx = ridx;
				FLAG = A[j].flag;
				if( FLAG == 0 ){
					//tapos = A[j].refPos;
					appPos = max( 0LL, (long long)(A[j].refPos - startPos[i]) );
				}
				else{
					
					appPos = max( 0LL, (long long)(A[j].refPos - ( reads[ii].length() - startPos[i] )) );
				}
				
				allP.push_back( appPos );
			}

		}

		for( int j = 0; j < A.size(); j++ ){
			curMapping.push_back( A[j] );
			//checkedPos.insert( A[j].refPos );
		}
	}
	if( appPos == -1 ){
		mappingInfo mi = mappingInfo( readsName[ii], reads[ii], 0, 0, 0, 4, qc[ii] );
		mi.order = ii;
		results[idx].push_back( mi );
	}
	else{
		//sort(
		string repairCigar = "";
		string realRead = ( FLAG == 0 ) ? mos : rev;
		pll alg = longReadAlignment( appPos, mos, FLAG, rev, segLen, repairLen, idx, repairCigar );
		int mpError = alg.first;
		long long realPos = alg.second;
		if( mpError <= 0.3 * realRead.length() ){
			//repairCigar = finCig;
			repairCigar = packCigar( repairCigar );
			mappingInfo mps = mappingInfo( readsName[ii] , (FLAG == 0 ) ? mos : rev, appPos, reads[ii].length(), -1, FLAG, qc[ii] );
			mps.refPos = realPos;
			mps.cigar.first = repairCigar;
			mps.order = ii;
			results[idx].push_back( mps );
			resPtr[idx] = 1;
			if( resPtr[idx] ){
				mch[idx]++;
			}
			return;
		}
		else{
			/*int lenH = Ridx - Lidx;
			Lidx = max( 0, Lidx - lenH / 2 );
			Ridx = min( (int)tmp.length(), Ridx + lenH / 2 );
			*/
			mos = mos.substr( Lidx, Ridx - Lidx );
			rev = rev.substr( rev.length() - Ridx, Ridx - Lidx );
			realRead = ( FLAG == 0 ) ? mos : rev;
			repairCigar = "";
			int LH = 0, RH = 0;
			if( FLAG == 0 ){
				qc[ii] = qc[ii].substr( Lidx, Ridx - Lidx );
				LH = Lidx, RH = tmp.length() - Ridx;
				appPos += Lidx;
			}
			else{
				qc[ii] = qc[ii].substr( qc[ii].length() - Ridx, Ridx - Lidx );
				reverse( qc[ii].begin(), qc[ii].end() );
				LH = tmp.length() - Ridx;
				RH = tmp.length() - LH - ( Ridx - Lidx );
				appPos += tmp.length() - Ridx;
			}
			alg = longReadAlignment( appPos, mos, FLAG, rev, segLen, repairLen, idx, repairCigar );
			mpError = alg.first;
			realPos = alg.second;
			//cerr << mpError << ' ' << mos.length() << ' ' << tmp.length() << endl;
			//cerr << "L and R " << Lidx << ' ' << Ridx << endl;
			if( mpError <= 0.3 * realRead.length() ){
				//repairCigar = finCig;
				repairCigar = packCigar( repairCigar );
				
				if( LH )
					repairCigar = std::to_string( LH ) + "H" + repairCigar;
				if( RH )
					repairCigar += ( std::to_string( RH ) + "H" );
				
				/*mappingInfo mps = mappingInfo( readsName[ii] , (FLAG == 0 ) ? tmp : reverseComplement( tmp ), appPos, tmp.length(), -1, FLAG, qc[ii] );
				mps.refPos = realPos - LH;
				*/
				
				mappingInfo mps = mappingInfo( readsName[ii] , (FLAG == 0 ) ? mos : rev, appPos, mos.length(), -1, FLAG, qc[ii] );
				mps.refPos = realPos ;
				
				mps.cigar.first = repairCigar;
				mps.order = ii;
				results[idx].push_back( mps );
				resPtr[idx] = 1;
				if( resPtr[idx] ){
					mch[idx]++;
				}
				return;
			}
			
		}
		
		mappingInfo mi = mappingInfo( readsName[ii], reads[ii], 0, 0, 0, 4, qc[ii] );
		mi.order = ii;
		results[idx].push_back( mi );
	}
	//cerr << "done\n";
}

void hardReadsMapping( int idx, readMappingInfo& Rm, int ii ){
	//cerr << idx << ' ' << read << endl;
	int stepBase = 5;
	int part = 20;
	int error = 1;
	int step = stepBase;
	int gapVal = 0;
	long long appPos = -1;
	string mos = Rm.read;
	string rev = reverseComplement(Rm.read );
	vector< mappingInfo > validMappings;
	mxFailed[idx] = 1000000;
	int totMapped = 0;
	int tryNow = 0;
	vector< long long > checked;
	if( globalMode == "fast" )
		stepBase = 6;
	else if( globalMode == "normal" )
		stepBase = 5;
	else if( globalMode == "sensitive" )
		stepBase = 3;
	else stepBase = 1;
	for( int i = 0; i + part <= Rm.read.length(); i += step ){
		string subRead = Rm.read.substr( i, part );
		string subQual = Rm.qual.substr( i, part );
		readMappingInfo rm( subRead, subQual, "0", 5, 0, 0, error, gapVal, 0, -1 );
		vector< mappingInfo > allMp = getAllMapping( idx, rm );
		for( int K = 0; K < allMp.size(); K++ ){
			mappingInfo curMp = allMp[K];
			//cerr << i << ' ' << curMp.refPos << ' ' << curMp.flag << endl;
			//break;
			if( curMp.flag == 4 ){
				//i += 2 * step;
				step += stepBase;
				//step *= 2;
				//step++;
				continue;
			}
			totMapped ++;
			bool cc = false;
			for( int k = 0; k < checked.size(); k++ )
				if( abs( (long long)(curMp.refPos) - checked[k] ) <= 2 * Rm.read.length() ){
					cc = true;
					break;
				}
			if( cc )
				continue;
			//if( getPosition( curMp.refPos ).first == "chr5" )
			//cerr << i << ' ' << getPosition( curMp.refPos ).first << ' ' << getPosition( curMp.refPos ).second << endl;
			step = stepBase;
			//cerr << i << ' ' << curMp.flag << ' ' << subRead << ' ' << curMp.refPos << endl;
			for( int k = 0; k < validMappings.size(); k++ ){
				if( validMappings[k].flag == curMp.flag && abs( (long long)curMp.refPos - (long long)validMappings[k].refPos ) <=  ( 12 * Rm.read.length() ) / 10 ){
					checked.push_back( curMp.refPos );
					tryNow = 1;
					int FLAG = curMp.flag;
					if( FLAG == 0 )
						appPos = curMp.refPos - i;
					else{
						appPos = curMp.refPos - ( Rm.read.length() - i );
					}
					psu cgr = getCigar( appPos, ( FLAG == 0 ) ? mos : rev, ( 20 * Rm.read.length() ) / 100, idx );
					//cerr << errorThread[idx] << ' ' << cgr.first << ' ' << cgr.second << endl;
					if( errorThread[idx] <= 0.15 * Rm.read.length() ){
						// readMappingInfo( string _read, string _qual, string _readName, int _maxReport, int _bestOne, int _maxDiffMismatch, 
						//int _maxDiffEdit,int _gap, int _uniqeOption, double _noicePercent ){
						
						resPtr[idx] = 1;
						if( resPtr[idx] ){
							mch[idx]++;
						}
						
						mxFailed[idx] = 3;
						readMappingInfo rm( mos, Rm.qual, Rm.readName, globalMaxReport, 0, 0, -1, 1, 0, 0.15 );
						vector< mappingInfo > mps = getAllMapping( idx, rm );
						for( int j = 0; j < mps.size(); j++ ){
							if( mps[j].flag != 4 ){
								mps[j].cigar = getCigar( mps[j].refPos, ( mps[j].flag == 0 ) ? mos : rev, mps[j].acceptedValue, idx  );
							}
							mps[j].order = ii;
							results[idx].push_back( mps[j] );
						}
						
						/*mappingInfo mps = mappingInfo( Rm.readName , (FLAG == 0 ) ? mos : rev, appPos, Rm.read.length(), -1, FLAG, Rm.qual );
						mps.cigar = cgr;
						mps.order = ii;
						results[idx].push_back( mps );
						//cerr << ii << ' ' << "YES " << totLoop << endl;
						*/
						return;
					}
					else break;
				}
			}
			validMappings.push_back( curMp );
		}
	}
	
	/*readMappingInfo rm( Rm.read, Rm.qual, Rm.readName, 2, globalBestOne, globalMaxDiffMismatch, globalMaxDiffEdit, 1, 
			1, 0.15 );
	shortReadsMapping( idx, rm, ii ); 
	cerr << ii << " NO " << tryNow << ' ' << totMapped << ' ' << resPtr[idx] << endl;
	*/
	mappingInfo mi = mappingInfo( Rm.readName, Rm.read, 0, 0, 0, 4, Rm.qual );
	mi.order = ii;
	results[idx].push_back( mi );
	
	
}

void threadAlignmentJob( int ii, int idx ){
	/*if( totReads - cntReads + ii <= 750000 )
		return;
   	if( totReads - cntReads + ii > 750000 )
		cerr << ii << endl << reads[ii] << ' ' << qc[ii] << endl << endl;
	*/
	cntBut = 0;
	if( running_mode == "hard" ){
		readMappingInfo rm( reads[ii], qc[ii], readsName[ii], globalMaxReport, globalBestOne, globalMaxDiffMismatch, globalMaxDiffEdit, globalGap, 
			globalUniqeOption, globalNoisePercent );
		hardReadsMapping( idx, rm, ii );
	}
    else if( running_mode == "short" ){
		readMappingInfo rm( reads[ii], qc[ii], readsName[ii], globalMaxReport, globalBestOne, globalMaxDiffMismatch, globalMaxDiffEdit, globalGap, 
			globalUniqeOption, globalNoisePercent );
		shortReadsMapping( idx, rm, ii ); 
    }
    else{
    	readMappingInfo rm( reads[ii], qc[ii], readsName[ii], globalMaxReport, globalBestOne, globalMaxDiffMismatch, globalMaxDiffEdit, globalGap, 
			globalUniqeOption, globalNoisePercent );
		longReadsMapping( idx, rm, ii ); 
    }
    
}

void threadLoop( int idx ){
	double avg = 0;
	double tt = 0;
    while( readPtr < cntReads ){
        mu.lock();
        int selected = readPtr++;
        mu.unlock();
        if( selected < cntReads ){
            rc = 0;
            threadAlignmentJob( selected, idx );
            avg += rc;
            tt += 1;
        }
    }
    avg /= tt;
}

void Alignment(){
    vector< std::thread > th;
    readPtr = 0;
    for( int i = 0; i < core; i++ )
        th.push_back( thread( threadLoop, i ) );
    for( int i = 0; i < th.size(); i++ )
        th[i].join();
}



void Aligner( string fastqAdr ){
    ifstream fin( fastqAdr.c_str() );
    resetControllers();
    string line;
    cntReads = 0;
    int tt = 0, tmpName = 0, lpNum = 0;
    auto T1 = std::chrono::high_resolution_clock::now();
    while( true ){
        for( int i = 0; i < seed2; i++ )
            writing[i] = reading[i] = false;
        cntReads = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        long long totbp = 0;
        int mxLen = 0;
        while( getline( fin, line ) ){
            if( line.length() == 0 )
                    break;
            if( line[0] == '@' && line.length() > 1 ){
                readsName[cntReads] = line.substr( 1 );
                istringstream sin( readsName[cntReads] );
                string tnp;
                sin >> tnp;
                readsName[cntReads] = tnp;
            }
            else{
                ostringstream sout;
                sout << tmpName++;
                readsName[cntReads] = sout.str();
            }
            getline( fin, line );
            for( int i = 0; i < line.length(); i++ ){
                if( line[i] >= 'A' && line[i] <= 'Z' )
                    line[i] = 'a' + ( line[i] - 'A' );
				if( line[i] != 'a' && line[i] != 'c' && line[i] != 'g' && line[i] != 't' )
					line[i] = 'a';
			}
            totbp += line.length();
            mxLen = max( mxLen, (int)line.length());
            reads[cntReads] = line;
            getline( fin, line );
            getline( fin, line );
			//if( totReads + cntReads >= 750000 )
			//	cerr << reads[cntReads] << endl;
            qc[cntReads++] = line;
            if( cntReads >= MAXREADSIZE || totbp > 100000000)
                break;
        }
		
        if( cntReads == 0 )
            break;
        totReads += cntReads;
        cerr << "\nProcessing " << cntReads << " reads  (" << totbp << "bp)" << "...\n";
        cerr << "max length = " << mxLen << endl;
        Alignment();
        long long matched = 0;
        for( int i = 0; i < core; i++ )
            matched += mch[i];
        tt++;
        auto t2 = std::chrono::high_resolution_clock::now();
        auto dd = std::chrono::duration_cast< std::chrono::seconds>(t2 - t1);
        cerr << "###### " << lpNum++ << ": mapping time = " << dd.count() << " seconds ######" << endl;
        double nesbat = (double)matched * 100. / (double)totReads; 
        cerr << "aligned reads : " << matched << "/" << totReads << " -- " << nesbat << "%\n\n";
        prepareForWriting();
		//break;
    }
    auto T2 = std::chrono::high_resolution_clock::now();
    auto DD = std::chrono::duration_cast< std::chrono::seconds>(T2 - T1);
    cerr << "Total time = " << DD.count() << "s" << endl;
    int matched = 0 ;
}

void loadRef( string pref ){
    string filename = pref + ".rinfo";
    FILE* fin2 = fopen( filename.c_str(), "rb" );
    if( fin2 == NULL ){
    	cerr << "can not open index files\n";
    	cerr << "Error: " << strerror(errno) << endl;
    	exit( 0 );
    }
    ifstream fin( pref + ".cinfo" );
    Ref.load( fin2 );
    refLen = Ref.sz;
    cerr << "ref length = " << Ref.sz << endl;
    fclose( fin2 );
    unsigned int pos;
    string tmp;
    while( fin >> tmp >> pos ){
        refNames.push_back( tmp );
        refOffSets.push_back( pos );
    }
    refOffSets.push_back( 4000000000ULL );
}

bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.length() >= suffix.length() &&
           str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

bool has_prefix(const std::string &str, const std::string &suffix)
{
    return str.length() >= suffix.length() &&
    	str.substr( 0, suffix.length() ) == suffix;
           
}

string getValidPrefix( string adr ){
	string folderPath;
	int ptr = adr.length() - 1;
	while( ptr >= 0 && adr[ptr] != '/' && adr[ptr] != '\\' )
		ptr--;
	folderPath = adr.substr( 0, ptr + 1 );
	string curPref = adr.substr( ptr + 1 );
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (folderPath.c_str())) != NULL) {
	  /* print all the files and directories within directory */
	  while ((ent = readdir (dir)) != NULL) {
	  	string curFilePath = ent->d_name;
    	string validSuf = ".rinfo";
        if( has_suffix( curFilePath, validSuf ) && has_prefix( curFilePath, curPref ) ) {
        	return folderPath + curFilePath.substr( 0, curFilePath.length() - 6 );
        }
	  }
	  closedir (dir);
	} else {
	  perror ("");
	}
    return "";
}

static void
print_usage(const char *program)
{
    fprintf(stderr, "\n" );
    fprintf(stderr, "USAGE: %s <Algorithm> [OPTIONS] <ref.fa> <query file>\n", program);
	fprintf(stderr, "  choose an algorithm for aignment: :[short, long, hard] [mandatory]\n");
    fprintf(stderr, "  ref.fa: prefix path to the directory containing the reference genome and its index files[mandatory]\n");
    fprintf(stderr, "  query file : file containing queries in fastq format[mandatory]\n");
    fprintf(stderr, "  OPTIONS: \n");
    fprintf(stderr, "\t -a : report all locations [default : disabled]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -b : report bestOne  [default : disabled]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -g : map with  indels [default : disabled]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -u : report only unique alignments [default : disabled]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -t val[int]  : number of threads [default : 1]\n");
    fprintf(stderr, "\t -k val[int]  : maximum number of alignments to be reported for each read [default : 1]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -e val[double]  : estimated error rate of reads [default : 0.05]{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -v val[int]  : maximum number of allowed mismatch (included gap if any) for read mapping{available for 'short' and 'hard' algorithms}\n");
    fprintf(stderr, "\t -o val[string]  : output sam file [defaul : report.sam]\n");
    fprintf(stderr, "\t -m choose sensitivity threshold {fast, normal, sensitive, very-sensitive}[default: normal] \n");
	fprintf(stderr, "\t -L fragment length in approximate global alignment{available for 'long' algorithm}\n");
    fprintf(stderr, "EXAMPLE: %s short -t 2 -u -e 0.1 -o output.sam /path/to/index/directory/hg19.fa query.fastq\n",program);
	fprintf(stderr, "EXAMPLE: %s hard -t 2 -u -e 0.1 -o output.sam /path/to/index/directory/hg19.fa query.fastq\n",program);
	fprintf(stderr, "EXAMPLE: %s long -t 2 -o output.sam /path/to/index/directory/hg19.fa query.fastq\n",program);
    fprintf(stderr, "\n");
    return;
}

int main(int argc, char** argv) {
	std::srand ( unsigned ( std::time(0) ) );
    vector< unsigned short > test;
    int32_t opt,nqrys,maxqry,i;
    char* idxname;char* qryname;
    FILE* f;
    uint8_t** queries;
    char buf[4096];
    uint32_t start,stop,cnt;
    string args[100];
	
    /* parse command line parameter */
    if (argc < 4) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    for( int i = 0; i < argc; i++ )
		args[i] = argv[i];
	running_mode = args[1];
	if( running_mode != "short" && running_mode != "long" && running_mode != "hard" ){
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	int argPtr = 2;
	bool isSetD = false;
	while( argPtr < argc && args[argPtr][0] == '-' ){
		if( args[argPtr].length() != 2 ){
			print_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
		char option = args[argPtr][1];
		argPtr++;
		if( argPtr >= argc ){
			print_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
		if( option == 'b' ){
			globalBestOne = true;
			globalMaxReport = MAXREPORTSIZE;
			continue;
		}
		else if( option == 'g' ){
			globalGap = 1;
			continue;
		}
		else if( option == 'h' ){
			hardToMap = 1;
			continue;
		}
		else if( option == 'a' ){
			globalMaxReport = MAXREPORTSIZE;
			continue;
		}
		else if( option == 'u' ){
			globalUniqeOption = 1;
			globalMaxReport = 2;
			continue;
		}
		string operand = args[argPtr++];
		if( option == 't' ){
			core = std::stoi( operand );
		}
		else if( option == 'k' ){
			globalMaxReport = std::stoi( operand );
		}
		else if( option == 'e' ){
			globalNoisePercent = std::stod( operand );
		}
		else if( option == 'v' ){
			globalMaxDiffEdit = std::stoi( operand );
			isSetD = true;
		}
		else if( option == 'L' ){
			fragLen = std::stoi( operand );
		}
		else if( option == 'o' ){
			outputAdr = operand;
		}
		else if( option == 'm' ){
			globalMode = ( operand );
			if( globalMode != "normal" && globalMode != "fast" && globalMode != "sensitive" && globalMode != "very-sensitive" ){
				print_usage(argv[0]);
				exit(EXIT_FAILURE);	
			}
		}
		
		else{
			print_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	if( globalNoisePercent < 0 && !isSetD )
		globalNoisePercent = 0.05;
	if( core < 0 )
		core = 1;
	if( core > MAXTHREADS )
		core = MAXTHREADS;
    idxname = qryname = NULL;
    /* read ref name */
	if( argPtr < argc )
		idxname = argv[argPtr++];
		
	/* read filenames */
    if(argPtr < argc) { 
        qryname = argv[argPtr++];
    }
	
    if(qryname==NULL || idxname == NULL ) {
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
    }
    {
        cerr << endl;
        cerr << "\t** Running options: **\n";
        cerr << "\tnumber of threads : " << core << endl;
		if( running_mode != "long" ){
			cerr << "\tmaximum number of alignment reports per read : " << globalMaxReport << endl;
			cerr << "\tbest one option = " << ( globalBestOne ? "enable" : "disable" ) << endl;
			if( globalNoisePercent >= 0 )
					cerr << "\terror rate = " << globalNoisePercent << endl;
			if( globalBestOne == 0 && globalNoisePercent < 0 )
					cerr << "\tnumber of allowed mismatch: " << ( globalMaxDiffEdit == -1 ? 1 : globalMaxDiffEdit ) << endl;
			cerr << "\tconsider indel in mapping: " << ( globalGap ? "Yes" : "No" ) << endl; 
		}
        cerr << "\toutput sam file adress = " << outputAdr << endl;
    }
	ofstream fout( outputAdr );
	fout.close();
    string refPrefix = idxname;
    refPrefix = getValidPrefix( refPrefix );
    //cerr << refPrefix << endl;
    string fastqAdr = qryname;
    string fwAdr = refPrefix + ".fm";
    string rvAdr = refPrefix + ".rev.fm";
    string rfInf = refPrefix;
    loadRef( rfInf );
    cerr << "--Ref loading is done!\n";
    /* load index */
    fmIdx.load(fwAdr);
    cerr << "--forward BWT loading is done!\n";
    revIdx.load( rvAdr );
    cerr << "--Reverse BWT loading is done!\n";
    totChar = fmIdx.sigma;
    Aligner( fastqAdr );
    return (EXIT_SUCCESS);
}

