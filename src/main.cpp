#include <iostream>
#include <map>
#include <algorithm> 
#include <vector>
#include <deque>
#include <string>
#include <google/sparse_hash_map>
#include <gatb/gatb_core.hpp>

//class nbSeq;
using namespace std;
using google::sparse_hash_map;
//static const size_t span = KSIZE_1;
int maxDepth = 10000;

sparse_hash_map <string,int> KmerHash;
sparse_hash_map < string, vector<int> > KmerLong;
vector <string> readsLong;
vector <string> vect;
Graph graph;
int kmerSize;
int seuil;
int cover;
int maxDiff;

/**creates complement DNA sequence**/
string complement(string dna){
	
	std::string rev = "";
	
	for(int i = 1;i <= dna.length();i++){
		
       char nucleotide = dna[dna.length()-i];
	   
       switch (nucleotide){
		   
		   case 'A': rev = rev+"T"; break;
		   case 'C': rev = rev+"G"; break;
		   case 'G': rev = rev+"C"; break;
		   case 'T': rev = rev+"A"; break;
       }
	   
     }
	 
	 return rev;
}

/** returns the first k-mer in lexicographic order
 * between the k-mer its self and its reverse complement**/
string complementCanonical(string k1mer){
	
	std::string revk1 = "";
	
	for(int i = 1;i <= k1mer.length();i++){
		
       char nucleotide = k1mer[k1mer.length()-i];
	   
       switch (nucleotide){
		   
		   case 'A': revk1 = revk1+"T"; break;
		   case 'C': revk1 = revk1+"G"; break;
		   case 'G': revk1 = revk1+"C"; break;
		   case 'T': revk1 = revk1+"A"; break;
       }
	   
     }
	 
	 if(k1mer.compare(revk1)<0){
	 
		 return k1mer;
		 
	 }else{
	 
		 return revk1;
	}

}

/** filpping pairs from the hash-map of k-mers into the multiset of frequencies**/
template<typename A, typename B>
pair<B,A> flip_pair(const pair<A,B> &p)
{

    return pair<B,A>(p.second, p.first);

}

/** construct the multiset of frequencies from the hash-map of k-mers**/
template<typename A, typename B>
multiset<pair<B,A> > flip_map(const sparse_hash_map<A,B> &src)
{

    multiset<pair<B,A> > dst;
    transform(src.begin(), src.end(), inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
    
}

/**computes the levenshtein distance between two sequences**/
template <class T> unsigned int levenshtein_distance(const T& s1, 
													 const T& s2)
{

	const size_t len1 = s1.size(), len2 = s2.size();
	vector<vector<unsigned int> > d(len1 + 1, vector<unsigned int>(len2 + 1));
 
	d[0][0] = 0;
	for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
	for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;
 
	for(unsigned int i = 1; i <= len1; ++i)
		for(unsigned int j = 1; j <= len2; ++j)
 
                      d[i][j] = std::min( std::min(d[i - 1][j] + 1,d[i][j - 1] + 1),
                                          d[i - 1][j - 1] + (s1[i - 1] ==  s2[j - 1] ? 0 : 1) );
                                          
	return d[len1][len2];
	
}

/**constructs the map arc frequency with the frequency of the canonical complement frequency**/
void mapFreqInsert(string k1mer, vector<pair <string,int> > & mapFreq){

	int freq = KmerHash[complementCanonical(k1mer)];
	
	if(freq != 0)
		mapFreq.push_back(pair<string,int> (k1mer,freq));
		
}

/** constructs the frequency map of output arcs**/
void outFreqConstr(string node, string letterA, vector<pair <string,int> > & mapFreq){

	if(letterA != "A") mapFreqInsert(node+"A", mapFreq);
	if(letterA != "C") mapFreqInsert(node+"C", mapFreq);
	if(letterA != "G") mapFreqInsert(node+"G", mapFreq);
	if(letterA != "T") mapFreqInsert(node+"T", mapFreq);
	
}

/** constructs the frequency map of input arcs**/
void inFreqConstr(string node, string letterB, vector<pair <string,int> > & mapFreq){

	if(letterB != "A") mapFreqInsert("A"+node, mapFreq);
	if(letterB != "C") mapFreqInsert("C"+node, mapFreq);
	if(letterB != "G") mapFreqInsert("G"+node, mapFreq);
	if(letterB != "T") mapFreqInsert("T"+node, mapFreq);

}

/** atomoical pattern research**/
std::pair<string,int> findAtomPattern(string pattern){
	
	string atomPattern = "";
	int length = 2;

	string smallPattern = "";
	string restPattern = "";
	int nbCopiesAtomic = 1;
	bool repAtom = false;
	
	while(!repAtom && length < pattern.length() ){

		smallPattern = pattern.substr(0,length); //le motif atomic possible
		nbCopiesAtomic = 1;
		bool repInt = true;
		restPattern = pattern.substr(length);

		while( repInt && restPattern.length() >= length ){;

			if(restPattern.find(smallPattern) != 0) repInt = false;
			nbCopiesAtomic++;
			restPattern = restPattern.substr(length);
		}
		if(restPattern.length() > 0 && smallPattern.find(restPattern)  != 0 ) repInt = false;
					
		if (repInt == true) {
				repAtom = true;
				atomPattern = smallPattern;
		}
		length++;
	}

		if(atomPattern.length() == 0) {
			atomPattern = pattern;
			nbCopiesAtomic = 1;
		}
		return std::pair<string,int>(atomPattern,nbCopiesAtomic);
		
}

/**computes the number of mismatch between two sequences**/
int difference(string & pattern, string & pattRead){

	int diff = 0;
	for(int i = 0;i<pattern.size();i++)
				if(pattern[i] != pattRead[i]) diff++;
	return diff;
	
}

/**computes the number of approximate copies in the long read at the left of the pattern**/
int nbCopExtLeft(string & readSeq, string & pattern, int & diffBefore, 
			int & firstLength){
	
	int i = 0;
	string pattRead;
	int pattLength = pattern.length();
	int readLength = readSeq.length();
	float realDiff = (diffBefore*100)/firstLength;
	
	while (readLength>=pattLength && realDiff<maxDiff ){
	
		pattRead = readSeq.substr(readLength-pattLength);
		int diff = diffBefore+levenshtein_distance(pattern,pattRead);
		diffBefore = diff;
		firstLength = firstLength+pattLength;
		realDiff = (diff*100)/firstLength;
		
		if(realDiff<maxDiff){
		
			i++;
			readSeq = readSeq.substr(0,readLength-pattLength);
			readLength = readSeq.length();
			
		}
		
	}
	
	return i;
	
}

/**computes the number of approximate copies in the long read at the right of the pattern**/
int nbCopExtRight(string & readSeq, string & pattern, int & diffBefore, 
				int & firstLength){
	
	int i = 0;
	string pattRead;
	int pattLength = pattern.length();
	int readLength = readSeq.length();
	float realDiff = (diffBefore*100)/firstLength;
	
	while (readLength>=pattLength && realDiff<maxDiff ){
	
		pattRead = readSeq.substr(0,pattLength);
		int diff = diffBefore+levenshtein_distance(pattern,pattRead);
		diffBefore = diff;
		firstLength = firstLength+pattLength;
		realDiff = (diff*100)/firstLength;
		
		if(realDiff<maxDiff){
		
			i++;
			readSeq = readSeq.substr(pattLength);
			readLength = readSeq.length();
			
		}
		
	}
	
	return i;
}

/**Checks if the pattern is present in a long read**/
sparse_hash_map<string, pair<int,string> > oneReadCheck(
				sparse_hash_map<string, pair <int,int> > patternsCycle, 
				string & read, string & k1start){
				
	//cout<<"read in one read check "<<read<<endl;
	//cout<<pattern<<" "<<k1start<<" "<<pos<<endl;
	sparse_hash_map<string, pair <int,int> >::iterator itPattCycle;
	sparse_hash_map<string, pair<int,string> > patternRead;
	sparse_hash_map<string, pair<int,string> >::iterator itPattReads;
	//int nbCopReel = 0;
	string pattern;
	int pos;
	int nbCop;
	int lengthP;
	int i; float realDiff = 0; int debut;
	string pattRead;
	string::size_type start = 0;
	
	while((start = read.find(k1start,start)) != string::npos){
	
	//for(int j=pos;j<read.size()-kmerSize;j++){
		//string kmer_read=read.substr(j,kmerSize);
			//if(kmer_read.compare(k1start) == 0){
		for(itPattCycle = patternsCycle.begin();itPattCycle != patternsCycle.end();++itPattCycle){
		
				pattern = itPattCycle->first;
				nbCop = itPattCycle->second.first;
				pos = itPattCycle->second.second;
				lengthP = pattern.length();
				//cout<<"One read search "<<pattern<<" "<<nbCop<<" "<<pos<<endl; 
				i = 0;
				debut = start-pos;
				
				if((debut>=0)&& ((debut+lengthP)<=read.length())){
				
					pattRead = read.substr(debut,lengthP);
					int diff = levenshtein_distance(pattern,pattRead);
					realDiff = (diff*100)/lengthP;
					//cout<<"Comp: "<<lengthP<<" "<<diff<<" "<<realDiff<<" ";
					
					if(realDiff<=maxDiff){
					
						std::pair<string,int> pattIndent=
										findAtomPattern(pattern);
						i = pattIndent.second;
						pattern = pattIndent.first;
						int firstLength = lengthP;
						lengthP = pattern.length();
						
						if(debut>=lengthP) {
						
							pattRead = read.substr(0,debut);
							i = i+nbCopExtLeft(pattRead,pattern,diff,firstLength);
							
						}
						
						if((debut+firstLength+lengthP)<=read.length() ) {
						
							int fin = firstLength - (firstLength%lengthP);
							pattRead = read.substr(debut+fin);
							i = i+nbCopExtRight(pattRead,pattern,diff,firstLength);
							
						}
					
						/*itPattReads = patternRead.find(pattern);
						if(itPattReads != patternRead.end()){
							if(itPattReads->second < i) itPattReads->second = i;
						}else patternRead.insert(pair<string,int> (pattern,i) );*/
					}
					
					itPattReads = patternRead.find(pattern);
					if(itPattReads != patternRead.end()){
					
						if(itPattReads->second.first < i) {
						
							itPattReads->second.first = i;							
							itPattReads->second.second = itPattCycle->first;
						}
						
					}else{
					
						pair<int,string> pattCycle = pair<int,string>(i,itPattCycle->first);
						patternRead.insert(pair<string,pair<int,string> > (pattern,pattCycle) );
						
					} 
				}	
			//}
		}		
		start++;	
		
	}
	//cout<<endl;
	return patternRead;
	
}

/**Computes the set of long reads in which the seed k-mer is present**/
sparse_hash_map<string, pair<int,string> > longReadCheck(
				sparse_hash_map<string, pair <int,int> > & patternsCycle, 
				string & k1start){
	
	sparse_hash_map<string, pair<int,string> > maxPatternReads;
	sparse_hash_map<string, pair<int,string> > patternsReads;
	sparse_hash_map<string, pair<int,string> >::iterator itPattReads;
	sparse_hash_map<string, pair<int,string> >::iterator itMaxPattReads;
	
	//cout<<"pattern: "<<pattern<<" "<<k1start<<" "<<pos<<endl;
	//int minDiff = pattern.length();
	
	sparse_hash_map< string, vector<int> >::iterator itMapLong;
	itMapLong = KmerLong.find(k1start);
	
	if(itMapLong != KmerLong.end()){
	
		int nbReads = itMapLong->second.size();
		vector<int> numeroReads = itMapLong->second;
		//cout<<"Nb Reads "<< nbReads<<endl;
		
		for(int i = 0;i<nbReads;i++){
		
			string read = readsLong[numeroReads[i]];
			
			patternsReads = oneReadCheck(patternsCycle,read,k1start);
			//cout<<"Long reads search "<<patternsReads.size()<<endl;
			
			for(itPattReads = patternsReads.begin(); itPattReads != patternsReads.end(); ++itPattReads){
			
				string patt = itPattReads->first;
				//cout <<" Pattern "<<patt<<endl;
				itMaxPattReads = maxPatternReads.find(patt);
				if(itMaxPattReads != maxPatternReads.end()){
				
					if(itPattReads->second.first > itMaxPattReads->second.first)
					
						itMaxPattReads->second = itPattReads->second;
					
				}else {
				
					maxPatternReads.insert(pair<string,pair<int,string> > (patt,itPattReads->second));
					//cout<<"insert "<<maxPatternReads.size()<<endl;
				}
			}
			//cout<<"Max long reads search "<<maxPatternReads.size()<<endl;
		}
		//cout<<" Comp: "<<pattern.length()<<" "<<minDiff<<endl;
	}//else cout<<"Not found"<<endl;
	
	return maxPatternReads;
	
}

/**Cleans the extra frequency in a cycle and 
*verifies if the remaining frequency corresponds to a Tandem Repeat Patter**/
void searchForTR(vector<string> cycle,	vector<pair <string,int> > cycleFreq,
				vector<pair <string,int> > inputFreq,
				vector<pair <string,int> > outputFreq ){
				
	//cout<<"Search for TR"<<endl;				
	int cycleSize = cycle.size();
	int freqArcSize = cycleFreq.size();
	if(cycleSize != freqArcSize) cout<<"ERROR!!"<<endl;

	std::vector<int> suppFreqPos(cycleSize);
	std::vector<int> suppFreqNeg(cycleSize);
	std::vector<int> negFreqNodes(cycleSize);
	
	for(int i = 0;i<cycleSize;i++){
	
		suppFreqPos[i] = 0;
		suppFreqNeg[i] = 0;
		negFreqNodes[i] = 0;
		
	}

	//compute for each node the in freq (suppFreqPos) and the out freq (suppFreqNeg) to discover the negative freq nodes
	for(int i = 0;i<cycleSize;i++){
	
		//cout<<cycle[i]<<" ";
		
		for(int j = 0;j<inputFreq.size();j++){
		
			if(inputFreq[j].first.substr(1).compare(cycle[i]) == 0) 
				suppFreqPos[i] = suppFreqPos[i]+inputFreq[j].second;
				
		}
		
		for(int j = 0;j<outputFreq.size();j++){
		
			if(outputFreq[j].first.substr(0,kmerSize-1).compare(cycle[i]) == 0) 
				suppFreqNeg[i] = suppFreqNeg[i]-outputFreq[j].second;
				
		}
		//cout<<suppFreqPos[i]<<" "<<suppFreqNeg[i]<<endl;
	}


	std::vector<int> freqReal(freqArcSize);
	int freqCumulate;

	/*****if indetifReads**/
	/*std::vector <string> possibleEntries;
	std::vector <string> possibleExits;
	std::vector <string> possiblePatterns;
	std::vector <int> possibleNbCopies;*/
	/*****end if indetifReads**/
	sparse_hash_map<string, pair <int,int> > patternsCycle; //patter,<nbCop,pos1stkmer>

	for(int i = 0;i<inputFreq.size();i++){
	
		for(int j = 0;j<outputFreq.size();j++){//for each pair of possible entries and exits

			freqCumulate = 0;
			int startNode = 0;
			int endNode = 0;

			//cout<<inputFreq[i].first<<" "<<inputFreq[i].second<<endl;
			//cout<<outputFreq[j].first<<" "<<outputFreq[i].second<<endl;
			for(int l = 0;l<cycleSize;l++){ //on calcule des point de depart
			
				if(inputFreq[i].first.substr(1).compare(cycle[l]) == 0) 
					startNode = l;
				if(outputFreq[j].first.substr(0,kmerSize-1).compare(cycle[l]) == 0) 
					endNode = l;
					
			}
			//cout<<startNode<<" "<<endNode<<endl;
			//remouve the min freq between the entry and exit=the number of times we go throught the repetiton
			int min;
			if(inputFreq[i].second>=outputFreq[j].second) {
				min = outputFreq[j].second;
			}else{
				min = inputFreq[i].second;
			}

			suppFreqPos[startNode] = suppFreqPos[startNode]-min;
			suppFreqNeg[endNode] = suppFreqNeg[endNode]+min;


			//remouve supp feq by the negative freq nodes!!à modifier
			for(int l = 0; l<freqArcSize ; l++){

						freqCumulate = freqCumulate + suppFreqPos[l] + suppFreqNeg[l];
						
						if (freqCumulate < 0) {
						
								negFreqNodes[l] = freqCumulate;
								freqCumulate = 0;
								
						}	
						//cout<<freqCumulate<<" ";					
						freqReal[l] = cycleFreq[l].second-freqCumulate;
			}
			//cout<<endl;
			int remainingExits = 0;
			for(int l = 0; l<cycleSize;l++){
			
				remainingExits = remainingExits + negFreqNodes[l];
				
			}
			
			if(freqCumulate>0 && remainingExits<0){//s'il reste des choses à enlever, 2nd pass
			
				for(int l = 0; l<cycleSize;l++){
				
							freqCumulate = freqCumulate+negFreqNodes[l];
							//cout<<freqCumulate<<" ";
							if (freqCumulate < 0) {
							
								negFreqNodes[l] = freqCumulate;
								freqCumulate = 0;
								
							}	
							freqReal[l] = freqReal[l]-freqCumulate;
				}
				//cout<<endl;
			
			}	
			//replace the freq of the entry and the exit on the nodes for the next possible cupple
			int freqIn = min;
			
			suppFreqPos[startNode] = suppFreqPos[startNode]+min;
			suppFreqNeg[endNode] = suppFreqNeg[endNode]-min;
			/*
			cout<<"Frequence clean: "<<cycle.size()<<endl;
			for(int l=0;l<freqReal.size();l++){
				cout<<freqReal[l]<<" ";
			}

			cout<<endl;*/
			//compute approximation of the freq by dividing it by the cover
			bool repetition = true;
			for(int l = 0;l<freqReal.size();l++){
			
				if(freqReal[l]<0) {
				
						freqReal[l] = 0;
						repetition = false; //if there is a node with freq 0 then there is no repetiton in the cycle
						
				}else {
				
						float fValeur = freqReal[l]/cover;
						float fDecimal = fValeur-floor(fValeur);
				   		if (fDecimal< 0.5)
				     		freqReal[l] = floor(fValeur);
				   		else
				     		freqReal[l] = ceil(fValeur);
				     		
				}
				
			}
			
			int nbCopies = 0;
			//cout<<"1. Division par couv "<<repetition<<endl;
			if(repetition) {
			
				/*cout<<"Frequence clean2: "<<endl;
				for(int l = 0;l<freqReal.size();l++){
					cout<<freqReal[l]<<" ";
				}
				cout<<endl;*/
				//test frequency pattern and create pattern
				nbCopies = freqReal[startNode];
				if(startNode == endNode){
						for(int l = 0;l<freqArcSize;l++) if (abs(freqReal[l]-nbCopies)>2) repetition = false;
				}else{
					if(startNode<endNode){
							for(int l = startNode+1;l<endNode;l++) if (abs(freqReal[l]-nbCopies)>2) repetition = false;
							for(int l = endNode;l<freqArcSize;l++)  if (abs(freqReal[l]-nbCopies)>3) repetition = false;		
							for(int l = 0;l<startNode;l++)  if (abs(freqReal[l]-nbCopies)>3) repetition = false;	
					}else{
							for(int l = startNode+1;l<freqArcSize;l++) if (abs(freqReal[l]-nbCopies)>2) repetition = false;						
							for(int l = 0;l<endNode;l++)  if (abs(freqReal[l]-nbCopies)>2) repetition = false;
							for(int l = endNode;l<startNode;l++)  if (abs(freqReal[l]-nbCopies)>3) repetition = false;		
					}
				}
		}
		
		if(nbCopies == 0) repetition = false;

		std::string pattern = "";
		//cout<<"2. Freq pattern "<<repetition<<endl;
		if(repetition) {

			//recuperation du pattern
			pattern = cycle[startNode];
			if(startNode<=endNode){
					for(int l = startNode+1;l<cycleSize;l++)
					pattern = pattern + cycle[l].substr(kmerSize-2);
								
					for(int l = 0;l<=endNode;l++)
					pattern = pattern + cycle[l].substr(kmerSize-2);
								
							
			}else{
				if(startNode<cycleSize-1)
					for(int l = startNode+1;l<cycleSize;l++)
						pattern = pattern + cycle[l].substr(kmerSize-2);
				for(int l=0;l<cycleSize;l++)
					pattern = pattern + cycle[l].substr(kmerSize-2);
				for(int l = 0;l<=endNode;l++)
					pattern = pattern + cycle[l].substr(kmerSize-2);
							
			}						

			//cout<<"REP FULL "<<pattern<<" ";
			sparse_hash_map< string, pair<int,int> >::iterator itPattCycle;
			itPattCycle = patternsCycle.find(pattern);
			
			if(itPattCycle != patternsCycle.end()){
			
				if(itPattCycle->second.first<nbCopies){
				
					itPattCycle->second.first = nbCopies;
					itPattCycle->second.second = startNode;
					
				}
				
			}else{
			
				pair<int,int> nbCopPos = pair<int,int>(nbCopies,startNode);
				patternsCycle.insert(pair<string, pair<int,int> >(pattern,nbCopPos));
				
			}
			
		}
			
		 // cout<<"3. Attomic pattern "<<repetition<<endl;
			/*if(repetition){
				string startSeq = inputFreq[i].first.substr(0,kmerSize-2);
				string endSeq = outputFreq[j].first.substr(1);*/


	/*****if indetifReads*/
				/*possibleEntries.push_back(startSeq);
				possibleExits.push_back(endSeq);
				possiblePatterns.push_back(pattern);
				possibleNbCopies.push_back(nbCopies*nbCopiesAtomic);*/
				
				
				/*if(nbCopiesAtomic*nbCopiesRead>1){
				 cout<<"Real Repeat REPETITION "<<" "<<pattern<<" "<<nbCopiesRead<<" "<<nbCopies<<" "<<nbCopiesAtomic<< " "<< freqIn <<" "<<startSeq<<" "<<endSeq<<endl;
				 }else{
				 	cout<<"REPETITION "<<" "<<pattern<<" "<<nbCopiesRead<<" "<<nbCopies<<" "<<nbCopiesAtomic<< " "<< freqIn <<" "<<startSeq<<" "<<endSeq<< endl;
				 	}*/
				
	/**/
				//initialization of frequency
				for(int l = 0;l<freqReal.size();l++){
					freqReal[l] = 0;
				}

			//}
			/*else{
						int p = 0;
						while(p<patternsInCycle.size()){
							if(	edge_in[i].first == ioArcs[p].first.first && 
									edge_in[i].second == ioArcs[p].first.second &&
									edge_out[j].first == ioArcs[p].second.first && 
									edge_out[j].second == ioArcs[p].second.second){
						cout<<"REPETITION LOST FREQUENCY " << patternsInCycle[p]
								<<" "<<ioArcs[p].first.first<<"->"<<ioArcs[p].first.second
 								<<" "<<ioArcs[p].second.first<<"->"<<ioArcs[p].second.second
								<<endl;
						patternsInCycle.erase(patternsInCycle.begin()+p);
						ioArcs.erase(ioArcs.begin()+p);
						}else p++;
						}

			}*/
		}
	}
	
	cout<<"nb patterns cycle " <<patternsCycle.size()<<endl;
	int nbCopiesRead;
	int nbCopiesAtomic;	
	int nbCopies;
	string pattern;
	string patternAtom;
	string patternCycle;
	
	if(patternsCycle.size()>0) {
	
			sparse_hash_map< string, pair<int,string> > patternReads;
			patternReads = longReadCheck(patternsCycle, cycleFreq[0].first);
			cout<<"nb patterns reads " <<patternReads.size()<<endl;
			
			sparse_hash_map< string, pair<int,string> >::iterator itPattReads;
			sparse_hash_map< string, pair<int,int> >::iterator itPattCycle;
			//int comp = 0;
			for(itPattReads = patternReads.begin();itPattReads != patternReads.end();++itPattReads){
			
				//comp++;
				//cout<<"OK"<<comp<<endl;
				pattern = itPattReads->first;
				nbCopiesRead = itPattReads->second.first;
				patternCycle = itPattReads->second.second;
				//itPattCycle = patternsCycle.find(pattern);
				//nbCopies = itPattCycle->second.first;
				//std::pair<string,int> pattCop = findAtomPattern(pattern);
				//patternAtom = pattCop.first;
				//nbCopiesAtomic = pattCop.second;
				//if(patternAtom.length() != 0)
				
					if(nbCopiesRead>1) {
					
						cout<<"Real REPETITION "<<" "<<pattern<<" "<<nbCopiesRead<<" "<<patternCycle<<endl;
						//<<patternAtom<<" "<<nbCopiesRead<<" "<<nbCopies<<" "<<nbCopiesAtomic<<endl;
						
					}else{
					
						cout<<"REPETITION "<<" "<<pattern<<" "<<nbCopiesRead<<" "<<patternCycle<<endl;
						//<<patternAtom<<" "<<nbCopiesRead<<" "<<nbCopies<<" "<<nbCopiesAtomic<< endl;
						
					}
					
			}
			
	}
					
}

/**Computes the extra frequency of a cycle once the bening and ending arcs were chosen**/
void computeExtraFreq(vector<string> cycle){
	
	/*for(int i = 0;i<cycle.size();i++)
		cout<<cycle[i]<<" ";
	cout<<endl;*/
	//cout<<"Compute Extra Freq"<<endl;
	vector<pair <string,int> > cycleFreq;
	vector<pair <string,int> > inputFreq;
	vector<pair <string,int> > outputFreq;
	string k1mer;
	int freq;
	int cycleSize=cycle.size();
	sparse_hash_map<string,int>::iterator k1_it;
	vector <pair<string,int> >::iterator it;
	
	k1mer = cycle[0]+cycle[1].substr(kmerSize-2,1);
	
	cycleFreq.push_back(std::pair<string,int> (k1mer, KmerHash[complementCanonical(k1mer)]));
	//output arcs
	outFreqConstr(cycle[0], cycle[1].substr(kmerSize-2,1), outputFreq);
	//input arcs
	inFreqConstr(cycle[0], cycle[cycleSize-1].substr(0,1), inputFreq);
	
	
	for (int i = 1;i<cycleSize-1;i++){
		
		//arcs of the cycle
		k1mer = cycle[i]+cycle[i+1].substr(kmerSize-2,1);
		cycleFreq.push_back(std::pair<string,int> (k1mer, KmerHash[complementCanonical(k1mer)]));
		
		//output arcs
		outFreqConstr(cycle[i], cycle[i+1].substr(kmerSize-2,1), outputFreq);
			
		//input arcs
		inFreqConstr(cycle[i], cycle[i-1].substr(0,1), inputFreq);
				
	}
	
	k1mer = cycle[cycleSize-1]+cycle[0].substr(kmerSize-2,1);
	cycleFreq.push_back(std::pair<string,int> (k1mer, KmerHash[complementCanonical(k1mer)]));
	//output arcs
	outFreqConstr(cycle[cycleSize-1], cycle[0].substr(kmerSize-2,1), outputFreq);
			
	//input arcs
	inFreqConstr(cycle[cycleSize-1], cycle[cycleSize-2].substr(0,1), inputFreq);
	
	/*cout<<"cycle frequencies"<<endl;
	for(it = cycleFreq.begin();it != cycleFreq.end();it++){
		cout<<it->first<<" "<<it->second<<endl;
	}
	cout<<"output frequencies"<<endl;*/
	int sumOut = 0;
	
	for(it = outputFreq.begin();it != outputFreq.end();it++){
	
		//cout<<it->first<<" "<<it->second<<endl;
		sumOut = sumOut+it->second;
		
	}
	//cout<<"output frequencies total "<<sumOut<<endl;
	//cout<<"input frequencies"<<endl;
	int sumIn = 0;
	
	for(it = inputFreq.begin();it != inputFreq.end();it++){
	
		//cout<<it->first<<" "<<it->second<<endl;
		sumIn = sumIn+it->second;
		
	}
	//cout<<"input frequencies total "<<sumIn<<endl;
	//cout<<"end frequencies"<<endl;
	
	
	searchForTR(cycle,cycleFreq, inputFreq, outputFreq);
	
}


/**search from cycles starting from a specific edge startNode->vertexNode**/
void findCycles(string & startNode, string & vertexNode,string & k1start, int freqStart, vector<string> & stack, 
				int & maxRepLenMin,int & maxRepLenMax,int & depth,bool & maxDepthReached){
	//cout<<"find cycle"<<endl;
	//cout<<graph.toString(vertexNode)<<endl;
	//put the current vertex in the stack
	stack.push_back(vertexNode);
	//cout<<vertexNode<<endl;
	//get its neighbors
	/*pour version de gatb-core 1.1.0*/
	Graph::Vector<Node> neighborsGraph = graph.successors<Node> (graph.buildNode ((char*)vertexNode.c_str()));
	
	/*pour version de gatb_core 1.2.1
	Node currentN=graph.buildNode ((char*)vertexNode.c_str());
	Graph::Vector<Node> neighborsGraph = graph.successors (currentN);*/
	
	//cout<<"stack: "<<stack.size()<<" max size"<<stack.max_size()<<endl;
	multimap<int,string> neighbors;
	string successor;
	string k1mer;
	int k1freq;
	
	for(int i = 0;i<neighborsGraph.size();i++){
	
		successor = graph.toString(neighborsGraph[i]);
		k1mer = vertexNode+successor.substr(kmerSize-2,1);
		sparse_hash_map<string,int>::iterator k1_it = 
						KmerHash.find(complementCanonical(k1mer));
						
		if(k1_it != KmerHash.end()){ //if the arc (k-mer) exists
				
			k1freq = k1_it->second;
				
			if(k1freq>=seuil/2){
			
				neighbors.insert(pair<int,string>(k1freq,k1mer));
				
			}
		}
	
	}
	
	//cout<<"neighbors "<<total<<" "<<neighborsGraph.size()<<endl;;
	if( ( !maxDepthReached && depth<maxDepth && stack.size()<=maxRepLenMax ) ||
		( maxDepthReached && stack.size()<=maxRepLenMin ) ){//the maximal length of cycles to be explored
		
		depth = depth+neighbors.size();
		
		for(multimap<int,string>::reverse_iterator it = neighbors.rbegin(); it != neighbors.rend(); ++it){
			
			
			//cout<<successor<<endl;
			k1mer = (*it).second;
			k1freq = (*it).first;
			successor = k1mer.substr(1);
			string k1merOrig = k1mer;
			k1mer = complementCanonical(k1mer);
			
			//if it has a frequency higher then the threshold and 
			//if it was already traversed in the cycle and
			//if it is below the start node in the stack 
			//we explore it
			
			if((std::find(stack.begin(), stack.end(), successor) == stack.end())&&
					((k1freq<freqStart) || 
						((k1freq == freqStart)&&(k1mer.compare(k1start)<=0)) ||
						KmerLong.find(k1merOrig) == KmerLong.end())){
								
					//cout<<"k+1 mer "<<k1mer<<" freq "<<k1freq<<endl;
			
					//if we reached the first node in the stack we have found a cycle
					if(startNode.compare(successor) == 0){
					
						std::vector<string>::iterator it = stack.begin();
						vector<string> cycle;
						cycle.push_back(startNode);
						
						while (it != stack.end()){
						
							cycle.push_back(*it++);
							
						}
						
						cout<<"Cycle found "<<cycle.size()<<" "<<depth<<endl;
						computeExtraFreq(cycle);
						/*for(int j = 0;j<cycle.size();j++){
							cout<<cycle[j]<<" ";
						}
						cout<<endl;*/
					}else{
					
						//the k-mer representing the arc we are about to take
						if(!maxDepthReached && depth>=maxDepth){
						
							maxDepthReached = true;
							
						}
						
						findCycles(startNode,successor,k1start,freqStart,stack,maxRepLenMin,maxRepLenMax,depth,maxDepthReached);
					}
			
			}
		}
	}
	
	stack.pop_back();
	
}

/**Initialisation of the data**/
int main(int argc, char **argv)
{

	string arg = "";
	vector<string> files;
	string	fileLongRead;
	int maxRepLenMax;
	int maxRepLenMin;
	int depth = 0;
	
	if(argc < 8){
	
		cout<< "Usage is -k <kmerLength> -t <frequency threshold> -s <short read files> -c <short read cover> -l <long read files> -d <difference threshold> -r <repetition max length>\n";
		cin.get();
		exit(0);
		
	}else {
		
		int i = 1;
		while( i < argc ){
		
			arg = string(argv[i]);
			
			if(arg == "-k") {			
				kmerSize = atoi(argv[i+1]); i++;
			}
			
			if(arg == "-t") {
				seuil = atoi(argv[i+1]); i++;
			}
			
			if(arg == "-s") {
			
				while (string(argv[i+1]) != "-c"){
					files.push_back(argv[i+1]);
					i++;
				}
				
			}
			
			if(arg == "-c") {
				cover = atoi(argv[i+1]); i++;
			}
			
			if(arg == "-l") {
				fileLongRead = argv[i+1]; i++;
			}	
					
			if(arg == "-d") {
				maxDiff = atoi(argv[i+1]); i++;
			}
			
			if(arg == "-r") {
				maxRepLenMin = atoi(argv[i+1]); i++;
			}
			
			if(arg == "-R") {
				maxRepLenMax = atoi(argv[i+1]); i++;
			}
			
			i++;
		}
			
	}

	cout<<maxRepLenMax<<" "<<maxRepLenMin<<endl;
	//maxRepLenMax = maxRepLen+2;
	/**construction of hash-map for short reads**/

	cout<<"Step 1: construction of hash-map for short reads"<<endl;
	
	sparse_hash_map<string,int>::iterator itMap;
	
	//shrot reads files iterator
	string fileT="../tests/pairedShortReads_1.fasta,../tests/pairedShortReads_2.fasta";
	const string& fileTest=fileT;
	//BankFasta b (fileT);
	//BankFasta::Iterator itSeq (b);
	IBank* bank = Bank::open (fileTest);
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	
	//cout<<"hash map bank size "<<b.getSize()<<endl;
	//short reads kmer model
	Kmer<64>::ModelDirect model (kmerSize);
	Kmer<64>::ModelDirect::Iterator itKmer (model);
	
	int nbKmers = 0;//nb of all kmers
	int nbSeq = 0;//nb of reads
	string readShort;
	//iterate the short reads
	for (itSeq.first(); !itSeq.isDone(); itSeq.next()){
			
		itKmer.setData (itSeq->getData());
		readShort = itSeq->toString();
		//iterate the kmers.
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{
		
			nbKmers++;
			//compute the minimum between the direct kmer and its reverse complement
			string k1merHash = complementCanonical(model.toString(itKmer->value()));
			
			itMap = KmerHash.find(k1merHash);
			
			if(itMap == KmerHash.end()){
				KmerHash.insert(std::pair<string,int>(k1merHash,1));
			}else{
				 KmerHash[k1merHash]++;
			} 
			
		}
		
		nbSeq++;
		if(nbSeq % 100000 == 0) cout<<nbSeq<<endl;
	}
	cout<<"last read "<<readShort<<endl;
	cout<<"number of unique kmers "<<KmerHash.size()<<endl;
	cout<<"number of all k-mers "<<nbKmers<<endl;
	cout<<"number of reads "<<nbSeq<<endl;

	/**construction of the vector with k-mers for short reads ordered by their frequency:
	 * from the most frequent to the less frequent
	 * if 2 k-mers have the same frequency then they are ordered in inverse lexicographic order
	 **/ 
 	cout<<"Step 2: construction of the vector with k-mers for short reads ordered by their frequency"<<endl;

    //multiset with key:frequency value:kmer
	multiset<pair<int, string> > dst = flip_map(KmerHash);
	multiset<pair<int,string> >::iterator itMulti;
	
	for (multiset<pair<int,string> >::reverse_iterator rit = dst.rbegin(); rit != dst.rend(); ++rit){
	
		if(rit->first>=seuil)
			vect.push_back(rit->second);
			
	}
	dst.clear();
	cout<<"vect size"<<vect.size()<<endl;
	/*
	for(int i = 0;i<vect.size();i++){
		cout<<vect[i]<<endl;
	}*/
	
	
	/**construction of the graph from short reads with k-1 length for k-mers**/
	
	cout<<"Step 3: construction of the graph"<<endl;
	
	
	//nb of threads for graph construction
	int threads = 2;
	
	//IBank* bank = new BankFasta (fileT);
	cout<<"Bank construction"<<endl;
	cout<<"graph bank size "<<bank->getSize()<<endl;
	//graph = Graph::create (bank, "-kmer-size %d -abundance-min %d -bloom cache -debloom original -nb-cores %d", 
					// kmerSize-1, 1, threads);
	graph = Graph::create (bank, "-kmer-size %d -abundance-min %d -nb-cores %d", 
					   kmerSize-1, 1, threads);
	/*
	// We get an iterator for all nodes of the graph.
    Graph::Iterator<Node> itGraph = graph.iterator<Node> ();
    // We loop each node. Note the structure of the for loop.
    for (itGraph.first(); !itGraph.isDone(); itGraph.next())
    {
        // The currently iterated node is available with it.item()
        // We dump an ascii representation of the current node.
        std::cout << graph.toString (itGraph.item()) << std::endl;
    }
	*/
	cout<<"k-1 mer size"<< graph.getKmerSize ()<<endl;
	cout<<"graph info "<<graph.getInfo()<<endl;
	
	

	
	
	
	/**construction of the hash-map for long reads**/

	cout<< "Step 4: construction of the hash-map for long reads" <<endl;
	sparse_hash_map< string, vector<int> >::iterator itMapLong;

	//long reads file iterator
	BankFasta bLong (fileLongRead.c_str());
	BankFasta::Iterator itSeqLong (bLong);
	
	nbKmers = 0;//nb of all k-mers
	nbSeq = -1;//nb of long reads
	
	vector<int> readPos;//vector containing the nb of reads where a k-mer appears
	
	//iterate the long reads
	for (itSeqLong.first(); !itSeqLong.isDone(); itSeqLong.next()){
		
		nbSeq++;
		readsLong.push_back(itSeqLong->toString());
		//cout<<itSeqLong->toString()<<endl;
		itKmer.setData (itSeqLong->getData());
		
		//iterate the kmers
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{			
			nbKmers++;
			//compute the minimum between the direct kmer and its reverse complement
			string k1merHash = model.toString(itKmer->value());
			//cout<<k1merHash<<endl;
			string k1merHashCanon = complementCanonical(model.toString(itKmer->value()));
			itMap = KmerHash.find(k1merHashCanon);
			if(itMap != KmerHash.end()){//if the k-mer exist in the short reads hash-map

				itMapLong = KmerLong.find(k1merHash);
			
				if(itMapLong == KmerLong.end()){//if the k-mer does not exist in the hash-map
			
					//create a vector with the index of the current read
					readPos.clear();
					readPos.push_back(nbSeq);
					KmerLong.insert(std::pair< string, vector<int> >(k1merHash,readPos));
					
				}else{
					if(itMapLong->second[itMapLong->second.size()-1] != nbSeq){
						itMapLong->second.push_back(nbSeq);
					}
				} 
			}
		}
		
		if(nbSeq % 1000 == 0) cout<<nbSeq<<" "<< KmerLong.size() <<" "<< nbKmers<<endl;
	}
	
	cout<<KmerLong.size()<<endl;
	cout<<nbKmers<<" "<<nbSeq;
	cout<<readsLong.at(readsLong.size()-1)<<endl;
	
/**cycle research in the short read graph**/

	cout<<"Step 5: cycle research"<<endl;

 //#pragma omp parallel num_threads(8)
 //{
	 	
	//#pragma omp parallel for ordered schedule(dynamic)
	//#pragma omp for
	for(int i = 0;i<vect.size()-2;i++){ 
		
		vector<string> stack;
		stack.reserve(maxRepLenMax);
		depth = 0;
		//each k-mer in the vector represent an arc in the graph of k-1 mers
		string k1Start = vect[i];
		string node1 = k1Start.substr(0,kmerSize-1);
		string node2 = k1Start.substr(1);
		
		itMap = KmerHash.find(k1Start);
		itMapLong = KmerLong.find(k1Start);
		
		if(itMap != KmerHash.end()&&itMapLong != KmerLong.end()){
		
			int freqStart = itMap->second;
			if(freqStart>=seuil){
			
				cout<<i<<" k1 start "<<k1Start<<" frequency start "<<freqStart<<endl;
				
				if(node1.compare(node2) == 0){
				
					vector<string> cycle;
					cycle.push_back(node1);
					cout<<"Cycle found "<<cycle.size()<<endl;
					//computeExtraFreq(cycle);
					
				}else {
				
					pair<string,int> atomPat = findAtomPattern(k1Start);
					int atomPatleg = atomPat.first.length();
					
				/*	if(atomPatleg>0&&atomPat.second>=2){
						
						atomPatleg=atomPatleg+2;
						findCycles(node1,node2,k1Start,freqStart,stack,atomPatleg,depth);

					}else{*/
					    bool maxDepthReached = false;
						findCycles(node1,node2,k1Start,freqStart,stack,maxRepLenMin,maxRepLenMax,depth,maxDepthReached);
				//  }
				}
				
			}
		}
		
		vector<string> stackCom;
		stackCom.reserve(maxRepLenMax);
		depth = 0;
		string k1StartCom = complement(k1Start);
		string node1Comp = k1StartCom.substr(0,kmerSize-1);
		string node2Comp = k1StartCom.substr(1);
		
		
		itMapLong = KmerLong.find(k1StartCom);
		
		if(itMap != KmerHash.end()&&itMapLong != KmerLong.end()){
		
			int freqStart = itMap->second;
			
			if(freqStart>=seuil){
			
				cout<<"k1 start Complement "<<k1StartCom<<" frequency start "<<freqStart<<endl;

				if(node1Comp.compare(node2Comp) == 0){
				
					vector<string> cycle;
					cycle.push_back(node1Comp);
					cout<<"Cycle found "<<cycle.size()<<endl;
					//computeExtraFreq(cycle);
					
				}else {
				
					pair<string,int> atomPat = findAtomPattern(k1Start);
					int atomPatleg = atomPat.first.length();
					/*if(atomPatleg>0&&atomPat.second>=2){
						
						atomPatleg=atomPatleg+2;
						findCycles(node1Comp,node2Comp,k1Start,freqStart,stackCom,atomPatleg,depth);

					}else{*/
						bool maxDepthReached = false;
						findCycles(node1Comp,node2Comp,k1Start,freqStart,stackCom,maxRepLenMin,maxRepLenMax,depth,maxDepthReached);
					//}
				}
			}
		}
	}
 //}
	
	return 0;
}
