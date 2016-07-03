#include <iostream>
#include <fstream>
#include <map>
#include <algorithm> 
#include <vector>
#include <deque>
#include <string>
#include <google/sparse_hash_map>
#include <gatb/gatb_core.hpp>
#include <seqan/align.h>

//using namespace seqan;
using namespace std;
using google::sparse_hash_map;

const char* tab = "\t";
const char* esp = " ";
string ssakePath="/home/andreea/src/LongReads/ssake_v3-8-2/SSAKE";
string mrepPath="/home/andreea/src/LongReads/mreps/mreps";
string resultPath="smallAssem";

	
string complement(string dna){
		
	std::string rev="";
	
	for(int i=1;i<=dna.length();i++){
		
       char nucleotide=dna[dna.length()-i];
	   
       switch (nucleotide){
		   
		   case 'A': rev=rev+"T"; break;
		   case 'C': rev=rev+"G"; break;
		   case 'G': rev=rev+"C"; break;
		   case 'T': rev=rev+"A"; break;
       }
	   
     }
	 
	 return rev;
}
void eraseSSAKEfile(string pattern){
	cout<<"----------Remouve temp files-----------"<<endl;
	string fileName;
	fileName=resultPath+pattern+".txt";
	remove( fileName.c_str() );
	fileName=resultPath+pattern+"Assam.contigs";
	remove( fileName.c_str() );
	fileName=resultPath+pattern+"Assam.log";
	remove( fileName.c_str() );
	fileName=resultPath+pattern+"Assam.short";	
	remove( fileName.c_str() );
	fileName=resultPath+pattern+"Assam.singlets";
	remove( fileName.c_str() );
	fileName=resultPath+pattern+"_ATR.txt";
	remove( fileName.c_str() );
	cout<<"----------Remouve temp files end-----------"<<endl;
}

void mrepDetect(string pattern, map<string,string> & contigs ){
	
	string fileS=resultPath+"Seeds.fasta";
	string fileC=resultPath+"Contigs.fasta";
	
	ofstream fileSeeds;
	fileSeeds.open (fileS.c_str(),std::ios_base::app);
	ofstream fileContigs;
	fileContigs.open (fileC.c_str(),std::ios_base::app);
	
	string fileName=resultPath+pattern+"Assam.contigs";
	string path=mrepPath;
	string command=path+" -res 1 -minperiod 2 -maxperiod 20 -fasta "+fileName+" > "+
			resultPath+pattern+"_ATR.txt";
	cout<<"-------"<<command<<"-------"<<endl;
	int retSys=std::system(command.c_str());
	cout<<"-------"<<retSys<<"-------"<<endl;
	
	fileName=resultPath+pattern+"_ATR.txt";
	ifstream mrepFile(fileName.c_str());
	string line;
	
	string firstW;
	int pos;
	string contigName;
	string contigString;
	string limitLine=" ";
	
	bool seed=false;
	bool contig=false;
	string posDeb;
	string posFin;
	
	for (int i=0;i<93;i++){
		limitLine=limitLine+"-";
	}
	for(int i=0;i<9;i++){
		getline(mrepFile,line);
	}
	
	while ( getline(mrepFile,line) ){
		seed=false;
		contig=false;
		pos=line.find(esp);
		firstW=line.substr(0,pos);
	

		if(firstW.compare("Processing")==0){

		  	line=line.substr(pos+1);
			pos=line.find(esp);
			line=line.substr(pos+1);
			contigName=line.substr(1,line.length()-2);
			contigName=">"+contigName;
			//cout<<"CONTIG NAME"<<contigName<<endl;
			
			contigString=contigs.find(contigName)->second;
			//cout<<"CONTIG String"<<contigString<<endl;
			
			for(int i=0;i<4;i++){
				getline(mrepFile,line);
			}
			
			pos=line.find(esp);
			firstW=line.substr(0,pos);
			
			while( firstW.compare("The")==0 ){
					getline(mrepFile,line);
					pos=line.find(esp);
					firstW=line.substr(0,pos);
			}
			
			if( firstW.compare("RESULTS:")==0 ){
				
				//cout<<"no repeat contig"<<endl;
			
			} else{
				
				for(int i=0;i<2;i++){
					getline(mrepFile,line);
				}
				while(line.compare(limitLine)!=0){
					
					pos=line.find("->");
					posDeb=line.substr(0,pos-2);
					line=line.substr(pos+1);
					pos=posDeb.find_last_of(esp);
					posDeb=posDeb.substr(pos+1);
					
					pos=line.find(":");					
					posFin=line.substr(0,pos-1);
					pos=posFin.find_last_of(esp);
					posFin=posFin.substr(pos+1);
					
					//cout<<posDeb<<" "<<posFin<<" "<< contigString.length()<<endl;
					if( (atoi(posDeb.c_str())==1) ||
						(atoi(posFin.c_str())==contigString.length()) )
						
						seed=true;
					else
						contig=true;
						
					getline(mrepFile,line);	
				}
				
			}
			
			
		}
		if( seed ){
							
			fileSeeds<<contigName<<endl;
			fileSeeds<<contigString<<endl;			
		}else if(contig){
			fileContigs<<contigName<<endl;
			fileContigs<<contigString<<endl;
		}
		
	}
	fileContigs.close();
	fileSeeds.close();
	eraseSSAKEfile(pattern);
	contigs.clear();
}

void ssake(string pattern) {
	
	string fileName=resultPath+pattern+".txt";
	string path=ssakePath;
	string command=path+" -f "+fileName+" -w 1 -b "+resultPath+pattern+"Assam";
	cout<<"-------"<<command<<"-------"<<endl;
	int retSys=std::system(command.c_str());
	cout<<"-------"<<retSys<<"-------"<<endl;
	
}

void ssakeExtendingSeeds(string pairedReads){
	
	string fileName=resultPath+"Seeds.fasta";
	
	string path=ssakePath;
	string command=path+" -f "+pairedReads+" -p 1 -w 1 -s "+fileName+" -b "+resultPath+"/ExtendedSeeds";
	cout<<"-------"<<command<<"-------"<<endl;
	int retSys=std::system(command.c_str());
	cout<<"-------"<<retSys<<"-------"<<endl;
	
	
}

void chargeContigs(string pattern,map<string,string> & contigs ){

	cout<<"-------Contigs charging-------"<<endl;
	string fileName=resultPath+pattern+"Assam.contigs";
	ifstream ssakeFile(fileName.c_str());
	string contigName;
	string contigString;
	while ( getline(ssakeFile,contigName) ){
		
		getline(ssakeFile,contigString);
		contigs.insert(pair<string,string>(contigName,contigString));
		
	}
	cout<<"-------Contigs charged: "<<contigs.size()<<"-------"<<endl;
	
}
/*
 int alignReads(string & pattern, string & read, bool initGapPattern){
	
	typedef seqan::String<char> TSequence;                             // sequence type
	typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;     // align type
    typedef seqan::Row<TAlign>::Type TRow;                 // gapped sequence type
	typedef seqan::Iterator<TRow>::Type TRowIterator;
	
	
	TSequence seq1;
    TSequence seq2;
	TAlign align;
  
	seqan::AlignmentStats stats;
	int alignLength;
	int matchingScore;
	
	
	seq1=pattern;
	seq2=read;
	
	int pattLength=pattern.length();
	
	resize(rows(align), 2);
	TRow & row1= row(align,0);	
	TRow & row2= row(align,1);
	
	assignSource(row1, seq1);
	assignSource(row2, seq2);
	int score;
	TRowIterator it;
	TRowIterator itEnd;
	int c=0;
	if(initGapPattern){
		score = seqan::globalAlignment(align, 
					seqan::Score<int, seqan::Simple>(1, -1, -1), 
					seqan::AlignConfig<true, false, true, false>());
		it = begin(row2);
		itEnd = end(row2);
		for (; it != itEnd; ++it)
		{
			if (isGap(it))
				c++;
			else
				break;
		}
	}else {
		score = seqan::globalAlignment(align, 
					seqan::Score<int, seqan::Simple>(1, -1, -1), 
					seqan::AlignConfig<false, true, false, true>());
		it = end(row1);
		itEnd = begin(row1);
		for (; it != itEnd; --it)
		{
			if (isGap(it))
				c++;
			else
				break;
		}
		
	}				
	
	
	
	alignLength=pattLength-c;
	if(alignLength>16) {
		
		int scoreVal = seqan::computeAlignmentStats(stats, align, 
					seqan::Score<int, seqan::Simple>(1, -1, -1));
		SEQAN_ASSERT_EQ(scoreVal, stats.alignmentScore);
	
		matchingScore=(stats.numMatches*100)/alignLength;
	}else matchingScore=0;
	
	return matchingScore;
	
}
*/
void patternReadFile(string pattern, vector<string> readsShort){
	

	cout<<"-------Searching for reads containing the pattern-------"<<endl;
	ofstream myfile;
	string fileName=resultPath+pattern+".txt";
	string comPattern=complement(pattern);
	myfile.open (fileName.c_str());
	int nbSeq=0;
	string readName;
	int pattLength=pattern.length();
	for (int i=0;i<readsShort.size();i++){	
		/*
		if(pattLength>17){
			
			int matchingScore=alignReads(pattern, readsShort[i],true);
			if(matchingScore<80) 
				matchingScore=alignReads(pattern, readsShort[i],false);
			if(matchingScore>80) {
				nbSeq++;
				ostringstream convert;
				convert << nbSeq;
				readName=">read"+convert.str();
				myfile <<readName<<endl;
				myfile << readsShort[i] << endl;
			}else {
				matchingScore=alignReads(comPattern, readsShort[i],true);
				if(matchingScore<80)
					matchingScore=alignReads(comPattern, readsShort[i],false);
				if(matchingScore>80) {
					nbSeq++;
					ostringstream convert;
					convert << nbSeq;
					readName=">read"+convert.str();
					myfile <<readName<<endl;
					myfile << complement(readsShort[i]) << endl;
				}
			}
		}else {*/
	
			if(readsShort[i].find(pattern) != string::npos ){
				nbSeq++;
				ostringstream convert;
				convert << nbSeq;
				readName=">read"+convert.str();
				//cout<<readName<<endl;
				myfile <<readName<<endl;
				myfile << readsShort[i] << endl;
			}else if(readsShort[i].find(comPattern) != string::npos){
				nbSeq++;
				ostringstream convert;
				convert << nbSeq;
				readName=">read"+convert.str();
				myfile <<readName<<endl;
				myfile << complement(readsShort[i]) << endl;
			}
		//}
			
	}
	cout<<"-------Number of reads containing the pattern "<<pattern<<" "<<nbSeq
	<<"-------"<<endl;
	myfile.close();
}

void recupMotifsAlgo(ifstream & file, set<string> & listeMotifs){
	
	cout<<"-------Patterns charging-------"<<endl;
	string line;
	string motif="";
	int maxTaille=0;
	string firstW;
	int pos;
	while ( getline(file,line) ){
		
		pos=line.find(esp);
		firstW=line.substr(0,pos);
		if(firstW.compare("Real")==0){
			line=line.substr(pos+1);
			pos=line.find(esp);
			line=line.substr(pos+2);
			pos=line.find(esp);
			motif=line.substr(0,pos);
			//cout<<"Motif: "<<motif<<endl;
			//motif=atomicCanonical(motif);
			//cout<<"Motif Compl: "<<motif<<endl; 
			//if(motif.length()>6){
				listeMotifs.insert(motif);
				if(maxTaille<motif.length()) maxTaille=motif.length();
			/*}else{
				line=line.substr(pos+1);
				pos=line.find(esp);
				line=line.substr(pos+1);
				motif=line;
				listeMotifs.insert(motif);
			}*/
			//cout<<"taile liste "<<listeMotifs.size()<<endl;
		}
		
		/*if(firstW.compare("REPETITION")==0){
			line=line.substr(pos+2);
			pos=line.find(esp);
			motif=line.substr(0,pos);
			pair<string,int> atom=findAtomPattern(motif);
			motif=atom.first;
			//cout<<"Motif: "<<motif<<endl;
			//motif=atomicCanonical(motif);
			//cout<<"Motif Compl: "<<motif<<endl; 
			if(motif.length()<21){
				listeMotifs.insert(atomicCanonical(motif));
				if(maxTaille<motif.length()) maxTaille=motif.length();
			}
			//cout<<"taile liste "<<listeMotifs.size()<<endl;
		}*/
		
	}
	cout<<"-------Motif le plus long "<<maxTaille<<"-------"<<endl;
}
int main( int argc, char **argv ){

	vector<string> readsShort;
	vector<string> fileShortRead;
	fileShortRead.push_back(argv[1]);
	fileShortRead.push_back(argv[2]);
	
	//short reads file iterator
	int nbSeq=0;
	BankFasta bShort (fileShortRead);
	BankFasta::Iterator itSeqShort (bShort);
	for (itSeqShort.first(); !itSeqShort.isDone(); itSeqShort.next()){
		
		nbSeq++;
		readsShort.push_back(itSeqShort->toString());
		//if(nbSeq % 1000==0) cout<<nbSeq<<endl;
	}
	cout <<"----------Number of short reads "<< nbSeq<<"-------"<<endl;
	

	ifstream algoFile(argv[3]);
	string pairedReads=string(argv[4]);
	
	string executionSuff=string(argv[5]);
	
	resultPath=resultPath+executionSuff+"/";
	
	set<string> ATRalgo;
	recupMotifsAlgo(algoFile,ATRalgo);
	cout<<"----------Total ATR algo : "<<ATRalgo.size()<<"-------"<<endl;
	
	set<string>::iterator itMotifAlgo;
	string pattern;
	int i=0;
	for(itMotifAlgo=ATRalgo.begin();itMotifAlgo!=ATRalgo.end();++itMotifAlgo){
		i++;
		
		pattern=*itMotifAlgo;
		cout<<"-------Pattern "<<i<<": "<<pattern<<"-------"<<endl;
		patternReadFile(pattern, readsShort);
		map<string,string> contigs;
		ssake(pattern);	
		chargeContigs(pattern, contigs );
		//cout<<contigs.size()<<endl;
		mrepDetect(pattern,contigs);
		
	}
	
	cout<<"-------Seeds extending-------"<<endl;
	ssakeExtendingSeeds(pairedReads);
	cout<<"-------Seeds extended-------"<<endl;
	return 0;
}
