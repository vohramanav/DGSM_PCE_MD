#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <boost/algorithm/string.hpp>

using namespace boost::algorithm;
using namespace std;
using std::ifstream;

const int MAX_CHARS_PER_LINE = 1024;

/* 
  Local functions 
*/

/* Parse a text file into a vector of strings, with each line as a
   separate string */
void getTextFile(char *filein, vector<string> &ftext) ;
/* Determine the start line and end line for a specific section in a
   kinetic model */
void getSection(vector<string> mech, string sname, int *sst, int *sen) ;

void getReacLines(vector<string> mech, int rst, int ren, vector<int> &reacIdx) ;

void checkSpecInList(string sline, vector<string> specDelList, vector<string> &validSpec) ;

void checkSpecInMech(vector<string> &mechfile, vector<int> reacIdx, vector<int> reacFR, 
		     int ispec, vector<string> &specList, vector<string> &delList) ;

void getReacString(string reacline, vector<string> &rcomp) ;

bool checkValidReac(vector<string> mechfile, int ireac, vector<string> specDelList) ;

string delSpecFromReac(string mline, vector<string> specDelList) ;

bool isSpecInVecStr(vector<string> rcomp, vector<string> specDelList) ;

bool isSpecNameInVecStr(string rcomp, vector<string> specDelList) ;

/* Returns the position of the first upper case letter in a string */
size_t firstCap(string strin);

void checkDup(vector<string> mechfile, vector<int> reacIdx, vector<int> reacFR, 
              vector<int> &reacDup);

#ifdef __cplusplus
extern "C" {
#endif

void TC_reduce (char *mechIn, char *slist, char *mechOut) ;
void TC_splitFR(char *minp, int *rlist) ;
int  TC_parseFG(char *minp);

#ifdef __cplusplus
}
#endif

#undef VERBOSE
#undef DEBUG


int TC_parseFG(char *minp) {

  int ierr=0;
  vector<string> mechfile, reacTxt, specDelList, arrhrev, listFixSpec;

  /* read mechanism from file and identify sections */
  getTextFile((char *)minp, mechfile) ;
#ifdef VERBOSE
  cout << "-> Done reading mechanism" << endl;
#endif
  
  vector<int> reacIdx;
  int ElemStart, ElemEnd, SpecStart, SpecEnd, TherStart, TherEnd, ReacStart, ReacEnd;
  getSection(mechfile,string("ELEM"),&ElemStart,&ElemEnd);
  getSection(mechfile,string("SPEC"),&SpecStart,&SpecEnd);
  getSection(mechfile,string("THER"),&TherStart,&TherEnd);
  getSection(mechfile,string("REAC"),&ReacStart,&ReacEnd);
  getReacLines(mechfile,ReacStart,ReacEnd,reacIdx);
  int Nreac = reacIdx.size()-1;

#ifdef VERBOSE
  cout<<"Elements  :"<<ElemStart<<"->"<<ElemEnd<<endl<<flush;
  cout<<"Species   :"<<SpecStart<<"->"<<SpecEnd<<endl<<flush;
  cout<<"Reactions :"<<ReacStart<<"->"<<ReacEnd<<endl<<flush;
  cout<<"Done itemizing reactions"<<endl<<flush;
#endif

  ofstream fileout;
  fileout.open("chem_mod.inp");
  ofstream fgout;
  fgout.open("fg.dat");
  ofstream fglowout;
  fglowout.open("fglow.dat");

  /* write elements */
  fileout << "ELEMENTS" <<endl;
  for ( int i = ElemStart; i<=ElemEnd; i++ )
    fileout << mechfile[i] << endl;
  fileout << "END" <<endl;
  
  /* write species */
  fileout << "SPECIES" <<endl;
  vector<string> delSpec;
  for ( int i = SpecStart; i<=SpecEnd; i++ ) 
    fileout << mechfile[i] << endl;
  fileout << "END" <<endl;

  /* write THERMO data */
  for ( int i = TherStart-1; i<=TherEnd; i++ )
    fileout << mechfile[i] << endl;
  fileout << "END" <<endl;
  
  /* write reactions */
  fileout << "REACTIONS" <<endl;
  for ( int i = 0; i< (int) reacIdx.size()-1; i++ ) {

    char *rline = new char[mechfile[reacIdx[i]].length()+1];
    strncpy(rline, mechfile[reacIdx[i]].c_str(), mechfile[reacIdx[i]].length());
    rline[mechfile[reacIdx[i]].length()] = '\0';
    string rlin(rline);
    delete [] rline;
    trim(rlin);

    vector< string > rsp;
    split( rsp, rlin, is_any_of(" "), token_compress_on );

    if ( ( (int) rsp.size() != 4 ) && ( (int) rsp.size() != 6 ) ) {
      ierr=1;
      cout<<"Reaction line " << i << " has an unknown number of strings: "<<(int) rsp.size()<<endl;
      cout<<mechfile[reacIdx[i]]<<endl;
      for ( int j = 0; j < rsp.size(); j++ )
        cout<<"\""<<rsp[j]<<"\""<<endl;
    }
    for ( int j = 0; j < 4; j++ ) fileout<<rsp[j]<<"   ";
    fileout<<endl;
    if ( (int) rsp.size() == 6 ) fgout<<i<<"   "<<rsp[4]<<"   "<<rsp[5]<<endl;

    /* Write auxiliary info */
    for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
      vector< string > auxline;
      trim(mechfile[j]);
      split( auxline, mechfile[j], is_any_of(" ,/"), token_compress_on );
      if (auxline.back()==string("")) auxline.pop_back();
      if ( auxline[0] == string("LOW") ) {
        if ( ( (int) auxline.size() != 4 ) && ( (int) auxline.size() != 6 )) {
          ierr=2;
          cout<<"Reaction " << i << " has an unknown number of LOW values "<<(int)auxline.size()-1<<endl;
          cout<<mechfile[j]<<endl;
          for ( int k = 0; k < auxline.size(); k++ )
            cout<<"\""<<auxline[k]<<"\""<<endl;
        }
        else {
          fileout<<auxline[0]<<" / "<<auxline[1]<<" "<<auxline[2]<<" "<<auxline[3]<<" / "<<endl;
          if ( (int) auxline.size() == 6 )
            fglowout<<i<<" "<<auxline[4]<<" "<<auxline[5]<<endl;
        }
      }
      else
        fileout << mechfile[j] << endl;
    }

  } 
  fileout << "END" <<endl;
  
  /* done */
  fileout.close() ;
  fgout.close() ;
  fglowout.close() ;

  return (ierr);
  
}

/* ----------------------------------------------------------------------------------- */
void TC_splitFR(char *minp,int *rlist) {

  vector<string> mechfile, reacTxt, specDelList, arrhrev, listFixSpec;

  /* read mechanism from file and identify sections */
  getTextFile((char *)minp, mechfile) ;
#ifdef VERBOSE
  cout << "-> Done reading mechanism" << endl;
#endif
  
  vector<int> reacIdx;
  int ElemStart, ElemEnd, SpecStart, SpecEnd, ReacStart, ReacEnd;
  getSection(mechfile,string("ELEM"),&ElemStart,&ElemEnd);
  getSection(mechfile,string("SPEC"),&SpecStart,&SpecEnd);
  getSection(mechfile,string("REAC"),&ReacStart,&ReacEnd);
  getReacLines(mechfile,ReacStart,ReacEnd,reacIdx);
  int Nreac = reacIdx.size()-1;

  ElemStart -= 1;  /* Adjust indices to account for instances  */
  ElemEnd   += 1;  /* when the keywords are on the same lines  */
  SpecStart -= 1;  /* as some of the elements or species names */
  SpecEnd   += 1;

#ifdef VERBOSE
  cout<<"Elements  :"<<ElemStart<<"->"<<ElemEnd<<endl<<flush;
  cout<<"Species   :"<<SpecStart<<"->"<<SpecEnd<<endl<<flush;
  cout<<"Reactions :"<<ReacStart<<"->"<<ReacEnd<<endl<<flush;
  cout<<"Done itemizing reactions"<<endl<<flush;
#endif

  /* read species that need to be present no matter what the reduced mechanism */
  getTextFile((char *)"setspec.txt", listFixSpec) ;
#ifdef VERBOSE
  cout << "-> Done reading list of species" << endl;
#endif

  /* read list of reactions */
  vector<int> reacFR(2*Nreac);
  for ( int i = 0; i<2*Nreac; i++ ) reacFR[i] = rlist[i];

  vector<int> reacDup;
  checkDup(mechfile,reacIdx,reacFR,reacDup);

  /* read Arrhenius factors */
  getTextFile((char *)"arrhrev.dat", arrhrev) ;

  ofstream fileout;
  fileout.open("chem_red.inp");

  /* ----------------------------------ELEMENTS---------------------------------- */
  //fileout << "ELEMENTS" <<endl;
  for ( int i = ElemStart; i<=ElemEnd; i++ )
    fileout << mechfile[i] << endl;
  //fileout << "END" <<endl;
  
  /* ----------------------------------SPECIES---------------------------------- */
  fileout << "SPECIES" <<endl;
  vector<string> delSpec;
  for ( int i = SpecStart; i<=SpecEnd; i++ ) {
    vector<string> validSpec;

    checkSpecInMech(mechfile,reacIdx,reacFR,i,validSpec,delSpec);

    for ( int j = 0; j < (int) validSpec.size(); j++ )
      fileout << validSpec[j] <<" " ;
    if (validSpec.size()>0) fileout << endl;

    for ( int j=0; j<validSpec.size(); j++ ) {
      vector<string>::iterator a1=std::find(listFixSpec.begin(),listFixSpec.end(),validSpec[j]);
      if ( a1 != listFixSpec.end() ) listFixSpec.erase(a1,a1+1);
    }

  }

  for ( int j = 0; j < (int) listFixSpec.size(); j++ ) {
    fileout << listFixSpec[j] <<" " ;
    vector<string>::iterator a1=std::find(delSpec.begin(),delSpec.end(),listFixSpec[j]);
    if ( a1 != delSpec.end() ) delSpec.erase(a1,a1+1);
  }
  if (listFixSpec.size()>0) fileout << endl;
  fileout << "END" <<endl;
  
  /* ----------------------------------REACTIONS---------------------------------- */
  fileout << "REACTIONS" <<endl;
  for ( int i = 0; i< (int) reacIdx.size()-1; i++ ) {

    if ( ( reacFR[i] == 0 ) && ( reacFR[Nreac+i] == 0 ) ) continue;
    
    vector<string> rcomp;
    getReacString(mechfile[reacIdx[i]], rcomp) ;
    
    if ( ( reacFR[i] == 1 ) && ( reacFR[Nreac+i] == 1 ) ){

      fileout << rcomp[0]<<"<=>"<<rcomp[1]<<" "<<rcomp[2]<<" "<<rcomp[3]<<" "<<rcomp[4]<<endl;

      /* write auxiliary info */
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
	      if ( ( mechfile[j].find("DUP") != string::npos) && (reacDup[i] == 0) ) continue;
	      string cleanLine = delSpecFromReac(mechfile[j],delSpec);
	      if ( cleanLine.size() > 0 )
	        fileout << cleanLine << endl ;
	    }
    }
    else if ( reacFR[i] == 1 ) {
      fileout << rcomp[0]<<"=>"<<rcomp[1]<<" "<<rcomp[2]<<" "<<rcomp[3]<<" "<<rcomp[4]<<endl;

      /* write auxiliary info */
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
	      if ( mechfile[j].find("REV") != string::npos) continue;
	      if ( ( mechfile[j].find("DUP") != string::npos) && (reacDup[i] == 0) ) continue;
	      string cleanLine = delSpecFromReac(mechfile[j],delSpec);
	      if ( cleanLine.size() > 0 )
	        fileout << cleanLine << endl ;
      }
    }
    else if ( reacFR[Nreac+i] == 1 ) {

      vector<string> asplit;
      split( asplit, arrhrev[i], is_any_of(" "), token_compress_on );
      fileout << rcomp[1]<<"=>"<<rcomp[0]<<" "<<asplit[1]<<" "<<asplit[2]<<" "<<asplit[3]<<endl;

      /* write auxiliary info */
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
	      if ( mechfile[j].find("REV") != string::npos) continue;
	      if ( ( mechfile[j].find("DUP") != string::npos) && (reacDup[i] == 0) ) continue;
	      string cleanLine = delSpecFromReac(mechfile[j],delSpec);
	      if ( cleanLine.size() > 0 )
	        fileout << cleanLine << endl ;
      }
    }
    //break;
  } 
  fileout << "END" <<endl;
  
  /* done */
  fileout.close() ;

  return ;
  
}

void checkDup(vector<string> mechfile, vector<int> reacIdx, vector<int> reacFR, 
              vector<int> &reacDup) {

  int Nreac = (int) reacIdx.size()-1 ;
#ifdef DEBUG
  cout << reacDup.size() <<endl;
#endif


  /* Check which reactions have the duplicate keyword */
  for ( int i = 0; i< Nreac; i++ ) {
    if ( reacIdx[i]+1 == reacIdx[i+1] ) {
      reacDup.push_back(0);
    }
    else {
      bool foundDUP = false;
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
        size_t sFound=mechfile[j].find("DUP") ;
        if ( sFound != string::npos) {
          reacDup.push_back(1);
          foundDUP = true;
	}
      }
      if ( !foundDUP ) reacDup.push_back(0);
    }
  }

#ifdef DEBUG
  cout << reacDup.size() <<endl;
  ofstream fileout;
  fileout.open("debug_dup1.dat");
  for ( int i = 0; i< Nreac; i++ )
    fileout<<i+1<<": "<<reacDup[i]<<endl;
  fileout.close();
#endif

  for ( int i = 0; i< Nreac; i++ ) {
    if ( reacDup[i] == 0 ) continue ;
    vector<string> ri;
    getReacString(mechfile[reacIdx[i]], ri);
    bool foundGoodDup = false ;
    for ( int j = 0; j< Nreac; j++ ) {
      if ( ( i==j ) || ( reacDup[j] == 0 ) ) continue;
      vector<string> rj;
      getReacString(mechfile[reacIdx[j]], rj);
      if ( ( ri[0] == rj[0] ) && ( ri[1] == rj[1] ) ) {
	//cout<< ri[0] <<"<->"<< rj[0]<<endl;
	//cout<< ri[1] <<"<->"<< rj[1]<<endl;
        /* Case 1 - both reactions are at least forward */
        if ( ( reacFR[i] == reacFR[j] ) && ( reacFR[i] == 1) ) {
          foundGoodDup = true;
          //cout << i+1 <<":" << reacFR[i] <<","<< reacFR[Nreac+i] <<endl;
          //cout << j+1 <<":" << reacFR[j] <<","<< reacFR[Nreac+j] <<endl;
          break ;
	}
        /* Case 2 - both reactions are at least reverse */
        if ( ( reacFR[Nreac+i] == reacFR[Nreac+j] ) && ( reacFR[Nreac+i] == 1) ) {
          foundGoodDup = true;
          //cout << i+1 <<":" << reacFR[i] <<","<< reacFR[Nreac+i] <<endl;
          //cout << j+1 <<":" << reacFR[j] <<","<< reacFR[Nreac+j] <<endl;
          break ;
	}
      } /* Found another reaction with the same reactants/products */
    }
    if ( !foundGoodDup ) {
      reacDup[i] = 0;
      //cout << i+1 << ": "<< reacFR[i] << " "<<  reacFR[Nreac+i] << endl;
    }
  }

#ifdef DEBUG
  fileout.open("debug_dup2.dat");
  for ( int i = 0; i< Nreac; i++ )
    fileout<<i+1<<": "<<reacDup[i]<<endl;
  fileout.close();
#endif

  return ;

}

void getReacString(string reacline, vector<string> &rcomp) {
 
  vector<string> rtmp1, rtmp2;

  split( rtmp1, reacline, is_any_of("<=>"), token_compress_on );
  rtmp1[0].erase(std::remove(rtmp1[0].begin(),rtmp1[0].end(),' '),rtmp1[0].end());

  split( rtmp2, rtmp1[1], is_any_of(" "), token_compress_on );
  int rcsize = rtmp2.size() ;
  for ( int j = rcsize-1; j>=0; j-- )  
    if ( rtmp2[j] == string("")) rtmp2.erase (rtmp2.begin()+j);
  rcsize = rtmp2.size() ;

  for ( int j = 1; j < rcsize-3; j++ ) rtmp2[0] += rtmp2[j];

  rcomp.resize(5);
  rcomp[0] = rtmp1[0];
  rcomp[1] = rtmp2[0];
  rcomp[2] = rtmp2[rcsize-3];
  rcomp[3] = rtmp2[rcsize-2];
  rcomp[4] = rtmp2[rcsize-1];

  return ;

}

void TC_reduce(char *mechIn, char *slist, char *mechOut) {

  vector<string> mechfile, specDelList;
  /* read mechanism from file */
  getTextFile(mechIn, mechfile) ;
  cout << "Done reading mechanism" << endl;
  /* read species list */
  getTextFile(slist,  specDelList) ;
  cout << "Done reading list of species" << endl;
  
  vector<int> reacIdx;
  int ElemStart, ElemEnd, SpecStart, SpecEnd, ReacStart, ReacEnd;
  getSection(mechfile,string("ELEM"),&ElemStart,&ElemEnd);
  getSection(mechfile,string("SPEC"),&SpecStart,&SpecEnd);
  getSection(mechfile,string("REAC"),&ReacStart,&ReacEnd);
  getReacLines(mechfile,ReacStart,ReacEnd,reacIdx);

#ifdef VERBOSE
  cout<<"Elements  :"<<ElemStart<<"->"<<ElemEnd<<endl<<flush;
  cout<<"Species   :"<<SpecStart<<"->"<<SpecEnd<<endl<<flush;
  cout<<"Reactions :"<<ReacStart<<"->"<<ReacEnd<<endl<<flush;
  cout<<"Done itemizing reactions"<<endl<<flush;
#endif

  ofstream fileout;
  fileout.open(mechOut);

  /* write elements */
  fileout << "ELEMENTS" <<endl;
  for ( int i = ElemStart; i<=ElemEnd; i++ )
    fileout << mechfile[i] << endl;
  fileout << "END" <<endl;
  
  /* write species */
  fileout << "SPECIES" <<endl;
  for ( int i = SpecStart; i<=SpecEnd; i++ ) {
    vector<string> validSpec;
    checkSpecInList(mechfile[i],specDelList,validSpec);
    for ( int j = 0; j < (int) validSpec.size(); j++ )
      fileout << validSpec[j] <<" " ;
    fileout << endl;
  }
  fileout << "END" <<endl;
  
  /* write reactions */
  fileout << "REACTIONS" <<endl;
  for ( int i = 0; i< (int) reacIdx.size()-1; i++ ) {

    vector<string> validSpec;
    if ( checkValidReac(mechfile,reacIdx[i],specDelList) ) {

      fileout << mechfile[reacIdx[i]] << endl;
      for ( int j = reacIdx[i]+1; j < reacIdx[i+1]; j++ ) {
	string cleanLine = delSpecFromReac(mechfile[j],specDelList);
	if ( cleanLine.size() > 0 )
	  fileout << cleanLine << endl ;
      }

    }

  } 
  fileout << "END" <<endl;
  
  /* done */
  fileout.close() ;

  return ;
  
}

void getTextFile(char *filein, vector<string> &ftext) {

  ifstream fin ;
  fin.open(filein); 
  if (!fin.good()) {
    cout << "getTextFile() Could not open "<<filein<<endl<<flush ;
    exit(0);
  }

  while (!fin.eof())
  {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    /* Replace tabs with spaces */
    for (int i=0; i<MAX_CHARS_PER_LINE; i++) if (buf[i]=='\t') buf[i]=32;
    /* Replace CR with spaces */
    for (int i=0; i<MAX_CHARS_PER_LINE; i++) if (buf[i]==13)   buf[i]=32;
    string bufstr(buf) ;
    
    /* Trim extra spaces */
    //trim(bufstr);
    size_t found=bufstr.find(string("!"));
    if ( found != string::npos ) {
      bufstr.erase (found, bufstr.length());
    }
    if ( bufstr.size() == 0 ) continue ;
    transform(bufstr.begin(), bufstr.end(),bufstr.begin(), ::toupper);
    ftext.push_back(bufstr);
  }
  fin.close();
  return ;

}

void getSection(vector<string> mech, string sname, int *sst, int *sen) {
  
  bool foundSection = false;
  int i = 0 ;
  *sst = *sen = 0;
  while ( ( i < (int) mech.size() ) && ( *sen == 0 ) ){
     
    size_t sFound=mech[i].find(sname) ;
    if ( sFound != string::npos) {
      foundSection = true ;
      *sst = i+1;
    }
    if ( foundSection ) {
      size_t eFound=mech[i].find(string("END")) ;
      if ( eFound != string::npos) {
	foundSection = false ;
	*sen = i-1;
      }
    }
    i++ ;
  }

  return ;
  
}
 
void getReacLines(vector<string> mech, int rst, int ren, vector<int> &reacIdx) {

  if ( ( rst == 0 ) && ( ren == 0) ) return ;

  for ( int i = rst; i <= ren; i++ ) {
    size_t rFound=mech[i].find(string("=")) ;
    if ( rFound != string::npos) reacIdx.push_back(i) ;
  }
  reacIdx.push_back(ren+1);

  return ;

}

void checkSpecInList(string sline, vector<string> specDelList, vector<string> &validSpec) {
  
  split( validSpec, sline, is_any_of(" ,/<=>"), token_compress_on );

  for ( int j = 0; j < (int) specDelList.size(); j++ )
    for ( int i = 0; i < (int) validSpec.size(); i++ )
      if ( validSpec[i] == specDelList[j] ) {
	validSpec.erase(validSpec.begin()+i);
        i--;
      }
    
  return ;

}

void checkSpecInMech(vector<string> &mechfile, vector<int> reacIdx, vector<int> reacFR, 
                      int ispec, vector<string> &specList, vector<string> &delList) {

  vector<string> slist;
  string tmpstr = mechfile[ispec];
  trim(tmpstr);
  split( slist, tmpstr, is_any_of(" "), token_compress_on );
  if ( slist.back() == string("END") ) slist.pop_back();

  for ( int i = 0; i < reacFR.size()/2; i++ ) {

    if ( ( reacFR[i] == 0 ) && ( reacFR[reacFR.size()/2+i]==0 ) ) continue ;

    vector<string> rcomp;
    getReacString(mechfile[reacIdx[i]], rcomp);
    string react=rcomp[0], prod=rcomp[1];
    rcomp.clear();

    string stmp;
    size_t findTrdB;

    /* Check third-body in reactants */
    findTrdB=react.find(string("(+"));
    if ( findTrdB != string::npos ) {
      stmp=react.substr(findTrdB) ;
      react.erase(react.begin()+findTrdB,react.end());
      stmp.erase(stmp.begin(),stmp.begin()+2);
      stmp.erase(stmp.end()-1,stmp.end());
      rcomp.push_back(stmp);
      if ( isSpecNameInVecStr(rcomp[0],slist) ) {
        if (find(specList.begin(),specList.end(),rcomp[0]) == specList.end())
	      specList.push_back(rcomp[0]);
      }
      stmp.clear();
      rcomp.clear();
    } 

    /* Check third-body in products */
    findTrdB=prod.find(string("(+"));
    if ( findTrdB != string::npos ) {
      stmp=prod.substr(findTrdB) ;
      prod.erase(prod.begin()+findTrdB,prod.end());
      stmp.erase(stmp.begin(),stmp.begin()+2);
      stmp.erase(stmp.end()-1,stmp.end());
      rcomp.push_back(stmp);
      if ( isSpecNameInVecStr(rcomp[0],slist) ) {
        if (find(specList.begin(),specList.end(),rcomp[0]) == specList.end())
	      specList.push_back(rcomp[0]);
      }
      stmp.clear();
      rcomp.clear();
    }

    split( rcomp, react, is_any_of("+"), token_compress_on );
    for ( int j = 0; j < rcomp.size(); j++ ) 
      if ( isSpecNameInVecStr(rcomp[j],slist) ) {
        string snm = rcomp[j].substr(firstCap(rcomp[j]));
        if (find(specList.begin(),specList.end(),snm) == specList.end())
	      specList.push_back(snm);
      }
    rcomp.clear();

    split( rcomp, prod, is_any_of("+"), token_compress_on );
    for ( int j = 0; j < rcomp.size(); j++ ) {
      if ( isSpecNameInVecStr(rcomp[j],slist) ) {
        string snm = rcomp[j].substr(firstCap(rcomp[j]));
        if (find(specList.begin(),specList.end(),snm) == specList.end())
	        specList.push_back(snm);
      }
    }
    rcomp.clear();

#ifdef INCLTHB
    for ( int k = reacIdx[i]+1; k < reacIdx[i+1]; k++ ) {
      split( rcomp, mechfile[k], is_any_of(" /"), token_compress_on );
      for ( int j = 0; j < rcomp.size()-1; j++ ) {
	      if ( isSpecNameInVecStr(rcomp[j],slist)) {
          string snm = rcomp[j].substr(firstCap(rcomp[j]));
          if (find(specList.begin(),specList.end(),snm) == specList.end())
	          specList.push_back(snm);
	      }
      }
    }
#endif

  }

  /* If a certain species in the original list was not selected, then
     add it to the list of species to be deleted */
  for ( int i = 0; i < slist.size(); i++ )
    if (find(specList.begin(),specList.end(),slist[i]) == specList.end())
      delList.push_back(slist[i]);

  return ;

}

bool checkValidReac(vector<string> mechfile, int ireac, vector<string> specDelList) {

  vector<string> rcomp;

  split( rcomp, mechfile[ireac], is_any_of(" <=>"), token_compress_on );
  string react=rcomp[0], prod=rcomp[1];
  rcomp.clear();


  string stmp;
  size_t findTrdB;

  /* Check third-body in reactants */
  findTrdB=react.find(string("(+"));
  if ( findTrdB != string::npos ) {
    stmp=react.substr(findTrdB) ;
    react.erase(react.begin()+findTrdB,react.end());
    stmp.erase(stmp.begin(),stmp.begin()+2);
    stmp.erase(stmp.end()-1,stmp.end());
    rcomp.push_back(stmp);
    if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
    stmp.clear();
    rcomp.clear();
  } 

  /* Check third-body in products */
  findTrdB=prod.find(string("(+"));
  if ( findTrdB != string::npos ) {
    stmp=prod.substr(findTrdB) ;
    prod.erase(prod.begin()+findTrdB,prod.end());
    stmp.erase(stmp.begin(),stmp.begin()+2);
    stmp.erase(stmp.end()-1,stmp.end());
    rcomp.push_back(stmp);
    if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
    stmp.clear();
    rcomp.clear();
  } 

  split( rcomp, react, is_any_of("+"), token_compress_on );
  if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
  rcomp.clear();

  split( rcomp, prod, is_any_of("+"), token_compress_on );
  if ( isSpecInVecStr(rcomp,specDelList) ) return ( false );
  rcomp.clear();

  return (true) ;

}

bool isSpecInVecStr(vector<string> rcomp, vector<string> specDelList) {

  for (int i = 0 ; i < (int) rcomp.size(); i++) {
    if (firstCap(rcomp[i]) > rcomp[i].size() ) {
      cout<<"error parsing reaction file for "<<endl<<rcomp[i]<<endl<<flush ;
    }
    string spec=rcomp[i].substr(firstCap(rcomp[i]));
    for  (int j = 0; j < (int) specDelList.size(); j++) {
      if (specDelList[j] == spec) return ( true ) ;
    }
  }

  return ( false );

}

bool isSpecNameInVecStr(string rcomp, vector<string> specDelList) {

  if (firstCap(rcomp) > rcomp.size() ) return (false) ;
  string spec=rcomp.substr(firstCap(rcomp));
  for  (int j = 0; j < (int) specDelList.size(); j++) {
    if (specDelList[j] == spec) return ( true ) ;
  }

  return ( false );

}

string delSpecFromReac(string mline, vector<string> specDelList) {

  vector<string> rcomp;
  string tmpstr = mline;
  trim(tmpstr);
  split( rcomp, tmpstr, is_any_of(" /"), token_compress_on );
  if (rcomp.back()==string("")) rcomp.pop_back();

  bool foundmatch=false;
  for (int i = 0 ; i < (int) rcomp.size()-1; i++) {
    for  (int j = 0; j < (int) specDelList.size(); j++ ) 
      if (specDelList[j] == rcomp[i]) {
        foundmatch=true;
        //cout<<"\""<<specDelList[j]<<"\"  \""<<rcomp[i]<<"\" "<<foundmatch<<endl ;
       }
  }


  if ( foundmatch ) {
    //for (int i = 0 ; i < (int) rcomp.size(); i++)
    //  cout<<i<<": "<<rcomp[i]<<endl;
    assert(rcomp.size()%2==0);
    string rstr;
    for (int i = 0 ; i < (int) rcomp.size()/2; i++) {
      foundmatch=false;
      for  (int j = 0; j < (int) specDelList.size(); j++ ) 
	      if (specDelList[j] == rcomp[2*i]) foundmatch=true;
      if ( !foundmatch ) {
        rstr.append(rcomp[2*i].c_str()) ;
        rstr.append(" /") ;
        rstr.append(rcomp[2*i+1].c_str()) ;
        rstr.append("/ ") ;
      }
    }
    //cout << rstr << endl ;
    return (rstr);
  }
  else
    return ( mline );

}    

size_t firstCap(string strin){

  size_t i=0;
  while (strin[i])
  {
    if (isupper(strin[i])) return (i);
    i++;
  }
  return (strin.size()+1);

}

