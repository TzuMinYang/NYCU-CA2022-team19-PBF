#include <Utils/Utils.h>

//assume that the gravity is in the y-direction (this can easily be changed if need be), and this value gives its magnitude. 
double Globals::g = -9.8;
//this is the direction of the up-vector
V3D Globals::worldUp = V3D(0, 1, 0);
//and the ground plane
Plane Globals::groundPlane = Plane(P3D(0,0,0), Globals::worldUp);
//this is the total ellapsed sim time
double Globals::currentSimulationTime = 0;

// given a list of keywords that map strings to integer values denoting keyword types, this method will determine the type of command that is passed in
int getLineType(char* &line, KeyWord* keywords, int nKeywords){
	for (int i = 0;i < nKeywords;i++) {
		if (strncmp(line, keywords[i].keyWord, strlen(keywords[i].keyWord)) == 0 && isWhiteSpace(line[strlen(keywords[i].keyWord)])) {
			line += strlen(keywords[i].keyWord);
			return keywords[i].retVal;
		}
	}
	return -1;
}

// given a list of keywords that map strings to integer values denoting keyword types, this method will determine the string corresponding to the token passed in
char* getKeyword(int lineType, KeyWord* keywords, int nKeywords) {
	for (int i = 0;i<nKeywords;i++) {
		if (lineType == keywords[i].retVal)
			return keywords[i].keyWord;
	}

	return NULL;
}

