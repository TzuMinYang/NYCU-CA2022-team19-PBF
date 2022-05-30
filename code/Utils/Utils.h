#pragma once

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Logger.h"

#include <MathLib/V3D.h>
#include <MathLib/Plane.h>

#include <vector>
#include <string>

#pragma warning( disable : 4996)

#define DynamicArray std::vector
typedef unsigned int uint;


class Globals {
public:
	// gravitational acceleration 
	static double g;
	//this is the direction of the up-vector
	static V3D worldUp;
	//and the ground plane
	static Plane groundPlane;
	//this is the total ellapsed sim time
	static double currentSimulationTime;
};

typedef struct key_word {
	char keyWord[50];
	int retVal;
}KeyWord;

// given a list of keywords that map strings to integer values denoting types of keywords, this method will determine the type of command that is passed in
// 'line' will be modified after this function call to point to the first character after the keyword. If no keywords are matched, the method returns -1
int getLineType(char* &line, KeyWord* keywords, int nKeywords);

// given a list of keywords that map strings to integer values denoting types of keywords, this method will determine the string corresponding to the token passed in
char* getKeyword(int lineType, KeyWord* keywords, int nKeywords);

/**
	This method throws an error with a specified text and arguments 
*/
inline void throwError(const char *fmt, ...){		// Custom error creation method
	char *pBuffer = NULL;
	GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

	Logger::print("Error Thrown: %s\n", pBuffer);
	
	throw pBuffer;

	RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

//interpret the array of characters as a sequence of strings...
inline void getCharSeparatedStringList(const char* cString, std::vector<std::string>& lines, char separator = '\n') {
	lines.clear();
	std::string line;
	for (const char *c = cString;*c;c++) {
		if (*c && *c != separator)
			line.append(1, *c);
		else {
			lines.push_back(line);
			line.clear();
		}
	}
	if (line.size() > 0)
		lines.push_back(line);
}

/**
	This method reads all the doubles from the given file and stores them in the array of doubles that is passed in
*/
inline void readDoublesFromFile(FILE* f, DynamicArray<double> *d){
	double temp;
	char line[200];
//	while (fscanf(f, "%lf\n", &temp) == 1)
//		d->push_back(temp);

	while (!feof(f)){
		if (fscanf(f, "%lf\n", &temp) == 1)
			d->push_back(temp);
		else
			fgets(line, sizeof(line)/sizeof(line[0])-1, f);
	}
}

/**
	This method returns a pointer to the first non-white space character location in the provided buffer
*/
inline char* lTrim(char* buffer){
	while (*buffer==' ' || *buffer=='\t' || *buffer=='\n' || *buffer=='\r')
		buffer++;
	return buffer;
}

inline char* rTrim(char* buffer){
	int index = (int)strlen(buffer) - 1;
	while (index>=0){
		if (buffer[index]==' ' || buffer[index]=='\t' || buffer[index]=='\n' || buffer[index]=='\r'){
			buffer[index] = '\0';
			index--;
		}
		else
			break;
	}
	return buffer;
}

inline char* trim(char* buffer){
	return rTrim(lTrim(buffer));
}

/**
	This method reads a line from a file. It does not return empty lines or ones that start with a pound key - those are assumed to be comments.
	This method returns true if a line is read, false otherwise (for instance the end of file is met).
*/
inline bool readValidLine(char* line, int nChars, FILE* fp){
	line[0] = '\0';
	while (!feof(fp)){
		fgets(line, nChars, fp);
		if ((int)strlen(line)>=nChars)
			throwError("The input file contains a line that is longer than the buffer - not allowed");
		char* tmp = trim(line);
		if (tmp[0]!='#' && tmp[0]!='\0')
			return true;
	}

	return false;
}


/**
	This method reads a line from a file. It does not return empty lines or ones that start with a pound key - those are assumed to be comments.
	This method returns true if a line is read, false otherwise (for instance the end of file is met).
*/
inline bool readValidLine(char* line, FILE* fp, int n){
	line[0] = '\0';
	while (!feof(fp)){
		fgets(line, n, fp);
		if ((int)strlen(line) > n-5)
			Logger::print("rearValidLine: Warning. The read file exceeds the buffer limits. The line was probably not read in its entirety.");
		char* tmp = trim(line);
		if (tmp[0]!='#' && tmp[0]!='\0')
			return true;
	}

	return false;
}

inline bool isWhiteSpace(char ch){
	return (ch==' ' || ch=='\t' || ch=='\n' || ch=='\r' || ch=='\0');
}


/**
	This method returns a DynamicArray of char pointers that correspond to the addressed
	of the tokens that are separated by white space in the string that is passed in as a pointer.
*/
inline DynamicArray<char*> getTokens(char* input){
	DynamicArray<char*> result;
	input = lTrim(input);
	//read in the strings one by one - assume that each tokens are less than 100 chars in length
	while (input[0]!='\0'){
		result.push_back(input);
		char tempStr[100];
		sscanf(input, "%s", tempStr);
		input = lTrim(input + strlen(tempStr));
	}
	return result;
}


