#include "Logger.h"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <Utils/Utils.h>

std::vector<std::string> Logger::consoleOutput;
int Logger::maxConsoleLineCount = 10;

void Logger::consolePrint(const char *fmt, ...) {
	char *pBuffer = NULL;
	GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);


	static FILE *fp = fopen("..\\out\\console.txt", "wt");
	fprintf(fp, "%s", pBuffer);
	fflush(fp);

	std::vector<std::string> newLines;
	getCharSeparatedStringList(pBuffer, newLines);

	for (uint i = 0; i < newLines.size();i++)
		consoleOutput.push_back(newLines[i]);

	while ((int)consoleOutput.size() > maxConsoleLineCount) consoleOutput.erase(consoleOutput.begin());

	RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

void Logger::print(const char *fmt, ...){
	char *pBuffer = NULL;
	GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

	printf(pBuffer);

	RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

void Logger::logPrint(const char *fmt, ...) {
	char *pBuffer = NULL;
	GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

	static FILE *fp = fopen("..\\out\\log.txt", "wt");

	fprintf(fp, "%s", pBuffer);
	fflush(fp);


	RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}


