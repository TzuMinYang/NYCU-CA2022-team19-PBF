#pragma once

#pragma warning( disable : 4996)

#include <stdarg.h>
#include <vector>
#include <string>
#include <cstdio>

#define GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer)						\
{																		\
	va_list args;														\
	pBuffer = NULL;														\
	int length = 1024;													\
	int result = -1;													\
	while (result == -1) {												\
		delete[] pBuffer;												\
		pBuffer = new char[length + 1];									\
		memset(pBuffer, 0, length + 1);									\
		va_start(args, fmt);											\
		result = vsnprintf(pBuffer, length, fmt, args);	\
		va_end(args);													\
		if (result >= length)											\
			result = -1;												\
		length *= 2;													\
	}																	\
}																		\

#define RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer)						\
{																		\
	delete[] pBuffer;													\
	pBuffer = NULL;														\
}

class Logger{
public:
	static void print(const char *fmt, ...);
	static void logPrint(const char *fmt, ...);
	static void consolePrint(const char*fmt, ...);

	static std::vector<std::string> consoleOutput;
	static int maxConsoleLineCount;

protected:
	Logger() {}
	virtual ~Logger() {}
};
