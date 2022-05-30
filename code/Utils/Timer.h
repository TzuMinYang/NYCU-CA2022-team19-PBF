#pragma once
#ifdef __linux__
#include <sys/time.h>
#endif

/*==================================================================================================================================================================*
 |	This is a simple timer class, that can be reset at any point, and can be used to get the ellapsed time since the last reset. In windows, it will be implemented |
 |  using the methods QueryPerformanceFrequency() and QueryPerformanceCounter(). In Linux, it can be implemented using the method gettimeofday(). The timeEllapsed  |
 |  method will return the time, in seconds but expressed as a double, since the timer was reset.                                                                   |
 *==================================================================================================================================================================*/

class Timer{
private:
#ifdef WIN32
	//this is the start time of the timer, in milliseconds.
	long long int startTime;
	//this is the frequency of the performance counter
	long long int frequency;
#else
	//! start time of reset
	struct timeval startTime;
#endif


public:
	/**
		Default constructor - resets the timer for the first time.
	*/
	Timer();
	/**
		Default destructor - doesn't do much.
	*/
	~Timer();
	/**
		This method resets the starting time. All the time Ellapsed function calls will use the reset start time for their time evaluations.
	*/
	void restart();
	
	/**
		This method returns the number of seconds that has ellapsed since the timer was reset.
	*/
	double timeEllapsed();

#ifdef WIN32
#else
protected:
	/*! Subtract two timevals.
	 *
	 * @param result	the result
	 * @param x			first timeval
	 * @param y			second timeval
	 * @return			1 if result is negative
	 */
	int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);

	void timeval_add(struct timeval *result, struct timeval *time1,struct timeval *time2);

	void timeval_create(struct timeval *result, long  seconds, long  microseconds);
#endif
};

