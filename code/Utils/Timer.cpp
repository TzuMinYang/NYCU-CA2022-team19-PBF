#include <Utils/Timer.h>

#ifdef WIN32
#include <time.h>
#include <windows.h>
#else
#include <sys/time.h>
#endif


/**
	This constructor initializes a timer.
*/
Timer::Timer(){
#ifdef WIN32
	//first, get an idea of the frequency...
	DWORD_PTR oldmask = SetThreadAffinityMask(GetCurrentThread(), (DWORD_PTR)1);
	QueryPerformanceFrequency((LARGE_INTEGER *)&this->frequency);
	SetThreadAffinityMask(GetCurrentThread(), oldmask);
#else
#endif
	restart();

}

/**
	Default timer destructor - doesn't do much.
*/
Timer::~Timer(){
}

/**
	This method resets the starting time.
*/
void Timer::restart(){
#ifdef WIN32
	DWORD_PTR oldmask = SetThreadAffinityMask(GetCurrentThread(), (DWORD_PTR)1);
	QueryPerformanceCounter((LARGE_INTEGER *)&this->startTime);
	SetThreadAffinityMask(GetCurrentThread(), oldmask);
#else
	gettimeofday(&startTime,0);

#endif
}

/**
	This method returns the number of milliseconds that has ellapsed since the timer was restarted.
*/
double Timer::timeEllapsed(){
#ifdef WIN32
	long long int tempTime;
	//force the thread to run on CPU 0 because the QPC method is buggy
	DWORD_PTR oldmask = SetThreadAffinityMask(GetCurrentThread(), (DWORD_PTR)1);
	QueryPerformanceCounter((LARGE_INTEGER *)&tempTime);
	//let it run wild and free again
	SetThreadAffinityMask(GetCurrentThread(), oldmask);
	if (tempTime<startTime)
		return 0;
	return (tempTime - startTime) / (double)this->frequency;
#else
	struct timeval endTime;
	struct timeval interval_elapsed;
	gettimeofday(&endTime,0);
	timeval_subtract(&interval_elapsed, &endTime, &startTime);
	return (double)interval_elapsed.tv_sec+(double)interval_elapsed.tv_usec/1000000;
#endif
}


#ifdef WIN32
#else

int Timer::timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

void Timer::timeval_add(struct timeval *result, struct timeval *time1,struct timeval *time2)
{

	/* Add the two times together. */
    result->tv_sec = time1->tv_sec + time2->tv_sec ;
    result->tv_usec = time1->tv_usec + time2->tv_usec ;
    if (result->tv_usec >= 1000000L) {			/* Carry? */
        result->tv_sec++ ;  result->tv_usec = result->tv_usec - 1000000L ;
    }

}

void Timer::timeval_create(struct timeval *result, long  seconds, long  microseconds)
{

    seconds = (seconds < 0) ? 0 : seconds ;
    microseconds = (microseconds < 0) ? 0 : microseconds ;

/* "Normalize" the time so that the microseconds field is less than a million. */

    while (microseconds >= 1000000L) {
        seconds++ ;  microseconds = microseconds - 1000000L ;
    }

/* Return the time in a TIMEVAL structure. */

    result->tv_sec = seconds ;
    result->tv_usec = microseconds ;
}
#endif

