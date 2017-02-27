// this files contains all the application independent little
// functions and macros used for the optimizer.
// In particular Peters debug macros and Dags stuff
// from dbasic.h cdefs, random,...

//////////////// stuff originally from debug.h ///////////////////////////////
// (c) 1997 Peter Sanders
// some little utilities for debugging adapted
// to the paros conventions


#ifndef UTIL
#define UTIL

// default debug level. will be overidden e.g. if debug.h is included
#ifndef DEBUGLEVEL
#define DEBUGLEVEL 1
#endif

#if DEBUGLEVEL >= 0
#define Debug0(A) A
#else
#define Debug0(A) 
#endif
#if DEBUGLEVEL >= 1
#define Debug1(A) A
#else
#define Debug1(A) 
#endif
#if DEBUGLEVEL >= 2
#define Debug2(A) A
#else
#define Debug2(A) 
#endif
#if DEBUGLEVEL >= 3
#define Debug3(A) A
#else
#define Debug3(A) 
#endif
#if DEBUGLEVEL >= 4
#define Debug4(A) A
#else
#define Debug4(A) 
#endif
#if DEBUGLEVEL >= 5
#define Debug5(A) A
#else
#define Debug5(A) 
#endif
#if DEBUGLEVEL >= 6
#define Debug6(A) A
#else
#define Debug6(A) 
#endif

#define Assert(c) if(!(c))\
  {cout << "\nAssertion violation " << __FILE__ << ":" << __LINE__ << endl;}
#define Assert0(C) Debug0(Assert(C))
#define Assert1(C) Debug1(Assert(C))
#define Assert2(C) Debug2(Assert(C))
#define Assert3(C) Debug3(Assert(C))
#define Assert4(C) Debug4(Assert(C))
#define Assert5(C) Debug5(Assert(C))

#define Error(s) {cout << "\nError:" << s << " " << __FILE__ << ":" << __LINE__ << endl;}

////////////// min, max etc. //////////////////////////////////////

#ifndef Max
#define Max(x,y) ((x)>=(y)?(x):(y))
#endif

#ifndef Min
#define Min(x,y) ((x)<=(y)?(x):(y))
#endif

#ifndef Abs
#define Abs(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef PI
#define PI 3.1415927
#endif

// is this the right definition of limit?
inline double limit(double x, double bound)
{
  if      (x >  bound) { return  bound; }
  else if (x < -bound) { return -bound; }
  else                   return x;
}

/////////////////////// timing /////////////////////
#include <time.h>


// elapsed CPU time see also /usr/include/sys/time.h
inline double cpuTime()
{ //struct timespec tp;

  return clock() * 1e-6;
//  clock_gettime(CLOCK_VIRTUAL, &tp);
//  return tp.tv_sec + tp.tv_nsec * 1e-9;
}

#endif

