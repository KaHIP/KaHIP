// File:   macros_common.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Common macros used in many places.
//
// Macros:
//
//   DEBUG    Iff 1 then debugging is enabled.  Default: 0.
//   VLOG(x)  x if DEBUG, the empty string otherwise.
//   VVLOG(x) x if DEBUG, the empty string otherwise.
//   STR(x)   Converts the symbol x to a string literal.
//   DISALLOW_COPY_AND_ASSIGN(TypeName)
//                              Declare copy constructor and assignment operator
//   DISALLOW_ASSIGN(TypeName)  Declare assignment operator.

#ifndef MACROS_COMMON_H 
#define MACROS_COMMON_H

// If DEBUG has not been defined yet then define it as false.
#ifndef DEBUG
# define DEBUG 0
#endif

#if DEBUG
# define VLOG(x) x
# define VVLOG(x) x
#else
# define VLOG(x)
# define VVLOG(x)
#endif


// A macro to define empty copy constructors and assignment operators.
// Use this after private: to hide them to the outside.
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)


// A macro to disallow the assignment operators.
#define DISALLOW_ASSIGN(TypeName) \
  void operator=(const TypeName&)


// Helper macro for STR.
#define ASSERT_H_XSTR(x) (#x)

// This macro allows to convert an expression to a string.
#define STR(x) ASSERT_H_XSTR(x)


#endif // ifndef MACROS_COMMON_H

