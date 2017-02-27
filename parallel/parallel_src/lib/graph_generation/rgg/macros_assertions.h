// File:   assert.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Assertions that do not kill the program but only print a warning.
//
// Macros:
//
//   ASSERT_TRUE(x)
//   ASSERT_NEQ(left, right)
//   ASSERT_EQ(left, right)
//   ASSERT_LT(left, right)
//   ASSERT_GT(left, right)
//   ASSERT_LEQ(left, right)
//   ASSERT_GEQ(left, right)
//   ASSERT_BETWEEN(left, x, right)
//   ASSERT_RANGE_GEQ(sequence, begin, end, x, i)

#ifndef ASSERT_H
#define ASSERT_H

#include <cassert>
#include <iostream> // cerr
#include <cstdio>   // fprintf(), stderr
#include <cstdlib>  // abort()

#include "macros_common.h"

// A custom assertion macro that does not kill the program but prints to
// stderr instead.  This is named ASSERT_TRUE as not to collide with METIS.
#ifdef NDEBUG
# define ASSERT_TRUE(x) do {} while (false);
#else
# define ASSERT_TRUE(expression) \
  do { \
    if (not (expression)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(expression) << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left != right.
#ifdef NDEBUG
# define ASSERT_NEQ(left, right) do {} while (false);
#else
# define ASSERT_NEQ(left, right) \
  do { \
    if ((left) == (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " != " << STR(right) << " but was " \
        << left << " == " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left == right.
#ifdef NDEBUG
# define ASSERT_EQ(left, right) do {} while (false);
#else
# define ASSERT_EQ(left, right) \
  do { \
    if ((left) != (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " == " << STR(right) << " but was " \
        << left << " != " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left < right.
#ifdef NDEBUG
# define ASSERT_LT(left, right) do {} while (false);
#else
# define ASSERT_LT(left, right) \
  do { \
    if ((left) >= (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " < " << STR(right) << " but was " \
        << left << " >= " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left > right.
#ifdef NDEBUG
# define ASSERT_GT(left, right) do {} while (false);
#else
# define ASSERT_GT(left, right) \
  do { \
    if ((left) <= (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " > " << STR(right) << " but was " \
        << left << " <= " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left <= right.
#ifdef NDEBUG
# define ASSERT_LEQ(left, right) do {} while (false);
#else
# define ASSERT_LEQ(left, right) \
  do { \
    if ((left) > (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " <= " << STR(right) << " but was " \
        << left << " > " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: left >= right.
#ifdef NDEBUG
# define ASSERT_GEQ(left, right) do {} while (false);
#else
# define ASSERT_GEQ(left, right) \
  do { \
    if ((left) < (right)) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(left) << " >= " << STR(right) << " but was " \
        << left << " < " << right << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: x <= y <= z.
#ifdef NDEBUG
# define ASSERT_BETWEEN(x, y, z) do {} while (false);
#else
# define ASSERT_BETWEEN(left, x, right) \
  do { \
    if (((left) > (x)) or ((right) < (x))) { \
      std::cerr << "ASSERTION FAILED [" << __FILE__ << ":" << __LINE__ << \
        "]. Asserted: " << STR(x) << " in {" << STR(left) << ", ..., " << \
        STR(right) << "} but was " << x << " not in {" << left << \
        ", ..., " << right << "}." << std::endl; \
      abort(); \
    } \
  } while (false)
#endif

// Assert: \forall begin <= i < end: sequence[i] > x.
#ifdef NDEBUG
# define ASSERT_RANGE_GT(sequence, begin, end, x, i) do {} while (false);
#else
# define ASSERT_RANGE_GT(sequence, begin, end, x, i) \
  for (int i = begin; i < end; ++i) { \
    ASSERT_GT( sequence[i], x ); \
  }
#endif

// Assert: \forall begin <= i < end: sequence[i] >= x.
#ifdef NDEBUG
# define ASSERT_RANGE_GEQ(sequence, begin, end, x, i) do {} while (false);
#else
# define ASSERT_RANGE_GEQ(sequence, begin, end, x, i) \
  for (int i = begin; i < end; ++i) { \
    ASSERT_GEQ( sequence[i], x ); \
  }
#endif

#endif // ifndef ASSERT_H

