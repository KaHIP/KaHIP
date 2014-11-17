/******************************************************************************
 * macros_common.h               
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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

