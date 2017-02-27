/******************************************************************************
 * error_codes.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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

#ifndef ERROR_CODES_H
#define ERROR_CODES_H

namespace distributed_graph {

// Generic error codes.
//

// Return code for "no error".
const int kOk = 0;
// Return code for generic errors.
const int kError = 1;
// Return value for invalid arguments.
const int kInvalidArguments = 2;


// I/O-Related errors codes.
//


// Code returned when the file could not be opened.
const int kCouldNotOpenFile = 3;
// Code returned when the file format was invalid.
const int kInvalidFormat = 4;
// Code returned when there was an error writing to the file.
const int kErrorWritingFile = 5;
// Code returned when there was an error reading the file.
const int kErrorReadingFile = 6;

}  // namespace distributed_graph

#endif // ifndef ERROR_CODES_H

