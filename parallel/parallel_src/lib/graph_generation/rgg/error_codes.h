/* File:   error_codes.h
 * Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
 *
 * This header defines the error codes used throughout the program.
 */
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

