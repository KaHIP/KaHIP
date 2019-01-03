/******************************************************************************
 * matrix.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MATRIX_BHHZ9T7P
#define MATRIX_BHHZ9T7P

class matrix {
public:
        matrix(unsigned int dim_x, unsigned int dim_y) {};
        matrix() {};
        virtual ~matrix() {};

        virtual int  get_xy(unsigned int x, unsigned int y)            = 0;
        virtual void set_xy(unsigned int x, unsigned int y, int value) = 0;

        virtual unsigned int get_x_dim() = 0;
        virtual unsigned int get_y_dim() = 0;
};


#endif /* end of include guard: MATRIX_BHHZ9T7P */
