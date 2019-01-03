/******************************************************************************
 * path.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PATH_X5LQS3DT
#define PATH_X5LQS3DT

#include "definitions.h"

class path {
        public:
                path( );
                path( const NodeID & v );
                virtual ~path();

                void init(const NodeID & v);

                NodeID get_tail() const;
                void set_tail(const NodeID & id);

                NodeID get_head() const;
                void set_head(const NodeID & id);

                void set_length(const EdgeID & length);
                EdgeID get_length() const;

                //returns wether the given node is an endpoint of the path
                bool is_endpoint(const NodeID & id) const;
                
                //returns wether the path is a cycle or not.  
                bool is_cycle() const;

                bool is_active() const;
                void set_active(const bool active);

        private:
                //Last vertex of the path. Cycles have head == tail
                NodeID head;

                //First vertex of the path. Cycles have head == tail
                NodeID tail;

                //Number of edges in the graph
                EdgeID length; 
               
                // True iff the parth is still in use. False iff it has been removed.  
                bool active;

};

inline void path::init(const NodeID & v) {
        head   = v;
        tail   = v;
        length = 0;
        active = true;
}

inline NodeID path::get_tail() const {
        return tail;        
}

inline void path::set_tail(const NodeID & v) {
        tail = v;        
}

inline NodeID path::get_head( ) const {
        return head;
}

inline void path::set_head(const NodeID & v) {
        head = v;        
}

inline EdgeID path::get_length( ) const {
        return length;        
}

inline void path::set_length(const EdgeID & len) {
        length = len;        
}

inline bool path::is_endpoint(const NodeID & v) const {
        return (v == tail) or (v == head);
}

inline bool path::is_cycle() const {
        return (head == tail) and (length > 0);       
}

inline bool path::is_active() const {
        return active;
}

inline void path::set_active(const bool act) {
        active = act;        
}


#endif /* end of include guard: PATH_X5LQS3DT */
