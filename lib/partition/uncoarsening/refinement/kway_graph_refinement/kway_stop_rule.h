/******************************************************************************
 * kway_stop_rule.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_STOP_RULE_ULPK0ZTF
#define KWAY_STOP_RULE_ULPK0ZTF

class kway_stop_rule {
public:
        kway_stop_rule(PartitionConfig & config) {};
        kway_stop_rule() {};
        virtual ~kway_stop_rule() {};

        virtual void push_statistics(Gain gain) = 0;
        virtual void reset_statistics() = 0;
        virtual bool search_should_stop(unsigned int min_cut_idx, 
                                        unsigned int cur_idx, 
                                        unsigned int search_limit) = 0;
};

class kway_simple_stop_rule : public kway_stop_rule {
public:
        kway_simple_stop_rule(PartitionConfig & config) {};
        virtual ~kway_simple_stop_rule() {};

        void push_statistics(Gain gain) {};
        void reset_statistics() {};
        bool search_should_stop(unsigned int min_cut_idx, 
                                unsigned int cur_idx, 
                                unsigned int search_limit);
};

inline bool kway_simple_stop_rule::search_should_stop(unsigned int min_cut_idx, 
                                                      unsigned int cur_idx, 
                                                      unsigned int search_limit) {
        return cur_idx - min_cut_idx > search_limit;
}


class kway_adaptive_stop_rule : public kway_stop_rule {
public:
        kway_adaptive_stop_rule(PartitionConfig & config) : m_steps(0), 
                                                            m_expected_gain(0.0), 
                                                            m_expected_variance2(0.0), 
                                                            pconfig(&config) {}
        virtual ~kway_adaptive_stop_rule() {};

        void push_statistics(Gain gain) {
                //erwartungstreue schätzer für varianz und erwartungswert
                m_expected_gain *= m_steps;
                m_expected_gain += gain;
                if(m_steps == 0) {
                        m_expected_variance2 = 0.0;
                } else {
                        m_expected_variance2 *= (m_steps-1);
                        //expected_variance2 += (gain - expected_gain)*(gain - expected_gain);
                        //real implementation in kaspar
                        m_expected_variance2 += (gain )*(gain );
                }
                m_steps++;
                
                m_expected_gain /= m_steps;
                if(m_steps > 1) 
                        m_expected_variance2 /= (m_steps-1);
        };

        void reset_statistics() {
                m_steps              = 0;
                m_expected_gain      = 0.0;
                m_expected_variance2 = 0.0;
        };

        bool search_should_stop(unsigned int min_cut_idx, 
                                unsigned int cur_idx, 
                                unsigned int search_limit);
private:
        unsigned m_steps;
        double   m_expected_gain;
        double   m_expected_variance2;
        PartitionConfig * pconfig;
};

inline bool kway_adaptive_stop_rule::search_should_stop(unsigned int min_cut_idx, 
                                                        unsigned int cur_idx, 
                                                        unsigned int search_limit) {

        return m_steps*m_expected_gain*m_expected_gain > 
                pconfig->kway_adaptive_limits_alpha * m_expected_variance2 + pconfig->kway_adaptive_limits_beta && (m_steps != 1);
}



#endif /* end of include guard: KWAY_STOP_RULE_ULPK0ZTF */
