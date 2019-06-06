#ifndef KAHIP_VIECLUS_ADAPTER_H
#define KAHIP_VIECLUS_ADAPTER_H

#include <vector>

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"

class VieClus_adapter {
public:
    VieClus_adapter();
    ~VieClus_adapter();

    VieClus_adapter(const VieClus_adapter&) = delete;
    VieClus_adapter(VieClus_adapter&&) = delete;
    VieClus_adapter& operator=(const VieClus_adapter&) = delete;
    VieClus_adapter& operator=(VieClus_adapter&&) = delete;

    std::vector<int> compute_modularity_clustering(graph_access &G, const PartitionConfig &partition_config);

    static VieClus_adapter &instance();
};

#endif // KAHIP_VIECLUS_ADAPTER_H
