#include <pybind11/pybind11.h>
#include "interface/kaHIP_interface.h"

pybind11::object wrap_kaffpa(
                const pybind11::object &vwgt,
                const pybind11::object &xadj,
                const pybind11::object &adjwgt,
                const pybind11::object &adjncy,
                int nparts,
                double imbalance,
                bool supress_output,
                int seed,
                int mode) {
        int n = pybind11::len(xadj) - 1;
        std::vector<int> xadjv, adjncyv, vwgtv, adjwgtv;
        for (auto it : xadj) 
                xadjv.push_back(pybind11::cast<int>(*it)); 

        for (auto it : adjncy) 
                adjncyv.push_back(pybind11::cast<int>(*it)); 

        for (auto it : vwgt) 
                vwgtv.push_back(pybind11::cast<int>(*it)); 

        for (auto it : adjwgt) 
                adjwgtv.push_back(pybind11::cast<int>(*it)); 

        int* part        = new int[n];
        int edge_cut     = 0;

        kaffpa(&n, &vwgtv[0], &xadjv[0], 
               &adjwgtv[0], &adjncyv[0], &nparts, 
               &imbalance, supress_output, 
               seed, mode, & edge_cut, part);

        pybind11::list returnblocksArray; 
        for (int i = 0; i<n; ++i) 
                returnblocksArray.append(part[i]); 

        delete[] part;

        return pybind11::make_tuple(edge_cut, returnblocksArray);
}

PYBIND11_MODULE(kahip, m) {
        m.doc() = "pybind11 example plugin"; // optional module docstring
        m.def("kaffpa", &wrap_kaffpa, "A function that partitions a graph.");
}
