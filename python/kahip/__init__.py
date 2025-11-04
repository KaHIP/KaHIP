"""
KaHIP - Karlsruhe High Quality Graph Partitioning Framework

Python bindings for the KaHIP graph partitioning library.
"""

try:
    from ._version import __version__
except ImportError:
    # Fallback version for development installs without setuptools-scm
    __version__ = "0.0.0+unknown"

try:
    from .kahip import kaffpa
except ImportError as e:
    raise ImportError(
        "Failed to import the KaHIP C++ extension. "
        "Make sure the package is properly installed."
    ) from e

class kahip_graph:
    """Graph class for constructing graphs and converting to CSR format for KaHIP."""
    
    def __init__(self):
        self.num_nodes = 0
        self.edges = {}
        self.node_weights = {}
    
    def set_num_nodes(self, n):
        self.num_nodes = n
        for i in range(n):
            if i not in self.node_weights:
                self.node_weights[i] = 1
    
    def add_undirected_edge(self, source, target, weight=1):
        if source >= self.num_nodes or target >= self.num_nodes:
            raise ValueError(f"Node index out of range. Graph has {self.num_nodes} nodes.")
        self.edges[(source, target)] = weight
        if source != target:
            self.edges[(target, source)] = weight
    
    def set_weight(self, node_id, weight):
        if node_id >= self.num_nodes:
            raise ValueError(f"Node index {node_id} out of range. Graph has {self.num_nodes} nodes.")
        self.node_weights[node_id] = weight
    
    def get_csr_arrays(self):
        if self.num_nodes == 0:
            return [], [0], [], []
        
        degrees = [0] * self.num_nodes
        for u, v in self.edges:
            degrees[u] += 1
        
        total_edges = sum(degrees)
        adjncy = [0] * total_edges
        adjcwgt = [0] * total_edges
        xadj = [0] * (self.num_nodes + 1)
        
        for i in range(self.num_nodes):
            xadj[i + 1] = xadj[i] + degrees[i]
        
        pos = xadj[:]
        for (u, v), weight in self.edges.items():
            adjncy[pos[u]] = v
            adjcwgt[pos[u]] = weight
            pos[u] += 1
        
        for i in range(self.num_nodes):
            start, end = xadj[i], xadj[i + 1]
            if end > start:
                pairs = list(zip(adjncy[start:end], adjcwgt[start:end]))
                pairs.sort()
                for j, (v, w) in enumerate(pairs):
                    adjncy[start + j] = v
                    adjcwgt[start + j] = w
        
        vwgt = [self.node_weights.get(i, 1) for i in range(self.num_nodes)]
        return vwgt, xadj, adjcwgt, adjncy

__all__ = ["kaffpa", "kahip_graph", "__version__"]

# Partitioning mode constants
FAST = 0
ECO = 1
STRONG = 2
FASTSOCIAL = 3
ECOSOCIAL = 4
STRONGSOCIAL = 5
