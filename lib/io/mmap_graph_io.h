#pragma once

#include <cctype>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace kahip {
namespace mmap_io {
namespace {
struct MappedFile {
  const int fd;
  std::size_t position;
  const std::size_t length;
  char *contents;

  inline bool valid_position() const { return position < length; }
  inline char current() const { return contents[position]; }
  inline void advance() { ++position; }
};

int open_file(const std::string &filename) {
  int fd = open(filename.c_str(), O_RDONLY);
  if (fd < 0) {
    std::cerr << "Error while opening " << filename;
    std::exit(-1);
  }
  return fd;
}

std::size_t file_size(const int fd) {
  struct stat file_info {};
  if (fstat(fd, &file_info) == -1) {
    close(fd);
    std::cerr << "Error while determining file size";
    std::exit(-1);
  }
  return static_cast<std::size_t>(file_info.st_size);
}

MappedFile mmap_file_from_disk(const std::string &filename) {
  const int fd = open_file(filename);
  const std::size_t length = file_size(fd);

  char *contents = static_cast<char *>(mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0));
  if (contents == MAP_FAILED) {
    close(fd);
    std::cerr << "Error while mapping file to memory";
    std::exit(-1);
  }

  return {
      .fd = fd,
      .position = 0,
      .length = length,
      .contents = contents,
  };
}

void munmap_file_from_disk(const MappedFile &mapped_file) {
  if (munmap(mapped_file.contents, mapped_file.length) == -1) {
    close(mapped_file.fd);
    std::cerr << "Error while unmapping file from memory";
    std::exit(-1);
  }
  close(mapped_file.fd);
}

inline void skip_spaces(MappedFile &mapped_file) {
  while (mapped_file.valid_position() && mapped_file.current() == ' ') {
    mapped_file.advance();
  }
}

inline void skip_comment(MappedFile &mapped_file) {
  while (mapped_file.valid_position() && mapped_file.current() != '\n') {
    mapped_file.advance();
  }
  if (mapped_file.valid_position()) {
    mapped_file.advance();
  }
}

inline void skip_nl(MappedFile &mapped_file) { mapped_file.advance(); }

inline std::uint64_t scan_uint(MappedFile &mapped_file) {
  std::uint64_t number = 0;
  while (mapped_file.valid_position() && std::isdigit(mapped_file.current())) {
    const int digit = mapped_file.current() - '0';
    number = number * 10 + digit;
    mapped_file.advance();
  }
  skip_spaces(mapped_file);
  return number;
}
} // namespace

struct GraphHeader {
  uint64_t number_of_nodes;
  uint64_t number_of_edges;
  bool has_node_weights;
  bool has_edge_weights;
};

GraphHeader read_graph_header(MappedFile &mapped_file) {
  skip_spaces(mapped_file);
  while (mapped_file.current() == '%') {
    skip_comment(mapped_file);
    skip_spaces(mapped_file);
  }

  const std::uint64_t number_of_nodes = scan_uint(mapped_file);
  const std::uint64_t number_of_edges = scan_uint(mapped_file);
  const std::uint64_t format =
      (mapped_file.current() != '\n') ? scan_uint(mapped_file) : 0;
  skip_nl(mapped_file);

  const bool has_node_weights = (format % 100) / 10; // == x1x
  const bool has_edge_weights = format % 10;         // == xx1

  return {
      .number_of_nodes = number_of_nodes,
      .number_of_edges = number_of_edges,
      .has_node_weights = has_node_weights,
      .has_edge_weights = has_edge_weights,
  };
}

void graph_from_metis_file(graph_access &G, const std::string &filename) {
  MappedFile mapped_file = mmap_file_from_disk(filename);
  const GraphHeader header = read_graph_header(mapped_file);
  G.start_construction(header.number_of_nodes, 2 * header.number_of_edges);

  for (NodeID u = 0; u < header.number_of_nodes; ++u) {
    G.new_node();
    G.setPartitionIndex(u, 0);

    skip_spaces(mapped_file);
    while (mapped_file.current() == '%') {
      skip_comment(mapped_file);
      skip_spaces(mapped_file);
    }

    G.setNodeWeight(u, header.has_node_weights ? scan_uint(mapped_file) : 1);
    while (std::isdigit(mapped_file.current())) {
      const EdgeWeight weight = (header.has_edge_weights) ? scan_uint(mapped_file) : 1;
      const EdgeID e = G.new_edge(u, scan_uint(mapped_file) - 1);
      G.setEdgeWeight(e, weight);
    }
    if (mapped_file.current() == '\n') {
      skip_nl(mapped_file);
    }
  }
  G.finish_construction();
  munmap_file_from_disk(mapped_file);
}
} // namespace mmap_io
} // namespace kahip
