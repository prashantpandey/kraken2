/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "taxonomy.h"
#include "kv_store.h"
#include "math.h"

using std::string;
using std::map;
using std::set;
using std::ifstream;
using std::istringstream;
using std::queue;
using std::ofstream;

namespace kraken2 {

NCBITaxonomy::NCBITaxonomy(string nodes_filename, string names_filename) {
  ifstream nodes_file(nodes_filename.c_str());
  if (! nodes_file.good())
    err(EX_NOINPUT, "error opening %s", nodes_filename.c_str());

  ifstream names_file(names_filename.c_str());
  if (! names_file.good())
    err(EX_NOINPUT, "error opening %s", names_filename.c_str());

  string line;
  uint64_t node_id, parent_id;
  string name, rank, junk;

  const string delim = "\t|\t";
  while (getline(nodes_file, line)) {
    line.pop_back();  // discard trailing
    line.pop_back();  //   "\t|"
    size_t pos1, pos2;
    pos1 = 0;
    int field_ct = 0;
    bool finished = false;
    while (field_ct++ < 10 && ! finished) {  // tokenizing loop
      pos2 = line.find(delim, pos1);
      string token;
      if (pos2 == string::npos) {
        token = line.substr(pos1);
        finished = true;
      }
      else {
        token = line.substr(pos1, pos2 - pos1);
        pos1 = pos2 + delim.size();
      }

      switch (field_ct) {  // 1-based counting
        case 1 :  // node ID
          node_id = (uint64_t) stoul(token);
          if (node_id == 0)
            errx(EX_DATAERR, "attempt to create taxonomy w/ node ID == 0");
          break;
        case 2 : // parent ID
          parent_id = (uint64_t) stoul(token);
          break;
        case 3 : // rank
          rank = token;
          finished = true;
          break;
      }
    }  // end tokenizing loop
    if (node_id == 1)
      parent_id = 0;
    parent_map_[node_id] = parent_id;
    child_map_[parent_id].insert(node_id);
    rank_map_[node_id] = rank;
    known_ranks_.insert(rank);
  }

  while (getline(names_file, line)) {
    line.pop_back();  // discard trailing
    line.pop_back();  //   "\t|"
    size_t pos1, pos2;
    pos1 = 0;
    int field_ct = 0;
    bool finished = false;
    while (field_ct++ < 10 && ! finished) {  // tokenizing loop
      pos2 = line.find(delim, pos1);
      string token;
      if (pos2 == string::npos) {
        token = line.substr(pos1);
        finished = true;
      }
      else {
        token = line.substr(pos1, pos2 - pos1);
        pos1 = pos2 + delim.size();
      }

      switch (field_ct) {  // 1-based counting
        case 1 :  // node ID
          node_id = (uint64_t) stoul(token);
          if (node_id == 0)
            errx(EX_DATAERR, "attempt to create taxonomy w/ node ID == 0");
          break;
        case 2 : // name
          name = token;
          break;
        case 4 : // type of name - determines whether or not to record in map
          if (token == "scientific name")
            name_map_[node_id] = name;
          finished = true;
          break;
      }
    }  // end tokenizing loop
  }

  marked_nodes_.insert(1);  // mark root node
}

// Mark the given taxonomy node and all its unmarked ancestors
void NCBITaxonomy::MarkNode(uint64_t taxid) {
  while (marked_nodes_.count(taxid) == 0) {
    marked_nodes_.insert(taxid);
    taxid = parent_map_[taxid];
  }
}

#define SEED 2038074761

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;
	uint64_t h = seed ^ (len * m);
	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);
	while(data != end)
	{
		uint64_t k = *data++;
		k *= m; 
		k ^= k >> r; 
		k *= m; 
		h ^= k;
		h *= m; 
	}
	const unsigned char * data2 = (const unsigned char*)data;
	switch(len & 7)
	{
		case 7: h ^= (uint64_t)data2[6] << 48;
		case 6: h ^= (uint64_t)data2[5] << 40;
		case 5: h ^= (uint64_t)data2[4] << 32;
		case 4: h ^= (uint64_t)data2[3] << 24;
		case 3: h ^= (uint64_t)data2[2] << 16;
		case 2: h ^= (uint64_t)data2[1] << 8;
		case 1: h ^= (uint64_t)data2[0];
						h *= m;
	};
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

uint32_t maxDepth(uint64_t node_id, const Taxonomy& taxo) {
	  const TaxonomyNode node = taxo.nodes()[node_id];
    if (node.child_count == 0)
        return 1;
    else
    {
			uint32_t max = 0;
			uint32_t start = node.first_child;
			uint32_t end = node.first_child + node.child_count;
			  for (uint32_t i = start; i < end; i++) {
					uint32_t depth = maxDepth(i, taxo);
					max = max > depth ? max : depth;
				}

        return(max + 1);
    }
}

void Taxonomy::generatePrefixEncoding(const std::string &filename) const {
	uint32_t num_hash_bits = log2 (node_count_);
	uint32_t num_hash_bytes = (num_hash_bits+8)/8;
	uint32_t max_depth = maxDepth(1, *this);
	
	uint32_t encoding_len_bytes = num_hash_bytes * max_depth;
	std::unordered_map<uint64_t, std::string> node_encoding_map;

  ofstream taxo_file(filename + ".encoding");
  queue<std::pair<uint64_t, uint32_t>> bfs_queue;
  bfs_queue.push({1, 0});
  while (! bfs_queue.empty()) {
		auto node = bfs_queue.front();
		bfs_queue.pop();

		char prefix_encoding[encoding_len_bytes];
		memset(prefix_encoding, 0, encoding_len_bytes);
		uint64_t hash = MurmurHash3(node.first);

		memcpy(prefix_encoding + node.second, (char*)&hash, num_hash_bits);

		uint64_t parent = nodes_[node.first].parent_id;
		uint32_t height = node.second - 1;
		while (parent != 0) {
			auto encoding = node_encoding_map[parent].c_str();
			memcpy(prefix_encoding + height, &encoding, num_hash_bits);
			parent = nodes_[parent].parent_id;
			height--;
		}

		node_encoding_map.insert({node.first, prefix_encoding});
		taxo_file << node.first << '\t' << prefix_encoding << '\n';

		uint32_t start = nodes_[node.first].first_child;
		uint32_t end = nodes_[node.first].first_child +
			nodes_[node.first].child_count;
		for (uint32_t i = start; i < end; i++) {
			bfs_queue.push({i, node.second + 1});
		}
	}

  taxo_file.close();
}

void NCBITaxonomy::ConvertToKrakenTaxonomy(const char *filename) {
  Taxonomy taxo;
  TaxonomyNode zeroes_node = { 0, 0, 0, 0, 0, 0, 0 };

  taxo.node_count_ = marked_nodes_.size() + 1;  // +1 because 0 is illegal value
  taxo.nodes_ = new TaxonomyNode[taxo.node_count_];
  taxo.nodes_[0] = zeroes_node;

 	string name_data;
  string rank_data;

  // Because so many of the node rank names are shared, we only store one copy
  // of each rank
  map<string, uint64_t> rank_offsets;
  for (string rank : known_ranks_) {
    rank_offsets[rank] = rank_data.size();
    rank_data.append(rank.c_str(), rank.size() + 1);
  }

  uint64_t internal_node_id = 0;
  map<uint64_t, uint64_t> external_id_map;  // keys: ext. ID, values: int. ID
  external_id_map[0] = 0;
  external_id_map[1] = 1;  // 1 is root in both NCBI and Kraken taxonomies

  // Breadth-first search through NCBI taxonomy, assigning internal IDs
  // in sequential order as nodes are encountered via BFS.
  queue<uint64_t> bfs_queue;
  bfs_queue.push(1);
  while (! bfs_queue.empty()) {
    ++internal_node_id;
    uint64_t external_node_id = bfs_queue.front();
    bfs_queue.pop();
    external_id_map[external_node_id] = internal_node_id;

    TaxonomyNode node = zeroes_node;  // just to initialize node to zeros
    node.parent_id = external_id_map[ parent_map_[external_node_id] ];
    node.external_id = external_node_id;
    node.rank_offset = rank_offsets[ rank_map_[external_node_id] ];
    node.name_offset = name_data.size();
    node.first_child = internal_node_id + 1 + bfs_queue.size();
    for (uint64_t child_node : child_map_[external_node_id]) {
      // Only add marked nodes to our internal tree
      if (marked_nodes_.count(child_node) > 0) {
        bfs_queue.push(child_node);
        node.child_count++;
      }
    }
    taxo.nodes_[internal_node_id] = node;

    string name = name_map_[external_node_id];
    name_data.append(name.c_str(), name.size() + 1);
  }  // end BFS while loop

  taxo.rank_data_ = new char[ rank_data.size() ];
  memcpy(taxo.rank_data_, rank_data.data(), rank_data.size());
  taxo.rank_data_len_ = rank_data.size();

  taxo.name_data_ = new char[ name_data.size() ];
  memcpy(taxo.name_data_, name_data.data(), name_data.size());
  taxo.name_data_len_ = name_data.size();

  taxo.WriteToDisk(filename);
}

Taxonomy::Taxonomy(const string &filename, bool memory_mapping) {
  Init(filename.c_str(), memory_mapping);
}

Taxonomy::Taxonomy(const char *filename, bool memory_mapping) {
  Init(filename, memory_mapping);
}

void Taxonomy::Init(const char *filename, bool memory_mapping) {
  if (memory_mapping) {
    taxonomy_data_file_.OpenFile(filename);
    file_backed_ = true;
    char *ptr = taxonomy_data_file_.fptr();
    if (strncmp(FILE_MAGIC, ptr, strlen(FILE_MAGIC)) != 0) {
      errx(EX_DATAERR,
          "attempt to load taxonomy from malformed file %s", filename);
    }
    ptr += strlen(FILE_MAGIC);
    memcpy((char *) &node_count_, ptr, sizeof(node_count_));
    ptr += sizeof(node_count_);
    memcpy((char *) &name_data_len_, ptr, sizeof(name_data_len_));
    ptr += sizeof(name_data_len_);
    memcpy((char *) &rank_data_len_, ptr, sizeof(rank_data_len_));
    ptr += sizeof(rank_data_len_);
    nodes_ = (TaxonomyNode *) ptr;
    ptr += sizeof(*nodes_) * node_count_;
    name_data_ = ptr;
    ptr += name_data_len_;
    rank_data_ = ptr;
  }
  else {
    std::ifstream ifs(filename);
    file_backed_ = false;
    char magic[strlen(FILE_MAGIC) + 1];
    memset(magic, 0, strlen(FILE_MAGIC) + 1);
    ifs.read(magic, strlen(FILE_MAGIC));
    if (strcmp(magic, FILE_MAGIC) != 0)
      errx(EX_DATAERR, "malformed taxonomy file %s", filename);
    ifs.read((char *) &node_count_, sizeof(node_count_));
    ifs.read((char *) &name_data_len_, sizeof(name_data_len_));
    ifs.read((char *) &rank_data_len_, sizeof(rank_data_len_));
    nodes_ = new TaxonomyNode[node_count_];
    ifs.read((char *) nodes_, sizeof(*nodes_) * node_count_);
    name_data_ = new char[name_data_len_];
    ifs.read((char *) name_data_, name_data_len_);
    rank_data_ = new char[rank_data_len_];
    ifs.read((char *) rank_data_, rank_data_len_);
    if (! ifs)
      errx(EX_DATAERR, "read exhausted taxonomy information in %s", filename);
  }
}

Taxonomy::~Taxonomy() {
  // If file backed, deleting would be... bad.
  if (! file_backed_) {
    delete[] nodes_;
    delete[] name_data_;
    delete[] rank_data_;
  }
}

// Logic here depends on higher nodes having smaller IDs
// Idea: advance B tracker up tree, A is ancestor iff B tracker hits A
bool Taxonomy::IsAAncestorOfB(uint64_t a, uint64_t b) {
  if (! a || ! b)
    return false;
  while (b > a) {
    b = nodes_[b].parent_id;
  }
  return b == a;
}

// Logic here depends on higher nodes having smaller IDs
// Idea: track two nodes, advance lower tracker up tree, trackers meet @ LCA
uint64_t Taxonomy::LowestCommonAncestor(uint64_t a, uint64_t b) {
  if (! a || ! b)  // LCA(x,0) = LCA(0,x) = x
    return a ? a : b;
  while (a != b) {
    if (a > b)
      a = nodes_[a].parent_id;
    else
      b = nodes_[b].parent_id;
  }
  return a;
}

// Dump binary data to file
void Taxonomy::WriteToDisk(const char *filename) {
  ofstream taxo_file(filename);
  taxo_file.write(FILE_MAGIC, strlen(FILE_MAGIC));
  taxo_file.write((char *) &node_count_, sizeof(node_count_));
  taxo_file.write((char *) &name_data_len_, sizeof(name_data_len_));
  taxo_file.write((char *) &rank_data_len_, sizeof(rank_data_len_));
  taxo_file.write((char *) nodes_, sizeof(*nodes_) * node_count_);
  taxo_file.write(name_data_, name_data_len_);
  taxo_file.write(rank_data_, rank_data_len_);
  if (! taxo_file.good())
    errx(EX_OSERR, "error writing taxonomy to %s", filename);
  taxo_file.close();
}

void Taxonomy::GenerateExternalToInternalIDMap() {
  external_to_internal_id_map_[0] = 0;
  for (size_t i = 1; i < node_count_; i++) {
    external_to_internal_id_map_[nodes_[i].external_id] = i;
  }
}

}
