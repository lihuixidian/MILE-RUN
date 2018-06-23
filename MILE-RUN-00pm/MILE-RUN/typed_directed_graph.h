#pragma once

#include <set> // std::set
#include <vector> // std::vector
#include <utility> // std::pair
#include <sstream> // std::stringstream
#include <fstream> // std::ofstream

#include <boost/tuple/tuple.hpp> // boost::tuple
#include <boost/unordered_map.hpp> // for consistency reasons, we use boost::unordered_map instead of std::unordered_map
#include <boost/functional/hash.hpp> // boost::hash_combine
#include <boost/heap/fibonacci_heap.hpp> // boost::heap::fibonacci_heap

#include "geo_defines.h" // geo_cand_rtree; the better implementation is to inherit class "typed_directed_graph" and put this function in a newly derived class


//================================================================================
// Notations:
//	Starting and ending vertices are the both endpoints of a directed edge.
//	Source and target vertices are the vertices of a path.
//	Each vertex consist of an index number and a type, which is one of the following (must be lower case):
//		v - ordinary vertex;
//		f - facility vertex;
//		c - candidate vertex;
//		r - reference location vertex.
//	The index number of a vertex is a non-negative and mutually exclusive integer value, and inconsecutive values is permitted.
// Implementation:
//	Parallel edge is not allow, which is guaranteed via std::set for out-edges.
//	Use adjacent edges for each vertex, with redundant vertices of in-edges aiming at edge split and vertex replacing.
//	Reference locations are created onto graph as virtual vertices, each of which has a zero-distance unidirection edge from it to the overlapped vertex.
// Remark:
//	Facilities are strictly forbidden to overlap with each other.
//	Candidates are strictly forbidden to overlap with facilities or each other.
//	The edges file format must be "edge_id vs_id ve_id dist".
//	The facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve".
//	The candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat".
//	The reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat".

//================================================================================
// vertex associated with type

struct vertex
{
	int id; // vertex index number
	char type; // vertex type: 'v', 'f', 'c', 'r'

	vertex()
		: id(-1) // default with invalid id (-1)
		, type('v') {} // default is ordinary vertex 'v'

	vertex(int id_in, char type_in)
		: id(id_in)
		, type(type_in) {}

	vertex(const vertex &rhs)
	{
		if (this == &rhs)
			return;
		id = rhs.id;
		type = rhs.type;
	}

	vertex& operator=(const vertex &rhs)
	{
		if (this == &rhs)
			return *this;
		id = rhs.id;
		type = rhs.type;
		return *this;
	}

	bool operator<(const vertex &rhs) const
	{
		if (type == rhs.type) // if with the same type
			return id < rhs.id; // first compare id
		else
			return type < rhs.type; // then compare type
	}

	bool operator==(const vertex &rhs) const
	{
		return id == rhs.id && type == rhs.type;
	}

	// boost::hash is implemented by calling the function hash_value
	friend std::size_t hash_value(const vertex& v)
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, v.id);
		boost::hash_combine(seed, v.type);
		return seed;
	}
};

//================================================================================
// directed out-edge with edge weight (distance) but without starting vertex

struct out_edge
{
	vertex ve; // ending vertex
	float dist; // edge weight (distance)

	out_edge()
		: ve() // ve.id is default -1
		, dist(1.0f) {} // default weight

	out_edge(int ve_id, char ve_type, float dist_in)
		: ve(ve_id, ve_type)
		, dist(dist_in) {}

	out_edge(const out_edge &rhs)
	{
		if (this == &rhs)
			return;
		ve = rhs.ve;
		dist = rhs.dist;
	}

	bool operator<(const out_edge &rhs) const
	{
		return ve < rhs.ve;
	}

	bool operator==(const out_edge &rhs) const
	{
		return ve.id == rhs.ve.id && ve.type == rhs.ve.type && dist == rhs.dist;
	}
};

//================================================================================
// a tuple that represents a directed edge

typedef boost::tuple < vertex, // starting vertex
	vertex, // ending vertex
	float > // edge weight (distance)
directed_edge;

enum directed_edge_element { VS = 0, VE, DIST };

//================================================================================
// hash map for <vertex, adjacent edges> is implemented based on unordered_map;
// adjacent edges is a pair std::sets of out-edges (out_edge) and in-edges (vertex)

typedef std::pair < std::set<out_edge>, // out-edges
	std::set<vertex> > // starting vertice of in-edges
vs_adjacent_edges;

typedef boost::unordered_map < vertex, // starting vertex
	vs_adjacent_edges > // adjacent edges of starting vertex
adjacent_edges_hash_map;

//================================================================================
// vertex that has been inserted to split an original (ordinary) edge of graph

struct inserted_vertex
	: public vertex
{
	float dist_vs; // the distance from vs to vertex "inserted"

	inserted_vertex(int id_in, char type_in, float dist_vs_in)
		: vertex(id_in, type_in)
		, dist_vs(dist_vs_in) {}
};

//--------------------------------------------------------------------------------
// hash map for split edges with the corresponding inserted vertices is implemented based on unordered_map

typedef boost::unordered_map < std::string, // starting and ending vertice string, which format must be "vs_id-ve_id"
	std::vector<inserted_vertex> > // vertices that have been inserted to split <vs, ve>
split_edges_hash_map;

//================================================================================
// target vertex associated with the shortest path distance from source vertex

struct target_vertex
	: public vertex
{
	float dist; // the shortest path distance from source vertex

	target_vertex(int id_in, char type_in, float dist_in)
		: vertex(id_in, type_in)
		, dist(dist_in) {}

	target_vertex(const target_vertex &rhs)
	{
		if (this == &rhs)
			return;
		id = rhs.id;
		type = rhs.type;
		dist = rhs.dist;
	}

	// remark: the same vertex validation is not implemented
	bool operator>(const target_vertex &rhs) const
	{
		return dist > rhs.dist; // only compare shortest distance
	};
};

//--------------------------------------------------------------------------------
// a mutable fibonacci min-heap for ordering target vertices without indication of source vertex, which is used in dijkstra algorithm;
// an assistant hash map for storing vertices and the corresponding min-heap handles & current shortest distance from source vertex;
// a hash map for storing vertices with definite shortest distance from source vertex

typedef boost::heap::fibonacci_heap < target_vertex, // target vertex associated with the (current) shortest path distance from source vertex to itself
	boost::heap::compare < std::greater < target_vertex > >, // shift default max-heap to min-heap by replacing "less" function with "greater"
	boost::heap::stable <false> > // mutable min-heap
vertex_min_fibonacci_heap;

typedef vertex_min_fibonacci_heap::handle_type vertex_handle; // handle of target vertex in the min-heap

typedef boost::unordered_map < vertex, // target vertex (base class is enough)
	std::pair < vertex_handle, // handle of target vertex in the min-heap
	float > > // current shortest distance from source vertex to target vertex
target_hash_map; // assistant hash map

typedef boost::unordered_map < vertex, // target vertex with definite shortest path distance from source vertex
	float > // the definite shortest path distance
definite_hash_map;

//================================================================================
// voronoi vertex

typedef boost::unordered_map < vertex,
	boost::tuple < float, // (current) distance to NN facility (or candidate)
	int, // NN facility (or candidate) id
	char, // NN type: facility 'f' (or candidate 'c')
	bool > > // optimal mark (true) or not (false)
voronoi_hash_map;

enum voronoi_vertex_element { CUR_DIST = 0, NN_ID, NN_TYPE, MARK };

//================================================================================
// candidate associated with the expected reduction of distance (ERD)

struct candidate
{
	int id; // candidate index number
	float ERD; // the expected reduction of distance (ERD)

	candidate()
		: id(-1)
		, ERD(0.0f) {}

	candidate(int id_in, float ERD_in)
		: id(id_in)
		, ERD(ERD_in) {}

	// remark: the same candidate validation is not implemented
	bool operator<(const candidate &rhs) const
	{
		return ERD < rhs.ERD; // only compare ERD
	};
};

//--------------------------------------------------------------------------------
// a mutable fibonacci max-heap for ordering candidates by ERD;
// an assistant hash map for storing max-heap handles and ERDs of candidates

typedef boost::heap::fibonacci_heap<candidate> cand_max_fibonacci_heap;

typedef cand_max_fibonacci_heap::handle_type cand_handle; // handle of candidate in the fibonacci max-heap

typedef boost::unordered_map < int, // candidate id
	std::pair < cand_handle, // handle of candidate in the fibonacci max-heap
	float > > // ERD of candidate
cand_hash_map; // assistant hash map

//================================================================================
// facilities associated with the endpoints

typedef boost::unordered_map < int, // facility id
	boost::tuple < int, // vs_id
	int, // ve_id
	float, // vs_dist
	float > > // ve_dist
fac_hash_map;

enum fac_hash_map_element { VS_ID = 0, VE_ID, VS_DIST, VE_DIST };

//--------------------------------------------------------------------------------
// reference locations associated with attractor distances

typedef boost::unordered_map < int, // reference location id
	boost::tuple < int, // vs_id
	int, // ve_id
	float, // vs_dist
	float, // ve_dist
	float, // prob
	int, // NN facility id
	float > > // distance to NN facility
refloc_hash_map;

enum refloc_hash_map_element { RF_PROB = 4, RF_NN_F, RF_NN_DIST }; // other use "fac_hash_map_element"

//================================================================================
// base classes of all event visitors, which are function objects and invoked during traversing graph by dijkstra algorithm;
// "top event" occurs when the shortest path distance of any target vertex (from the source vertex) is determined;
// "extend event" occurs when new out-vertices of the current target vertex are being extended

class top_event_visitor
{
public:
	top_event_visitor() : no_need_extend_out_edges(false) {} // default Dijkstra mode
	virtual bool // return "true" to terminate graph traversal (by dijkstra algorithm); otherwise (false), continue the traversal
		operator()(const target_vertex &v, // top target vertex
		const definite_hash_map &hash_map) = 0; // the hash map for storing vertices with definite shortest distance from source vertex
	bool no_need_extend_out_edges; // indicate to extend out edges of the current visited vertex (false, default Dijkstra mode) or not (true)
};

class extend_event_visitor
{
public:
	virtual bool // return "true" to terminate graph traversal (by dijkstra algorithm); otherwise (false), continue the traversal
		operator()(const target_vertex &v, // the vertex has definite shortest path distance
		const out_edge &e) = 0; // the extended out-edge from the target vertex
};

//================================================================================
// typed directed graph

class typed_directed_graph
{
public:
	typed_directed_graph()
		: is_to_bidirectional(false)
		, is_graph_init(false)
		, is_facilities_deployed(false)
		, is_candidates_deployed(false) {}
	virtual ~typed_directed_graph(){}

public:
	std::ofstream *out_file; // the output file for testing status texts; default value is NULL, which means no need to output status texts
	bool is_to_bidirectional; // indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
	bool is_graph_init; // indicate whether graph is succussfully initialized (true) or not (false)
	bool is_facilities_deployed; // indicate whether facilities are succussfully deployed (true) or not (false)
	bool is_candidates_deployed; // indicate whether candidates are succussfully deployed (true) or not (false)
	adjacent_edges_hash_map adjacent_edges; // vertices (key) each of which binds the corresponding out-edges and in-edges (value)
	split_edges_hash_map split_edges; // split edges with the corresponding inserted vertices

public:
	static bool init_graphs(const char *edges_file_path, typed_directed_graph &graph, typed_directed_graph &reverse_graph);

public:
	void set_testing(std::ofstream *out_file_in = NULL) { out_file = out_file_in; };
	void set_to_bidirectional(bool to_bidirectional) { is_to_bidirectional = to_bidirectional; };
	bool get_is_facilities_deployed() { return is_facilities_deployed; };
	bool get_is_candidates_deployed() { return is_candidates_deployed; };
	bool init_graph(const char *edges_file_path, bool reverse_graph = false);
	bool deploy_facilities(const char *facs_file_path, fac_hash_map *hash_map = NULL);
	bool deploy_reflocs(const char *reflocs_file_path, refloc_hash_map *hash_map = NULL);
	bool deploy_candidates(const char *cands_file_path, cand_max_fibonacci_heap *max_heap = NULL, cand_hash_map *hash_map = NULL);
	void insert_vertex(const vertex &vs, const vertex &ve, const vertex &inserted, float dist_vs_inserted, float dist_inserted_ve,
		std::vector<directed_edge> *removed_edges = NULL, std::vector<directed_edge> *inserted_edges = NULL, bool is_replace_overlap = false);
	void restore_graph(const vertex &inserted, const std::vector<directed_edge> &removed_edges, const std::vector<directed_edge> &inserted_edges);
	void dijkstra(const vertex &source, top_event_visitor *top_visitor = NULL, extend_event_visitor *extend_visitor = NULL);
	void voronoi_diagram(voronoi_hash_map &voronoi_vertices, refloc_hash_map &attr_dists);

	// the better implementation is to inherit class "typed_directed_graph" and put this function in a newly derived class
	bool deploy_candidates_rtree(const char *cands_file_path, geo_cand_rtree &cand_rtree);

private:
	void overlap_vertex(const vertex &overlapped, const vertex &inserted, std::vector<directed_edge> *removed_edges = NULL,
		std::vector<directed_edge> *inserted_edges = NULL, bool is_replace_overlap = false);
	bool split_edge(const vertex &vs, const vertex &ve, const vertex &inserted, float dist_vs_inserted, float dist_inserted_ve,
		std::vector<directed_edge> *removed_edges = NULL, std::vector<directed_edge> *inserted_edges = NULL);
	void split_split_edge(const split_edges_hash_map::iterator &iter, int vs_id, int ve_id, const vertex &inserted,
		float dist_vs, float dist_ve, std::vector<directed_edge> *removed_edges, std::vector<directed_edge> *inserted_edges);
	void replace_vertex(const vertex &replaced, const vertex &inserted, std::vector<directed_edge> *removed_edges = NULL,
		std::vector<directed_edge> *inserted_edges = NULL);
};