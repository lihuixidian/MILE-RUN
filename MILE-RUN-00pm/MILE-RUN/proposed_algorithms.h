#pragma once

#include <string> // std::string
#include <fstream> // std::ifstream
#include <chrono> // std::chrono::high_resolution_clock

#include "geo_defines.h" // R*-tree, EARTH_RADIUS, PI
#include "typed_directed_graph.h"

//================================================================================
// REMARK:
//	if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp), we view its ERD as 0


//================================================================================
// testing for output edges in graph

//#define ALGO_TESTING
#define ALGO_TESTING_NO_TRACE_GRAPH
//#define ALGO_TESTING_NO_EARLY_STOPPING // early stopping "ol.ERD > next upper bound" is ignored for testing

static std::ofstream algo_out_file("D:\\Experiment\\MILE-RUN\\datasets\\algo_testing.txt", std::ofstream::out | std::ofstream::trunc);

void print_graph(const typed_directed_graph &graph)
{
	std::ofstream out_file("D:\\Experiment\\MILE-RUN\\datasets\\graph.txt", std::ofstream::out | std::ofstream::trunc);
	if (!out_file.fail())
	{
		out_file << "vs ve dist\n---------------\n";
		int edges_count = 0;
		adjacent_edges_hash_map::const_iterator citer_v = graph.adjacent_edges.cbegin();
		for (; citer_v != graph.adjacent_edges.cend(); ++citer_v)
		{
			std::set<out_edge>::const_iterator citer_out_edge = citer_v->second.first.cbegin();
			for (; citer_out_edge != citer_v->second.first.cend(); ++citer_out_edge)
			{
				out_file << citer_v->first.type << citer_v->first.id << ' ' << citer_out_edge->ve.type << citer_out_edge->ve.id << ' ' << citer_out_edge->dist << '\n';
				++edges_count;
			}
		}
		out_file << "---------------\nedges count: " << edges_count << '\n';
		out_file.flush();
	}
}


#pragma region ALGO_EN
//================================================================================
// top event visitor for EN algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			max-fibonacci-heap will never record any ERD information of the reference location, and we view its ERD as 0

class EN_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'c')
		{
			dists[v.id] = v.dist; // record the shortest path distance from reference location to the candidate
#ifdef ALGO_TESTING
			algo_out_file << "top: c" << v.id << ", d(r,c) = " << v.dist << '\n';
#endif
		}
		else if (v.type == 'f') // nearest facility is found
		{
#ifdef ALGO_TESTING
			algo_out_file << "NND = " << v.dist << '\n';
#endif
			boost::unordered_map<int, float>::iterator iter = dists.begin();
			for (; iter != dists.end(); ++iter) // iterate each encountered candidate before the nearest facility
			{
#ifdef ALGO_TESTING
				float old_ERD = (*ptr_hash_map)[iter->first].second;
#endif
				cand_handle handle = (*ptr_hash_map)[iter->first].first; // retrieve the max-heap handle of candidate
				float new_ERD = (*ptr_hash_map)[iter->first].second += (v.dist - iter->second) * prob; // left to "+=" is ERD, right to "+=" is ERd
				ptr_max_heap->increase(handle, candidate(iter->first, new_ERD)); // update ordered ERD value of the candidate in max-heap
#ifdef ALGO_TESTING
				algo_out_file << "ERD: c" << iter->first << ", current: " << old_ERD << ", new: " << new_ERD << '\n';
#endif
			}

			// this sentence is removed with reset_dists() function, which is invoked outside before traversing each reference location
			//dists.clear(); // reset temporary hash map for next reference location

			return true; // a flag to terminate graph traversal
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(cand_max_fibonacci_heap *max_heap, cand_hash_map *hash_map){ ptr_max_heap = max_heap; ptr_hash_map = hash_map; }
	void set_prob(float prob_in) { prob = prob_in; }
	void reset_dists() { dists.clear(); }

private:
	cand_max_fibonacci_heap *ptr_max_heap; // a pointer to a mutable fibonacci max-heap ordered candidates by ERD
	cand_hash_map *ptr_hash_map; // a pointer to an assistant hash map for storing max-heap handles and ERDs of candidates
	float prob; // present probability of the reference location
	boost::unordered_map < int, // candidate id
		float > // shortest path distance from reference location to the candidate
		dists; // hash map for temporarily storing shortest distances from a reference location to encountered candidates before nearest facility is found
};

//--------------------------------------------------------------------------------
// EN algorithm
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [re] candidate: the optimal candidate; default candidate object (id is -1) means candidates are not deployed successfully or opening reference locations file fails
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [in] cands_file_path: file path of candidates
// [in] k: return top-k candidates, defualt is 1; as k must >= 1, 0 is viewed as 1
// [in] topk: the vector for recording top-k candidates, default is NULL

candidate EN(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path, const char *cands_file_path,
			 unsigned k = 1, std::vector<candidate> *topk = NULL)
{
	// init graph, and deploy facilities
	typed_directed_graph graph;
#ifndef ALGO_TESTING_NO_TRACE_GRAPH
	graph.set_testing(&algo_out_file);
#endif
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);
	graph.deploy_facilities(facs_file_path);
	
	// deploy candidates and construct related data structures
	cand_max_fibonacci_heap max_heap;
	cand_hash_map hash_map;
	graph.deploy_candidates(cands_file_path, &max_heap, &hash_map);

	if (!graph.get_is_candidates_deployed()) // check graph, facilities and candidates status
		return candidate();

	// create an top event visitor for EN algorithm
	EN_visitor visitor;
	visitor.set_data_structures(&max_heap, &hash_map); // set data structures

	// deal with each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			// lon and lat are useless for graph

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

#ifdef ALGO_TESTING
			algo_out_file << "Ref Loc: " << refloc_id << ", prob = " << prob << '\n';
#endif
			// traverse graph from the reference location by dijkstra algorithm
			visitor.set_prob(prob); // set the present probability of the reference location
			visitor.reset_dists(); // reset the hash map that temporarily stores shortest distances from a reference location to encountered candidates before nearest facility is found
			graph.dijkstra(refloc_v, &visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		candidate ol = max_heap.top();
#ifdef ALGO_TESTING
		algo_out_file << "optimal: " << ol.id << ", ERD = " << ol.ERD << '\n';
		algo_out_file << std::endl;
#endif

		// return k-top candidates
		if (topk != NULL && k > 1)
		{
			unsigned n = 0;
			while (n < k)
			{
				topk->push_back(max_heap.top());
				++n;
				max_heap.pop(); // pop the top vertex
			}
		}

		return ol; // the optimal candidate
	}
	return candidate(); // opening reference locations file fails
}
#pragma endregion ALGO_EN


//================================================================================
// Local Network Table for LNB algorithm;
// an assistant hash map for local network edges, which consists of shortest distances from a reference location to a vertex, and that distance plus out-edge weight (distance)
// an assistant hash map for information of candidates
// an assistant hash map for reference locations that overlap some vertices

typedef boost::unordered_map < std::pair < int, int >, // ids of starting and ending vertices of a directed edge, whose format is <vs_id, ve_id>
	boost::tuple < float, // prob
	float, // dist(r, vs), use DDU for alternative
	float, // dist(r, vs) + dist(vs, ve), use DDL for alternative
	int > > // ve_id, use OFFSET for alternative
loc_edges_hash_map; // assistant hash map

typedef boost::tuple < float, // prob
	float, // Dd+ (Dd upper bound)
	float, // Dd- (Dd lower bound)
	float > // offset
ref_loc_entry;

enum ref_loc_entry_element { PROB = 0, DDU, DDL, OFFSET }; // OFFSET also represent "ve_id" when for "loc_edges_hash_map", as "VE_ID" has been defined in "fac_hash_map_element"

typedef boost::unordered_map < std::pair < int, int >, // ids of starting and ending vertices of a directed edge, whose format is <vs_id, ve_id>
	std::pair < float, // ERD+ (ERD upper bound on a directed edge)
	std::vector<ref_loc_entry> > >
LNT_hash_map;

typedef boost::unordered_map < int, // candidate id
	boost::tuple < LNT_hash_map::const_iterator, // const iterator pointing to edge "vs_id-ve_id" in LNT
	LNT_hash_map::const_iterator, // const iterator pointing to edge "ve_id-vs_id" in LNT
	float, // dist_vs
	float > > // dist_ve
cand_info_hash_map; // assistant hash map

enum cand_info_hash_map_element { ITER_VS_VE = 0, ITER_VE_VS }; // the value of VS_DIST, VE_DIST are defined and consistent with "fac_hash_map_element"

typedef boost::unordered_map < int, // vertex id
	float > // accumulated nnd
vertex_nnd_hash_map; // assistant hash map

//--------------------------------------------------------------------------------
// top event visitor for LNB algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			LNT and vertex_nnd hash maps will never record any local network of the reference location, and we view its ERD as 0

class LNB_top_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map)
	{
		if (v.type == 'f') // nearest facility is found
		{
			float nnd = v.dist; // nearset neighbor distance

			if (r_v_id != -1) // the reference location overlaps a vertex
			{
#ifdef ALGO_TESTING
				algo_out_file << "r overlaps v" << r_v_id << ", nnd = " << v.dist << '\n';
#endif
				vertex_nnd_hash_map::iterator iter_nnd = ptr_nnd_hash_map->find(r_v_id);
				if (iter_nnd != ptr_nnd_hash_map->end()) // the vertex has been overlapped by some reference location
					iter_nnd->second += nnd; // accumulate nnd
				else // the vertex is overlapped by this reference location first time
					(*ptr_nnd_hash_map)[r_v_id] = nnd; // initial nnd
			}

			// deal with each edge in the local network of a reference location
			loc_edges_hash_map::iterator iter_edge = ptr_loc_edges_hash_map->begin();
			for (; iter_edge != ptr_loc_edges_hash_map->end(); ++iter_edge)
			{
				// for clarity, we define two variables below
				std::pair<int, int> vs_ve = iter_edge->first; // the edge vs->ve, <vs_id, ve_id>
				int ve_id = iter_edge->second.get<OFFSET>(); // use OFFSET for alternative, actually is "ve_id"

				// compute offset
				float offset = 0.0f;
				definite_hash_map::const_iterator iter_definite = hash_map.find(vertex(ve_id, 'v'));
				if (iter_definite != hash_map.cend() // ve has definite shortest path distance d(r, ve)
					&& iter_definite->second < iter_edge->second.get<DDL>()) // d(r, ve) < d(r, vs) + d(vs, ve)
					offset = (iter_edge->second.get<DDL>() - iter_definite->second) / 2.0f; // (d(r, vs) + d(vs, ve) - d(r, ve)) / 2

				// create new entry in LNT
				float Dd_u = nnd - iter_edge->second.get<DDU>(); // Dd+ = nnd - dist(r, vs)
				LNT_hash_map::iterator iter_LNT = ptr_LNT_hash_map->find(vs_ve);
				if (iter_LNT != ptr_LNT_hash_map->end()) // edge "vs_id-ve_id" has already been in local network of other reference locations
				{
					iter_LNT->second.first += (Dd_u * iter_edge->second.get<PROB>()); // new_ERD+ = old_ERD+ + (Dd+ * prob)
					iter_LNT->second.second.push_back(boost::make_tuple(iter_edge->second.get<PROB>(), // prob
						Dd_u, // Dd+
						nnd - iter_edge->second.get<DDL>(), // Dd- = nnd - (dist(r, vs) + dist(vs, ve))
						offset)); // offset
				}
				else // "vs_id-ve_id" first time presents
				{
					(*ptr_LNT_hash_map)[vs_ve] = std::make_pair(Dd_u, // Dd+ as initial ERD+
						std::vector<ref_loc_entry>());
					(*ptr_LNT_hash_map)[vs_ve].second.push_back(boost::make_tuple(iter_edge->second.get<PROB>(), // prob
						Dd_u, // Dd+
						nnd - iter_edge->second.get<DDL>(), // Dd- = nnd - (dist(r, vs) + dist(vs, ve))
						offset)); // offset
				}
			}
			return true; // a flag to terminate graph traversal
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(LNT_hash_map *LNT_hash_map_in, loc_edges_hash_map *loc_edges_hash_map_in, vertex_nnd_hash_map *nnd_hash_map)
		{ ptr_LNT_hash_map = LNT_hash_map_in; ptr_loc_edges_hash_map = loc_edges_hash_map_in; ptr_nnd_hash_map = nnd_hash_map; }
	void set_overlap(int r_v_id_in) { r_v_id = r_v_id_in; };

private:
	LNT_hash_map *ptr_LNT_hash_map; // a pointer to a hash map as Local Network Table
	loc_edges_hash_map *ptr_loc_edges_hash_map; // a pointer to an assistant hash map for local network edges
	vertex_nnd_hash_map *ptr_nnd_hash_map; // a pointer to an assistant hash map for reference locations that overlap some vertices
	int r_v_id; // indicate whether this reference location overlaps a vertex (v_id) or not (-1)
};

//--------------------------------------------------------------------------------
// extend event visitor for LNB algorithm

class LNB_extend_visitor
	: public extend_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const out_edge &e)
	{
		if (v.type != 'r') // impossible to be 'f', then must be 'v'
		{
			if (e.ve.type != 'f') // must be 'v' (e.dist is impossible to be zero), then v->e.ve
			{
				(*ptr_loc_edges_hash_map)[std::pair<int, int>(v.id, e.ve.id)] = boost::make_tuple(prob, v.dist, v.dist + e.dist, e.ve.id);
#ifdef ALGO_TESTING
				algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + e.dist << '\n';
#endif
				return false; // a flag to continue graph traversal
			}
			else // e.ve.type == 'f'
			{
				if (e.dist == 0) // f overlaps v, then no need to consider this virtual edge
					return false; // a flag to continue graph traversal
				else // e.dist != 0, e.ve is a facility that must split an edge
				{
					fac_hash_map::iterator iter_fac = ptr_fac_hash_map->find(e.ve.id); // the facility must exist
					if (v.id == iter_fac->second.get<VS_ID>()) // f.vs(v)->f.ve
					{
						(*ptr_loc_edges_hash_map)[std::pair<int, int>(v.id, iter_fac->second.get<VE_ID>())] =
							boost::make_tuple(prob, v.dist, v.dist + iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>(), iter_fac->second.get<VE_ID>());
#ifdef ALGO_TESTING
						algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>() << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
					else // v.id == iter_fac->second.get<VE_ID>(), f.ve(v)->f.vs
					{
						(*ptr_loc_edges_hash_map)[std::pair<int, int>(v.id, iter_fac->second.get<VS_ID>())] =
							boost::make_tuple(prob, v.dist, v.dist + iter_fac->second.get<VE_DIST>() + iter_fac->second.get<VS_DIST>(), iter_fac->second.get<VS_ID>());
#ifdef ALGO_TESTING
						algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + iter_fac->second.get<VE_DIST>() + iter_fac->second.get<VS_DIST>() << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
				}
			}
		}
		else // v.type == 'r'
		{
			if (e.dist == 0) // v(r) overlaps e.ve, whatever e.ve.type is 'v' or 'f', no need to consider the virtual edge
				return false; // a flag to continue graph traversal
			else // e.dist != 0, v(r) must split an edge
			{
				if (e.ve.type != 'f') // e.ve is ordinary vertex
				{
					if (e.ve.id == ve_id) // r.vs->r.ve
					{
						(*ptr_loc_edges_hash_map)[std::pair<int, int>(vs_id, ve_id)] = boost::make_tuple(prob, 0.0f, e.dist, ve_id);
#ifdef ALGO_TESTING
						algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << e.dist << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
					else // e.ve.id == vs_id, r.ve->r.vs
					{
						(*ptr_loc_edges_hash_map)[std::pair<int, int>(ve_id, vs_id)] = boost::make_tuple(prob, 0.0f, e.dist, vs_id);
#ifdef ALGO_TESTING
						algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << e.dist << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
				}
				else // e.ve.type == 'f', r and f locate on and split the same edge
				{
					fac_hash_map::iterator iter_fac = ptr_fac_hash_map->find(e.ve.id); // the facility must exist
					if (iter_fac->second.get<VS_ID>() == vs_id)
					{
						if (iter_fac->second.get<VS_DIST>() < dist_vs) // v(r).ve->v(r).vs
						{
							(*ptr_loc_edges_hash_map)[std::pair<int, int>(ve_id, vs_id)] = boost::make_tuple(prob, 0.0f, dist_vs, vs_id);
#ifdef ALGO_TESTING
							algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_vs << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
						else // iter_fac->second.get<VS_DIST>() > dist_vs (== never holds), v(r).vs->v(r).ve
						{
							(*ptr_loc_edges_hash_map)[std::pair<int, int>(vs_id, ve_id)] = boost::make_tuple(prob, 0.0f, dist_ve, ve_id);
#ifdef ALGO_TESTING
							algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_ve << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
					}
					else // iter_fac->second.get<VS_ID>() == ve_id
					{
						if (iter_fac->second.get<VS_DIST>() < dist_ve) // v(r).vs->v(r).ve
						{
							(*ptr_loc_edges_hash_map)[std::pair<int, int>(vs_id, ve_id)] = boost::make_tuple(prob, 0.0f, dist_ve, ve_id);
#ifdef ALGO_TESTING
							algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_ve << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
						else // iter_fac->second.get<VS_DIST>() > dist_ve (== never holds), v(r).ve->v(r).vs
						{
							(*ptr_loc_edges_hash_map)[std::pair<int, int>(ve_id, vs_id)] = boost::make_tuple(prob, 0.0f, dist_vs, vs_id);
#ifdef ALGO_TESTING
							algo_out_file << ss_vs_ve.str() << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_vs << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
					}
				}
			}
		}
	}

	void set_data_structure(fac_hash_map *fac_hash_map_in, loc_edges_hash_map *loc_edges_hash_map_in)
		{ ptr_fac_hash_map = fac_hash_map_in; ptr_loc_edges_hash_map = loc_edges_hash_map_in; }
	void set_vs_ve(int vs_in, int ve_in, float dist_vs_in, float dist_ve_in) { vs_id = vs_in; ve_id = ve_in; dist_vs = dist_vs_in; dist_ve = dist_ve_in; }
	void set_prob(float prob_in) { prob = prob_in; }

private:
	fac_hash_map *ptr_fac_hash_map; // a pointer to a hash map of facilities associated with endpoints
	loc_edges_hash_map *ptr_loc_edges_hash_map; // a pointer to an assistant hash map for local network edges
	int vs_id, ve_id; // the original starting and ending ordinary vertices of the edge which the reference location splits
	float dist_vs, dist_ve; // the distances from reference location to vs and ve
	float prob; // present probability of the reference location
};

//--------------------------------------------------------------------------------
// LNB algorithm: construct LNT
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] LNT; a hash map as Local Network Table, which needs to be constructed in this function
// [out] nnd_map; a hash map for reference locations that overlap some vertices, which needs to be constructed in this function
// [out] graph: the typed directed graph, which needs to be constructed in this function
// [out] facs_map: a hash map for storing facilities associated with the endpoints, which needs to be constructed in this function

void LNB_construct_LNT(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
					   LNT_hash_map &LNT, vertex_nnd_hash_map &nnd_map, typed_directed_graph &graph, fac_hash_map &facs_map)
{
	// init graph, and deploy facilities
#ifndef ALGO_TESTING_NO_TRACE_GRAPH
	graph.set_testing(&algo_out_file);
#endif
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);

	// deploy facilities and construct related data structure
	graph.deploy_facilities(facs_file_path, &facs_map);

	if (!graph.get_is_facilities_deployed()) // check graph and facilities status
		return;

	// create an extend event visitor for LNB algorithm
	LNB_extend_visitor extend_visitor;
	loc_edges_hash_map refloc_loc_edges_hash_map; // an assistant hash map for local network edges
	extend_visitor.set_data_structure(&facs_map, &refloc_loc_edges_hash_map); // set data structures

	// create a top event visitor for LNB algorithm
	LNB_top_visitor top_visitor;
	top_visitor.set_data_structures(&LNT, &refloc_loc_edges_hash_map, &nnd_map); // set data structures

	// deal with each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			// lon and lat are useless for network traversal

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			extend_visitor.set_vs_ve(vs_id, ve_id, dist_vs, dist_ve); // the original starting and ending vertices of the reference location
			extend_visitor.set_prob(prob); // set the present probability of the reference location
			refloc_loc_edges_hash_map.clear(); // must reset the assistant hash map for each new reference location
			if (vs_id == ve_id)
				top_visitor.set_overlap(vs_id); // reference location overlaps a vertex
			else
				top_visitor.set_overlap(-1); // reference location doesn't overlap a vertex
			graph.dijkstra(refloc_v, &top_visitor, &extend_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}
}

//--------------------------------------------------------------------------------
// LNB algorithm: query optimal candidate based on LNT
// remark: the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat";
//		   for simplicity, we use LNT as a non-const parameter, hence, in this function, be careful to use it, and DON'T change its value anytime
// note: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp), we view its ERD as 0
// [re] candidate: the optimal candidate; default candidate object means opening candidates file fails
// [out] checked_num: the number of candidates actually checked
// [in] LNT: the constructed hash map as LNT
// [in] nnd_map; the constructed hash map for reference locations that overlap some vertices
// [in] graph: the constructed graph
// [in] facs_map: a hash map for storing facilities associated with the endpoints
// [in] cands_file_path: file path of candidates

candidate LNB_query(int &checked_num, LNT_hash_map &LNT, const vertex_nnd_hash_map &nnd_map, const typed_directed_graph &graph,
					const fac_hash_map &facs_map, const char *cands_file_path)
{
	cand_max_fibonacci_heap max_heap; // max-heap ordered by ERD+ in LNT
	cand_info_hash_map hash_map; // hash map for information of candidates

	// construct a max-heap for candidates
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail()) // candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, lon_pos - dist_ve_pos - 1).c_str()));
			// lon and lat are useless for query

			if (vs_id != ve_id) // candidate splits an edge
			{
				// here, we explicitly indicate "vs" and "ve"
				std::pair<int, int> vs_ve(vs_id, ve_id), // vs->ve
					ve_vs(ve_id, vs_id); // ve->vs, reverse direction

				// retrieve ERD+ of edge
				float ERD_u = 0.0f; // ERD+ for ordering candidates in max-heap
				LNT_hash_map::const_iterator iter_vs_ve = LNT.find(vs_ve); // vs->ve
				LNT_hash_map::const_iterator iter_ve_vs = LNT.find(ve_vs); // ve->vs
				if (iter_vs_ve != LNT.cend() && iter_ve_vs != LNT.cend()) // vs->ve & ve->vs are both valid
					ERD_u = iter_vs_ve->second.first + iter_ve_vs->second.first; // larger than upper
				else if (iter_vs_ve != LNT.cend()) // only vs->ve is valid
					ERD_u = iter_vs_ve->second.first;
				else if (iter_ve_vs != LNT.cend()) // only ve->vs is valid
					ERD_u = iter_ve_vs->second.first;
				else // neither are valid, ERD+ should be 0.0f
					continue; // no need to consider this candidate

				max_heap.push(candidate(cand_id, ERD_u)); // push each candidate
				hash_map[cand_id] = boost::make_tuple(iter_vs_ve, iter_ve_vs, dist_vs, dist_ve); // record information of each candidate
			}
			else // vs_id == ve_id, candidate overlaps a vertex
				 // NOTE: 1) facility cannot overlap the vertex, otherwise, candidate overlaps facility, which is strictly forbidden;
				 //			 hence, any network traversal cannot end at the vertex, thus considering only out-edges can avoid double counting,
				 //			 where candidate is viewed as an intermediate vertex in a path, namely as vs and ve for two different edges, for computing ERd;
				 //		  2) also, if traversal ends immediately before the vertex (i.e., r->f->c), candidate is outside local network and has no benefit;
				 //		  3) moreover, if a reference location splits an edge with candidate (the vertex) as one of the two endpoints,
				 //			 namely out-edge c->r->ve, then Rd of c must > nnd (Dd+), thus it has to consider the reverse direction
				 //		  4) conversely, only considering in-edges is not enough, the counter example is reference location overlaps the vertex,
				 //			 and for this reference location, there is only out-edge
				 //		  5) in the special case, where c overlaps v and some r(s) also overlap v, we compute ERD with the help of vertex_nnd_hash_map
			{
				float ERD = 0.0f; // for a candidate that overlaps a vertex, we directly compute its ERD, because ERd can be obtained via Dd+ or Dd-

				adjacent_edges_hash_map::const_iterator iter_vertex = graph.adjacent_edges.find(vertex(vs_id, 'v')); // the overlapped vertex
				
				// considering only out-edges is enough; hence, no need to consider in-edges
				std::set<out_edge>::const_iterator iter_out_edge = iter_vertex->second.first.cbegin();
				for (; iter_out_edge != iter_vertex->second.first.cend(); ++iter_out_edge) // compute ERd for each out-edge
				{
					// fault-tolerant check in case out-vertex is a facility
					int out_v_id = iter_out_edge->ve.id;
					float dist_vs_ve = iter_out_edge->dist;
					if (iter_out_edge->ve.type != 'v') // must be 'f'
					{
						fac_hash_map::const_iterator iter_fac = facs_map.find(out_v_id); // as facility exists, this iterator must not be end()
						if (iter_fac->second.get<VS_ID>() == vs_id)
							out_v_id = iter_fac->second.get<VE_ID>();
						else
							out_v_id = iter_fac->second.get<VS_ID>();
						dist_vs_ve = iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>(); // edge distance
					}

					LNT_hash_map::const_iterator iter_vs_ve = LNT.find(std::make_pair(vs_id, out_v_id)); // vs->ve
					if (iter_vs_ve != LNT.cend()) // vs->ve takes effect on ERd
					{
						std::vector<ref_loc_entry>::const_iterator iter_entry = iter_vs_ve->second.second.cbegin();
						for (; iter_entry != iter_vs_ve->second.second.cend(); ++iter_entry) // each reference location related to this out-edge
						{
							// 1) we view edge distance as d(c, ve), then d(c, ve) must > offset, thus offset has no effect
							// 2) c (say the overlapped vertex) has definite distance from r, hence d(r, vs) <= nnd, thus Rd is impossible < 0 (locates outside local network)
							float Rd = iter_entry->get<DDL>() + dist_vs_ve; // Dd- + d(c, ve) >= 0
							if (Rd > iter_entry->get<DDU>()) // > Dd+, must consider reverse direction
								continue;
							else // c takes effect
								ERD += Rd * iter_entry->get<PROB>();
						}
					}
				}

				// fault-tolerant check in case the candidate overlaps a vertex, which is also overlapped by some reference location(s)
				unsigned out_edges_size = static_cast<unsigned>(iter_vertex->second.first.size());
				vertex_nnd_hash_map::const_iterator iter_nnd = nnd_map.find(vs_id);
				if (iter_nnd != nnd_map.cend()) // some reference location(s) overlaps the vertex
					ERD -= (out_edges_size - 1) * iter_nnd->second; // ERD = ERD - size * accumulated nnd (minus all out-edges) + accumulated nnd (reserve one out-edge)

				max_heap.push(candidate(cand_id, ERD)); // push candidate and actual ERD
				hash_map[cand_id] = boost::make_tuple(LNT.cend(), LNT.cend(), 0.0f, 0.0f); // set overlap flags (i.e., LNT.cend() and 0.0f) for the candidate
			}
		}
	}
	else // opening candidates file fails
		return candidate();

	candidate ol(-1, 0.0f); // initialize optimal location
	checked_num = 0; // intialize
	while (!max_heap.empty())
	{
		// retrieve the top candidate in max-heap
		candidate top_c = max_heap.top();
#ifdef ALGO_TESTING
		algo_out_file << "ERD+: c" << top_c.id << ": " << top_c.ERD << '\n';
#endif

#ifndef ALGO_TESTING_NO_EARLY_STOPPING
		if (top_c.ERD < ol.ERD) // current ol candidate is the optimal location
			break;
#endif
		++checked_num; // need to check more candidates

		if (hash_map[top_c.id].get<VE_DIST>() != 0 && hash_map[top_c.id].get<VS_DIST>() != 0) // c splits an edge, but not overlaps a vertex, as overlapping has definite ERD
		{
			top_c.ERD = 0.0f; // prepare to accumulate actual ERD of this candidate

			// iteratively deal with two possible edges vs->ve and ve->vs
			LNT_hash_map::const_iterator iters[2] = { hash_map[top_c.id].get<ITER_VS_VE>(), hash_map[top_c.id].get<ITER_VE_VS>() };
			float dist_c_ve[2] = { hash_map[top_c.id].get<VE_DIST>(), hash_map[top_c.id].get<VS_DIST>() }; // for vs->ve, VE_DIST; for ve->vs, VS_DIST
			for (int i = 0; i < 2; ++i)
			{
				std::vector<ref_loc_entry>::const_iterator iter_entry, iter_entry_end; // use "const", as ITER_VS_VE and ITER_VE_VS are both const_iterator
				if (iters[i] != LNT.cend()) // vs->ve or ve->vs exists
				{
					// iterate all related reference locations
					iter_entry_end = iters[i]->second.second.cend();
					for (iter_entry = iters[i]->second.second.cbegin(); iter_entry != iter_entry_end; ++iter_entry)
					{
						// the case that offset > 0
						if (dist_c_ve[i] < iter_entry->get<OFFSET>()) // the condition that offset takes effect
						{
#ifdef ALGO_TESTING
							algo_out_file << iters[i]->first << ", d(ve) = " << dist_c_ve[i] << " < offset = " << iter_entry->get<OFFSET>() << '\n';
#endif

							// utilizing virtual candidate c', which is the mapping of c with respect to lc, to compute ERd
							float dist_v_c_ve = 2.0f * iter_entry->get<OFFSET>() - dist_c_ve[i]; // d(c', ve) = 2 * offset - d(c, ve)
							top_c.ERD += (iter_entry->get<DDL>() + dist_v_c_ve) * iter_entry->get<PROB>();
							continue;
						}

						float Rd = iter_entry->get<DDL>() + dist_c_ve[i]; // Dd- + d(c, ve)
						if (Rd <= 0) // c locates outside local network (<), or has no benefit (==)
							continue;
						else // Rd > 0
						{
							if (Rd > iter_entry->get<DDU>()) // > Dd+, must consider reverse direction
								continue;
							else // c takes effect
								top_c.ERD += Rd * iter_entry->get<PROB>();
						}
					}
				}
			}
		}

#ifdef ALGO_TESTING
		algo_out_file << "ERD: c" << top_c.id << ": " << top_c.ERD << '\n';
#endif
		// update ol, if necessary
		if (top_c.ERD > ol.ERD)
		{
			ol.id = top_c.id;
			ol.ERD = top_c.ERD;
		}

		max_heap.pop(); // pop the top vertex
	}

	return ol; // the optimal candidate
}


#pragma region ALGO_NSJ
//================================================================================
// NNFC hash map for NSJ algorithm;
// an assistant hash map for sets of reference locations which cover some specific candidate;
// an assistant hash map for Euclidean distance dE(r,c)

typedef boost::unordered_map < int, // reference location id
	boost::tuple < float, // prob
	geo_point, // coordinate
	float, // nnd
	geo_box, // NNFC MBR
	int > > // number of covered candidates
nnfc_hash_map; // NNFC hash map

enum nnfc_hash_map_element { // PROB = 0, use this enumeration value in "ref_loc_entry_element"
	COORDINATE = 1,	NND, NNFC_MBR, COVER_NUM };

enum direction { NORTHWARD = 0, SOUTHWARD, EASTWARD, WESTWARD };

typedef boost::unordered_map < int, // candidate id
	std::pair < float, // UB(c)
	std::vector <int> > > // ids of reference locations covering the candidate
cover_hash_map;

typedef boost::unordered_map < int, // candidate id
	float > // ERD(c)
ERD_hash_map;

//--------------------------------------------------------------------------------
// top event visitor for calculate NNFC
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//		   its NNFC radius is 0, and we view its effect on ERD of any candidate as 0, then NNFC hash map never records this reference location

class NNFC_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'f') // nearest facility is found
		{
			float nnd = v.dist; // nearset neighbor distance, namely the NNFC radius
#ifdef ALGO_TESTING
			algo_out_file << "NND(r" << refloc_id << ") = " << nnd << '\n';
#endif

			// calculate offsets coordinate for NNFC MBR
			float min_lon, min_lat, max_lon, max_lat;
			geo_offset(lon, lat, nnd, WESTWARD, min_lon);
			geo_offset(lon, lat, nnd, SOUTHWARD, min_lat);
			geo_offset(lon, lat, nnd, EASTWARD, max_lon);
			geo_offset(lon, lat, nnd, NORTHWARD, max_lat);

			(*ptr_NNFC_hash_map)[refloc_id] = // reference location id
				boost::make_tuple(prob, // present probability
				geo_point(lon, lat), // coordinate
				nnd, // NND
				geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)), // NNFC
				0); // candidates being covered

			return true; // a flag to terminate graph traversal
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structure(nnfc_hash_map *ptr_NNFC_hash_map_in) { ptr_NNFC_hash_map = ptr_NNFC_hash_map_in; }
	void set_refloc(int refloc_id_in, float lon_in, float lat_in, float prob_in) { refloc_id = refloc_id_in; lon = lon_in; lat = lat_in; prob = prob_in;	}

private:
	// calculate offset coordinate after moving to a direction
	// [in] lon_c: longitude of source coordinate
	// [in] lat_c: latitude of source coordinate
	// [in] distance: distance to be moved
	// [in] dir: the direction
	// [out] lon_or_lat: resulting longitude or latitude
	void geo_offset(float lon_c, float lat_c, float distance, direction dir, float &lon_or_lat)
	{
		if (dir == NORTHWARD || dir == SOUTHWARD) // the same longitude
		{
			float degree = 360.0f * distance / (2.0f * PI * EARTH_RADIUS);

			// lon_or_lat now is lat
			if (dir == NORTHWARD)
				lon_or_lat = lat_c + degree;
			else // dir == SOUTHWARD
				lon_or_lat = lat_c - degree;
		}
		else if (dir == EASTWARD || dir == WESTWARD) // the same latitude
		{
			float radian_latitude = lat_c * PI / 180.0f;
			float degree = 360.0f * distance / (2.0f * PI * EARTH_RADIUS * cos(radian_latitude));

			// lon_or_lat now is lon
			if (dir == EASTWARD)
				lon_or_lat = lon_c + degree;
			else // dir == WESTWARD
				lon_or_lat = lon_c - degree;
		}
	}

private:
	nnfc_hash_map *ptr_NNFC_hash_map; // a pointer to a hash map for NNFCs
	int refloc_id; // 
	float lon; // longitude of the reference location
	float lat; // latitude of the reference location
	float prob; // present probability of the reference location
};

//--------------------------------------------------------------------------------
// top event visitor for NSJ algorithm

class NSJ_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'c')
		{
#ifdef ALGO_TESTING
			float old_ERD = (*ptr_ERDs)[v.id];
#endif
			// update ERD(c)
			float ERd = ((*ptr_NNFCs)[refloc_id].get<NND>() - v.dist) * prob;
			float new_ERD = ((*ptr_ERDs)[v.id] += ERd); // update ERD(c)

			// update ol, if necessary
			if (new_ERD > ptr_ol->ERD)
			{
				ptr_ol->id = v.id;
				ptr_ol->ERD = new_ERD;
			}

#ifdef ALGO_TESTING
			algo_out_file << "ERD: c" << v.id << ", current: " << old_ERD << ", new: " << new_ERD << '\n';
#endif

			if (--((*ptr_NNFCs)[refloc_id].get<COVER_NUM>()) // a candidate is encountered, then the number of remnant candidates is decreases by 1
				== 0) // all d(r,c) <= nnd, which means all candidates have benefits, hence, no need to traverse until facility
				return true; // a flag to terminate graph traversal
		}
		else if (v.type == 'f') // nearest facility is found
		{
			(*ptr_NNFCs)[refloc_id].get<COVER_NUM>() = 0; // already traversed from this reference location, then reset this value as a flag

			return true; // a flag to terminate graph traversal
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(nnfc_hash_map *NNFCs, ERD_hash_map *ERDs, candidate *ol) { ptr_NNFCs = NNFCs; ptr_ERDs = ERDs; ptr_ol = ol; }
	void set_refloc(int refloc_id_in, float prob_in) { refloc_id = refloc_id_in; prob = prob_in; }

private:
	int refloc_id; // reference location id
	float prob; // present probability of reference location
	nnfc_hash_map *ptr_NNFCs; // a pointer to a hash map for NNFCs
	ERD_hash_map *ptr_ERDs; // a pointer to a hash map for ERDs
	candidate *ptr_ol; // a pointer to current optimal candidate
};

//--------------------------------------------------------------------------------
// NSJ algorithm: construct NNFCs R-tree
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as data structure is the same
// [out] NNFCs; a hash map for NNFCs, which needs to be constructed in this function
// [out] graph: the typed directed graph, which needs to be constructed in this function
// [reserved/out] NNFC_rtree; this parameter is now reserved; an R-tree for NNFCs, which needs to be constructed in this function

void NSJ_construct_NNFC(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
						fac_hash_map &reflocs, nnfc_hash_map &NNFCs, typed_directed_graph &graph)
{
	// init graph, and deploy facilities
#ifndef ALGO_TESTING_NO_TRACE_GRAPH
	graph.set_testing(&algo_out_file);
#endif
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);

	// deploy facilities and construct related data structure
	graph.deploy_facilities(facs_file_path);

	if (!graph.get_is_facilities_deployed()) // check graph and facilities status
		return;

	// create a top event visitor for NNFCs
	NNFC_visitor top_visitor;
	top_visitor.set_data_structure(&NNFCs); // set data structure

	// compute NNFC of each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// record network location information for each reference location
			reflocs[refloc_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve);

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id, lon, lat, prob); // set id, geo-information and present probability of a reference location
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}
}


//--------------------------------------------------------------------------------
// NSJ algorithm: query optimal candidate based on NNFCs
// remark: the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
// [re] candidate: the id of optimal candidate
// [out] checked_num: the number of candidates actually checked
// [out] rtree_time: the time (ms) of querying R*-tree
// [in] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as data structure is the same
// [in/out] NNFCs; a hash map for NNFCs, COVER_NUM element will be changed
// [in/out] graph: the constructed graph, candidates will be further deployed
// [in] cands_file_path: file path of candidates

candidate NSJ_query(int &checked_num, __int64 &rtree_time, fac_hash_map &reflocs, nnfc_hash_map &NNFCs, typed_directed_graph &graph, const char *cands_file_path)
{
	auto begin_time = std::chrono::high_resolution_clock::now();

	// deploy candidates onto the directed graph and construct candidates R*-tree
	geo_cand_rtree cand_rtree;
	graph.deploy_candidates_rtree(cands_file_path, cand_rtree);

	// construct hash map for <c, R(c)> entries
	cover_hash_map c_rs_map; // sets <c, R(c)>
	nnfc_hash_map::iterator iter_NNFC = NNFCs.begin();
	for (; iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		// range query for candidates which are covered by a NNFC
		std::vector<geo_cand> covered_cands;
		cand_rtree.query(boost::geometry::index::within(iter_NNFC->second.get<NNFC_MBR>()), std::back_inserter(covered_cands));
		iter_NNFC->second.get<COVER_NUM>() = static_cast<int>(covered_cands.size()); // record the number of covered candidates

		// record the set of ids of reference locations covering the candidate
		std::vector<geo_cand>::iterator iter_cand = covered_cands.begin();
		for (; iter_cand != covered_cands.end(); ++iter_cand)
		{
			// compute ERd+ for pair <r,c>
			float distE = boost::geometry::distance(iter_NNFC->second.get<COORDINATE>(), iter_cand->first) * EARTH_RADIUS; // dE(r,c)
			float ERd_u = (iter_NNFC->second.get<NND>() - distE) * iter_NNFC->second.get<PROB>(); // ERd+(F,r,c) = (nnd - dE(r,c)) * Pr(r)
			if (ERd_u <= 0.0f) // impossible for benefit
							   // <: candidate locates inside the 4 corners between NNFC and its MBR
							   // ==: exactly on the NNFC bound
			{
				continue;
			}

			// record for each candidate which NNFCs cover itself and its UB(c)
			cover_hash_map::iterator iter_c_rs_find = c_rs_map.find(iter_cand->second); // find candidate
			if (iter_c_rs_find != c_rs_map.end()) // candidate has already been created
			{
				iter_c_rs_find->second.first += ERd_u; // UB(c) += ERd+(F,r,c)
				iter_c_rs_find->second.second.push_back(iter_NNFC->first); // record id of reference location covering the candidate
			}
			else // need to create a set for the candidate
			{
				c_rs_map[iter_cand->second] = std::make_pair(ERd_u, std::vector<int>()); // ERd_u as the initial UB value
				c_rs_map[iter_cand->second].second.push_back(iter_NNFC->first); // record id of reference location covering the candidate
			}
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	rtree_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();

	cand_max_fibonacci_heap max_heap; // max-heap ordered by UB(c)
	ERD_hash_map ERDs; // actual ERDs for all candidates

	// construct max-heap ordered by UB(c)
	cover_hash_map::iterator iter_c_rs = c_rs_map.begin();
	for (; iter_c_rs != c_rs_map.end(); ++iter_c_rs) // for each candidate
	{
		max_heap.push(candidate(iter_c_rs->first, iter_c_rs->second.first)); // push each candidate into max-heap
		ERDs[iter_c_rs->first] = 0.0f; // initial ERD value is 0

#ifdef ALGO_TESTING
		algo_out_file << "c" << iter_c_rs->first << " is covered by " << iter_c_rs->second.second.size() << " reflocs:" << '\n';
		for (std::vector<int>::iterator iter_r = iter_c_rs->second.second.begin(); iter_r != iter_c_rs->second.second.end(); ++iter_r)
			algo_out_file << "    NNFC(r" << *iter_r << ")" << '\n';
#endif
	}

	// create a top event visitor for NSJ algorithm
	NSJ_visitor top_visitor;
	candidate ol(-1, 0.0f); // initialize optimal location
	top_visitor.set_data_structures(&NNFCs, &ERDs, &ol); // set data structures

	checked_num = 0; // initialize
	while (!max_heap.empty())
	{
		// retrieve the top candidate in max-heap
		candidate top_c = max_heap.top();
#ifdef ALGO_TESTING
		algo_out_file << "ERD(UB): c" << top_c.id << ": " << top_c.ERD << '\n';
#endif

#ifndef ALGO_TESTING_NO_EARLY_STOPPING
		if (top_c.ERD < ol.ERD) // current ol candidate is the optimal location
			break;
#endif
		++checked_num; // need to check more candidates

		// traverse each reference location whose NNFC covers top_c
		std::vector<int>::iterator iter_refloc = c_rs_map[top_c.id].second.begin(),
			iter_refloc_end = c_rs_map[top_c.id].second.end();
		for (; iter_refloc != iter_refloc_end; ++iter_refloc)
		{
			int refloc_id = *iter_refloc;
			if (NNFCs[refloc_id].get<COVER_NUM>() == 0) // all candidates that can affect the reference location have been traversed by a certain candidate prior to this top_c
				continue; // traverse from next reference location

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(reflocs[refloc_id].get<VS_ID>(), 'v'), vertex(reflocs[refloc_id].get<VE_ID>(), 'v'), refloc_v,
				reflocs[refloc_id].get<VS_DIST>(), reflocs[refloc_id].get<VE_DIST>(), &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id, NNFCs[refloc_id].get<PROB>()); // set id and present probability of a reference location
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

#ifdef ALGO_TESTING
		algo_out_file << "current optima: c" << ol.id << ": " << ol.ERD << '\n';
#endif

		max_heap.pop(); // pop the top vertex
	}

	return ol; // the optimal candidate
}
#pragma endregion ALGO_NSJ