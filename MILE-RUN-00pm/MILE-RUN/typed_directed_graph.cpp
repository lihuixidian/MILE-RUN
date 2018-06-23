#include "stdafx.h"
#include "typed_directed_graph.h"

#include <fstream> // std::ifstream
#include <map> // ordered std::map

#include "float.h" // FLT_MAX constant


//--------------------------------------------------------------------------------
// initialize the directed graph with typed vertices for both normal and reverse versions
// redundant std::set of vertices of in-edges is for edge split and vertex replacing
// remark: the edges file format must be "edge_id vs_id ve_id dist"
// [re] bool: return true if initialization is successful, otherwise, return false
// [in] edges_file_path: file path of edges
// [out] normal_graph: normal graph
// [out] reverse_graph: reverse graph
bool typed_directed_graph::init_graphs(const char *edges_file_path, typed_directed_graph &normal_graph, typed_directed_graph &reverse_graph)
{
	// fault-tolerant check in case re-initialize existing graphs
	if (normal_graph.is_graph_init)
	{
		normal_graph.adjacent_edges.clear();
		normal_graph.is_graph_init = false;
	}
	if (reverse_graph.is_graph_init)
	{
		reverse_graph.adjacent_edges.clear();
		reverse_graph.is_graph_init = false;
	}

	std::ifstream ifs_edges(edges_file_path); // read edges file
	if (!ifs_edges.fail())	// edges file exists
	{
		while (!ifs_edges.bad() && ifs_edges.good())
		{
			char buf[1024];
			ifs_edges.getline(buf, 1024);
			std::string str_buf(buf);

			// the edges file format must be "edge_id vs_id ve_id dist"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_pos = str_buf.find(' ', ve_pos + 1);

			// edge index number is useless
			//int edge_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());

			// init vertex for normal graph or reverse graph
			int vs_id[2], ve_id[2]; // for normal graph [0], and reverse graph [1]
			vs_id[0] = ve_id[1] = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			ve_id[0] = vs_id[1] = atoi(str_buf.substr(ve_pos + 1, dist_pos - ve_pos - 1).c_str());
			float dist = static_cast<float>(atof(str_buf.substr(dist_pos + 1, str_buf.size() - dist_pos - 1).c_str()));

			typed_directed_graph graph[2] = { normal_graph, reverse_graph };
			for (int i = 0; i < 2; ++i)
			{
				// insert out-edge <vs, ve> for vs
				vertex vs(vs_id[i], 'v'); // starting vertex of edge <vs, ve>
				adjacent_edges_hash_map::iterator iter_vs = graph[i].adjacent_edges.find(vs);
				if (iter_vs != graph[i].adjacent_edges.end()) // vs had been inserted
					iter_vs->second.first.insert(out_edge(ve_id[i], 'v', dist)); // an out-edge <vs, ve> of vs
				else // vs is not found
				{
					graph[i].adjacent_edges[vs] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for vs
					graph[i].adjacent_edges[vs].first.insert(out_edge(ve_id[i], 'v', dist)); // an out-edge <vs, ve> of vs
				}

				// insert in-edge <vs, ve> for ve
				vertex ve(ve_id[i], 'v'); // starting vertex of reverse direction edge <ve, vs>
				adjacent_edges_hash_map::iterator iter_ve = graph[i].adjacent_edges.find(ve);
				if (iter_ve != graph[i].adjacent_edges.end()) // reverse starting vertex (ve) exists
					iter_ve->second.second.insert(vertex(vs_id[i], 'v')); // an in-edge <vs, ve> of reverse starting vertex (ve)
				else
				{
					graph[i].adjacent_edges[ve] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for ve
					graph[i].adjacent_edges[ve].second.insert(vertex(vs_id[i], 'v')); // an in-edge <vs, ve> of reverse starting vertex (ve)
				}

				if (graph[i].is_to_bidirectional) // convert undirected graph into bidirectional graph
				{
					// vertices vs and ve must exist now
					graph[i].adjacent_edges[ve].first.insert(out_edge(vs_id[i], 'v', dist)); // an out-edge <ve, vs> of ve
					graph[i].adjacent_edges[vs].second.insert(vertex(ve_id[i], 'v')); // an in-edge <ve, vs> of vs
				}
			}
		}
		return normal_graph.is_graph_init = reverse_graph.is_graph_init = true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// initialize the directed graph with typed vertices;
// redundant std::set of vertices of in-edges is for edge split and vertex replacing
// remark: the edges file format must be "edge_id vs_id ve_id dist"
// [re] bool: return true if initialization is successful, otherwise, return false
// [in] edges_file_path: file path of edges
// [in] reverse_graph: indicate whether constuct a reverse graph (true) or normal graph (false)

bool typed_directed_graph::init_graph(const char *edges_file_path, bool reverse_graph /* false */)
{
	// fault-tolerant check in case re-initialize existing graph
	if (is_graph_init)
	{
		adjacent_edges.clear();
		is_graph_init = false;
	}

	std::ifstream ifs_edges(edges_file_path); // read edges file
	if (!ifs_edges.fail())	// edges file exists
	{
		while (!ifs_edges.bad() && ifs_edges.good())
		{
			char buf[1024];
			ifs_edges.getline(buf, 1024);
			std::string str_buf(buf);

			// the edges file format must be "edge_id vs_id ve_id dist"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_pos = str_buf.find(' ', ve_pos + 1);

			// edge index number is useless
			//int edge_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());

			// init vertex for normal graph or reverse graph
			int vs_id, ve_id;
			if (!reverse_graph) // normal graph
			{
				vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
				ve_id = atoi(str_buf.substr(ve_pos + 1, dist_pos - ve_pos - 1).c_str());
			}
			else // reverse graph
			{
				ve_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
				vs_id = atoi(str_buf.substr(ve_pos + 1, dist_pos - ve_pos - 1).c_str());
			}
			float dist = static_cast<float>(atof(str_buf.substr(dist_pos + 1, str_buf.size() - dist_pos - 1).c_str()));

			// insert out-edge <vs, ve> for vs
			vertex vs(vs_id, 'v'); // starting vertex of edge <vs, ve>
			adjacent_edges_hash_map::iterator iter_vs = adjacent_edges.find(vs);
			if (iter_vs != adjacent_edges.end()) // vs had been inserted
				iter_vs->second.first.insert(out_edge(ve_id, 'v', dist)); // an out-edge <vs, ve> of vs
			else // vs is not found
			{
				adjacent_edges[vs] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for vs
				adjacent_edges[vs].first.insert(out_edge(ve_id, 'v', dist)); // an out-edge <vs, ve> of vs
			}

			// insert in-edge <vs, ve> for ve
			vertex ve(ve_id, 'v'); // starting vertex of reverse direction edge <ve, vs>
			adjacent_edges_hash_map::iterator iter_ve = adjacent_edges.find(ve);
			if (iter_ve != adjacent_edges.end()) // reverse starting vertex (ve) exists
				iter_ve->second.second.insert(vertex(vs_id, 'v')); // an in-edge <vs, ve> of reverse starting vertex (ve)
			else
			{
				adjacent_edges[ve] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for ve
				adjacent_edges[ve].second.insert(vertex(vs_id, 'v')); // an in-edge <vs, ve> of reverse starting vertex (ve)
			}

			if (is_to_bidirectional) // convert undirected graph into bidirectional graph
			{
				// vertices vs and ve must exist now
				adjacent_edges[ve].first.insert(out_edge(vs_id, 'v', dist)); // an out-edge <ve, vs> of ve
				adjacent_edges[vs].second.insert(vertex(ve_id, 'v')); // an in-edge <ve, vs> of vs
			}
		}
		return is_graph_init = true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// deploy facilities onto the directed graph
// remark: facilities are reasonably and strictly forbidden to overlap with each other;
//         the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve"
// [re] bool: return true if deployment is successful, otherwise, return false
// [in] facs_file_path: file path of facilities
// [out/def] hash_map: facilities associated with the endpoints; default is NULL

bool typed_directed_graph::deploy_facilities(const char *facs_file_path, fac_hash_map *hash_map /* = NULL */)
{
	if (!is_graph_init) // must deploy on existing graph
		return false;

	// insert facilities as vertices into graph
	std::ifstream ifs_facs(facs_file_path); // read facilities file
	if (!ifs_facs.fail())	// facilities file exists
	{
		while (!ifs_facs.bad() && ifs_facs.good())
		{
			char buf[1024];
			ifs_facs.getline(buf, 1024);
			std::string str_buf(buf);

			// the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);

			int fac_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, str_buf.size() - dist_ve_pos - 1).c_str()));

			// construct data structure for facilities
			if (hash_map) // used for storing endpoints of facilities
				(*hash_map)[fac_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve);

			// normally insert a vertex to <vs, ve>, all the details are transparent
			insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), vertex(fac_id, 'f'), dist_vs, dist_ve);
		}
		return is_facilities_deployed = true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// deploy reference locations onto the directed graph
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_r dist_r_ve prob lon lat"
// [re] bool: return true if deployment is successful, otherwise, return false
// [in] reflocs_file_path: file path of reference locations
// [out/def] hash_map: reference locations associated with with attractor distances; default is NULL

bool typed_directed_graph::deploy_reflocs(const char *reflocs_file_path, refloc_hash_map *hash_map /* NULL */)
{
	if (!is_graph_init) // must deploy on existing graph
		return false;

	// insert reference locations as vertices into graph
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail())	// reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

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

			// construct data structure for reference locations
			if (hash_map) // used for storing endpoints of reference locations
				(*hash_map)[refloc_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve, prob, -1, FLT_MAX);

			// normally insert a vertex to <vs, ve>, all the details are transparent
			insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), vertex(refloc_id, 'r'), dist_vs, dist_ve);
		}
		return true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// deploy candidates onto the directed graph where facilities have already been deployed
// remark: candidates are reasonably and strictly forbidden to overlap with facilities or each other;
//         the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat";
//	       the max_heap and hash_map must be initialized outside this function;
//		   if max_heap or hash_map is NULL, it means we do not need to construct data structures for candidates
// [re] bool: return true if deployment is successful, otherwise, return false
// [in] cands_file_path: file path of candidates
// [out/def] max_heap: the candidates max-heap ordered by ERD, from which the optimal result will be retrieved; default is NULL
// [out/def] hash_map: the assistant hash map for storing handles of candidates in max-heap and corresponding ERDs; default is NULL

bool typed_directed_graph::deploy_candidates(const char *cands_file_path, cand_max_fibonacci_heap *max_heap /* = NULL */, cand_hash_map *hash_map /* = NULL */)
{
	if (!is_graph_init || !is_facilities_deployed) // must deploy on existing graph with facilities
		return false;

	// insert candidates as vertices into graph
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
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
			// lon and lat are useless for graph

			// construct data structures for candidates
			if (max_heap && hash_map) // the former is used for the optimal result, the latter is an assistant
			{
				// insert a candidate to max-heap ordered by ERD, and initialize the assistant hash-map, where the initial ERD must be 0.0
				(*hash_map)[cand_id] = std::make_pair(max_heap->push(candidate(cand_id, 0.0f)), 0.0f);
			}

			// normally insert a vertex to <vs, ve>, all the details are transparent
			insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), vertex(cand_id, 'c'), dist_vs, dist_ve);
		}
		return is_candidates_deployed = true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// insert vertex to a directed edge
// remark: if removed_edges or inserted_edges is NULL, it means we do not need to record and return related edges;
//         these two output containers must be initialized outside this function;
//		   the format of KEY of split_edges_hash_map must be "vs_id-ve_id"
// [in] vs: the starting vertex of the operated edge
// [in] ve: the ending vertex of the operated edge
// [in] inserted: new vertex to be inserted
// [in] dist_vs_inserted: distance from starting vertex to the vertex to be inserted; 0 if the inserted vertex overlaps vs
// [in] dist_inserted_ve: distance from the vertex to be inserted to ending vertex; 0 if the inserted vertex overlaps ve
// [out/def] removed_edges: record and return the removed edge(s), which would be used for restoring the graph later; default is NULL
// [out/def] inserted_edges: record and return the inserted edge(s), which would be used for restoring the graph later; default is NULL
// [in/def] is_replace_overlap: indicate whether to create a virtual vertex (false/default) or replace (true) overlapped vertex with "inserted';
//								virtual vertex for "inserted" is a (two) zero-distance unidirectional/bidirectional edge(s) between it and the overlapped vertex

void typed_directed_graph::insert_vertex(const vertex &vs, const vertex &ve, const vertex &inserted, float dist_vs_inserted, float dist_inserted_ve,
										 std::vector<directed_edge> *removed_edges /* = NULL */, std::vector<directed_edge> *inserted_edges /* = NULL */,
										 bool is_replace_overlap /* = false */)
{
	// check whether inserted vertex overlaps a certain existing vertex
	if (dist_vs_inserted == 0) // overlap vs
	{
		overlap_vertex(vs, inserted, removed_edges, inserted_edges, is_replace_overlap);
		return;
	}
	else if (dist_inserted_ve == 0) // overlap ve
	{
		overlap_vertex(ve, inserted, removed_edges, inserted_edges, is_replace_overlap);
		return;
	}

	// forward direction, i.e., <vs, ve>
	std::stringstream ss_vs_ve;
	ss_vs_ve << vs.id << "-" << ve.id; // format must be "vs_id-ve_id"
	split_edges_hash_map::iterator iter_vs_ve = split_edges.find(ss_vs_ve.str());
	if (iter_vs_ve != split_edges.end()) // split edge KEY "vs_id-ve_id" is found
	{
		// insert to further split an edge which has been split before
		split_split_edge(iter_vs_ve, vs.id, ve.id, inserted, dist_vs_inserted, dist_inserted_ve, removed_edges, inserted_edges);
		if (inserted.type != 'r') // no need to record the splitting for reference location, as it will be sooner removed
			iter_vs_ve->second.push_back(inserted_vertex(inserted.id, inserted.type, dist_vs_inserted)); // record <vs, ve> is split by only facility and candidate vertex
	}
	else // edge <vs, ve> has not been split yet
	{
		if (split_edge(vs, ve, inserted, dist_vs_inserted, dist_inserted_ve, removed_edges, inserted_edges) // normally insert a vertex to <vs, ve>
			&& inserted.type != 'r') // no need to record the splitting for reference location, as it will be sooner removed
		{
			split_edges[ss_vs_ve.str()] = std::vector<inserted_vertex>(); // bind split vertices to edge <vs, ve>
			split_edges[ss_vs_ve.str()].push_back(inserted_vertex(inserted.id, inserted.type, dist_vs_inserted)); // record <vs, ve> is split by a vertex
		}
	}

	// reverse direction, i.e., <ve, vs>
	std::stringstream ss_ve_vs;
	ss_ve_vs << ve.id << "-" << vs.id; // format must be "ve_id-vs_id"
	split_edges_hash_map::iterator iter_ve_vs = split_edges.find(ss_ve_vs.str());
	if (iter_ve_vs != split_edges.end()) // split edge KEY "ve_id-vs_id" is found
	{
		// insert to further split an edge which has been split before
		split_split_edge(iter_ve_vs, ve.id, vs.id, inserted, dist_inserted_ve, dist_vs_inserted, removed_edges, inserted_edges);
		if (inserted.type != 'r') // no need to record the splitting for reference location, as it will be sooner removed
			iter_ve_vs->second.push_back(inserted_vertex(inserted.id, inserted.type, dist_inserted_ve)); // record <vs, ve> is split by only facility and candidate vertex
	}
	else // edge <ve, vs> has not been split yet
	{
		if (split_edge(ve, vs, inserted, dist_inserted_ve, dist_vs_inserted, removed_edges, inserted_edges) // normally insert a vertex to <ve, vs>
			&& inserted.type != 'r') // no need to record the splitting for reference location, as it will be sooner removed
		{
			split_edges[ss_ve_vs.str()] = std::vector<inserted_vertex>();
			split_edges[ss_ve_vs.str()].push_back(inserted_vertex(inserted.id, inserted.type, dist_inserted_ve)); // record <ve, vs> is split by a facility
		}
	}
}


//--------------------------------------------------------------------------------
// restore graph by inverse operating on "removed edges" and "inserted edges" which are previously recorded
// remark: this function only restore for reference location now, as "split_edges" member is not considered here, which records no reference location splitting;
//		   these two containers "removed edges" and "inserted edges" are not set to empty in this function
// [in] inserted: the inserted vertex
// [in] removed_edges: the removed edge(s) that are previously recorded
// [in] inserted_edges: the inserted edge(s) that are previously recorded

void typed_directed_graph::restore_graph(const vertex &inserted, const std::vector<directed_edge> &removed_edges, const std::vector<directed_edge> &inserted_edges)
{
	// restore the removed edges by re-inserting
	std::vector<directed_edge>::const_iterator iter_removed = removed_edges.cbegin();
	for (; iter_removed != removed_edges.cend(); ++iter_removed)
	{
		// restore out-edge of vs
		adjacent_edges_hash_map::iterator iter_vs = adjacent_edges.find(iter_removed->get<VS>());
		if (iter_vs != adjacent_edges.end()) // starting vertex vs is found
			iter_vs->second.first.insert(out_edge(iter_removed->get<VE>().id, iter_removed->get<VE>().type, iter_removed->get<DIST>())); // restore <vs, ve> for vs
		else // starting vertex vs is not found
		{
			adjacent_edges[iter_removed->get<VS>()] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for vs
			adjacent_edges[iter_removed->get<VS>()].first.insert(out_edge(iter_removed->get<VE>().id, iter_removed->get<VE>().type, iter_removed->get<DIST>())); // restore <vs, ve> for vs
		}

		// restore in-edge of ve
		adjacent_edges_hash_map::iterator iter_ve = adjacent_edges.find(iter_removed->get<VE>());
		if (iter_ve != adjacent_edges.end()) // ending vertex ve is found
			iter_ve->second.second.insert(iter_removed->get<VS>()); // restore <vs, ve> for ve
		else // ending vertex ve is not found
		{
			adjacent_edges[iter_removed->get<VE>()] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for ve
			adjacent_edges[iter_removed->get<VE>()].second.insert(iter_removed->get<VS>()); // restore <vs, ve> for ve
		}
	}

	// remove the inserted vertex, so as the corresponding out-edges and in-edges
	adjacent_edges.erase(inserted);

	// remove the inserted edges except the out-edges and in-edges of "inserted";
	// note that any inserted edge must have an endpoint of vertex "inserted"
	std::vector<directed_edge>::const_iterator iter_inserted = inserted_edges.cbegin();
	for (; iter_inserted != inserted_edges.cend(); ++iter_inserted)
	{
		// remove out-edge
		if (iter_inserted->get<VE>() == inserted) // vs is not "inserted"
		{
			adjacent_edges[iter_inserted->get<VS>()].first.erase(out_edge(inserted.id, inserted.type, iter_inserted->get<DIST>()));
			continue;
		}

		// remove in-edge
		if (iter_inserted->get<VS>() == inserted) // ve is not "inserted"
		{
			adjacent_edges[iter_inserted->get<VE>()].second.erase(inserted);
			continue;
		}
	}
}


//--------------------------------------------------------------------------------
// traverse graph from source vertex by dijkstra algorithm, and many event visitors can be available
// note: if a reference location locates inside a sub-graph without facility, eventually the min-fibonacci-heap will be empty and this function end
// [in] source: the source vertex
// [in/out/def] top_visitor: a pointer to a top event visitor; default is NULL, which means top event won't be visited
// [in/out/def] extend_visitor: a pointer to a extend event visitor; default is NULL, which means extend event won't be visited

void typed_directed_graph::dijkstra(const vertex &source, top_event_visitor *top_visitor /* = NULL */, extend_event_visitor *extend_visitor /* = NULL */)
{
	vertex_min_fibonacci_heap min_heap; // min-heap for target vertices ordered by current shortest path distance from source vertex
	target_hash_map v_h_d_hash_map; // assistant hash map for storing vertices and the corresponding min-heap handles & current distance from source vertex
	definite_hash_map v_d_hash_map; // hash map for storing vertices with definite shortest distance from source vertex

	min_heap.push(target_vertex(source.id, source.type, 0.0f)); // as distance from source vertex to itself is zero, then no need to record its handle
	//if (out_file)
	//	(*out_file) << "source: " << source.type << source.id << '\n';

	while (!min_heap.empty())
	{
		// retrieve the top vertex in min-heap
		target_vertex top_v = min_heap.top();
		v_d_hash_map[top_v] = top_v.dist; // the shortest distance to top vertex from source vertex is definite
		
		//if (out_file)
		//	(*out_file) << "top: " << top_v.type << top_v.id << ", shortest dist = " << top_v.dist << '\n';

		/***** top event *****/
		if (top_visitor) // top event visitor is valid
			if ((*top_visitor)(top_v, v_d_hash_map)) // terminate flag is true
				return;

		if (top_visitor->no_need_extend_out_edges) // don't extend out edges of v
		{
			min_heap.pop(); // pop the top vertex
			continue; // go next top vertex in min-heap
		}

		// push out-vertices of top target vertex into min-heap
		std::set<out_edge>::iterator iter_out_edge = adjacent_edges[top_v].first.begin();
		for (; iter_out_edge != adjacent_edges[top_v].first.end(); ++iter_out_edge) // for each out-vertex of top vertex
		{
			// if the out-vertex is already with definite distance, probe next out-vertex
			if (v_d_hash_map.find(iter_out_edge->ve) != v_d_hash_map.end())
			{
				//if (out_file)
				//	(*out_file) << "definite: " << top_v.type << top_v.id << " -> " << iter_out_edge->ve.type << iter_out_edge->ve.id << '\n';

				continue;
			}

			/***** extend event *****/
			if (extend_visitor) // extend event visitor is valid
				if ((*extend_visitor)(top_v, *iter_out_edge)) // terminate flag is true
					return;

			// compute the distance from source to the out-vertex, considering whether it has been traversed via another path
			float new_dist = top_v.dist + iter_out_edge->dist; // dist(s, top_v) + dist(top_v, v)
			target_hash_map::iterator iter_vertex = v_h_d_hash_map.find(iter_out_edge->ve);
			if (iter_vertex != v_h_d_hash_map.end()) // the out-vertex has already been traversed via another path
			{
				if (new_dist < iter_vertex->second.second) // if dist(s, top_v) + dist(top_v, v) < dist(s, v) (current distance), update the current shorter distance
				{
					// here "increase" means to raise the rank (its original meaning is to increase ordered value)
					min_heap.increase(iter_vertex->second.first, target_vertex(iter_out_edge->ve.id, iter_out_edge->ve.type, new_dist));
					iter_vertex->second.second = new_dist; // update the current shorter distance in assistant hash map
					
					//if (out_file)
					//	(*out_file) << "update adjacent: " << top_v.type << top_v.id << " -> " << iter_out_edge->ve.type << iter_out_edge->ve.id
					//		<< ", dist(" << source.type << source.id << ", " << iter_out_edge->ve.type << iter_out_edge->ve.id << ") = " << new_dist << '\n';
				}
				//else
				//{
				//	if (out_file)
				//		(*out_file) << "no update adjacent: " << top_v.type << top_v.id << " -> " << iter_out_edge->ve.type << iter_out_edge->ve.id
				//			<< ", dist(" << source.type << source.id << ", " << iter_out_edge->ve.type << iter_out_edge->ve.id << ") = " << new_dist << '\n';
				//}
			}
			else // first time encounter the out-vertex
			{
				v_h_d_hash_map[iter_out_edge->ve] = std::make_pair( // record min-heap handle and current distance from source vertex
					min_heap.push(target_vertex(iter_out_edge->ve.id, iter_out_edge->ve.type, new_dist)), // push a out-vertex of top vertex into min-heap
					new_dist); // current shortest path distance
				
				//if (out_file)
				//	(*out_file) << "new adjacent: " << top_v.type << top_v.id << " -> " << iter_out_edge->ve.type << iter_out_edge->ve.id
				//		<< ", dist(" << source.type << source.id << ", " << iter_out_edge->ve.type << iter_out_edge->ve.id << ") = " << new_dist << '\n';
			}
		}

		min_heap.pop(); // pop the top vertex
	}
}


//--------------------------------------------------------------------------------
// explore reverse graph, and obtain voronoi vertices
// [in/out] voronoi_vertices: voronoi vertices
// [in/out] attr_dists: reference locations' information, including attractor distances; this hash map must have been initialized already

void typed_directed_graph::voronoi_diagram(voronoi_hash_map &voronoi_vertices, refloc_hash_map &attr_dists)
{
	vertex_min_fibonacci_heap min_heap; // min-heap for voronoi vertices ordered by current optimal distance to NN facility
	target_hash_map v_h_d_hash_map; // assistant hash map for storing vertices and the corresponding min-heap handles & current distance to NN facility

	// init voronoi vertices and min-heap for them
	voronoi_vertices.clear();
	for (adjacent_edges_hash_map::iterator iter_v_adj = adjacent_edges.begin(); iter_v_adj != adjacent_edges.end(); ++iter_v_adj)
	{
		if (iter_v_adj->first.type == 'f')
		{
			voronoi_vertices[iter_v_adj->first] = boost::make_tuple(0.0f, // CUR_DIST
																	iter_v_adj->first.id, // NN_ID
																	'f', // NN_TYPE
																	false); // MARK
			v_h_d_hash_map[iter_v_adj->first] = std::make_pair( // record min-heap handle and current distance to NN facility
				min_heap.push(target_vertex(iter_v_adj->first.id, 'f', 0.0f)), // push vertex into min-heap
				0.0f); // current NN distance
		}
		else // normal vertices and reference locations
			voronoi_vertices[iter_v_adj->first] = boost::make_tuple(FLT_MAX, -1, ' ', false); // <CUR_DIST, NN_ID, NN_TYPE, MARK>
	}

	while (!min_heap.empty())
	{
		// retrieve the top vertex in min-heap
		target_vertex top_v = min_heap.top();
		voronoi_vertices[top_v].get<MARK>() = true; // now, top_v is optimal

		//if (out_file)
		//	(*out_file) << "top: " << top_v.type << top_v.id << ", NN dist = " << top_v.dist << '\n';

		// record attractor distances (without NN type, then facility is default type)
		if (top_v.type == 'r')
		{
			attr_dists[top_v.id].get<RF_NN_F>() = voronoi_vertices[top_v].get<NN_ID>(); // record the NN facility
			attr_dists[top_v.id].get<RF_NN_DIST>() = voronoi_vertices[top_v].get<CUR_DIST>(); // record the distance
		}

		// handle out-edges of top_v
		std::set<out_edge>::iterator iter_out_edge = adjacent_edges[top_v].first.begin(); // ".first" indicates the set of out-edges, i.e., std::set<out_edge>
		for (; iter_out_edge != adjacent_edges[top_v].first.end(); ++iter_out_edge) // for each out-vertex of top vertex
		{
			if (voronoi_vertices[iter_out_edge->ve].get<MARK>() == false) // ve of out-edge is not marked yet
			{
				float delta = voronoi_vertices[top_v].get<CUR_DIST>() // distance from top_v
					+ iter_out_edge->dist; // out_edge length
				if (voronoi_vertices[iter_out_edge->ve].get<CUR_DIST>() == FLT_MAX) // out vertex is first time encountered
				{
					voronoi_vertices[iter_out_edge->ve].get<CUR_DIST>() = delta;
					voronoi_vertices[iter_out_edge->ve].get<NN_ID>() = voronoi_vertices[top_v].get<NN_ID>();
					voronoi_vertices[iter_out_edge->ve].get<NN_TYPE>() = voronoi_vertices[top_v].get<NN_TYPE>();

					v_h_d_hash_map[iter_out_edge->ve] = std::make_pair( // record min-heap handle and current distance to NN facility
						min_heap.push(target_vertex(iter_out_edge->ve.id, iter_out_edge->ve.type, delta)), // push a out-vertex of top vertex into min-heap
						delta); // current NN distance
				}
				if (voronoi_vertices[iter_out_edge->ve].get<CUR_DIST>() < FLT_MAX // not first time
					&& delta < voronoi_vertices[iter_out_edge->ve].get<CUR_DIST>()) // shorter path is found
				{
					voronoi_vertices[iter_out_edge->ve].get<CUR_DIST>() = delta;
					voronoi_vertices[iter_out_edge->ve].get<NN_ID>() = voronoi_vertices[top_v].get<NN_ID>();
					voronoi_vertices[iter_out_edge->ve].get<NN_TYPE>() = voronoi_vertices[top_v].get<NN_TYPE>();

					target_hash_map::iterator iter_vertex = v_h_d_hash_map.find(iter_out_edge->ve);
					if (iter_vertex != v_h_d_hash_map.end()) // the out-vertex has already been explored
					{
						// here "increase" means to raise the rank (its original meaning is to increase ordered value)
						min_heap.increase(iter_vertex->second.first, target_vertex(iter_out_edge->ve.id, iter_out_edge->ve.type, delta));
						iter_vertex->second.second = delta; // update the current distance in assistant hash map
					}
				}
			}
		}
		
		min_heap.pop(); // pop the top vertex
	}
}


//--------------------------------------------------------------------------------
// deploy candidates onto the directed graph where facilities have already been deployed, also construct a R*-tree for candidates
// remark: similar to "deploy_candidates" function;
//		   the better implementation is to inherit class "typed_directed_graph" and put this function in a newly derived class
// [re] bool: return true if deployment is successful, otherwise, return false
// [in] cands_file_path: file path of candidates
// [out] cand_rtree: the candidates R*-tree

bool typed_directed_graph::deploy_candidates_rtree(const char *cands_file_path, geo_cand_rtree &cand_rtree)
{
	if (!is_graph_init || !is_facilities_deployed) // must deploy on existing graph with facilities
		return false;

	// insert candidates as vertices into graph, and construct R*-tree at the same time
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, lon_pos - dist_ve_pos - 1).c_str()));
			// differing from function "deploy_candidates", lon and lat are useful now
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct candidates R*-tree
			cand_rtree.insert(std::make_pair(geo_point(lon, lat), cand_id));

			// normally insert a vertex to <vs, ve>, all the details are transparent
			insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), vertex(cand_id, 'c'), dist_vs, dist_ve);
		}
		return is_candidates_deployed = true;
	}
	return false;
}


//--------------------------------------------------------------------------------
// overlap an existing vertex with a new inserted vertex, in an overlap way or replace way; the latter is not used now
// remark: if removed_edges or inserted_edges is NULL, it means we do not need to record and return related edges;
//         these two output containers must be initialized outside this function;
//		   if in the overlap way, removed_edges won't be changed, as no edge is removed
// [in] overlapped: the existing vertex to be overlapped
// [in] inserted: new vertex to be inserted
// [out/def] removed_edges: record and return the removed edge(s), which would be used for restoring the graph later; default is NULL
// [out/def] inserted_edges: record and return the inserted edge(s), which would be used for restoring the graph later; default is NULL
// [in/def] is_replace_overlap: indicate whether to create a virtual vertex (false/default) or replace (true) overlapped vertex with "inserted';
//								as the replace way is not used now, this parameter must be "false";
//								virtual vertex for "inserted" is a (two) zero-distance unidirectional/bidirectional edge(s) between it and the overlapped vertex;
//								specifically, unidirectional edge is for reference location ('r'), while bidirectional edge is for facility and candidate ('f' & 'c')

void typed_directed_graph::overlap_vertex(const vertex &overlapped, const vertex &inserted, std::vector<directed_edge> *removed_edges /* = NULL */,
										  std::vector<directed_edge> *inserted_edges /* = NULL */, bool is_replace_overlap /* = false */)
{
	if (is_replace_overlap) // in a replace way, which is not used now
		replace_vertex(overlapped, inserted, removed_edges, inserted_edges);
	else // in an overlap way, where no existing edge is removed
	{
		// first create a virtual vertex for "inserted"
		adjacent_edges[inserted] = std::make_pair(std::set<out_edge>(), std::set<vertex>());

		// edge <inserted, overlapped> is inserted for all kinds of vertex ('f', 'c' and 'r')
		adjacent_edges[inserted].first.insert(out_edge(overlapped.id, overlapped.type, 0.0f)); // new zero-distance out-edge <inserted, overlapped> for virtual vertex
		if (inserted_edges)
			inserted_edges->push_back(boost::make_tuple(inserted, overlapped, 0.0f)); // record the newly inserted out-edge <inserted, overlapped>
		adjacent_edges[overlapped].second.insert(inserted); // insert an in-edge <inserted, overlapped> for overlapped

		// an <inserted, overlapped> out-edge is enough for reference location, while reverse edge <overlapped, inserted> is meaningless
		if (inserted.type == 'r')
			return;

		// below, reverse edge <overlapped, inserted> is inserted for 'f' and 'c'
		adjacent_edges[overlapped].first.insert(out_edge(inserted.id, inserted.type, 0.0f)); // new zero-distance out-edge <overlapped, inserted> for "overlapped"
		if (inserted_edges)
			inserted_edges->push_back(boost::make_tuple(overlapped, inserted, 0.0f)); // record the newly inserted out-edge <overlapped, inserted>
		adjacent_edges[inserted].second.insert(overlapped); // insert an in-edge <overlapped, inserted> for virtual vertex
	}
}


//--------------------------------------------------------------------------------
// split a directed edge with an inserted vertex
// remark: if removed_edges or inserted_edges is NULL, it means we do not need to record and return related edges;
//         these two output containers must be initialized outside this function
// [re] bool: return true if <vs, ve> exists such that can be split, otherwise, return false
// [in] vs: the starting vertex of the edge to be split
// [in] ve: the ending vertex of the edge to be split
// [in] inserted: new vertex to be inserted
// [in] dist_vs_inserted: non-zero distance from starting vertex to the vertex to be inserted
// [in] dist_inserted_ve: non-zero distance from the vertex to be inserted to ending vertex
// [out/def] removed_edges: record and return the removed edge(s), which would be used for restoring the graph later; default is NULL
// [out/def] inserted_edges: record and return the inserted edge(s), which would be used for restoring the graph later; default is NULL

bool typed_directed_graph::split_edge(const vertex &vs, const vertex &ve, const vertex &inserted, float dist_vs_inserted, float dist_inserted_ve,
									  std::vector<directed_edge> *removed_edges /* = NULL */, std::vector<directed_edge> *inserted_edges /* = NULL */)
{
	adjacent_edges_hash_map::iterator iter_vs = adjacent_edges.find(vs);
	if (iter_vs != adjacent_edges.end()) // starting vertex vs is found
	{
		// record & remove the split edge
		bool vs_ve_exist = false; // indicate whether <vs, ve> exists (true) or not (false)
		std::set<out_edge>::iterator iter_out_edge = iter_vs->second.first.begin(); // function "find" is unavailable as "out-edge.dist" is unknown
		for (; iter_out_edge != iter_vs->second.first.end(); ++iter_out_edge) // scan out-edges of vs until <vs, ve> is found or to the end
		{
			if (iter_out_edge->ve == ve) // ve is found, thus directed edge <vs, ve> is found
			{
				vs_ve_exist = true; // <vs, ve> exists

				// remove out-edge <vs, ve> of vs
				if (removed_edges)
					removed_edges->push_back(boost::make_tuple(vs, ve, iter_out_edge->dist)); // first record the removed (split) edge
				iter_vs->second.first.erase(iter_out_edge); // then remove the edge from out-edges

				// in-edge <vs, ve> of ve must exist, then remove it accordingly
				adjacent_edges[ve].second.erase(vs);

				break;
			}
		}

		// insert two successive edges <vs, inserted> and <inserted, ve> to replace the split edge <vs, ve>
		if (vs_ve_exist)
		{
			// insert out-edge <vs, inserted> of vs
			iter_vs->second.first.insert(out_edge(inserted.id, inserted.type, dist_vs_inserted)); // insert new edge <vs, inserted>
			if (inserted_edges)
				inserted_edges->push_back(boost::make_tuple(vs, inserted, dist_vs_inserted)); // record the newly inserted edge

			// insert in-edge <vs, inserted> of vertex "inserted"
			adjacent_edges_hash_map::iterator iter_inserted = adjacent_edges.find(inserted); // maybe "inserted" exists resulted from reverse direction
			if (iter_inserted != adjacent_edges.end()) // vertex "inserted" is found (from the reverese direction)
				iter_inserted->second.second.insert(vs); // in-edge <vs, inserted> of vertex "inserted"
			else
			{
				adjacent_edges[inserted] = std::make_pair(std::set<out_edge>(), std::set<vertex>()); // bind adjacent edges for "inserted"
				adjacent_edges[inserted].second.insert(vs); // in-edge <vs, inserted> of vertex "inserted"
			}

			// insert out-edge <inserted, ve> of "inserted"
			adjacent_edges[inserted].first.insert(out_edge(ve.id, ve.type, dist_inserted_ve)); // insert new edge <inserted, ve>
			if (inserted_edges)
				inserted_edges->push_back(boost::make_tuple(inserted, ve, dist_inserted_ve)); // record the newly inserted edge

			// insert in-edge <inserted, ve> of ve
			adjacent_edges[ve].second.insert(inserted); // in-edge <inserted, ve> of ve
		}

		return vs_ve_exist;
	}
	return false;
}


//--------------------------------------------------------------------------------
// insert a new vertex to further split an edge which has been split before
// [in] iter: the iterator pointing to an edge that has been split before
// [in] vs_id: the starting vertex id of the operated edge
// [in] ve_id: the ending vertex id of the operated edge
// [in] inserted: new vertex to be inserted
// [in] dist_vs: non-zero distance from starting vertex to the vertex to be inserted
// [in] dist_ve: non-zero distance from the vertex to be inserted to ending vertex
// [out/def] removed_edges: record and return the removed edge(s), which would be used for restoring the graph later; default is NULL
// [out/def] inserted_edges: record and return the inserted edge(s), which would be used for restoring the graph later; default is NULL

void typed_directed_graph::split_split_edge(const split_edges_hash_map::iterator &iter, int vs_id, int ve_id, const vertex &inserted,
											float dist_vs, float dist_ve, std::vector<directed_edge> *removed_edges /* = NULL */,
											std::vector<directed_edge> *inserted_edges /* = NULL */)
{
	// order all vertices as scale marks by the distance from vs to themselves;
	// as facilities and candidates cannot overlap, scale marks must not overlap, either
	std::map<float, vertex> ordered; // <dist(vs, vertex), vertex>
	ordered[0.0f] = vertex(vs_id, 'v'); // vs scale mark
	ordered[dist_vs + dist_ve] = vertex(ve_id, 'v'); // ve scale mark
	std::vector<inserted_vertex>::iterator iter_v = iter->second.begin();
	for (; iter_v != iter->second.end(); ++iter_v) // iterate all inserted vertices on the edge
		ordered[iter_v->dist_vs] = vertex(iter_v->id, iter_v->type);

	// obtain the interval for new vertex "inserted"
	vertex pre_v, next_v; // the previous and next scale marks of "inserted"
	float dist_pre, dist_next; // distance from vs to pre_v and next_v
	for (std::map<float, vertex>::iterator iter_order = ordered.begin(); iter_order != ordered.end(); ++iter_order)
	{
		if (dist_vs == iter_order->first) // overlap occurs (must be reference location vertex)
		{
			// only reference location can overlap facility or candidate, hence to create virtual vertex (is_replace_overlap = false)
			overlap_vertex(iter_order->second, inserted, removed_edges, inserted_edges);

			return;
		}
		else if (iter_order->first < dist_vs)
		{
			pre_v = iter_order->second;
			dist_pre = iter_order->first;
		}
		else if (dist_vs < iter_order->first)
		{
			next_v = iter_order->second;
			dist_next = iter_order->first;
			break;
		}
	}

	// normally split edge <pre, next> by vertex "inserted"
	split_edge(pre_v, next_v, inserted, dist_vs - dist_pre, dist_next - dist_vs, removed_edges, inserted_edges);
}


//--------------------------------------------------------------------------------
// replace an existing vertex with a new inserted vertex and the corresponding edge(s); this function is not used now
// remark: if removed_edges or inserted_edges is NULL, it means we do not need to record and return related edges;
//         these two output containers must be initialized outside this function
// [in] replaced: the existing vertex to be replaced
// [in] inserted: new vertex to be inserted
// [out/def] removed_edges: record and return the removed edge(s), which would be used for restoring the graph later; default is NULL
// [out/def] inserted_edges: record and return the inserted edge(s), which would be used for restoring the graph later; default is NULL

void typed_directed_graph::replace_vertex(const vertex &replaced, const vertex &inserted, std::vector<directed_edge> *removed_edges /* = NULL */,
										  std::vector<directed_edge> *inserted_edges /* = NULL */)
{
	// first insert the new vertex "inserted", which will replace the vertex "replaced"
	adjacent_edges[inserted] = std::make_pair(std::set<out_edge>(), std::set<vertex>());

	// replace out-edges and in-edges of vertex "replaced" with "inserted"
	adjacent_edges_hash_map::iterator iter_replaced = adjacent_edges.find(replaced);
	if (iter_replaced != adjacent_edges.end()) // vertex "replaced" exists
	{
		// replace out-edges
		std::set<out_edge>::iterator iter_out_edge = iter_replaced->second.first.begin();
		for (; iter_out_edge != iter_replaced->second.first.end(); ++iter_out_edge) // deal with each out-edge of vertex "replaced"
		{
			// record the removed out-edge <replaced, out-vertex> of vertex "replaced"
			if (removed_edges)
				removed_edges->push_back(boost::make_tuple(replaced, iter_out_edge->ve, iter_out_edge->dist)); // record an out-edge to be removed
			adjacent_edges[inserted].first.insert(out_edge(iter_out_edge->ve.id, iter_out_edge->ve.type, iter_out_edge->dist)); // new out-edge <inserted, out-vertex> for replacing
			if (inserted_edges)
				inserted_edges->push_back(boost::make_tuple(inserted, iter_out_edge->ve, iter_out_edge->dist)); // record the newly inserted edge

			// replace in-edge of out-vertex
			adjacent_edges[iter_out_edge->ve].second.erase(replaced);
			adjacent_edges[iter_out_edge->ve].second.insert(inserted);
		}

		// replace in-edges
		std::set<vertex>::iterator iter_in_edge = iter_replaced->second.second.begin();
		for (; iter_in_edge != iter_replaced->second.second.end(); ++iter_in_edge) // deal with each in-edge of vertex "replaced"
		{
			// insert in-edge of "inserted" for replacing
			vertex vs_in(iter_in_edge->id, iter_in_edge->type); // out-vertex of "replaced", i.e., starting vertex of in-edge
			adjacent_edges[inserted].second.insert(vs_in);

			// record the removed out-edge <vs_in, replaced> of out-vertex of "replaced" (i.e., vs_in)
			adjacent_edges_hash_map::iterator iter_vs_in = adjacent_edges.find(vs_in);
			if (iter_vs_in != adjacent_edges.end()) // out-vertex of "replaced" (vs_in) exists
			{
				std::set<out_edge>::iterator iter_vs_in_out_edge = iter_vs_in->second.first.begin(); // function "find" is unavailable as "out-edge.dist" is unknown
				for (; iter_vs_in_out_edge != iter_vs_in->second.first.end(); ++iter_vs_in_out_edge) // scan out-edges of vs_in until <vs_in, replaced> is found or to the end
				{
					if (iter_vs_in_out_edge->ve == replaced) // "replaced" is found, thus directed edge <vs_in, replaced> is found
					{
						// remove out-edge <vs_in, replaced> of vs_in
						if (removed_edges)
							removed_edges->push_back(boost::make_tuple(vs_in, replaced, iter_vs_in_out_edge->dist)); // first record the removed edge
						iter_vs_in->second.first.insert(out_edge(inserted.id, inserted.type, iter_vs_in_out_edge->dist)); // insert <vs_in, inserted> for replacing
						if (inserted_edges)
							inserted_edges->push_back(boost::make_tuple(vs_in, inserted, iter_vs_in_out_edge->dist)); // record the newly inserted edge
						iter_vs_in->second.first.erase(iter_vs_in_out_edge); // then remove the edge from out-edges of vs_in

						break;
					}
				}
			}
		}

		// finally, remove the vertex "replaced", and so are the corresponding out-edges and in-edges
		adjacent_edges.erase(iter_replaced);
	}
}