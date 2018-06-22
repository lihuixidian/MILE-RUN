#pragma once

#include <string> // std::string
#include <fstream> // std::ifstream
#include <chrono> // std::chrono::high_resolution_clock

#include "typed_directed_graph.h"


#pragma region ALGO_Blossom
//================================================================================
// attr_set hash map for Blossom and OTF algorithms

typedef boost::unordered_map < int, // candidate id
	boost::unordered_map < int, // reference location id
	float > > // dist(r, c)
attr_set_hash_map; // attraction set hash map

//--------------------------------------------------------------------------------
// compute merits and return optimal candidate, for both Blossom and OTF
// [re] candidate: the optimal candidate
// [in] attr_dists: reference locations' information, including attractor distances
// [in] attr_set: attraction set

candidate merit(refloc_hash_map &attr_dists, attr_set_hash_map &attr_set)
{
	candidate op;

	// compute ERD for each candidate
	for (attr_set_hash_map::iterator iter_set = attr_set.begin(); iter_set != attr_set.end(); ++iter_set)
	{
		float ERD = 0.0f;
		boost::unordered_map<int, float>::iterator iter_ref = iter_set->second.begin(); // ".second" for each reference locations
		for (; iter_ref != iter_set->second.end(); ++iter_ref)
			ERD += (attr_dists[iter_ref->first].get<RF_NN_DIST>() - iter_ref->second) * attr_dists[iter_ref->first].get<RF_PROB>(); // (a(c) - d(c, v)) * prob
		if (ERD > op.ERD)
		{
			op.id = iter_set->first;
			op.ERD = ERD;
		}
	}
	
	return op;
}

//--------------------------------------------------------------------------------
// top event visitor for Blossom algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			max-fibonacci-heap will never record any ERD information of the reference location, and we view its ERD as 0

class Blossom_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.dist >= iter_r->second.get<RF_NN_DIST>()) // d(c, v) >= a(c), i.e., d(r, v) >= nnd(r)
			return true; // a flag to terminate graph traversal
		else if (v.type == 'c') // d(c, v) < a(c)
		{
			// update attraction set
			attr_set_hash_map::iterator iter_c = ptr_attr_set->find(v.id);
			if (iter_c != ptr_attr_set->end())
				iter_c->second.insert(std::make_pair(iter_r->first, v.dist));
			else
			{
				ptr_attr_set->insert(std::make_pair(v.id, boost::unordered_map<int, float>()));
				iter_c = ptr_attr_set->find(v.id);
				iter_c->second.insert(std::make_pair(iter_r->first, v.dist));
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(attr_set_hash_map *attr_set){ ptr_attr_set = attr_set; }
	void set_refloc(refloc_hash_map::iterator &iter_ref) { iter_r = iter_ref; }

private:
	attr_set_hash_map *ptr_attr_set; // attraction set pointer
	refloc_hash_map::iterator iter_r;
};

//--------------------------------------------------------------------------------
// Blossom algorithm
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [re] candidate: the optimal candidate; default candidate object (id is -1) means candidates are not deployed successfully or opening reference locations file fails
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [in] cands_file_path: file path of candidates

candidate Blossom(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path, const char *cands_file_path)
{
	// init normal and reverse graphs
	typed_directed_graph graph, re_graph;
	graph.set_to_bidirectional(is_to_bidirectional);
	re_graph.set_to_bidirectional(is_to_bidirectional);
	//typed_directed_graph::init_graphs(edges_file_path, graph, re_graph);
	//graph.init_graph(edges_file_path);
	re_graph.init_graph(edges_file_path, true);

	// deploy facilities and reference locations to re_graph
	refloc_hash_map attr_dists;
	re_graph.deploy_facilities(facs_file_path);
	re_graph.deploy_reflocs(reflocs_file_path, &attr_dists);

	// use voronoi diagram to compute "a(c)"s based on reverse graph
	voronoi_hash_map voronoi_vertices;
	re_graph.voronoi_diagram(voronoi_vertices, attr_dists);

	// deploy facilities and candidates to graph
	//graph.deploy_facilities(facs_file_path);
	//graph.deploy_candidates(cands_file_path);

	// create an top event visitor for Blossom algorithm
	//Blossom_visitor visitor;
	attr_set_hash_map attr_set;
/*	visitor.set_data_structures(&attr_set); // set data structures

	// deal with each reference location
	for (refloc_hash_map::iterator iter_ref = attr_dists.begin(); iter_ref != attr_dists.end(); ++iter_ref)
	{
		// insert a reference location vertex to <vs, ve>
		vertex refloc_v(iter_ref->first, 'r');
		std::vector<directed_edge> removed_edges, inserted_edges;
		graph.insert_vertex(vertex(iter_ref->second.get<VS_ID>(), 'v'), vertex(iter_ref->second.get<VE_ID>(), 'v'), refloc_v,
			iter_ref->second.get<VS_DIST>(), iter_ref->second.get<VE_DIST>(), &removed_edges, &inserted_edges,
			false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

		// traverse graph from the reference location by dijkstra algorithm
		visitor.set_refloc(iter_ref); // set data structures
		graph.dijkstra(refloc_v, &visitor);

		// restore the graph by removing the inserted reference location vertex
		graph.restore_graph(refloc_v, removed_edges, inserted_edges);
	}
	*/

	return merit(attr_dists, attr_set); // compute merits and return optimal candidate
}
#pragma endregion ALGO_Blossom


#pragma region ALGO_OTF
//================================================================================
// top event visitor for OTF algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			max-fibonacci-heap will never record any ERD information of the reference location, and we view its ERD as 0

class OTF_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		float lambda = (*ptr_voronoi_vertices)[v].get<CUR_DIST>(); // distance to v's NN facility

		if (v.type == 'r' && v.dist <= lambda)
		{
			// update attraction set
			attr_set_hash_map::iterator iter_c = ptr_attr_set->find(cand_id);
			if (iter_c != ptr_attr_set->end())
				iter_c->second.insert(std::make_pair(v.id, v.dist));
			else
			{
				ptr_attr_set->insert(std::make_pair(cand_id, boost::unordered_map<int, float>()));
				iter_c = ptr_attr_set->find(cand_id);
				iter_c->second.insert(std::make_pair(v.id, v.dist));
			}
		}
		if (v.dist > lambda)
			no_need_extend_out_edges = true; // don't extend out edges of v

		return false; // a flat to continue graph traversal
	}

	void set_data_structures(attr_set_hash_map *attr_set, voronoi_hash_map *voronoi_vertices)
		{ ptr_attr_set = attr_set; ptr_voronoi_vertices = voronoi_vertices; }
	void set_cand_id(int cand_id_in) { cand_id = cand_id_in; }

private:
	attr_set_hash_map *ptr_attr_set; // attraction set pointer
	voronoi_hash_map *ptr_voronoi_vertices; // voronoi vertices pointer
	int cand_id;
};

//--------------------------------------------------------------------------------
// OTF algorithm
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [re] candidate: the optimal candidate; default candidate object (id is -1) means candidates are not deployed successfully or opening reference locations file fails
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [in] cands_file_path: file path of candidates

candidate OTF(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path, const char *cands_file_path)
{
	// init reverse graphs
	typed_directed_graph re_graph;
	re_graph.set_to_bidirectional(is_to_bidirectional);
	re_graph.init_graph(edges_file_path, true);

	// deploy facilities and reference locations to re_graph
	refloc_hash_map attr_dists;
	re_graph.deploy_facilities(facs_file_path);
	re_graph.deploy_reflocs(reflocs_file_path, &attr_dists);

	// use voronoi diagram to compute "a(c)"s based on reverse graph
	voronoi_hash_map voronoi_vertices;
	re_graph.voronoi_diagram(voronoi_vertices, attr_dists);

	// create an top event visitor for Blossom algorithm
	OTF_visitor visitor;
	attr_set_hash_map attr_set;
	visitor.set_data_structures(&attr_set, &voronoi_vertices); // set data structures

	// deal with each candidate
	std::ifstream ifs_cand(cands_file_path); // read candidates file
	if (!ifs_cand.fail()) // candidates file exists
	{
		while (!ifs_cand.bad() && ifs_cand.good())
		{
			char buf[1024];
			ifs_cand.getline(buf, 1024);
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

			// insert a candidate vertex to <vs, ve>
			vertex cand_v(cand_id, 'c');
			std::vector<directed_edge> removed_edges, inserted_edges;
			re_graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), cand_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse re_graph from the candidate by dijkstra algorithm
			visitor.set_cand_id(cand_id); // set the present probability of the reference location
			re_graph.dijkstra(cand_v, &visitor);

			// restore the graph by removing the inserted candidate vertex
			re_graph.restore_graph(cand_v, removed_edges, inserted_edges);
		}
	}

	return merit(attr_dists, attr_set); // compute merits and return optimal candidate
}
#pragma endregion ALGO_OTF