#pragma once

#include <utility>
#include <map>
#include <set>
#include <fstream>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>

//================================================================================
// Notations:
// Starting and ending vertices are the both endpoints of a directed edge.
// Source and target vertices are the vertices of a path.
// Each vertex consist of an index number and a type, which can be set as follows:
//		v - ordinary vertex;
//		f - facility vertex;
//		c - candidate vertex;
//		r - reference location vertex.

//--------------------------------------------------------------------------------
// vertex associated with type

struct vertex
{
	int id; // vertex index number
	char type; // vertex type

	vertex()
		: id(0)
		, type('v') {}

	vertex(int id_in, char type_in)
		: id(id_in)
		, type(type_in) {}

	vertex(const vertex &rhs)
		: vertex(rhs.id, rhs.type) {}

	friend bool operator<(const vertex &lhs, const vertex &rhs)
	{
		if (lhs.type == rhs.type) // first consider type
			return lhs.id < rhs.id; // compare index number
		else
			return lhs.type < rhs.type; // then compare type
	}

	friend bool operator==(const vertex &lhs, const vertex &rhs)
	{
		return lhs.id == rhs.id && lhs.type < rhs.type;
	}
};

//--------------------------------------------------------------------------------
// directed out-edge with ending vertex and edge weight (distance) but without starting vertex

struct out_edge
{
	vertex ve; // ending vertex
	float dist; // edge weight (distance)

	out_edge()
		: ve()
		, dist(0.0f) {}

	out_edge(int id_in, char type_in, float dist_in)
		: ve(id_in, type_in)
		, dist(dist_in) {}

	friend bool operator<(const out_edge &lhs, const out_edge &rhs)
	{
		return lhs.ve < rhs.ve;
	}
};

//--------------------------------------------------------------------------------
// a tuple that represents a directed edge

typedef boost::tuple < vertex, // starting vertex
	vertex, // ending vertex
	float > // edge weight (distance)
	directed_edge;

//--------------------------------------------------------------------------------
// target vertex associated with the distance from source vertex

struct target_vertex
	: public vertex
{
	float dist; // the distance from source vertex

	target_vertex(int id_in, char type_in, float dist_in)
		: vertex(id_in, type_in)
		, dist(dist_in) {}

	// remark: the same vertex validation is not implemented
	friend bool operator>(const target_vertex &lhs, const target_vertex &rhs)
	{
		return lhs.dist > rhs.dist;
	};
};

//--------------------------------------------------------------------------------
// a mutable fibonacci min-heap for ordering target vertices without indication of source vertex

typedef boost::heap::fibonacci_heap < target_vertex, // target vertex associated with the distance from source vertex to it
									  boost::heap::compare < std::greater < target_vertex > >, // shift max-heap to min-heap
									  boost::heap::stable <false> > // mutable
									  min_fibonacci_heap;

typedef min_fibonacci_heap::handle_type vertex_handle;

//--------------------------------------------------------------------------------
// typed directed graph

class typed_directed_graph
{
public:
	typed_directed_graph(){}
	virtual ~typed_directed_graph(){}

public:
	std::map<vertex, std::set<out_edge>> out_edges; // vertices (key) each of which binds the corresponding out-edges (value)
	min_fibonacci_heap min_heap;

public:
	// remark: the edges file format must be "edge_id vs_id ve_id dist"
	void init_graph(const char *edges_file_path, // file path of edges
					bool is_bidirectional = false) // must explicitly indicate to convert undirected graph into bidirectional graph
	{
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

				// int edge_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str()); // edge index number is useless
				int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
				int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_pos - ve_pos - 1).c_str());
				float dist = static_cast<float>(atof(str_buf.substr(dist_pos + 1, str_buf.size() - dist_pos - 1).c_str()));

				vertex vs(vs_id, 'v');
				std::map<vertex, std::set<out_edge>>::iterator iter = out_edges.find(vs);
				if (iter != out_edges.end()) // vs had been inserted
					(*iter).second.insert(out_edge(ve_id, 'v', dist));
				else
				{
					out_edges[vs] = std::set<out_edge>(); // bind a set for out-edges
					out_edges[vs].insert(out_edge(ve_id, 'v', dist));
				}

				if (is_bidirectional) // convert undirected graph into bidirectional graph
				{
					vertex re_vs(ve_id, 'v'); // starting vertex of reverse direction edge <ve, vs>
					std::map<vertex, std::set<out_edge>>::iterator iter = out_edges.find(re_vs);
					if (iter != out_edges.end()) // reverse vs (original ve) had been inserted
						(*iter).second.insert(out_edge(vs_id, 'v', dist)); // original vs as reverse ve
					else
					{
						out_edges[re_vs] = std::set<out_edge>();
						out_edges[re_vs].insert(out_edge(vs_id, 'v', dist)); // original vs as reverse ve
					}
				}
			}
		}

		/* test
		std::map<vertex, std::set<out_edge>>::const_iterator citer = out_edges.cbegin();
		for (; citer != out_edges.cend(); ++citer)
		{
			std::set<out_edge>::const_iterator cit = citer->second.cbegin();
			for (; cit != citer->second.cend(); ++cit)
			{
				std::cout << citer->first.id << ' ' << cit->ve.id << ' ' << cit->dist << '\n';
			}
		}
		std::cout << std::endl;*/
	}

	// return the number of split edge(s)
	int split_edge(std::vector<directed_edge> &split_edges, // record the edge(s) being split
				   vertex vs, // the starting vertex of the edge being split
				   vertex ve, // the ending vertex of the edge being split
				   vertex inserted, // vertex being inserted
				   float dist_vs_inserted, // distance from starting vertex to the vertex being inserted
				   float dist_inserted_ve) // distance from the vertex being inserted to ending vertex 
	{
		int re_num = 0; // init the number of edge(s) being split

		// forward direction
		bool forward_direction_exist = false; // check whether forward direction edge exists or not
		std::map<vertex, std::set<out_edge>>::iterator iter_vs = out_edges.find(vs);
		if (iter_vs != out_edges.end()) // vs is found
		{
			// record & remove the split edge
			std::set<out_edge>::iterator iter_edge = iter_vs->second.begin();
			for (; iter_edge != iter_vs->second.end(); ++iter_edge)
			{
				if (iter_edge->ve == ve) // ve is found, thus edge <vs, ve> is found
				{
					split_edges.push_back(boost::make_tuple(vs, ve, iter_edge->dist)); // record the forward direction edge being split
					iter_vs->second.erase(iter_edge); // remove the edge being split

					++re_num; // increase the number of split edge(s) by 1
					forward_direction_exist = true;
					break;
				}
			}
			
			// insert two successive edges to replace the split edge
			iter_vs->second.insert(out_edge(inserted.id, inserted.type, dist_vs_inserted)); // new edge <vs, inserted>
			out_edges[inserted] = std::set<out_edge>();
			out_edges[inserted].insert(out_edge(ve.id, ve.type, dist_inserted_ve)); // new edge <inserted, ve>
		}

		// reverse direction
		std::map<vertex, std::set<out_edge>>::iterator iter_ve = out_edges.find(ve);
		if (iter_ve != out_edges.end()) // ve is found
		{
			// record & remove the split edge
			std::set<out_edge>::iterator iter_edge = iter_ve->second.begin();
			for (; iter_edge != iter_ve->second.end(); ++iter_edge)
			{
				if (iter_edge->ve == vs) // vs is found, thus edge <ve, vs> is found
				{
					split_edges.push_back(boost::make_tuple(ve, vs, iter_edge->dist)); // record the reverse direction edge being split
					iter_ve->second.erase(iter_edge); // remove the edge being split

					++re_num; // increase the number of split edge(s) by 1
					break;
				}
			}

			// insert two successive edges to replace the split edge
			iter_ve->second.insert(out_edge(inserted.id, inserted.type, dist_inserted_ve)); // new edge <ve, inserted>
			if (!forward_direction_exist)
				out_edges[inserted] = std::set<out_edge>();
			out_edges[inserted].insert(out_edge(vs.id, vs.type, dist_vs_inserted)); // new edge <inserted, vs>
		}

		return re_num;
	}

	void test()
	{
		vertex_handle v1 = min_heap.push(target_vertex(1, 'v', 5.5f));
		vertex_handle v2 = min_heap.push(target_vertex(2, 'v', 3.5f));
		vertex_handle v3 = min_heap.push(target_vertex(3, 'v', 1.5f));
		vertex_handle v4 = min_heap.push(target_vertex(4, 'v', 4.5f));
		vertex_handle v5 = min_heap.push(target_vertex(5, 'v', 2.5f));

		bool flag = true;
		while (!min_heap.empty())
		{
			std::cout << min_heap.top().type << min_heap.top().id << ' ' << min_heap.top().dist << '\n';
			min_heap.pop();
			if (flag)
			{
				min_heap.increase(v4, target_vertex(4, 'v', 2.0f));
				flag = false;
			}
			std::cout << min_heap.size() << '\n';
		}
		std::cout << std::endl;
	};
};

