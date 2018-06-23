#pragma once

#define EARTH_RADIUS 6378.137
#define EARTH_PERIMETER 2.0 * 3.14 * EARTH_RADIUS

#include <string>
#include <fstream>
#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

typedef boost::property<boost::vertex_name_t, char> VertexTypeProperty;
typedef boost::property<boost::vertex_distance_t, double, VertexTypeProperty> VertexDistanceTypeProperties;
typedef boost::property<boost::vertex_index2_t, int, VertexDistanceTypeProperties> VertexProperties;
typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::property<boost::edge_index_t, int> EdgeProperties;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProperties, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

void init_graph(const std::string &str_vertices_path, std::string &str_edges_path)
{
	static Graph s_g; // Global graph instance.

	static boost::property_map<Graph, boost::vertex_name_t>::type s_vertex_name_map = boost::get(boost::vertex_name, s_g);
	static boost::property_map<Graph, boost::vertex_distance_t>::type s_vertex_distance_map = boost::get(boost::vertex_distance, s_g);
	static boost::property_map<Graph, boost::vertex_index2_t>::type s_vertex_index_map = boost::get(boost::vertex_index2, s_g);
	static boost::property_map<Graph, boost::edge_weight_t>::type s_edge_weight_map = boost::get(boost::edge_weight, s_g);
	static boost::property_map<Graph, boost::edge_index_t>::type s_edge_index_map = boost::get(boost::edge_index, s_g);

	std::set<int> vertices_index;
	std::ifstream ifs_vertices(str_vertices_path.c_str());
	if (!ifs_vertices.fail())	// Vertices data file exists.
	{
		while (!ifs_vertices.bad() && ifs_vertices.good())
		{
			char buf[1024];
			ifs_vertices.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, lon_pos, lat_pos;
			lon_pos = str_buf.find(' ', begin_pos);
			//lat_pos = str_buf.find(' ', lon_pos + 1);

			int iVertex = atoi(str_buf.substr(begin_pos, lon_pos - begin_pos).c_str());
			//double dLon = atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str());
			//double dLat = atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str());

			vertices_index.insert(iVertex);

			/*Vertex v = boost::add_vertex(s_g);
			s_vertex_name_map[v] = 'v';
			s_vertex_distance_map[v] = EARTH_PERIMETER;
			s_vertex_index_map[v] = iVertex;*/
		}
	}

	std::ifstream ifs_edges(str_edges_path.c_str());
	if (!ifs_edges.fail())	// Edges data file exists.
	{
		while (!ifs_edges.bad() && ifs_edges.good())
		{
			char buf[1024];
			ifs_edges.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, vi_pos, vj_pos, dist_pos;
			vi_pos = str_buf.find(' ', begin_pos);
			vj_pos = str_buf.find(' ', vi_pos + 1);
			dist_pos = str_buf.find(' ', vj_pos + 1);

			int iEdge = atoi(str_buf.substr(begin_pos, vi_pos - begin_pos).c_str());
			int iVi = atoi(str_buf.substr(vi_pos + 1, vj_pos - vi_pos - 1).c_str());
			int iVj = atoi(str_buf.substr(vj_pos + 1, dist_pos - vj_pos - 1).c_str());
			double dDist = atof(str_buf.substr(dist_pos + 1, str_buf.size() - dist_pos - 1).c_str());



			Edge e = boost::add_edge();
		}
	}
}

