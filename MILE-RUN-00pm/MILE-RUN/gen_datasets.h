#pragma once

#include <string> // str::string
#include <utility> // std::pair
#include <vector> // std::vector
#include <fstream> // std::ifstream, std::ofstream
#include <iostream> // std::cout
#include <sstream> // std::stringstream
#include <ctime> // std::time

#include <boost/unordered_map.hpp> // for consistency reasons, we use boost::unordered_map instead of std::unordered_map
#include <boost/unordered_set.hpp> // for consistency reasons, we use boost::unordered_set instead of std::unordered_set
#include <boost/functional/hash.hpp> // boost::hash_combine
#include <boost/tuple/tuple.hpp> // boost::tuple
#include <boost/random/uniform_int_distribution.hpp> // boost::random::uniform_int_distribution
#include <boost/random/discrete_distribution.hpp> // boost::random::discrete_distribution
#include <boost/random/mersenne_twister.hpp> // boost::random::mt11213b

#include "geo_defines.h" // R*-tree, EARTH_RADIUS


//================================================================================
// Notations:
//	This header file implements a generator which integrates facilities onto graph edges according to coordinates,
//  and uniformly generates candidates and clients with reference locations.
//  To be specific, facilities cannot overlap one another, and candidates cannot overlap facilities and each other.
//  For simplicity, we DON'T parse the probabilities of different counts of reference locations, instead, we use constant probability array.

//--------------------------------------------------------------------------------
// testing for output procedures

//#define GEN_TESTING // only for "gen_datasets.h"
#ifdef _WIN32
	#include <windows.h> // for (_WIN32) MessageBoxA function
#endif

//--------------------------------------------------------------------------------
// facility or candidate location

struct fac_or_cand_loc
{
	int vs_id; // vs id
	int ve_id; // ve id
	float vs_dist; // distance from vs

	fac_or_cand_loc()
		: vs_id(-1)
		, ve_id(-1)
		, vs_dist(0.0f) {}

	fac_or_cand_loc(int vs_in, int ve_in, float dist_in)
		: vs_id(vs_in)
		, ve_id(ve_in)
		, vs_dist(dist_in) {}

	fac_or_cand_loc(const fac_or_cand_loc &rhs)
	{
		if (this == &rhs)
			return;
		vs_id = rhs.vs_id;
		ve_id = rhs.ve_id;
		vs_dist = rhs.vs_dist;
	}

	bool operator==(const fac_or_cand_loc &rhs) const
	{
		return vs_id == rhs.vs_id && ve_id == rhs.ve_id && vs_dist == rhs.vs_dist;
	}

	// boost::hash is implemented by calling the function hash_value
	friend std::size_t hash_value(const fac_or_cand_loc& loc)
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, loc.vs_id);
		boost::hash_combine(seed, loc.ve_id);
		boost::hash_combine(seed, loc.vs_dist);
		return seed;
	}
};

//--------------------------------------------------------------------------------
// a hash map that stores every vertex with longitude and latitude
// a hash map that stores every edge, which KEY is edge id (i.e., index for random chosen)

typedef boost::unordered_map < int, // vertex id
	geo_point > // <longitude, latitude>
raw_vertices;

typedef boost::unordered_map < int, // edge id (index)
	boost::tuple < int, // vs id
	int, // ve id
	float > > // edge distance, which is evalulated by invoking "boost::geometry::distance"
raw_edges;

enum raw_edges_element { RE_VS_ID = 0, RE_VE_ID, RE_DIST };

//--------------------------------------------------------------------------------
// poi category
// BAR, HOSPITAL, PO, PARK, SCHOOL for "gen_datasets.h"
// STATION, CAFE, LOGISTICS, BANK, SCHOOL for "convert_datasets.h"

enum poi_category { BAR = 0, STATION = 0, HOSPITAL = 1, CAFE = 1, PO = 2, LOGISTICS = 2, PARK = 3, BANK = 3, SCHOOL, CATEGORY };

//--------------------------------------------------------------------------------
// random generator

static boost::random::mt11213b random_gen(static_cast<uint32_t>(std::time(0)));

//--------------------------------------------------------------------------------
// load vertices and edges datasets; also construct an R*-tree to organize edges
// remark: the vertices file format must be "vertex_id lon lat"; the edges file format must be "edge_id vs_id ve_id dist_vs";
//		   the NULL-checking of "geo_range_file_path" is only for integrating real dataset
// [re] bool: return true if loading is successful, otherwise, return false
// [out] vertices: a hash map for vertices
// [out] edges: a hash map for edges
// [out] rtree: an R*-tree for edges
// [in] geo_range_file_path: file path of geo-range
// [in] vertices_file_path: file path of vertices
// [in] edges_file_path: file path of edges
// [in] check_directionality: indicate whether check edge bidirectionality (true) or not (false)

bool load_vertices_and_edges(raw_vertices &vertices, raw_edges &edges, geo_edge_rtree &rtree, const char *geo_range_file_path, const char *vertices_file_path,
	const char *edges_file_path, bool check_directionality)
{
	float min_lon = 180.0f, max_lon = -180.0f, min_lat = 90.0f, max_lat = -90.0f; // value for geo-range

	// first load vertices
	std::ifstream ifs_vertices(vertices_file_path); // read vertices file
	if (!ifs_vertices.fail())	// vertices file exists
	{
		while (!ifs_vertices.bad() && ifs_vertices.good())
		{
			char buf[1024];
			ifs_vertices.getline(buf, 1024);
			std::string str_buf(buf);

			// the vertices file format must be "vertex_id lon lat"
			std::string::size_type begin_pos = 0, lon_pos, lat_pos;
			lon_pos = str_buf.find(' ', begin_pos);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int vertex_id = atoi(str_buf.substr(begin_pos, lon_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// comparation for geo-range
			if (lon < min_lon)
				min_lon = lon;
			if (lon > max_lon)
				max_lon = lon;
			if (lat < min_lat)
				min_lat = lat;
			if (lat > max_lat)
				max_lat = lat;

			vertices[vertex_id] = geo_point(lon, lat); // insert a new vertex
#ifdef GEN_TESTING
			std::cout << "v: " << vertex_id << '\n';
#endif
		}
	}
	else
		return false; // loading fails, as reading vertices fails

	// output geo-range
	if (geo_range_file_path != NULL) // this checking is only for integrating real dataset
	{
		std::ofstream ofs_range(geo_range_file_path, std::ofstream::out | std::ofstream::trunc); // create geo_range file
		if (!ofs_range.fail())	// create geo_range file successfully
		{
			ofs_range << "min lon: " << min_lon << '\n'
				<< "max lon: " << max_lon << '\n'
				<< "min lat: " << min_lat << '\n'
				<< "max lat: " << max_lat << std::endl;
		}
		else
			return false; // creating fails
	}

	boost::unordered_set<std::string> reverse_edges; // used to check edge bidirectionality

	// then load edges
	std::ifstream ifs_edges(edges_file_path); // read edges file
	if (!ifs_edges.fail())	// edges file exists
	{
		while (!ifs_edges.bad() && ifs_edges.good())
		{
			char buf[1024];
			ifs_edges.getline(buf, 1024);
			std::string str_buf(buf);

			// the edges file format must be "edge_id vs_id ve_id dist_vs"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_pos = str_buf.find(' ', ve_pos + 1);

			int edge_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_pos - ve_pos - 1).c_str());

			// we use Euclidean distance between vertices instead of given road distance,
			// as candidates and reference locations cannot have a relatively accurate location on the unknown road
			//float dist = static_cast<float>(atof(str_buf.substr(dist_pos + 1, str_buf.size() - dist_pos - 1).c_str()));
			float dist = boost::geometry::distance(vertices[vs_id], vertices[ve_id]) * EARTH_RADIUS;
			
			edges[edge_id] = boost::make_tuple(vs_id, ve_id, dist); // insert a new edge
#ifdef GEN_TESTING
			std::cout << "e: " << edge_id << '\n';
#endif

			// check edge directionality
			if (check_directionality)
			{
				std::stringstream ss_vs_ve; // current edge
				ss_vs_ve << vs_id << "-" << ve_id; // format must be "vs_id-ve_id"
				if (reverse_edges.find(ss_vs_ve.str()) != reverse_edges.end()) // ve->vs has already been loaded, hence, now vs->ve makes vs-ve a bidirectional edge
				{
#ifdef GEN_TESTING					
					std::cout << "edge \"" << ss_vs_ve.str() << "\" is bidirectional" << '\n';
#ifdef _WIN32
					MessageBoxA(NULL, ss_vs_ve.str().c_str(), "check directionality", MB_OK);
#endif
#endif
				}

				std::stringstream ss_ve_vs; // reverse edge
				ss_ve_vs << ve_id << "-" << vs_id;
				reverse_edges.insert(ss_ve_vs.str());
			}

			// construct edge R*-tree
			rtree.insert(std::make_pair(geo_segment(vertices[vs_id], vertices[ve_id]), edge_id));
		}
	}
	else
		return false; // loading fails, as reading edges fails

	return true;
}

//--------------------------------------------------------------------------------
// load raw pois dataset, and pick pois on demand
// remark: the pois file format must be "poi_category lon lat";
//		   array pois_count must has CATEGORY elements, each of which must be initialized to 0 before invoking this function
// [re] bool: return true if loading is successful, otherwise, return false
// [in] pois_file_path: file path of pois
// [out] pois_count: the counts of each category of pois
// [out] bars: a vector for bars raw data
// [out] hospitals: a vector for hospitals raw data
// [out] pos: a vector for hospitals raw data
// [out] parks: a vector for hospitals raw data
// [out] schools: a vector for hospitals raw data

bool load_raw_pois(const char *pois_file_path, int pois_count[], std::vector<geo_point> &bars, std::vector<geo_point> &hospitals, std::vector<geo_point> &pos,
				   std::vector<geo_point> &parks, std::vector<geo_point> &schools)
{
	std::ifstream ifs_pois(pois_file_path); // read pois file
	if (!ifs_pois.fail())	// pois file exists
	{
		while (!ifs_pois.bad() && ifs_pois.good())
		{
			char buf[1024];
			ifs_pois.getline(buf, 1024);
			std::string str_buf(buf);

			// the pois file format must be "poi_category lon lat"
			std::string::size_type begin_pos = 0, lon_pos, lat_pos;
			lon_pos = str_buf.find(' ', begin_pos);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			std::string category = str_buf.substr(begin_pos, lon_pos - begin_pos);
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			if (category == "bar") // category id is 4
			{
				bars.push_back(geo_point(lon, lat));
				++pois_count[BAR];
			}
			else if (category == "hospital") // category id is 30
			{
				hospitals.push_back(geo_point(lon, lat));
				++pois_count[HOSPITAL];
			}
			else if (category == "po") // category id is 43
			{
				pos.push_back(geo_point(lon, lat));
				++pois_count[PO];
			}
			else if (category == "park") // category id is 40
			{
				parks.push_back(geo_point(lon, lat));
				++pois_count[PARK];
			}
			else if (category == "school") // category id is 50
			{
				schools.push_back(geo_point(lon, lat));
				++pois_count[SCHOOL];
			}
			else if (category == "sea") // out of useful categories
				break;
#ifdef GEN_TESTING
			std::cout << category << '\n';
#endif
		}
	}
	else
		return false; // loading fails, as reading pois fails

#ifdef GEN_TESTING
	for (int i = 0; i < CATEGORY; ++i)
		std::cout << "category " << i << ": " << pois_count[i] << '\n';
#endif

	return true;
}

//--------------------------------------------------------------------------------
// given a point, calculate its vertical point on a straight line, which is represented by two points
// [in] vs: vs point
// [in] ve: ve point
// [in] poi: a given poi point
// [in] edge_dist: the distance of the edge, which has been mutilplied EARTH_RADIUS
// [out] vs_dist: the distance from vs to vertical point of poi
// [out] ve_dist: the distance from vertical point of poi to ve
// [out/def] poi_on_road: the poi point on road; default is NULL, which means no on road poi point is needed

void calculate_vertical_point(const geo_point &vs, const geo_point &ve, const geo_point &poi, float edge_dist, float &vs_dist, float &ve_dist, geo_point *poi_on_road = NULL)
{
	if (vs.get<LON>() == ve.get<LON>() && vs.get<LAT>() == ve.get<LAT>()) // vs overlaps ve
	{
		vs_dist = ve_dist = 0.0f;
		if (poi_on_road != NULL)
		{
			poi_on_road->set<LON>(vs.get<LON>());
			poi_on_road->set<LAT>(vs.get<LAT>());
		}
	}
	else if (vs.get<LON>() == ve.get<LON>()) // the straight line is parallel with Y axis
	{
		if (vs.get<LAT>() < ve.get<LAT>()) // vs.lat < ve.lat
		{
			if (poi.get<LAT>() <= vs.get<LAT>()) // poi.lat <= vs.lat < ve.lat, poi_vp must locate vs
			{
				vs_dist = 0.0f;
				ve_dist = edge_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>());
					poi_on_road->set<LAT>(vs.get<LAT>());
				}
			}
			else if (ve.get<LAT>() <= poi.get<LAT>()) // vs.lat < ve.lat <= poi.lat, poi_vp must locate ve
			{
				vs_dist = edge_dist;
				ve_dist = 0.0f;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(ve.get<LON>());
					poi_on_road->set<LAT>(ve.get<LAT>());
				}
			}
			else // vs.lat < poi.lat < ve.lat, need to calculate d(vs, poi_vp)
			{
				vs_dist = (poi.get<LAT>() - vs.get<LAT>()) / (ve.get<LAT>() - vs.get<LAT>()) * edge_dist;
				ve_dist = edge_dist - vs_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>()); // vs.lon == ve.lon
					poi_on_road->set<LAT>(poi.get<LAT>()); // use poi.lat
				}
			}
		}
		else // ve.lat < vs.lat
		{
			if (poi.get<LAT>() <= ve.get<LAT>()) // poi.lat <= ve.lat < vs.lat, poi_vp must locate ve
			{
				vs_dist = edge_dist;
				ve_dist = 0.0f;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(ve.get<LON>());
					poi_on_road->set<LAT>(ve.get<LAT>());
				}
			}
			else if (vs.get<LAT>() <= poi.get<LAT>()) // ve.lat < vs.lat <= poi.lat, poi_vp must locate vs
			{
				vs_dist = 0.0f;
				ve_dist = edge_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>());
					poi_on_road->set<LAT>(vs.get<LAT>());
				}
			}
			else // ve.lat < poi.lat < vs.lat, need to calculate d(vs, poi_vp)
			{
				vs_dist = (vs.get<LAT>() - poi.get<LAT>()) / (vs.get<LAT>() - ve.get<LAT>()) * edge_dist;
				ve_dist = edge_dist - vs_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>()); // vs.lon == ve.lon
					poi_on_road->set<LAT>(poi.get<LAT>()); // use poi.lat
				}
			}
		}
	}
	else if (vs.get<LAT>() == ve.get<LAT>()) // the straight line is parallel with X axis
	{
		if (vs.get<LON>() < ve.get<LON>()) // vs.lon < ve.lon
		{
			if (poi.get<LON>() <= vs.get<LON>()) // poi.lon <= vs.lon < ve.lon, poi_vp must locate vs
			{
				vs_dist = 0.0f;
				ve_dist = edge_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>());
					poi_on_road->set<LAT>(vs.get<LAT>());
				}
			}
			else if (ve.get<LON>() <= poi.get<LON>()) // vs.lon < ve.lon <= poi.lon, poi_vp must locate ve
			{
				vs_dist = edge_dist;
				ve_dist = 0.0f;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(ve.get<LON>());
					poi_on_road->set<LAT>(ve.get<LAT>());
				}
			}
			else // vs.lon < poi.lon < ve.lon, need to calculate d(vs, poi_vp)
			{
				vs_dist = (poi.get<LON>() - vs.get<LON>()) / (ve.get<LON>() - vs.get<LON>()) * edge_dist;
				ve_dist = edge_dist - vs_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(poi.get<LON>()); // use poi.lon
					poi_on_road->set<LAT>(vs.get<LAT>()); // vs.lat == ve.lat
				}
			}
		}
		else // ve.lon < vs.lon
		{
			if (poi.get<LON>() <= ve.get<LON>()) // poi.lon <= ve.lon < vs.lon, poi_vp must locate ve
			{
				vs_dist = edge_dist;
				ve_dist = 0.0f;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(ve.get<LON>());
					poi_on_road->set<LAT>(ve.get<LAT>());
				}
			}
			else if (vs.get<LON>() <= poi.get<LON>()) // ve.lon < vs.lon <= poi.lon, poi_vp must locate vs
			{
				vs_dist = 0.0f;
				ve_dist = edge_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(vs.get<LON>());
					poi_on_road->set<LAT>(vs.get<LAT>());
				}
			}
			else // ve.lon < poi.lon < vs.lon, need to calculate d(vs, poi_vp)
			{
				vs_dist = (vs.get<LON>() - poi.get<LON>()) / (vs.get<LON>() - ve.get<LON>()) * edge_dist;
				ve_dist = edge_dist - vs_dist;
				if (poi_on_road != NULL)
				{
					poi_on_road->set<LON>(poi.get<LON>()); // use poi.lon
					poi_on_road->set<LAT>(vs.get<LAT>()); // vs.lat == ve.lat
				}
			}
		}
	}
	else // an arbitrary point
	{
		float A, B, C; // parameters of standard straight line formula Ax + By + C = 0, i.e., (y2-y1)x + (x1-x2)y + [(x2-x1)y1 - x1(y2-y1)] = 0, where p1 = vs, p2 = ve
		A = ve.get<LAT>() - vs.get<LAT>(); // y2 - y1
		B = vs.get<LON>() - ve.get<LON>(); // x1 - x2
		C = (ve.get<LON>() - vs.get<LON>()) * vs.get<LAT>() - vs.get<LON>() * (ve.get<LAT>() - vs.get<LAT>()); // (x2 - x1)y1 - x1(y2 - y1)
		float square_sum = A * A + B * B; // A^2 + B^2
		
		float dist_vp_vs, dist_vp_ve; // distances from vertical point to vs and ve, respectively
		if (poi_on_road != NULL)
		{
			poi_on_road->set<LON>((B * B * poi.get<LON>() - A * B * poi.get<LAT>() - A * C) / square_sum); // (B * B * px - A * B * py - A * C) / (A^2 + B^2)
			poi_on_road->set<LAT>((A * A * poi.get<LAT>() - A * B * poi.get<LON>() - B * C) / square_sum); // (A * A * py - A * B * px - B * C) / (A^2 + B^2)

			// calculate temporary distances from vertical point to vs and ve
			dist_vp_vs = boost::geometry::distance(*poi_on_road, vs) * EARTH_RADIUS;
			dist_vp_ve = boost::geometry::distance(*poi_on_road, ve) * EARTH_RADIUS;
		}
		else // use an alternative of poi_on_road
		{
			geo_point poi_vp; // the vertical point of poi on the line, which is an alternative if poi_on_road is NULL

			poi_vp.set<LON>((B * B * poi.get<LON>() - A * B * poi.get<LAT>() - A * C) / square_sum); // (B * B * px - A * B * py - A * C) / (A^2 + B^2)
			poi_vp.set<LAT>((A * A * poi.get<LAT>() - A * B * poi.get<LON>() - B * C) / square_sum); // (A * A * py - A * B * px - B * C) / (A^2 + B^2)

			// calculate temporary distances from vertical point to vs and ve
			dist_vp_vs = boost::geometry::distance(poi_vp, vs) * EARTH_RADIUS;
			dist_vp_ve = boost::geometry::distance(poi_vp, ve) * EARTH_RADIUS;
		}
		
		if (dist_vp_vs > dist_vp_ve) // point_on_road is near ve
		{
			if (dist_vp_vs < edge_dist) // point_on_road locate on edge vs-ve
			{
				vs_dist = dist_vp_vs;
				ve_dist = edge_dist - vs_dist;
			}
			else // dist_vp_vs >= edge_dist, vs - ve - point_on_road, point_on_road must locate ve
			{
				vs_dist = edge_dist;
				ve_dist = 0.0f;
			}
		}
		else if (dist_vp_vs < dist_vp_ve) // point_on_road is near vs
		{
			if (dist_vp_ve < edge_dist) // point_on_road locate on edge vs-ve
			{
				vs_dist = dist_vp_vs;
				ve_dist = edge_dist - vs_dist;
			}
			else // dist_vp_ve >= edge_dist, point_on_road - vs - ve, point_on_road must locate vs
			{
				vs_dist = 0.0f;
				ve_dist = edge_dist;
			}
		}
		else // dist_vp_vs == dist_vp_ve, point_on_road at the midpoint of vs-ve
			vs_dist = ve_dist = edge_dist / 2.0f;
	}
}

//--------------------------------------------------------------------------------
// integrate raw pois dataset into the data we need, and output results to files
// remark: the pois file format must be "poi_category lon lat";
//		   array pois_count must has CATEGORY elements, each of which must be initialized to 0 before invoking this function
// [re] bool: return true if integration is successful, otherwise, return false
// [in] edges: a hash map for edges, MUST be used as a const peremeter
// [in] rtree: an R*-tree for edges
// [in] pois: a vector for pois raw data
// [in] pois_file_path: output file path of pois
// [out] gen_facs: append each new generated facility to this hash set, used for uniqueness check
// [in/def] is_with_coordinates: whether the coordinates of facilities are appended at the end of each line

bool integrate_pois(raw_edges &edges, const geo_edge_rtree &rtree, const std::vector<geo_point> &pois, const char *pois_file_path,
					boost::unordered_set<fac_or_cand_loc> &gen_facs, bool is_with_coordinates = false)
{
	std::ofstream ofs_pois(pois_file_path, std::ofstream::out | std::ofstream::trunc); // create pois file
	if (!ofs_pois.fail())	// create pois file successfully
	{
		std::vector<geo_point>::const_iterator iter_poi = pois.cbegin();
		for (int index = 0; iter_poi != pois.cend(); ++iter_poi, ++index)
		{
			// query the nn edge of each poi
			std::vector<geo_edge> nn_edge;
			rtree.query(boost::geometry::index::nearest(*iter_poi, 1), std::back_inserter(nn_edge)); // 1 means top nn (i.e., for knn, k = 1)

			std::vector<geo_edge>::iterator iter_nn = nn_edge.begin();
			if (iter_nn != nn_edge.end()) // nn edge is found
			{
				float edge_dist = edges[iter_nn->second].get<RE_DIST>();
				int vs_id = edges[iter_nn->second].get<RE_VS_ID>();
				int ve_id = edges[iter_nn->second].get<RE_VE_ID>();
				float vs_dist, ve_dist;
				calculate_vertical_point(iter_nn->first.first, iter_nn->first.second, *iter_poi, edge_dist, vs_dist, ve_dist); // deal with arbitrary point
								
				// output the integrated poi, using hash set "gen_facs" for uniqueness check
				std::pair<boost::unordered_set<fac_or_cand_loc>::iterator, bool> return_value; // the 2nd component of return value of "insert" shows whether insertion takes place or not
				if (vs_dist == 0.0f) // overlaps vs
				{
					return_value = gen_facs.insert(fac_or_cand_loc(vs_id, vs_id, 0.0f));
					if (return_value.second) // insertion takes place
					{
						if (is_with_coordinates)
							ofs_pois << index << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << iter_poi->get<LON>() << ' ' << iter_poi->get<LAT>() << '\n';
						else
							ofs_pois << index << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << '\n';
					}
				}
				else if (ve_dist == 0.0f) // overlaps ve
				{
					return_value = gen_facs.insert(fac_or_cand_loc(ve_id, ve_id, 0.0f));
					if (return_value.second) // insertion takes place
					{
						if (is_with_coordinates)
							ofs_pois << index << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << iter_poi->get<LON>() << ' ' << iter_poi->get<LAT>() << '\n';
						else
							ofs_pois << index << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << '\n';
					}
				}
				else // between vs and ve
				{
					return_value = gen_facs.insert(fac_or_cand_loc(vs_id, ve_id, vs_dist));
					if (return_value.second) // insertion takes place
					{
						if (is_with_coordinates)
							ofs_pois << index << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << ve_dist << ' ' << iter_poi->get<LON>() << ' ' << iter_poi->get<LAT>() << '\n';
						else
							ofs_pois << index << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << ve_dist << '\n';
					}
				}
#ifdef GEN_TESTING
				std::cout << "gen poi: " << index << '\n';
#endif
			}			
		}

		ofs_pois.flush();
	}
	else
		return false; // creating fails

	return true;
}

//--------------------------------------------------------------------------------
// pick pois from integrated pois dataset, and output results to files
// [re] bool: return true if picking is successful, otherwise, return false
// [in] pois_file_path: input file path of integrated pois
// [in] picked_pois_file_path: output file path of picked pois
// [in] count: the count of picked pois

bool pick_pois(const char *pois_file_path, const char *picked_pois_file_path, const int count)
{
	// read integrated pois
	std::vector<std::string> pois;
	std::ifstream ifs(pois_file_path);
	if (!ifs.fail())	// integrated pois file exists
	{
		while (!ifs.bad() && ifs.good())
		{
			char buf[1024];
			ifs.getline(buf, 1024);
			if (std::string(buf).size() != 0) // avoid last line, i.e., "\n"
				pois.push_back(buf);
		}
	}

	// randomly pick pois
	std::ofstream ofs(picked_pois_file_path, std::ofstream::out | std::ofstream::trunc); // create picked pois file
	if (!ofs.fail())	// create pois file successfully
	{
		int selected = 0; // selected count
		while (selected < count)
		{
			boost::random::uniform_int_distribution<> random_index(0, static_cast<int>(pois.size()) - 1); // poi index begins from 0
			int selected_index = random_index(random_gen);
			
			// output picked poi
			ofs << pois[selected_index] << '\n';
#ifdef GEN_TESTING
			std::cout << "pick poi: " << pois[selected_index] << '\n';
#endif

			// delete the picked poi from all pois
			pois[selected_index] = pois.back();
			pois.pop_back();

			++selected;
		}
		ofs.flush();
	}
	else
		return false; // creating fails

	return true;
}

//--------------------------------------------------------------------------------
// output generated candidates
// [re] bool: return true if generation is successful, otherwise, return false
// [in] vertices: a hash map for vertices, MUST be used as a const peremeter
// [in] edges: a hash map for edges, MUST be used as a const peremeter
// [in] gen_cands_file_path: output file path of generated candidates
// [in] cands_count: the count of candidates to be generated
// [in/out] gen_facs: append each new generated candidate to this hash set, used for uniqueness check;
//					  NOTE facilities cannot overlap each other, and candidates cannot overlap facilities and each other

bool gen_cands(raw_vertices &vertices, raw_edges &edges, boost::unordered_set<fac_or_cand_loc> &gen_facs, const char *gen_cands_file_path, int cands_count)
{
	boost::random::uniform_int_distribution<> random_index(0, static_cast<int>(edges.size()) - 1); // edge index begins from 0
	boost::random::uniform_int_distribution<> random_percent(0, 1000); // distance percent from vs to candidate, accuracy is xx.x%

	std::ofstream ofs_cands(gen_cands_file_path, std::ofstream::out | std::ofstream::trunc); // create gen_cands file
	if (!ofs_cands.fail())	// create gen_cands file successfully
	{
		for (int i = 0; i < cands_count; ++i)
		{
			// pick an edge according to random edge index
			unsigned index = static_cast<unsigned>(random_index(random_gen));
			int vs_id = edges[index].get<RE_VS_ID>();
			int ve_id = edges[index].get<RE_VE_ID>();
			float dist = edges[index].get<RE_DIST>();

			// generate a candidate, and insert as a facility such that to avoid overlap
			float percent = -1.0f, vs_dist = -1.0f; // initial value is impossible
			std::pair<boost::unordered_set<fac_or_cand_loc>::iterator, bool> return_value; // the 2nd component of return value of "insert" shows whether insertion takes place or not
			do
			{
				percent = random_percent(random_gen) / 1000.0f; // vs->c / c->ve, accuracy is xx.x%
				vs_dist = dist * percent;
				if (vs_dist == 0.0f) // candidate overlaps vs
					return_value = gen_facs.insert(fac_or_cand_loc(vs_id, vs_id, 0.0f));
				else if (vs_dist == dist) // candidate overlaps ve
					return_value = gen_facs.insert(fac_or_cand_loc(ve_id, ve_id, 0.0f));
				else
					return_value = gen_facs.insert(fac_or_cand_loc(vs_id, ve_id, vs_dist));
			} while (!return_value.second); // insertion fails

			// output the generated candidate
			if (vs_dist == 0.0f) // candidate overlaps vs
				ofs_cands << i << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << vertices[vs_id].get<LON>() << ' ' << vertices[vs_id].get<LAT>() << '\n';
			else if (vs_dist == dist) // candidate overlaps ve
				ofs_cands << i << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << vertices[ve_id].get<LON>() << ' ' << vertices[ve_id].get<LAT>() << '\n';
			else
			{
				// c.x = "edge.length.x" * percent(vs->c) + vs.x = (ve.x - vs.x) * percent(vs->c) + vs.x; c.y is similar
				float lon = (vertices[ve_id].get<LON>() - vertices[vs_id].get<LON>()) * percent + vertices[vs_id].get<LON>();
				float lat = (vertices[ve_id].get<LAT>() - vertices[vs_id].get<LAT>()) * percent + vertices[vs_id].get<LAT>();
				ofs_cands << i << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << dist - vs_dist << ' ' << lon << ' ' << lat << '\n';
			}
#ifdef GEN_TESTING
			std::cout << "gen cand: " << i << '\n';
#endif
		}

		ofs_cands.flush();
	}
	else
		return false; // creating fails

	return true;
}

//--------------------------------------------------------------------------------
// output generated reference locations
// remark: for simplicity, we DON'T parse the input "ref_loc_points", instead, we use constant arrays in this function
// [re] bool: return true if generation is successful, otherwise, return false
// [in] vertices: a hash map for vertices, MUST be used as a const peremeter
// [in] edges: a hash map for edges, MUST be used as a const peremeter
// [in] gen_random_ref_locs_file_path: output file path of generated reference locations with randomly chosen present probabilities
// [in] gen_avg_ref_locs_file_path: output file path of generated reference locations with average present probabilities
// [in] clients_count: the count of clients to be generated
// [in] ref_loc_points: the probabilities of different counts of reference locations

bool gen_ref_locs(raw_vertices &vertices, raw_edges &edges, const char *gen_random_ref_locs_file_path, const char *gen_avg_ref_locs_file_path,
				  int clients_count, const char *ref_loc_points)
{
	// parse the probabilities of different counts of reference locations; the probabilities string format must be "count count#1 prob#1 ..."
	// here, we DON'T parse the input "ref_loc_points"
	// SHOULD parse "ref_loc_points" to obtain the following two arrays
	int counts[] = { 1, 2, 3, 4, 5, 6 };
	float probabilities[] = { 0.05f, 0.2f, 0.35f, 0.25f, 0.1f, 0.05f };

	boost::random::discrete_distribution<> random_count(probabilities); // reference locations count with different probabilities
	boost::random::uniform_int_distribution<> random_index(0, static_cast<int>(edges.size()) - 1); // edge index begins from 0
	boost::random::uniform_int_distribution<> random_percent(0, 1000); // distance percent from vs to candidate, accuracy is xx.x%

	std::ofstream ofs_random(gen_random_ref_locs_file_path, std::ofstream::out | std::ofstream::trunc); // create reference locations file with random present probabilities
	std::ofstream ofs_avg(gen_avg_ref_locs_file_path, std::ofstream::out | std::ofstream::trunc); // create reference locations file with average present probabilities
	if (!ofs_random.fail() && !ofs_avg.fail())	// create both files successfully
	{
		int id = 0; // global index of reference locations of all clients
		for (int i = 0; i < clients_count; ++i)
		{
			int count = counts[random_count(random_gen)]; // random count of reference locations of the client
			int remnant_prob = 1000; // the remnant present probability, accuracy is xx.x%
			for (int j = 0; j < count; ++j)
			{
				// compute random present probability
				boost::random::uniform_int_distribution<> random_present_prob(0, remnant_prob); // present probability, accuracy is xx.x%
				int gen_prob = random_present_prob(random_gen);
				float present_prob = gen_prob / 1000.0f; // accuracy is xx.x%
				remnant_prob -= gen_prob;

				// pick an edge according to random edge index
				unsigned index = static_cast<unsigned>(random_index(random_gen));
				int vs_id = edges[index].get<RE_VS_ID>();
				int ve_id = edges[index].get<RE_VE_ID>();
				float dist = edges[index].get<RE_DIST>();

				// generate a reference location
				float percent = random_percent(random_gen) / 1000.0f;
				float vs_dist = dist * percent;

				// output the generated reference location
				if (vs_dist == 0.0f) // reference location overlaps vs
				{
					ofs_random << id << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << present_prob << ' '
						<< vertices[vs_id].get<LON>() << ' ' << vertices[vs_id].get<LAT>() << '\n';
					ofs_avg << id << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << 1.0f / count << ' '
						<< vertices[vs_id].get<LON>() << ' ' << vertices[vs_id].get<LAT>() << '\n';
				}
				else if (vs_dist == dist) // reference location overlaps ve
				{
					ofs_random << id << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << present_prob << ' '
						<< vertices[ve_id].get<LON>() << ' ' << vertices[ve_id].get<LAT>() << '\n';
					ofs_avg << id << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << 1.0f / count << ' '
						<< vertices[ve_id].get<LON>() << ' ' << vertices[ve_id].get<LAT>() << '\n';
				}
				else
				{
					// c.x = "edge.length.x" * percent(vs->c) + vs.x = (ve.x - vs.x) * percent(vs->c) + vs.x; c.y is similar
					float lon = (vertices[ve_id].get<LON>() - vertices[vs_id].get<LON>()) * percent + vertices[vs_id].get<LON>();
					float lat = (vertices[ve_id].get<LAT>() - vertices[vs_id].get<LAT>()) * percent + vertices[vs_id].get<LAT>();
					ofs_random << id << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << dist - vs_dist << ' ' << present_prob << ' '
						<< lon << ' ' << lat << '\n';
					ofs_avg << id << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << dist - vs_dist << ' ' << 1.0f / count << ' '
						<< lon << ' ' << lat << '\n';
				}
#ifdef GEN_TESTING
				std::cout << "gen ref loc: " << id << '\n';
#endif

				++id; // increment global index of reference locations of all clients
			}
		}

		ofs_random.flush();
		ofs_avg.flush();
	}
	else
		return false; // creating fails

	return true;
}

//--------------------------------------------------------------------------------
// output generated edges based on raw edges dataset, where only distances are replaced with Euclidean distance between vertices
// [re] bool: return true if generation is successful, otherwise, return false
// [in] edges: a hash map for edges, whose distances are already replaced
// [in] gen_edges_file_path: output file path of generated edges

bool gen_edges(const raw_edges &edges, const char *gen_edges_file_path)
{
	std::ofstream ofs_edges(gen_edges_file_path, std::ofstream::out | std::ofstream::trunc); // create gen_edges file
	if (!ofs_edges.fail())	// create gen_edges file successfully
	{
		raw_edges::const_iterator iter_edge = edges.cbegin();
		for (; iter_edge != edges.cend(); ++iter_edge)
		{
			ofs_edges << iter_edge->first << ' '
				<< iter_edge->second.get<RE_VS_ID>() << ' '
				<< iter_edge->second.get<RE_VE_ID>() << ' '
				<< iter_edge->second.get<RE_DIST>() << '\n';
#ifdef GEN_TESTING
			std::cout << "gen edge: " << iter_edge->first << '\n';
#endif
		}

		ofs_edges.flush();
	}
	else
		return false; // creating fails

	return true;
}