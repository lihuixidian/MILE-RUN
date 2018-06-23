#pragma once

#include "gen_datasets.h"

//--------------------------------------------------------------------------------
// testing for output procedures

//#define COV_TESTING // only for "convert_datasets.h"

#ifdef COV_TESTING
static std::ofstream cov_out_file("D:\\Experiment\\MILE-RUN\\datasets\\BJ_datasets\\cov_testing.txt", std::ofstream::out | std::ofstream::trunc);
#endif

//--------------------------------------------------------------------------------
// edge and sub-segment information

struct edge_info
{
	int vs_id;
	int ve_id;
	float edge_dist;
	int seg_count;

	edge_info()
		: vs_id(-1)
		, ve_id(-1)
		, edge_dist(0.0f)
		, seg_count(0) {}

	edge_info(int vs_id_in, int ve_id_in, float edge_dist_in, int seg_count_in)
		: vs_id(vs_id_in)
		, ve_id(ve_id_in)
		, edge_dist(edge_dist_in)
		, seg_count(seg_count_in) {}

	edge_info(const edge_info &rhs)
	{
		if (this == &rhs)
			return;
		vs_id = rhs.vs_id;
		ve_id = rhs.ve_id;
		edge_dist = rhs.edge_dist;
		seg_count = rhs.seg_count;
	}
};

struct edge_subsegment_ids
{
	int edge_id; // edge id
	int seg_index; // the zero-based index from vs of the edge

	edge_subsegment_ids()
		: edge_id(-1)
		, seg_index(-1) {}

	edge_subsegment_ids(int edge_id_in, int seg_index_in)
		: edge_id(edge_id_in)
		, seg_index(seg_index_in) {}

	edge_subsegment_ids(const edge_subsegment_ids &rhs)
	{
		if (this == &rhs)
			return;
		edge_id = rhs.edge_id;
		seg_index = rhs.seg_index;
	}

	bool operator==(const edge_subsegment_ids &rhs) const
	{
		return edge_id == rhs.edge_id && seg_index == rhs.seg_index; // enough for unique segment
	}

	// boost::hash is implemented by calling the function hash_value
	friend std::size_t hash_value(const edge_subsegment_ids& ids)
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, ids.edge_id);
		boost::hash_combine(seed, ids.seg_index);
		return seed;
	}
};

struct edge_subsegment
{
	geo_segment seg; // the geo road segment with coordinates
	
	// note that the sum (dist_vs + dist + dist_ve) of the following three distances is the distance of the edge
	float dist; // distance of itself
	float dist_vs; // distance from vs to sub-vs
	float dist_ve; // distance from sub-ve to ve

	edge_subsegment()
		: seg() // no initialization
		, dist(0.0f)
		, dist_vs(0.0f)
		, dist_ve(0.0f) {}

	edge_subsegment(const geo_point &sub_vs, const geo_point &sub_ve, float dist_in, float dist_vs_in, float dist_ve_in)
		: seg(sub_vs, sub_ve)
		, dist(dist_in)
		, dist_vs(dist_vs_in)
		, dist_ve(dist_ve_in) {}

	edge_subsegment(const edge_subsegment &rhs)
	{
		if (this == &rhs)
			return;
		boost::geometry::set<0, 0>(seg, boost::geometry::get<0, 0>(rhs.seg));
		boost::geometry::set<0, 1>(seg, boost::geometry::get<0, 1>(rhs.seg));
		boost::geometry::set<1, 0>(seg, boost::geometry::get<1, 0>(rhs.seg));
		boost::geometry::set<1, 1>(seg, boost::geometry::get<1, 1>(rhs.seg));
		dist = rhs.dist;
		dist_vs = rhs.dist_vs;
		dist_ve = rhs.dist_ve;
	}
};

//--------------------------------------------------------------------------------
// a hash map for each sub-segment of all edges
// a hash map that stores every edge, which KEY is edge id (i.e., index for random chosen)

typedef boost::unordered_map < edge_subsegment_ids, // edge id & sub-segment index
	edge_subsegment > // segment information
seg_of_edge;

typedef boost::unordered_map < int, // edge id (index)
	edge_info > // edge information(excluding edge id)
raw_edges_ext;

//--------------------------------------------------------------------------------
// an R*-tree for all edge sub-segments

typedef std::pair < geo_segment, // edge sub-segment
	edge_subsegment_ids > // edge id & sub-segment index
geo_subsegment;
typedef boost::geometry::index::rtree<geo_subsegment, boost::geometry::index::rstar<16>> geo_subsegment_rtree; // a sub-segment R*-tree

//--------------------------------------------------------------------------------
// load edges and sub-segments (geos) datasets; also construct an R*-tree to organize sub-segments
// remark: because we use float type for longitude and latitude values, precision lost is permitted;
//		   the edges file format must be "edge_id\tvs_id\tve_id"; the geos file format must be "edge_id\tlat#1\tlon#1\t...";
//		   the NULL-checking of "geo_range_file_path" is only for integrating real dataset
// [re] bool: return true if loading is successful, otherwise, return false
// [out] edges: a hash map for edges
// [out] subsegments: a hash map for sub-segments
// [out] rtree: an R*-tree for sub-segments
// [in] geo_range_file_path: output file path of geo-range
// [in] edges_file_path: file path of edges
// [in] geos_file_path: file path of geos (sub-segments)

bool load_edges_and_geos(raw_edges_ext &edges, seg_of_edge &subsegments, geo_subsegment_rtree &rtree, const char *geo_range_file_path,
						 const char *edges_file_path, const char *geos_file_path)
{
	// first load edges
	std::ifstream ifs_edges(edges_file_path); // read edges file
	if (!ifs_edges.fail())	// edges file exists
	{
		while (!ifs_edges.bad() && ifs_edges.good())
		{
			char buf[1024];
			ifs_edges.getline(buf, 1024);
			std::string str_buf(buf);

			// the edges file format must be "edge_id\tvs_id\tve_id"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos;
			vs_pos = str_buf.find('\t', begin_pos);
			ve_pos = str_buf.find('\t', vs_pos + 1);

			int edge_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, str_buf.size() - ve_pos - 1).c_str());

			// we will accumulate Euclidean distances between sub-segment vertices for edge distance, and now we initialize it as 0.0f;
			// and the count of segments is also initialized as 0
			edges[edge_id] = edge_info(vs_id, ve_id, 0.0f, 0); // insert a new edge
#ifdef COV_TESTING
			std::cout << "e: " << edge_id << '\n';
#endif
		}
	}
	else
		return false; // loading fails, as reading edges fails

	float min_lon = 180.0f, max_lon = -180.0f, min_lat = 90.0f, max_lat = -90.0f; // value for geo-range

	// then load sub-segments (geos)
	std::ifstream ifs_geos(geos_file_path); // read geos file
	if (!ifs_geos.fail())	// sub-segments (geos) file exists
	{
		while (!ifs_geos.bad() && ifs_geos.good())
		{
			char buf[8192]; // use larger buffer as the number of segments is larger
			ifs_geos.getline(buf, 8192);
			std::string str_buf(buf);

			// the geos file format must be "edge_id\tlat#1\tlon#1\t..."
			std::string::size_type begin_pos = 0, coordinate_pos;
			coordinate_pos = str_buf.find('\t', begin_pos);
			int edge_id = atoi(str_buf.substr(begin_pos, coordinate_pos - begin_pos).c_str());

			// iteratively read each coordinate "lat\tlon"
			std::vector<geo_point> sub_vertices; // for all sub-vertices of sub-segments of the edge
			std::string::size_type lat_pos = coordinate_pos + 1, lon_pos, next_lat_pos;
			while (true) // this iteration will be terminated by "break"
			{
				lon_pos = str_buf.find('\t', lat_pos);
				next_lat_pos = str_buf.find('\t', lon_pos + 1);

				float lat = static_cast<float>(atof(str_buf.substr(lat_pos, lon_pos - lat_pos).c_str()));
				if (next_lat_pos != std::string::npos) // still more vertices
				{
					float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, next_lat_pos - lon_pos - 1).c_str()));
					sub_vertices.push_back(geo_point(lon, lat)); // record each sub-vertex
					lat_pos = next_lat_pos + 1; // for next sub-vertex

					// comparation for geo-range
					if (lon < min_lon)
						min_lon = lon;
					if (lon > max_lon)
						max_lon = lon;
					if (lat < min_lat)
						min_lat = lat;
					if (lat > max_lat)
						max_lat = lat;
				}
				else // next_lat_pos == std::string::npos, it's the last sub-vertex
				{
					float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, str_buf.size() - lon_pos - 1).c_str()));
					sub_vertices.push_back(geo_point(lon, lat)); // record each sub-vertex

					// comparation for geo-range
					if (lon < min_lon)
						min_lon = lon;
					if (lon > max_lon)
						max_lon = lon;
					if (lat < min_lat)
						min_lat = lat;
					if (lat > max_lat)
						max_lat = lat;

					break;
				}
			}

			// accumulate distances of all sub-segments of the edge
			float edge_dist = 0.0f;
			std::vector<float> segments_dist; // record distance of each segment
			for (unsigned i = 0; i < sub_vertices.size() - 1; ++i)
			{
				float dist = boost::geometry::distance(sub_vertices[i], sub_vertices[i + 1]) * EARTH_RADIUS;
				edge_dist += dist; // accumulation
				segments_dist.push_back(dist);
			}
			edges[edge_id].edge_dist = edge_dist; // set actual road distance
			edges[edge_id].seg_count = static_cast<int>(segments_dist.size()); // set actual count of segments
#ifdef COV_TESTING
			cov_out_file << "e: " << edge_id << ", dist: " << edges[edge_id].edge_dist << ", count: " << edges[edge_id].seg_count << std::endl;
#endif

			// construct sub-segments and R*-tree
			float dist_vs = 0.0f; // distance from vs of the edge to vs of the sub-segment
			for (unsigned i = 0; i < sub_vertices.size() - 1; ++i)
			{
				subsegments[edge_subsegment_ids(edge_id, i)] = edge_subsegment(sub_vertices[i], sub_vertices[i + 1], segments_dist[i], dist_vs, edge_dist - dist_vs - segments_dist[i]);
				dist_vs += segments_dist[i]; // accumulate distance from vs

				// construct sub-segments R*-tree
				rtree.insert(std::make_pair(geo_segment(sub_vertices[i], sub_vertices[i + 1]), edge_subsegment_ids(edge_id, i)));
#ifdef COV_TESTING
				std::cout << "e: " << edge_id << ", sub-seg: " << i << '\n';
#endif
			}
		}
	}
	else
		return false; // loading fails, as reading edges fails

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

	return true;
}

//--------------------------------------------------------------------------------
// load and pick pois, also including integration facilities
// remark: because we use float type for longitude and latitude values, precision lost is permitted;
//		   the facilities file format must be "lon\tlat"
// [re] bool: return true if loading and picking are successful, otherwise, return false
// [in] edges: a hash map for edges
// [in] subsegments: a hash map for each sub-segment of all edges, MUST be used as a const peremeter
// [in] rtree: an R*-tree for all edge sub-segments
// [in] raw_pois_file_path: file path of raw pois
// [in] out_pois_file_path: file path of output picked pois
// [in] poi_count: count of picked pois
// [out] picked_facs: append each picked facility to this hash set, used for uniqueness check
//					  NOTE facilities cannot overlap each other
// [in/def] is_with_coordinates: whether the coordinates of facilities are appended at the end of each line

bool load_and_pick_pois(const raw_edges_ext &edges, seg_of_edge &subsegments, const geo_subsegment_rtree &rtree, const char *raw_pois_file_path,
						const char *out_pois_file_path, int poi_count, boost::unordered_set<fac_or_cand_loc> &picked_facs, bool is_with_coordinates = false)
{
	// load raw pois
	std::vector<std::pair<float, float> > pois;
	std::ifstream ifs_pois(raw_pois_file_path); // read raw pois file
	if (!ifs_pois.fail())	// raw pois file exists
	{
		while (!ifs_pois.bad() && ifs_pois.good())
		{
			char buf[1024];
			ifs_pois.getline(buf, 1024);
			std::string str_buf(buf);

			// the raw pois file format must be "lon\tlat"
			std::string::size_type begin_pos = 0, lat_pos;
			lat_pos = str_buf.find('\t', begin_pos);

			float lon = static_cast<float>(atof(str_buf.substr(begin_pos, lat_pos - begin_pos).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			pois.push_back(std::make_pair(lon, lat)); // record each poi
		}
	}
	else
		return false; // loading fails

	// randomly pick pois
	std::ofstream ofs_pois(out_pois_file_path, std::ofstream::out | std::ofstream::trunc); // create picked pois file
	if (!ofs_pois.fail())	// create pois file successfully
	{
		int poi_index = 0; // index of picked poi, which is also the count of currently picked pois
		while (poi_index < poi_count)
		{
#ifdef COV_TESTING
			std::cout << "pick poi: " << poi_index << '\n';
#endif

			boost::random::uniform_int_distribution<> random_index(0, static_cast<int>(pois.size()) - 1); // poi index begins from 0
			int selected_index = random_index(random_gen);

			// query the nn edge (via sub-segment) of each poi
			std::vector<geo_subsegment> nn_subsegment;
			rtree.query(boost::geometry::index::nearest(geo_point(pois[selected_index].first, pois[selected_index].second), 32), // 32 means 32-nn (i.e., for knn, k = 32)
				std::back_inserter(nn_subsegment));																				 // use 32-nn to avoid exhaustion of "pois"
			std::vector<geo_subsegment>::iterator iter_nn = nn_subsegment.begin();
			for (; iter_nn != nn_subsegment.end(); ++iter_nn) // check for knn sub-segment, if necessary
			{
				int vs_id = edges.at(iter_nn->second.edge_id).vs_id; // vs id of the edge
				int ve_id = edges.at(iter_nn->second.edge_id).ve_id; // ve id of the edge
				float sub_seg_dist = subsegments[iter_nn->second].dist; // dist of the sub-segment
				float dist_vs = subsegments[iter_nn->second].dist_vs; // distance from edge vs to sub-segment vs
				float dist_ve = subsegments[iter_nn->second].dist_ve; // distance from sub-segment ve to edge ve
				float dist_vs_seg, dist_ve_seg; // distance from sub-segment vs to poi, and from poi to sub-segment ve
				calculate_vertical_point(iter_nn->first.first, iter_nn->first.second, geo_point(pois[selected_index].first, pois[selected_index].second),
					sub_seg_dist, dist_vs_seg, dist_ve_seg); // snap poi to road

				// output the poi, using hash set "picked_facs" for uniqueness check
				std::pair<boost::unordered_set<fac_or_cand_loc>::iterator, bool> return_value; // the 2nd component of return value of "insert" shows whether insertion takes place or not
				if (dist_vs + dist_vs_seg < 0.001f) // overlaps edge vs, where we view < 1m (0.001km) is also overlap
				{
					return_value = picked_facs.insert(fac_or_cand_loc(vs_id, vs_id, 0.0f));
					if (return_value.second) // insertion takes place
					{
						ofs_pois << poi_index << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< pois[selected_index].first << ' ' << pois[selected_index].second << '\n';
						++poi_index;
						break;
					}
				}
				else if (dist_ve_seg + dist_ve < 0.001f) // overlaps edge ve, where we view < 1m (0.001km) is also overlap
				{
					return_value = picked_facs.insert(fac_or_cand_loc(ve_id, ve_id, 0.0f));
					if (return_value.second) // insertion takes place
					{
						ofs_pois << poi_index << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< pois[selected_index].first << ' ' << pois[selected_index].second << '\n';
						++poi_index;
						break;
					}
				}
				else // between sub-segment vs and ve
				{
					return_value = picked_facs.insert(fac_or_cand_loc(vs_id, ve_id, dist_vs + dist_vs_seg));
					if (return_value.second) // insertion takes place
					{
						ofs_pois << poi_index << ' ' << vs_id << ' ' << ve_id << ' ' << dist_vs + dist_vs_seg << ' '
							<< dist_ve_seg + dist_ve << ' ' << pois[selected_index].first << ' ' << pois[selected_index].second << '\n';
						++poi_index;
						break;
					}
				}
			}

			// delete the picked poi from all pois
			pois[selected_index] = pois.back();
			pois.pop_back();
		}
		ofs_pois.flush();
	}
	else
		return false; // creating fails

	return true;
}

//--------------------------------------------------------------------------------
// output generated candidates
// [re] bool: return true if generation is successful, otherwise, return false
// [in] edges: a hash map for edges
// [in] subsegments: a hash map for each sub-segment of all edges, MUST be used as a const peremeter
// [in] gen_cands_file_path: output file path of generated candidates
// [in] cands_count: the count of candidates to be generated
// [in/out] picked_facs: append each new generated candidate to this hash set, used for uniqueness check;
//						 NOTE facilities cannot overlap each other, and candidates cannot overlap facilities and each other

bool gen_cands(const raw_edges_ext &edges, seg_of_edge &subsegments, boost::unordered_set<fac_or_cand_loc> &picked_facs, const char *gen_cands_file_path, int cands_count)
{
	boost::random::uniform_int_distribution<> random_edge_index(1, static_cast<int>(edges.size())); // from all edges, edge index begins from 1
	boost::random::uniform_int_distribution<> random_percent(0, 1000); // distance percent from sub-segment vs to candidate, accuracy is xx.x%

	std::ofstream ofs_cands(gen_cands_file_path, std::ofstream::out | std::ofstream::trunc); // create gen_cands file
	if (!ofs_cands.fail())	// create gen_cands file successfully
	{
		for (int i = 0; i < cands_count; ++i)
		{
			// pick an edge according to random edge index
			int edge_index = random_edge_index(random_gen);
			int edge_vs_id = edges.at(edge_index).vs_id;
			int edge_ve_id = edges.at(edge_index).ve_id;
			float edge_dist = edges.at(edge_index).edge_dist;
			int seg_count = edges.at(edge_index).seg_count;
			
			// pick a segment according to random segment index
			boost::random::uniform_int_distribution<> random_seg_index(0, seg_count - 1); // from the edge's sub-segments
			int seg_index = random_seg_index(random_gen);
			edge_subsegment_ids seg(edge_index, seg_index); // the picked segment

			// generate a candidate, and insert as a facility such that to avoid overlap
			float dist_percent = -1.0f, dist_sub_vs = -1.0f, dist_sub_ve = 0.0f; // initial values
			std::pair<boost::unordered_set<fac_or_cand_loc>::iterator, bool> return_value; // the 2nd component of return value of "insert" shows whether insertion takes place or not
			do
			{
				dist_percent = random_percent(random_gen) / 1000.0f; // sub_vs->c / c->sub_ve, accuracy is xx.x%
				dist_sub_vs = subsegments[seg].dist * dist_percent;
				dist_sub_ve = subsegments[seg].dist - dist_sub_vs;
				if (subsegments[seg].dist_vs + dist_sub_vs < 0.001f) // candidate overlaps edge vs, where we view < 1m (0.001km) is also overlap
					return_value = picked_facs.insert(fac_or_cand_loc(edge_vs_id, edge_vs_id, 0.0f));
				else if (dist_sub_ve + subsegments[seg].dist_ve < 0.001f) // candidate overlaps edge ve, where we view < 1m (0.001km) is also overlap
					return_value = picked_facs.insert(fac_or_cand_loc(edge_ve_id, edge_ve_id, 0.0f));
				else
					return_value = picked_facs.insert(fac_or_cand_loc(edge_vs_id, edge_ve_id, subsegments[seg].dist_vs + dist_sub_vs));
			} while (!return_value.second); // insertion fails

			// output the generated candidate
			if (subsegments[seg].dist_vs + dist_sub_vs < 0.001f) // candidate overlaps edge vs, where we view < 1m (0.001km) is also overlap
				ofs_cands << i << ' ' << edge_vs_id << ' ' << edge_vs_id << ' ' << 0.0f << ' ' << 0.0f << ' '
				<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LON>() << ' '
				<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LAT>() << '\n';
			else if (dist_sub_ve + subsegments[seg].dist_ve < 0.001f) // candidate overlaps edge ve, where we view < 1m (0.001km) is also overlap
				ofs_cands << i << ' ' << edge_ve_id << ' ' << edge_ve_id << ' ' << 0.0f << ' ' << 0.0f << ' '
				<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LON>() << ' '
				<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LAT>() << '\n';
			else
			{
				// c.x = "edge.length.x" * percent(vs->c) + vs.x = (ve.x - vs.x) * percent(vs->c) + vs.x; c.y is similar
				float lon = (subsegments[seg].seg.second.get<LON>() - subsegments[seg].seg.first.get<LON>()) * dist_percent + subsegments[seg].seg.first.get<LON>();
				float lat = (subsegments[seg].seg.second.get<LAT>() - subsegments[seg].seg.first.get<LAT>()) * dist_percent + subsegments[seg].seg.first.get<LAT>();
				ofs_cands << i << ' ' << edge_vs_id << ' ' << edge_ve_id << ' '
					<< subsegments[seg].dist_vs + dist_sub_vs << ' ' << edge_dist - subsegments[seg].dist_vs - dist_sub_vs << ' ' << lon << ' ' << lat << '\n';
			}
#ifdef COV_TESTING
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
// [in] edges: a hash map for edges
// [in] subsegments: a hash map for each sub-segment of all edges, MUST be used as a const peremeter
// [in] gen_random_ref_locs_file_path: output file path of generated reference locations with randomly chosen present probabilities
// [in] gen_avg_ref_locs_file_path: output file path of generated reference locations with average present probabilities
// [in] clients_count: the count of clients to be generated
// [in] ref_loc_points: the probabilities of different counts of reference locations

bool gen_ref_locs(const raw_edges_ext &edges, seg_of_edge &subsegments, const char *gen_random_ref_locs_file_path, const char *gen_avg_ref_locs_file_path,
				  int clients_count, const char *ref_loc_points)
{
	// parse the probabilities of different counts of reference locations; the probabilities string format must be "count count#1 prob#1 ..."
	// here, we DON'T parse the input "ref_loc_points"
	// SHOULD parse "ref_loc_points" to obtain the following two arrays
	int counts[] = { 1, 2, 3, 4, 5, 6 };
	float probabilities[] = { 0.05f, 0.2f, 0.35f, 0.25f, 0.1f, 0.05f };

	boost::random::discrete_distribution<> random_count(probabilities); // reference locations count with different probabilities
	boost::random::uniform_int_distribution<> random_edge_index(1, static_cast<int>(edges.size())); // from all edges, edge index begins from 1
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
				int edge_index = random_edge_index(random_gen);
				int edge_vs_id = edges.at(edge_index).vs_id;
				int edge_ve_id = edges.at(edge_index).ve_id;
				float edge_dist = edges.at(edge_index).edge_dist;
				int seg_count = edges.at(edge_index).seg_count;

				// pick a segment according to random segment index
				boost::random::uniform_int_distribution<> random_seg_index(0, seg_count - 1); // from the edge's sub-segments
				int seg_index = random_seg_index(random_gen);
				edge_subsegment_ids seg(edge_index, seg_index); // the picked segment

				// generate a reference location
				float dist_percent = random_percent(random_gen) / 1000.0f;
				float dist_sub_vs = subsegments[seg].dist * dist_percent;
				float dist_sub_ve = subsegments[seg].dist - dist_sub_vs;

				// output the generated reference location
				if (subsegments[seg].dist_vs + dist_sub_vs < 0.001f) // reference location overlaps edge vs, where we view < 1m (0.001km) is also overlap
				{
					ofs_random << id << ' ' << edge_vs_id << ' ' << edge_vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << present_prob << ' '
						<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LON>() << ' '
						<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LAT>() << '\n';
					ofs_avg << id << ' ' << edge_vs_id << ' ' << edge_vs_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << 1.0f / count << ' '
						<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LON>() << ' '
						<< subsegments[edge_subsegment_ids(edge_index, 0)].seg.first.get<LAT>() << '\n';
				}
				else if (dist_sub_ve + subsegments[seg].dist_ve < 0.001f) // reference location overlaps edge ve, where we view < 1m (0.001km) is also overlap
				{
					ofs_random << id << ' ' << edge_ve_id << ' ' << edge_ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << present_prob << ' '
						<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LON>() << ' '
						<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LAT>() << '\n';
					ofs_avg << id << ' ' << edge_ve_id << ' ' << edge_ve_id << ' ' << 0.0f << ' ' << 0.0f << ' ' << 1.0f / count << ' '
						<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LON>() << ' '
						<< subsegments[edge_subsegment_ids(edge_index, seg_count - 1)].seg.second.get<LAT>() << '\n';
				}
				else
				{
					// c.x = "edge.length.x" * percent(vs->c) + vs.x = (ve.x - vs.x) * percent(vs->c) + vs.x; c.y is similar
					float lon = (subsegments[seg].seg.second.get<LON>() - subsegments[seg].seg.first.get<LON>()) * dist_percent + subsegments[seg].seg.first.get<LON>();
					float lat = (subsegments[seg].seg.second.get<LAT>() - subsegments[seg].seg.first.get<LAT>()) * dist_percent + subsegments[seg].seg.first.get<LAT>();
					ofs_random << id << ' ' << edge_vs_id << ' ' << edge_ve_id << ' ' << subsegments[seg].dist_vs + dist_sub_vs << ' '
						<< edge_dist - subsegments[seg].dist_vs - dist_sub_vs << ' ' << present_prob << ' ' << lon << ' ' << lat << '\n';
					ofs_avg << id << ' ' << edge_vs_id << ' ' << edge_ve_id << ' ' << subsegments[seg].dist_vs + dist_sub_vs << ' '
						<< edge_dist - subsegments[seg].dist_vs - dist_sub_vs << ' ' << 1.0f / count << ' ' << lon << ' ' << lat << '\n';
				}
#ifdef COV_TESTING
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
// [in] edges: a hash map for edges, whose distances are already accumulated
// [in] gen_edges_file_path: output file path of generated edges

bool gen_edges(const raw_edges_ext &edges, const char *gen_edges_file_path)
{
	std::ofstream ofs_edges(gen_edges_file_path, std::ofstream::out | std::ofstream::trunc); // create gen_edges file
	if (!ofs_edges.fail())	// create gen_edges file successfully
	{
		raw_edges_ext::const_iterator iter_edge = edges.cbegin();
		for (; iter_edge != edges.cend(); ++iter_edge)
		{
			ofs_edges << iter_edge->first << ' '
				<< iter_edge->second.vs_id << ' '
				<< iter_edge->second.ve_id << ' '
				<< iter_edge->second.edge_dist << '\n';
#ifdef COV_TESTING
			std::cout << "gen edge: " << iter_edge->first << '\n';
#endif
		}

		ofs_edges.flush();
	}
	else
		return false; // creating fails

	return true;
}