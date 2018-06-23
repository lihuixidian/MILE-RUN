#pragma once

#include "convert_datasets.h" // also include "gen_datasets.h"


//--------------------------------------------------------------------------------
// testing for output procedures

//#define REAL_TESTING // only for "real_datasets.h"

//--------------------------------------------------------------------------------
// load and snap reference locations of BJ dataset
// remark: the reference location filename and file formats must be "url#.txt" and "lon\tlat\tprob"
// [re] bool: return true if loading and snaping are successful, otherwise, return false
// [in] edges: a hash map for edges
// [in] subsegments: a hash map for each sub-segment of all edges
// [in] rtree: an R*-tree for all edge sub-segments
// [in] user_count: user count of BJ dataset
// [in] pre_file_path: file path prefix of calculated reference locations and probabilities
// [in] real_file_path: output file path of real reference locations

bool load_and_snap_reflocs_BJ(const raw_edges_ext &edges, const seg_of_edge &subsegments, const geo_subsegment_rtree &rtree, int user_count,
							  const char *pre_file_path, const char *real_file_path)
{
	unsigned global_refloc_index = 0;
	std::ofstream ofs_real(real_file_path, std::ofstream::out | std::ofstream::trunc); // create real reference locations file
	if (!ofs_real.fail())	// create real reference locations file successfully
	{
		for (int u = 0; u < user_count; ++u) // user_count of BJ is 182, part files of users are absence as they are not in Beijing
		{
			// Get file path for each user.
			std::string user_id;
			std::stringstream ss;
			ss << u;
			ss >> user_id;
			std::string each_user = pre_file_path + user_id + ".txt";

#ifdef REAL_TESTING
			std::cout << "BJ: url" << user_id << '\n';
#endif

			std::ifstream ifs_user(each_user); // read user file
			if (!ifs_user.fail())	// user file exists
			{
				while (!ifs_user.bad() && ifs_user.good())
				{
					char buf[1024];
					ifs_user.getline(buf, 1024);
					std::string str_buf(buf);

					// the user file format must be "lon\tlat\tprob"
					std::string::size_type begin_pos = 0, lat_pos, prob_pos;
					lat_pos = str_buf.find('\t', begin_pos);
					prob_pos = str_buf.find('\t', lat_pos + 1);

					float lon = static_cast<float>(atof(str_buf.substr(begin_pos, lat_pos - begin_pos).c_str()));
					float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, prob_pos - lat_pos - 1).c_str()));
					float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, str_buf.size() - prob_pos - 1).c_str()));

					// query the nn edge (via sub-segment) of each reference location
					std::vector<geo_subsegment> nn_subsegment;
					rtree.query(boost::geometry::index::nearest(geo_point(lon, lat), 1), // 1 means 1-nn (i.e., for knn, k = 1)
						std::back_inserter(nn_subsegment));

					std::vector<geo_subsegment>::iterator iter_nn = nn_subsegment.begin();
					int vs_id = edges.at(iter_nn->second.edge_id).vs_id; // vs id of the edge
					int ve_id = edges.at(iter_nn->second.edge_id).ve_id; // ve id of the edge
					float sub_seg_dist = subsegments.at(iter_nn->second).dist; // dist of the sub-segment
					float dist_vs = subsegments.at(iter_nn->second).dist_vs; // distance from edge vs to sub-segment vs
					float dist_ve = subsegments.at(iter_nn->second).dist_ve; // distance from sub-segment ve to edge ve
					float dist_vs_seg, dist_ve_seg; // distance from sub-segment vs to reference location, and from reference location to sub-segment ve
					geo_point point_on_road; // the point-on-road, namely the vertical point of real poi
					calculate_vertical_point(iter_nn->first.first, iter_nn->first.second, geo_point(lon, lat),
						sub_seg_dist, dist_vs_seg, dist_ve_seg, &point_on_road); // snap reference location to road

#ifdef REAL_TESTING
					std::cout << "BJ real: " << global_refloc_index << '\n';
#endif

					if (dist_vs + dist_vs_seg < 0.001f) // overlaps edge vs, where we view < 1m (0.001km) is also overlap
						ofs_real << global_refloc_index << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< prob << ' ' << point_on_road.get<LON>() << ' ' << point_on_road.get<LAT>() << '\n';
					else if (dist_ve_seg + dist_ve < 0.001f) // overlaps edge ve, where we view < 1m (0.001km) is also overlap
						ofs_real << global_refloc_index << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< prob << ' ' << point_on_road.get<LON>() << ' ' << point_on_road.get<LAT>() << '\n';
					else // between sub-segment vs and ve
						ofs_real << global_refloc_index << ' ' << vs_id << ' ' << ve_id << ' ' << dist_vs + dist_vs_seg << ' ' << dist_ve_seg + dist_ve << ' '
							<< prob << ' ' << point_on_road.get<LON>() << ' ' << point_on_road.get<LAT>() << '\n';
					++global_refloc_index;
				}
			}
		}
	}
	else
		return false;
	ofs_real.flush();

	return true;
}


//--------------------------------------------------------------------------------
// load and snap reference locations of CA dataset, very similar to function load_and_snap_reflocs_BJ
// remark: the reference location filename and file formats must be "url#.txt" and "lon\tlat\tprob"
// [re] bool: return true if loading and snaping are successful, otherwise, return false
// [in] edges: a hash map for edges
// [in] rtree: an R*-tree for edges
// [in] user_count: user count of CA dataset
// [in] pre_file_path: file path prefix of calculated reference locations and probabilities
// [in] real_file_path: output file path of real reference locations

bool load_and_snap_reflocs_CA(const raw_edges &edges, const geo_edge_rtree &rtree, int user_count, const char *pre_file_path, const char *real_file_path)
{
	unsigned global_refloc_index = 0;
	std::ofstream ofs_real(real_file_path, std::ofstream::out | std::ofstream::trunc); // create real reference locations file
	if (!ofs_real.fail())	// create real reference locations file successfully
	{
		for (int u = 0; u < user_count; ++u) // user_count of CA is 26619
		{
			// Get file path for each user.
			std::string user_id;
			std::stringstream ss;
			ss << u;
			ss >> user_id;
			std::string each_user = pre_file_path + user_id + ".txt";

#ifdef REAL_TESTING
			std::cout << "CA: url" << user_id << '\n';
#endif

			std::ifstream ifs_user(each_user); // read user file
			if (!ifs_user.fail())	// user file exists
			{
				while (!ifs_user.bad() && ifs_user.good())
				{
					char buf[1024];
					ifs_user.getline(buf, 1024);
					std::string str_buf(buf);

					// the user file format must be "lon\tlat\tprob"
					std::string::size_type begin_pos = 0, lat_pos, prob_pos;
					lat_pos = str_buf.find('\t', begin_pos);
					prob_pos = str_buf.find('\t', lat_pos + 1);

					float lon = static_cast<float>(atof(str_buf.substr(begin_pos, lat_pos - begin_pos).c_str()));
					float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, prob_pos - lat_pos - 1).c_str()));
					float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, str_buf.size() - prob_pos - 1).c_str()));

					// query the nn edge of each reference location
					std::vector<geo_edge> nn_edge;
					rtree.query(boost::geometry::index::nearest(geo_point(lon, lat), 1), // 1 means 1-nn (i.e., for knn, k = 1)
						std::back_inserter(nn_edge));

					std::vector<geo_edge>::iterator iter_nn = nn_edge.begin();
					float edge_dist = edges.at(iter_nn->second).get<RE_DIST>();
					int vs_id = edges.at(iter_nn->second).get<RE_VS_ID>();
					int ve_id = edges.at(iter_nn->second).get<RE_VE_ID>();
					float vs_dist, ve_dist;
					calculate_vertical_point(iter_nn->first.first, iter_nn->first.second, geo_point(lon, lat), edge_dist, vs_dist, ve_dist); // deal with arbitrary point

#ifdef REAL_TESTING
					std::cout << "CA real: " << global_refloc_index << '\n';
#endif

					if (vs_dist < 0.001f) // overlaps edge vs, where we view < 1m (0.001km) is also overlap
						ofs_real << global_refloc_index << ' ' << vs_id << ' ' << vs_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< prob << ' ' << lon << ' ' << lat << '\n';
					else if (ve_dist < 0.001f) // overlaps edge ve, where we view < 1m (0.001km) is also overlap
						ofs_real << global_refloc_index << ' ' << ve_id << ' ' << ve_id << ' ' << 0.0f << ' ' << 0.0f << ' '
							<< prob << ' ' << lon << ' ' << lat << '\n';
					else // between vs and ve
						ofs_real << global_refloc_index << ' ' << vs_id << ' ' << ve_id << ' ' << vs_dist << ' ' << ve_dist << ' '
							<< prob << ' ' << lon << ' ' << lat << '\n';
					++global_refloc_index;
				}
			}
		}
	}
	else
		return false;
	ofs_real.flush();

	return true;
}