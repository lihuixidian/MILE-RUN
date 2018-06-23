#include "stdafx.h"

#include <chrono> // std::chrono::high_resolution_clock
#include <string> // std::string
#include <fstream> // std::ifstream, std::ofstream

#include "gen_datasets.h"
#include "convert_datasets.h"
#include "real_datasets.h"
#include "proposed_algorithms.h"
#include "ir_utility.h"
#include "compared_algorithms.h"


bool gen_datasets(_TCHAR* argv); // forward declaration
bool convert_datasets(_TCHAR* argv); // forward declaration
bool real_datasets(_TCHAR* argv); // forward declaration


int _tmain(int argc, _TCHAR* argv[])
{
	// parameter "-ap" means calculate AP for two ranked sequences
	// argv[2]: real data file path
	// argv[3]: query data file path
	// argv[4]: AP result file path
	if (std::wstring(argv[1]) == _T("-ap"))
	{
		std::ifstream in_real_file(argv[2]);
		std::map<unsigned, std::string> real_map;
		load_ranked_data(in_real_file, real_map, 0, '\t');

		std::ifstream in_query_file(argv[3]);
		std::map<unsigned, std::string> query_map;
		load_ranked_data(in_query_file, query_map, 0, '\t');

		double ap = 0.0;
		APatK(real_map, query_map, ap);

		std::ofstream ofs_ap(argv[4], std::ofstream::out | std::ofstream::trunc); // create AP result file
		if (!ofs_ap.fail())	// creating AP result file successfully
			ofs_ap << "AP@" << real_map.size() << ": " << ap << std::endl;
		else
			return -1;

		return 0;
	}

	// read config file for datasets files
	bool need_gen, need_cov, need_real;
	bool is_directed_graph;
	std::string //vertices_file_path, now this variable is useless
		edges_file_path, facs_file_path, reflocs_file_path, cands_file_path;
	std::string results_file_path;
	bool is_EN_50;
	bool is_execute_EN, is_execute_LNB, is_execute_NSJ;
	bool is_execute_Blossom, is_execute_OTF;
	std::ifstream ifs_config(argv[1]);
	if (!ifs_config.fail())	// Config file exists.
	{
		while (!ifs_config.bad() && ifs_config.good())
		{
			char buf[1024];
			ifs_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);
			
			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "need_gen")
				need_gen = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "need_cov")
				need_cov = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "need_real")
				need_real = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_directed_graph")
				is_directed_graph = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			// vertices <id lon lat> file is useless now
			//else if (str_param == "vertices")
			//	vertices_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "edges")
				edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "facs")
				facs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "reflocs")
				reflocs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "cands")
				cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "results")
				results_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "is_EN_50")
				is_EN_50 = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_EN")
				is_execute_EN = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_LNB")
				is_execute_LNB = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_NSJ")
				is_execute_NSJ = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_Blossom")
				is_execute_Blossom = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_OTF")
				is_execute_OTF = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
		}
	}

	if (need_gen)
		if (!gen_datasets(argv[2])) // generate candidates and reference locations
			return -1;

	if (need_cov)
		if (!convert_datasets(argv[3])) // convert facilities, then generate candidates and reference locations
			return -1;

	if (need_real)
		if (!real_datasets(argv[4])) // snap real reference locations to road network
			return -1;

	// output results
	std::ofstream ofs_re(results_file_path, std::ofstream::out | std::ofstream::trunc); // create results file
	if (ofs_re.fail())	// creating results file fails
		return -1;

	// time recorders
	auto begin_time = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();

	// EN for top-50 candidates
	if (is_EN_50)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		std::vector<candidate> topk;
		EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str(),
			50, &topk); // for top-50 candidates
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "EN took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n";
		for (int i = 0; i < 50; ++i)
			ofs_re << topk[i].id << '\t' << topk[i].ERD << '\n';
		ofs_re.flush();

		return 0;
	}

	candidate ol_EN, ol_LNB, ol_NSJ; // respective results

	if (is_execute_EN)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		ol_EN = EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "EN took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			   << "   result is c" << ol_EN.id << ", ERD = " << ol_EN.ERD << std::endl;
	}

	if (is_execute_LNB)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_LNB;
		fac_hash_map facs_map;
		LNT_hash_map LNT;
		vertex_nnd_hash_map nnd_map;
		LNB_construct_LNT(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), LNT, nnd_map, graph_for_LNB, facs_map);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB construct LNT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		begin_time = std::chrono::high_resolution_clock::now();
		int checked_num;
		ol_LNB = LNB_query(checked_num, LNT, nnd_map, graph_for_LNB, facs_map, cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			   << "    check num = " << checked_num << '\n'
			   << "    result is c" << ol_LNB.id << ", ERD = " << ol_LNB.ERD << std::endl;
	}

	if (is_execute_NSJ)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_NSJ;
		fac_hash_map reflocs;
		nnfc_hash_map NNFCs;
		NSJ_construct_NNFC(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), reflocs, NNFCs, graph_for_NSJ);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "NSJ construct NNFCs took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		begin_time = std::chrono::high_resolution_clock::now();
		int checked_num;
		__int64 rtree_time;
		ol_NSJ = NSJ_query(checked_num, rtree_time, reflocs, NNFCs, graph_for_NSJ, cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "NSJ took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			   << "    R*-tree query took " << rtree_time << "ms, check num = " << checked_num << '\n'
			   << "    result is c" << ol_NSJ.id << ", ERD = " << ol_NSJ.ERD << std::endl;
	}

	candidate ol_Blossom, ol_OTF; // respective results

	if (is_execute_Blossom)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		ol_Blossom = Blossom(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "Blossom took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "   result is c" << ol_Blossom.id << ", ERD = " << ol_Blossom.ERD << std::endl;
	}

	if (is_execute_OTF)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		ol_OTF = OTF(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "OTF took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "   result is c" << ol_OTF.id << ", ERD = " << ol_OTF.ERD << std::endl;
	}

	return 0;
}

// generate datasets based on CA datasets
bool gen_datasets(_TCHAR* argv)
{
	// read gen_config file for generating datasets files
	std::string geo_range_file_path;
	bool need_integrate_facs = false;
	bool fac_with_coordinates = false;
	int bar, hospital, po, park, school;
	std::string bars_file_path, hospitals_file_path, pos_file_path, parks_file_path, schools_file_path;
	std::string bars_200_file_path, hospitals_800_file_path, pos_800_file_path, parks_3200_file_path, schools_3200_file_path;
	bool need_gen_cands = false;
	std::string str_cands_count_1, str_cands_count_2, str_cands_count_3, str_cands_count_4, str_cands_count_5;
	int cands_count_1 = 0, cands_count_2 = 0, cands_count_3 = 0, cands_count_4 = 0, cands_count_5 = 0;
	std::string gen_cands_file_path;
	bool need_gen_ref_loc = false;
	std::string str_clients_count_1, str_clients_count_2, str_clients_count_3, str_clients_count_4, str_clients_count_5;
	int clients_count_1 = 0, clients_count_2 = 0, clients_count_3 = 0, clients_count_4 = 0, clients_count_5 = 0;
	std::string ref_loc_points, gen_random_ref_locs_file_path, gen_avg_ref_locs_file_path;
	bool need_check_directionality = false;
	bool need_gen_edges = false;
	std::string gen_edges_file_path;
	std::string vertices_file_path, edges_file_path, pois_file_path;

	std::ifstream ifs_gen_config(argv);
	if (!ifs_gen_config.fail())	// gen_config file exists.
	{
		while (!ifs_gen_config.bad() && ifs_gen_config.good())
		{
			char buf[1024];
			ifs_gen_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "geo_range")
				geo_range_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for facilities
			else if(str_param == "need_integrate_facs")
				need_integrate_facs = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "fac_with_coordinates")
				fac_with_coordinates = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "bar" && need_integrate_facs)
				bar = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "bars" && need_integrate_facs)
				bars_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "bars_200" && need_integrate_facs)
				bars_200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "hospital" && need_integrate_facs)
				hospital = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "hospitals" && need_integrate_facs)
				hospitals_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "hospitals_800" && need_integrate_facs)
				hospitals_800_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "po" && need_integrate_facs)
				po = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "pos" && need_integrate_facs)
				pos_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "pos_800" && need_integrate_facs)
				pos_800_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "park" && need_integrate_facs)
				park = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "parks" && need_integrate_facs)
				parks_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "parks_3200" && need_integrate_facs)
				parks_3200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "school" && need_integrate_facs)
				school = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "schools" && need_integrate_facs)
				schools_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "schools_3200" && need_integrate_facs)
				schools_3200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for candidates
			else if (str_param == "need_gen_cands")
				need_gen_cands = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "cands_count_1" && need_gen_cands)
			{
				str_cands_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_1 = atoi(str_cands_count_1.c_str());
			}
			else if (str_param == "cands_count_2" && need_gen_cands)
			{
				str_cands_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_2 = atoi(str_cands_count_2.c_str());
			}
			else if (str_param == "cands_count_3" && need_gen_cands)
			{
				str_cands_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_3 = atoi(str_cands_count_3.c_str());
			}
			else if (str_param == "cands_count_4" && need_gen_cands)
			{
				str_cands_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_4 = atoi(str_cands_count_4.c_str());
			}
			else if (str_param == "cands_count_5" && need_gen_cands)
			{
				str_cands_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_5 = atoi(str_cands_count_5.c_str());
			}
			else if (str_param == "gen_cands" && need_gen_cands)
				gen_cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for reference locations (clients)
			else if (str_param == "need_gen_ref_loc")
				need_gen_ref_loc = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "clients_count_1" && need_gen_ref_loc)
			{
				str_clients_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_1 = atoi(str_clients_count_1.c_str());
			}
			else if (str_param == "clients_count_2" && need_gen_ref_loc)
			{
				str_clients_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_2 = atoi(str_clients_count_2.c_str());
			}
			else if (str_param == "clients_count_3" && need_gen_ref_loc)
			{
				str_clients_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_3 = atoi(str_clients_count_3.c_str());
			}
			else if (str_param == "clients_count_4" && need_gen_ref_loc)
			{
				str_clients_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_4 = atoi(str_clients_count_4.c_str());
			}
			else if (str_param == "clients_count_5" && need_gen_ref_loc)
			{
				str_clients_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_5 = atoi(str_clients_count_5.c_str());
			}
			else if (str_param == "ref_loc_points" && need_gen_ref_loc)
				ref_loc_points = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "gen_random_ref_locs" && need_gen_ref_loc)
				gen_random_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "gen_avg_ref_locs" && need_gen_ref_loc)
				gen_avg_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for edges
			else if (str_param == "need_gen_edges")
				need_gen_edges = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "gen_edges")
				gen_edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for background datasets
			else if (str_param == "need_check_directionality")
				need_check_directionality = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "vertices")
				vertices_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "edges")
				edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "pois")
				pois_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
		}
	}
	else
		return false;

	if (need_integrate_facs || need_gen_cands || need_gen_ref_loc || need_gen_edges)
	{
		raw_vertices vertices;
		raw_edges edges;
		geo_edge_rtree rtree;
		if (!load_vertices_and_edges(vertices, edges, rtree, geo_range_file_path.c_str(), vertices_file_path.c_str(), edges_file_path.c_str(), need_check_directionality))
			return false;

		boost::unordered_set<fac_or_cand_loc> gen_facs; // locations of generated facilities
		if (need_integrate_facs)
		{
			// initialize counts array
			int pois_count[CATEGORY];
			for (int i = 0; i < CATEGORY; ++i)
				pois_count[i] = 0;

			// load raw pois
			std::vector<geo_point> bars, hospitals, pos, parks, schools; // raw poi containers
			if (!load_raw_pois(pois_file_path.c_str(), pois_count, bars, hospitals, pos, parks, schools))
				return false;

			// integrate pois
			/*if (!integrate_pois(edges, rtree, bars, bars_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, hospitals, hospitals_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, pos, pos_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, parks, parks_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, schools, schools_file_path.c_str(), gen_facs))*/
			if (!integrate_pois(edges, rtree, hospitals, hospitals_file_path.c_str(), gen_facs, fac_with_coordinates)
				|| !integrate_pois(edges, rtree, parks, parks_file_path.c_str(), gen_facs, fac_with_coordinates)
				|| !integrate_pois(edges, rtree, schools, schools_file_path.c_str(), gen_facs, fac_with_coordinates))
				return false;

			// pick pois
			/*if (!pick_pois(bars_file_path.c_str(), bars_200_file_path.c_str(), bar)
				|| !pick_pois(hospitals_file_path.c_str(), hospitals_800_file_path.c_str(), hospital)
				|| !pick_pois(pos_file_path.c_str(), pos_800_file_path.c_str(), po)
				|| !pick_pois(parks_file_path.c_str(), parks_3200_file_path.c_str(), park)
				|| !pick_pois(schools_file_path.c_str(), schools_3200_file_path.c_str(), school))*/
			if (!pick_pois(hospitals_file_path.c_str(), hospitals_800_file_path.c_str(), hospital)
				|| !pick_pois(parks_file_path.c_str(), parks_3200_file_path.c_str(), park)
				|| !pick_pois(schools_file_path.c_str(), schools_3200_file_path.c_str(), school))
				return false;
		}

		// generate candidates
		if (need_gen_cands)
		{
			if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_1 + ".txt").c_str(), cands_count_1))
				return false;
			if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_2 + ".txt").c_str(), cands_count_2))
				return false;
			if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_3 + ".txt").c_str(), cands_count_3))
				return false;
			if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_4 + ".txt").c_str(), cands_count_4))
				return false;
			if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_5 + ".txt").c_str(), cands_count_5))
				return false;
		}

		// generate reference locations
		if (need_gen_ref_loc)
		{
			if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(), clients_count_1, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(), clients_count_2, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(), clients_count_3, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(), clients_count_4, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(), clients_count_5, ref_loc_points.c_str()))
				return false;
		}

		// just output edges, whose distances are already replaced with Euclidean distance between vertices
		if (need_gen_edges
			&& !gen_edges(edges, gen_edges_file_path.c_str()))
			return false;
	}

	return true;
}

// convert and generate datasets based on Beijing datasets
bool convert_datasets(_TCHAR* argv)
{
	// read cov_config file for generating datasets files
	std::string cov_range_file_path;
	bool need_integrate_facs = false;
	bool fac_with_coordinates = false;
	int station, cafe, logistic, bank, school;
	std::string stations_file_path, cafes_file_path, logistics_file_path, banks_file_path, schools_file_path;
	std::string stations_1500_file_path, cafes_1500_file_path, logistics_2500_file_path, banks_2500_file_path, schools_2500_file_path;
	bool need_gen_cands = false;
	std::string str_cands_count_1, str_cands_count_2, str_cands_count_3, str_cands_count_4, str_cands_count_5;
	int cands_count_1 = 0, cands_count_2 = 0, cands_count_3 = 0, cands_count_4 = 0, cands_count_5 = 0;
	std::string gen_cands_file_path;
	bool need_gen_ref_loc = false;
	std::string str_clients_count_1, str_clients_count_2, str_clients_count_3, str_clients_count_4, str_clients_count_5;
	int clients_count_1 = 0, clients_count_2 = 0, clients_count_3 = 0, clients_count_4 = 0, clients_count_5 = 0;
	std::string ref_loc_points, gen_random_ref_locs_file_path, gen_avg_ref_locs_file_path;
	bool need_gen_edges = false;
	std::string gen_edges_file_path;
	std::string edges_file_path, geos_file_path;

	std::ifstream ifs_cov_config(argv);
	if (!ifs_cov_config.fail())	// cov_config file exists.
	{
		while (!ifs_cov_config.bad() && ifs_cov_config.good())
		{
			char buf[1024];
			ifs_cov_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "cov_range")
				cov_range_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for facilities
			else if (str_param == "need_integrate_facs")
				need_integrate_facs = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "fac_with_coordinates")
				fac_with_coordinates = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "station" && need_integrate_facs)
				station = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "stations" && need_integrate_facs)
				stations_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "stations_1500" && need_integrate_facs)
				stations_1500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "cafe" && need_integrate_facs)
				cafe = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "cafes" && need_integrate_facs)
				cafes_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "cafes_1500" && need_integrate_facs)
				cafes_1500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "logistic" && need_integrate_facs)
				logistic = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "logistics" && need_integrate_facs)
				logistics_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "logistics_2500" && need_integrate_facs)
				logistics_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "bank" && need_integrate_facs)
				bank = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "banks" && need_integrate_facs)
				banks_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "banks_2500" && need_integrate_facs)
				banks_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "school" && need_integrate_facs)
				school = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "schools" && need_integrate_facs)
				schools_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "schools_2500" && need_integrate_facs)
				schools_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for candidates
			else if (str_param == "need_gen_cands")
				need_gen_cands = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "cands_count_1" && need_gen_cands)
			{
				str_cands_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_1 = atoi(str_cands_count_1.c_str());
			}
			else if (str_param == "cands_count_2" && need_gen_cands)
			{
				str_cands_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_2 = atoi(str_cands_count_2.c_str());
			}
			else if (str_param == "cands_count_3" && need_gen_cands)
			{
				str_cands_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_3 = atoi(str_cands_count_3.c_str());
			}
			else if (str_param == "cands_count_4" && need_gen_cands)
			{
				str_cands_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_4 = atoi(str_cands_count_4.c_str());
			}
			else if (str_param == "cands_count_5" && need_gen_cands)
			{
				str_cands_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				cands_count_5 = atoi(str_cands_count_5.c_str());
			}
			else if (str_param == "gen_cands" && need_gen_cands)
				gen_cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for reference locations (clients)
			else if (str_param == "need_gen_ref_loc")
				need_gen_ref_loc = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "clients_count_1" && need_gen_ref_loc)
			{
				str_clients_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_1 = atoi(str_clients_count_1.c_str());
			}
			else if (str_param == "clients_count_2" && need_gen_ref_loc)
			{
				str_clients_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_2 = atoi(str_clients_count_2.c_str());
			}
			else if (str_param == "clients_count_3" && need_gen_ref_loc)
			{
				str_clients_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_3 = atoi(str_clients_count_3.c_str());
			}
			else if (str_param == "clients_count_4" && need_gen_ref_loc)
			{
				str_clients_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_4 = atoi(str_clients_count_4.c_str());
			}
			else if (str_param == "clients_count_5" && need_gen_ref_loc)
			{
				str_clients_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				clients_count_5 = atoi(str_clients_count_5.c_str());
			}
			else if (str_param == "ref_loc_points" && need_gen_ref_loc)
				ref_loc_points = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "gen_random_ref_locs" && need_gen_ref_loc)
				gen_random_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "gen_avg_ref_locs" && need_gen_ref_loc)
				gen_avg_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for edges
			else if (str_param == "need_gen_edges")
				need_gen_edges = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "gen_edges")
				gen_edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for background datasets
			else if (str_param == "edges")
				edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "geos")
				geos_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
		}
	}
	else
		return false;

	if (need_integrate_facs || need_gen_cands || need_gen_ref_loc || need_gen_edges)
	{
		raw_edges_ext edges;
		seg_of_edge subsegments; // sub-segments of edges
		geo_subsegment_rtree rtree; // sub-segments with edge id and sub-segment id
		if (!load_edges_and_geos(edges, subsegments, rtree, cov_range_file_path.c_str(), edges_file_path.c_str(), geos_file_path.c_str()))
			return false;

		// just output edges, whose distances are already accumulated by sub-segments
		if (need_gen_edges
			&& !gen_edges(edges, gen_edges_file_path.c_str()))
			return false;

		// load and pick pois
		boost::unordered_set<fac_or_cand_loc> picked_facs; // locations of picked facilities
		if (need_integrate_facs)
		{
			//boost::unordered_set<fac_or_cand_loc> picked_stations; // locations of picked stations
			//if (!load_and_pick_pois(edges, subsegments, rtree, stations_file_path.c_str(), stations_1500_file_path.c_str(), station, picked_stations, fac_with_coordinates))
			//	return false;
			boost::unordered_set<fac_or_cand_loc> picked_cafes; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, cafes_file_path.c_str(), cafes_1500_file_path.c_str(), cafe, picked_cafes, fac_with_coordinates))
				return false;
			boost::unordered_set<fac_or_cand_loc> picked_logistics; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, logistics_file_path.c_str(), logistics_2500_file_path.c_str(), logistic, picked_logistics, fac_with_coordinates))
				return false;
			boost::unordered_set<fac_or_cand_loc> picked_banks; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, banks_file_path.c_str(), banks_2500_file_path.c_str(), bank, picked_banks, fac_with_coordinates))
				return false;

			/*if (!load_and_pick_pois(edges, subsegments, rtree, stations_file_path.c_str(), stations_1500_file_path.c_str(), station, picked_stations))
				return false;

			boost::unordered_set<fac_or_cand_loc> picked_cafes; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, cafes_file_path.c_str(), cafes_1500_file_path.c_str(), cafe, picked_cafes))
				return false;

			boost::unordered_set<fac_or_cand_loc> picked_logistics; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, logistics_file_path.c_str(), logistics_2500_file_path.c_str(), logistic, picked_logistics))
				return false;

			boost::unordered_set<fac_or_cand_loc> picked_banks; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, banks_file_path.c_str(), banks_2500_file_path.c_str(), bank, picked_banks))
				return false;

			boost::unordered_set<fac_or_cand_loc> picked_schools; // locations of picked stations
			if (!load_and_pick_pois(edges, subsegments, rtree, schools_file_path.c_str(), schools_2500_file_path.c_str(), school, picked_schools))
				return false;

			// construct for all types of facilities
			picked_facs.insert(picked_stations.begin(), picked_stations.end());
			picked_facs.insert(picked_cafes.begin(), picked_cafes.end());
			picked_facs.insert(picked_logistics.begin(), picked_logistics.end());
			picked_facs.insert(picked_banks.begin(), picked_banks.end());
			picked_facs.insert(picked_schools.begin(), picked_schools.end());*/
		}

		// generate candidates
		if (need_gen_cands)
		{
			if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_1 + ".txt").c_str(), cands_count_1))
				return false;
			if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_2 + ".txt").c_str(), cands_count_2))
				return false;
			if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_3 + ".txt").c_str(), cands_count_3))
				return false;
			if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_4 + ".txt").c_str(), cands_count_4))
				return false;
			if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_5 + ".txt").c_str(), cands_count_5))
				return false;
		}

		// generate reference locations
		if (need_gen_ref_loc)
		{
			if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(), clients_count_1, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(), clients_count_2, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(), clients_count_3, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(), clients_count_4, ref_loc_points.c_str()))
				return false;
			if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(),
				std::string(gen_avg_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(), clients_count_5, ref_loc_points.c_str()))
				return false;
		}
	}

	return true;
}

// integrate real BJ and CA datasets
bool real_datasets(_TCHAR* argv)
{
	// read real_config file for integrate datasets files
	int BJ_user_count, CA_user_count;
	std::string BJ_edges, BJ_geos, BJ_refloc_pre, BJ_reflocs;
	std::string CA_edges, CA_vertices, CA_refloc_pre, CA_reflocs;

	std::ifstream ifs_real_config(argv);
	if (!ifs_real_config.fail())	// real_config file exists.
	{
		while (!ifs_real_config.bad() && ifs_real_config.good())
		{
			char buf[1024];
			ifs_real_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			//--------------------------------------------------------------------------------
			// for BJ
			if (str_param == "BJ_user_count")
				BJ_user_count = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "BJ_edges")
				BJ_edges = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "BJ_geos")
				BJ_geos = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "BJ_refloc_pre")
				BJ_refloc_pre = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "BJ_reflocs")
				BJ_reflocs = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			//--------------------------------------------------------------------------------
			// for CA
			else if (str_param == "CA_user_count")
				CA_user_count = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "CA_edges")
				CA_edges = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "CA_vertices")
				CA_vertices = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "CA_refloc_pre")
				CA_refloc_pre = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "CA_reflocs")
				CA_reflocs = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
		}
	}
	else
		return false;

	//--------------------------------------------------------------------------------
	// for CA real dataset
	raw_vertices vertices_CA;
	raw_edges edges_CA;
	geo_edge_rtree rtree_CA;
	std::cout << "load_vertices_and_edges is running ..." << std::endl;
	if (!load_vertices_and_edges(vertices_CA, edges_CA, rtree_CA, NULL, CA_vertices.c_str(), CA_edges.c_str(), false)) // NULL means geo-range file is ignored
		return false;
	if (!load_and_snap_reflocs_CA(edges_CA, rtree_CA, CA_user_count, CA_refloc_pre.c_str(), CA_reflocs.c_str()))
		return false;

	//--------------------------------------------------------------------------------
	// for BJ real dataset
	raw_edges_ext edges_BJ;
	seg_of_edge subsegments_BJ; // sub-segments of edges
	geo_subsegment_rtree rtree_BJ; // sub-segments with edge id and sub-segment id
	std::cout << "load_vertices_and_edges is running ..." << std::endl;
	if (!load_edges_and_geos(edges_BJ, subsegments_BJ, rtree_BJ, NULL, BJ_edges.c_str(), BJ_geos.c_str())) // NULL means geo-range file is ignored
		return false;
	if (!load_and_snap_reflocs_BJ(edges_BJ, subsegments_BJ, rtree_BJ, BJ_user_count, BJ_refloc_pre.c_str(), BJ_reflocs.c_str()))
		return false;

	return true;
}