#pragma once

#include <map>
#include <set>
#include <fstream>


// Load real or query data
// remark: (currently for simplicity) the useful object must be the first word in a line
// <T> type of data object, which must can be constructed using std::string and have operator==
// [in] file_path: file path name of data
// [out] data_map: ranked map of data
// [in] obj_index: the index of object in each line, default is 0, namely the first string
// [in] seperator: the seperator between every word in a line, default is '\n', which means there is only one word in a line

template <typename T>
void load_ranked_data(std::ifstream &in_file, std::map<unsigned, T> &data_map, unsigned obj_index = 0, char seperator = '\n')
{
	unsigned rank = 1;
	if (!in_file.fail()) // file exists
	{
		while (!in_file.bad() && in_file.good())
		{
			char buf[1024];
			in_file.getline(buf, 1024); // read a line
			std::string str_buf(buf);

			if (seperator == '\n')
				data_map[rank++] = str_buf;
			else // need to extract information from buffer string
			{
				// (currently for simplicity) the useful object must be the first word in a line
				std::string::size_type begin_pos = 0;
				std::string::size_type end_pos = str_buf.find(seperator, begin_pos);
				data_map[rank++] = str_buf.substr(begin_pos, end_pos);
			}
		}
	}
}

// Calculate Precision (P), Recall (R) and Average Precision (AP)
// <T> type of data object, which must can be constructed using std::string and have operator==
// [in] real_map: ranked map of real data
// [in] query_map: ranked map of query data
// [out] ap: AP value
// [out] p: precision value, default is NULL, which means no need for precision
// [out] r: recall value, default is NULL, which means no need for recall

template <typename T>
void APatK(const std::map<unsigned, T> &real_map, const std::map<unsigned, T> &query_map, double &ap, double *p = NULL, double *r = NULL)
{
	std::map<unsigned, unsigned> ap_map;
	unsigned hit_sum = 0;
	std::map<unsigned, T>::const_iterator iter_query = query_map.cbegin();
	for (; iter_query != query_map.cend(); ++iter_query)
	{
		std::map<unsigned, T>::const_iterator iter_real = real_map.cbegin();
		for (; iter_real != real_map.cend(); ++iter_real)
		{
			ap_map[iter_query->first] = 0;
			if (iter_real->second == iter_query->second)
			{
				hit_sum++;
				ap_map[iter_query->first] = 1;
				break;
			}
		}
	}

	if (p != NULL) // need to compute precision
		*p = static_cast<double>(hit_sum) / static_cast<double>(query_map.size());
	if (r != NULL) // need to compute recall
		*r = static_cast<double>(hit_sum) / static_cast<double>(real_map.size());

	ap = 0;
	unsigned hit = 0;
	std::map<unsigned, unsigned>::const_iterator iter_ap = ap_map.cbegin();
	for (; iter_ap != ap_map.cend(); ++iter_ap)
	{
		if (iter_ap->second == 1)
		{
			++hit;
			ap += static_cast<double>(hit) / static_cast<double>(iter_ap->first);
		}
	}
	ap /= static_cast<double>(ap_map.size());
}

