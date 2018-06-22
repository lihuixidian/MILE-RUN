#pragma once

#include <utility> // std::pair

#include <boost/geometry.hpp> // boost::geometry

//--------------------------------------------------------------------------------
// geo data types

enum geo { LON = 0, LAT };

typedef boost::geometry::model::point<float, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree> > geo_point;
typedef std::pair < geo_point, // candidate
	int > // candidate id
geo_cand;
typedef boost::geometry::index::rtree<geo_cand, boost::geometry::index::rstar<16>> geo_cand_rtree; // an candidates R*-tree
typedef boost::geometry::model::segment<geo_point> geo_segment;
typedef std::pair < geo_segment, // edge segment
	int > // edge id
geo_edge;
typedef boost::geometry::index::rtree<geo_edge, boost::geometry::index::rstar<16>> geo_edge_rtree; // an edges R*-tree
typedef boost::geometry::model::box<geo_point> geo_box;

//--------------------------------------------------------------------------------
// define static const values

static const float EARTH_RADIUS = 6378.137f; // km
static const float PI = 3.14159f; // pi