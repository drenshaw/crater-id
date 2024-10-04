#pragma once

#include <string>

enum columns {
  CRATER_ID,
  LAT_CIRC_IMG,
  LON_CIRC_IMG,
  LAT_ELLI_IMG,
  LON_ELLI_IMG,
  DIAM_CIRC_IMG,
  DIAM_CIRC_SD_IMG,
  DIAM_ELLI_MAJOR_IMG,
  DIAM_ELLI_MINOR_IMG,
  DIAM_ELLI_ECCEN_IMG,
  DIAM_ELLI_ELLIP_IMG,
  DIAM_ELLI_ANGLE_IMG,
  LAT_ELLI_SD_IMG,
  LON_ELLI_SD_IMG,
  DIAM_ELLI_MAJOR_SD_IMG,
  DIAM_ELLI_MINOR_SD_IMG,
  DIAM_ELLI_ANGLE_SD_IMG,
  DIAM_ELLI_ECCEN_SD_IMG,
  DIAM_ELLI_ELLIP_SD_IMG,
  ARC_IMG,
  PTS_RIM_IMG
};

typedef struct martian_crater {
  std::string crater_id;
	double  lat;
	double  lon;
	double  diam;
  double  ecc;
  uint    npts;
  std::string crater_name;
} martian_crater;

typedef struct lunar_crater {
  std::string crater_id;
	double  lat;
	double  lon;
	double  diam;
	double  ecc;
	double  ell;
	double  arc;
	uint    npts;
} lunar_crater;

typedef struct Point {
  double x;
  double y;
  double z;
} Point;
