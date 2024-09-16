#include <string>
#include <iomanip>

#include "io.h"
#include "structs.h"

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater) {
  // return os << "Crater: " << crater.crater_id;
  std::string latlon = io::stringify_latlon(crater);
  return os 
        << "ID: " << crater.crater_id
        << "  (" << latlon << ")"
        << "\tdiam: " << crater.diam << "km "
        << "\tell: " << crater.ell
        ;
}

namespace io {

bool readLunarCraterEntry(lunar_crater& crater,
                          const std::string entry, 
                          const char sep,
                          const double max_ell,
                          const double min_arc,
                          const double min_diam,
                          const double max_diam) {
  // lunar_crater crater;
  std::string lat, lon, diam, ecc, n_pts;
  std::istringstream str(entry);
  std::string token;
  uint count = 0;
  while(getline(str, token, sep)) {
    if(token.empty())
        return false;
    switch(count) {
      case CRATER_ID:
        crater.crater_id = token;
        break;
      case LAT_CIRC_IMG:
        crater.lat = stof(token);
        break;
      case LON_CIRC_IMG:
        crater.lon = stof(token);
        break;
      case DIAM_CIRC_IMG:
        crater.diam = stof(token);
        break;
      case DIAM_ELLI_ECCEN_IMG:
        crater.ecc = stof(token);
        break;
      case DIAM_ELLI_ELLIP_IMG:
        crater.ell = stof(token);
        break;
      case ARC_IMG:
        crater.arc = stof(token);
        break;
      default:
        break;
    }
    count++;
  }
  return crater.ell < max_ell 
      && crater.arc > min_arc
      && crater.diam < max_diam
      && crater.diam > min_diam;
}

bool readMartianCraterEntry(const std::string entry) {
  std::cout << "Reading Mars crater database: " << entry << std::endl;
  std::cerr << "Martian crater database parsing is not complete!\n";
  return false;
}

std::string stringify_lat(const double lat) {
  if(lat > 90. || lat < -90) {
    std::cerr << "Latitude falls outside [-90°,90°]: " << lat << std::endl;
    return "";
  }
  std::string n_or_s;
  double abs_lat;
  n_or_s = (lat > 0)?"N":"S";
  abs_lat = std::abs(lat);
  std::ostringstream outStream;
  outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lat << "°" << n_or_s;
  std::string out_str = outStream.str();
  return out_str;
}

std::string stringify_lon(const double lon) {
  double re_lon = lon;
  if(re_lon > 180)
    re_lon -= 360;
  if(re_lon > 180. || re_lon < -180) {
    std::cerr << "Longitude falls outside [-180°,180°]: " << lon << std::endl;
    return "";
  }
  double abs_lon = std::abs(re_lon);
  std::string e_or_w;
  e_or_w = (re_lon < 0)?"W":"E";
  std::ostringstream outStream;
  outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lon << "°" << e_or_w;
  std::string out_str = outStream.str();
  return out_str;
}

std::string stringify_latlon(const double lat, const double lon) {
  std::string str_lat, str_lon;
  str_lat = stringify_lat(lat);
  str_lon = stringify_lon(lon);  
  return (str_lat + " " + str_lon);
}
std::string stringify_lat(const lunar_crater crater) {
  return stringify_lat(crater.lat);
}

std::string stringify_lon(const lunar_crater crater) {
  return stringify_lon(crater.lon);
}

std::string stringify_latlon(const lunar_crater crater) {
  return stringify_latlon(crater.lat, crater.lon);
}

void printLunarCratersInfo(const lunar_crater crater) {
  // const char cout_sep = '\n';
  const std::string cout_sep = "\n\t";

  std::string latlon = stringify_latlon(crater);
  std::cout 
        << "CRATER_ID: " << crater.crater_id << cout_sep
        << "LAT_CIRC_IMG/LON_CIRC_IMG: " << latlon << cout_sep 
        << "DIAM_CIRC_IMG:       " << crater.diam << "km" << cout_sep 
        << "DIAM_ELLI_ECCEN_IMG: " << crater.ecc << cout_sep 
        << "DIAM_ELLI_ELLIP_IMG: " << crater.ell << cout_sep 
        << "ARC_IMG:             " << crater.arc << cout_sep
        << std::endl;
}

} // namespace
