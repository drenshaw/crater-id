#include <string>
#include <vector>
#include <tuple>

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

struct martian_crater {
    std::string crater_id;
	float  lat;
	float  lon;
	float  diam;
    float  ecc;
    uint   npts;
    std::string crater_name;
};

struct lunar_crater {
    std::string crater_id;
	float  lat;
	float  lon;
	float  diam;
	float  ecc;
	float  ell;
	float  arc;
	uint   npts;
};

bool readLunarCraterEntry(  const std::string, 
                            lunar_crater&,
                            const char, 
                            const float = 1.2, 
                            const float = 0.9,
                            const float = 60.0,
                            const float = 1000.0);

std::vector<lunar_crater> runCraterReader(const std::string,
                                          const char =',',
                                          const float =1.2,
                                          const float =0.9,
                                          const uint =10);

std::tuple<float, char> reinterpret_lat(const float);
std::tuple<float, char> reinterpret_lon(const float);
std::tuple<std::string, std::string> reinterpret_latlon(const float, const float);
void printLunarCratersInfo(const lunar_crater);
bool readLunarCraterEntry(const std::string, const char);
template <typename T>
void Permutation(std::vector<T>);
template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T>, const int);
std::ostream& operator<<(std::ostream&, const lunar_crater&);
