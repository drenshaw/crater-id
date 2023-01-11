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

typedef struct martian_crater {
    std::string crater_id;
	float  lat;
	float  lon;
	float  diam;
    float  ecc;
    uint   npts;
    std::string crater_name;
} martian_crater;

typedef struct lunar_crater {
    std::string crater_id;
	float  lat;
	float  lon;
	float  diam;
	float  ecc;
	float  ell;
	float  arc;
	uint   npts;
} lunar_crater;

typedef struct Point {
    float x;
    float y;
    float z;
} Point;

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
                                          const uint =100);
template <typename T>
std::tuple<float, char> stringify_lat(const T crater);
std::tuple<float, char> stringify_lat(const float);
template <typename T>
std::tuple<float, char> stringify_lon(const T crater);
std::tuple<float, char> stringify_lon(const float);
std::string stringify_latlon(const float, const float);
template <typename T>
std::string stringify_latlon(const T crater);
void printLunarCratersInfo(const lunar_crater);
bool readLunarCraterEntry(const std::string, const char);
template <typename T>
void Permutation(std::vector<T>);
template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T>, const uint);
template <typename T, typename N>
std::vector<std::tuple<T, T>> specialCombination(const std::vector<T>, const N =30.);
template <typename T, typename N>
void specialCombination(const std::vector<T>, 
                        const N,
                        std::vector<std::tuple<uint, uint>>&);
std::ostream& operator<<(std::ostream&, const lunar_crater&);
std::ostream& operator<<(std::ostream&, const Point&);
template <typename T>
T deg2rad(const T);
template <typename T>
T rad2deg(const T);
template <typename T>
Point latlon2unitVector(const T lat, const T lon);
template <typename T>
Point latlon2unitVector(const T crater);
template <typename R, typename T>
Point LLHtoECEF(const R, const T);
float vectorNorm(const Point pt);
void normalizeVector(Point&);
float dot(const Point, const Point);
float angularDistance(const Point, const Point);
float angularPseudodistance(const Point, const Point);
template <typename T>
float latlon_dist(const T, const T, const T, const T);
template <typename T>
float latlon_dist(const T, const T);
