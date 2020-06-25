#include <string>
#define ZERO_VERSION "1.0.0"
#define ZERO_VERSION_MAJOR "1"
#define ZERO_VERSION_MINOR "0"
#define ZERO_VERSION_PATCH "0"

class ZEROVersion {
public:
  ZEROVersion(std::string &major, std::string &minor, std::string &patch) {
    major = ZERO_VERSION_MAJOR;
    minor = ZERO_VERSION_MINOR;
    patch = ZERO_VERSION_PATCH;
  }
  ZEROVersion(std::string &version) { version = ZERO_VERSION; }
};
