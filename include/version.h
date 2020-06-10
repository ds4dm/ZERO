#include <string>
#define EPECSOLVE_VERSION "2.0.1"
#define EPECSOLVE_VERSION_MAJOR "2"
#define EPECSOLVE_VERSION_MINOR "0"
#define EPECSOLVE_VERSION_PATCH "1"

class EPECVersion {
public:
  EPECVersion(std::string &major, std::string &minor, std::string &patch) {
    major = EPECSOLVE_VERSION_MAJOR;
    minor = EPECSOLVE_VERSION_MINOR;
    patch = EPECSOLVE_VERSION_PATCH;
  }
  EPECVersion(std::string &version) { version = EPECSOLVE_VERSION; }
};
