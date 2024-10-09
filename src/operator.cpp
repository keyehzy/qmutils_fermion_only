#include "operator.h"

#include <sstream>

namespace qmutils {

std::string Operator::to_string() const {
  std::ostringstream oss;
  oss << (data_.type == static_cast<uint8_t>(Type::Creation) ? "c+" : "c")
      << "(" << (data_.spin == static_cast<uint8_t>(Spin::Up) ? "↑" : "↓")
      << "," << static_cast<int>(data_.orbital) << ")";
  return oss.str();
}

}  // namespace qmutils
