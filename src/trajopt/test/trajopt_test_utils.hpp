#include <json/json.h>
#include <fstream>
#include <stdexcept>

Json::Value readJsonFile(const std::string& fname) {
  Json::Value root;
  JSONCPP_STRING errs;
  Json::CharReaderBuilder builder;
  Json::CharReader* reader(builder.newCharReader());
  std::ifstream fh(fname.c_str());
  std::string fc((std::istreambuf_iterator<char>(fh)), std::istreambuf_iterator<char>());
  bool parse_success = reader->parse(fc.c_str(), fc.c_str() + fc.size(), &root, &errs);
  delete reader;
  if (!parse_success) throw std::runtime_error("failed to parse " + fname + "\n" + errs);
  return root;
}
