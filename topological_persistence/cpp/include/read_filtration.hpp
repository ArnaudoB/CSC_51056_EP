#pragma once
#include <string>
#include <vector>
#include <set>

typedef std::set<int> vertices;

struct simplex{
  int dim;
  float val;
  vertices vert;
};

std::vector<simplex> read_filtration(std::string name);




