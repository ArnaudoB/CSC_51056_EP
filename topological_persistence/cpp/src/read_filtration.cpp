#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <limits>

using namespace std;

typedef set<int> vertices;


struct simplex{
  int dim;
  float val;
  vertices vert;
};


vector<simplex> read_filtration(string name){

  vector<simplex> F; F.clear();
  char input_file[100];
  sprintf(input_file, "%s", (char*) name.c_str());
  ifstream input(input_file);

  if (input){
    string line;
    while(getline(input,line)){
      stringstream stream(line);
      simplex s; s.vert.clear();
      s.dim = -1; stream >> s.val; stream >> s.dim; int i = 0;
      while(i <= s.dim){
        int f; stream >> f;
	s.vert.insert(f); i++;
      }
      if(s.dim != -1)
        F.push_back(s);
    }
  }
  else{cout << "Failed to read file " << name << endl;}
  std::sort(F.begin(), F.end(), [](const simplex& a, const simplex& b) {
    return std::tie(a.val, a.vert, a.dim) < std::tie(b.val, b.vert, b.dim);});
  return F;

};


// int main(int argc, char** argv) {

//   if (argc != 2) {
//     cout << "Syntax: " << argv[0] << " <filtration_file>" << endl;
//     return 0;
//   }


//   string name = argv[1];

//   cout << "Reading filtration..." << endl;
//   vector<simplex> F = read_filtration(name);
//   cout << "Done." << endl;

//   for (vector<simplex>::iterator sit = F.begin(); sit != F.end(); sit++) {
//     cout << "{val=" << sit->val << "; dim=" 
// 	 << sit->dim << "; [";
//     for (vertices::iterator vit = sit->vert.begin();
// 	 vit != sit->vert.end(); vit++){
//       cout << *vit;
//       vertices::iterator vvit = vit; vvit++; 
//       if (vvit != sit->vert.end())
// 	cout << ", ";
//     }
//     cout << "]}" << endl;
//   }
// return 0;
// }