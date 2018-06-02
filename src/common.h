#pragma once
#include <string>
#include <vector>
#include <set>

using namespace std;

typedef unsigned char uchar;

string get_env_var(const string& key);
string get_basename(const string& path);
vector<string> split_path(const string& str);
void log(const char* fmt, ...);

// Constante utilizada para hacer igualdad aproximada por error de precisión
const double eps = 1e-5;
