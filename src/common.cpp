#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>
#include "common.h"

using namespace std;

const set<char> delims{'/'};

string get_env_var(const string& key) {
    char* val = getenv(key.c_str());
    return val == NULL ? string("") : string(val);
}

bool is_not_empty(string& s) {
    return !s.empty();
}

string get_basename(const string& str) {
    auto parts = split_path(str);
    auto basename = find_if(parts.rbegin(), parts.rend(), is_not_empty);
    return (*basename);
}

vector<string> split_path(const string& str) {
    vector<string> result;

    char const* pch = str.c_str();
    char const* start = pch;

    for(; *pch; ++pch) {
        if (delims.find(*pch) != delims.end()) {
            if (start != pch) {
                string str(start, pch);
                result.push_back(str);
            } else {
                result.push_back("");
            }
            start = pch + 1;
        }
    }

    result.push_back(start);

    return result;
}
