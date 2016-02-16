#pragma once

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

namespace myext
{
    // trim from start
    static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
    }

    static inline void tolower(const std::string &s, std::string &r) {
        std::transform(s.begin(), s.end(), r.begin(), std::tolower);
    }

    static inline int strcmp_nocase(const std::string &a, const std::string &b){
        std::string a1, b1;
        tolower(a, a1);
        tolower(b, b1);
        return std::strcmp(a1.c_str(), b1.c_str());
    }
}