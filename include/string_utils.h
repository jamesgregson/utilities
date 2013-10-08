#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include<string>
#include<sstream>
#include<cstdio>
#include<cstdarg>

#define STR_FORMAT_BUFFER_SIZE 2048

std::string str_format( const char *fmt, ... ){
    char buffer[STR_FORMAT_BUFFER_SIZE];
    va_list arg;
    va_start( arg, fmt );
    vsprintf( buffer, fmt, arg );
    va_end(arg);
    return std::string(buffer);
}

template<typename T>
std::string str( const T &in ){
    std::stringstream ss;
    ss << in;
    return ss.str();
}

#endif
