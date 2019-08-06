#pragma once

#include <cstddef>

#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)>(b)?(b):(a)

#define DEBUG 0

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long uint32;
typedef unsigned long long uint64;
typedef char int8;
typedef short int16;
typedef long int32;
typedef long long int64;

typedef int64 trieidx;
typedef std::size_t size_t;
