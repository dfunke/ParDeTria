#pragma once

// fix bug in TBB where emplaced_back is only defined when LIBCPP is defined
#ifndef _LIBCPP_VERSION
#define _LIBCPP_VERSION 1
#define _UNDEF_LIBCPP 1
#endif

#include <tbb/concurrent_vector.h>

#ifdef _UNDEF_LIBCPP
#undef _LIBCPP_VERSION
#undef _UNDEF_LIBCPP
#endif

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>