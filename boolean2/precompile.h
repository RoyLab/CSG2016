#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#undef min
#undef max

#include <fstream>
#include <cassert>

//#define XR_DEBUG
//#define USE_CGAL_PREDICATES
//#define USE_CGAL_PREDICATES_CHECK
//#define PREP_DEBUG_INFO

//#define NDEBUG

#ifdef XR_PROFILE

#ifdef XR_DEBUG
#undef XR_DEBUG
#endif

#ifdef USE_CGAL_PREDICATES
#undef USE_CGAL_PREDICATES
#endif

#ifdef USE_CGAL_PREDICATES_CHECK
#undef USE_CGAL_PREDICATES_CHECK
#endif

#ifdef PREP_DEBUG_INFO
#undef PREP_DEBUG_INFO
#endif

#ifndef NDEBUG
#define NDEBUG
#endif

#endif

#define MAX_MESH_COUNT 10000


#include "global.h"
#include "csgdefs.h"
#include "xgeometry.h"
#include "CGALext.h"
