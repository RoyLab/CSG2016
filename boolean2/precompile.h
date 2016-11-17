#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#undef min
#undef max

#include <fstream>
#include <cassert>

#define XR_DEBUG
//#define USE_CGAL_PREDICATES
#define USE_CGAL_PREDICATES_CHECK
#define PREP_DEBUG_INFO
#define MAX_MESH_COUNT 10000

//#define NDEBUG

#include "global.h"
#include "csgdefs.h"
#include "xgeometry.h"
#include "CGALext.h"
