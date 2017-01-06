#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#undef min
#undef max

#include <fstream>
#include <cassert>

#ifdef _DEBUG
#define XR_DEBUG
#define PREP_DEBUG_INFO
#endif

//#define NDEBUG

#ifdef XR_PROFILE

#ifdef XR_DEBUG
#undef XR_DEBUG
#endif

#ifdef PREP_DEBUG_INFO
#undef PREP_DEBUG_INFO
#endif

#ifndef NDEBUG
#define NDEBUG
#endif

#endif

#define MAX_MESH_COUNT 10000


//#include <vld.h>


#include "global.h"
#include "csgdefs.h"
#include "CGALext.h"
