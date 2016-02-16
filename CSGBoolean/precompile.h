#pragma once

#pragma warning(disable: 4005)

#define SAFE_RELEASE(ptr) if (ptr) {delete ptr; ptr = 0;}
#define SAFE_RELEASE_ARRAY(ptr) if (ptr) {delete    [] ptr; ptr = 0;}
#define CSG_EXPORTS

#ifndef _DEBUG
#define _USE_MATH_DEFINES
#define NDEBUG
#else
//#define NDEBUG
extern int countd1, countd2, countd3, countd4, countd5;
#endif

