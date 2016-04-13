#pragma once
#include <assert.h>
#include <string.h>
#include <string>


#ifdef NOT
#undef NOT
#endif

#define NOT(h) !(h)

#if defined(DEBUG) || defined(_DEBUG)
#ifndef V
#define V(x)           { hr = (x); if( NOT(hr) ) { assert( 0 || __FILE__); } }
#endif
#ifndef V_RETURN
#define V_RETURN(x)    { hr = (x); if( NOT(hr) ) { assert( 0 || __FILE__); return hr; } }
#endif
#else
#ifndef V
#define V(x)           { hr = (x); }
#endif
#ifndef V_RETURN
#define V_RETURN(x)    { hr = (x); if( NOT(hr) ) { return hr; } }
#endif
#endif

#ifndef SAFE_DELETE
#define SAFE_DELETE(p)       { if (p) { delete (p);     (p) = nullptr; } }
#endif
#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p) { if (p) { delete[] (p);   (p) = nullptr; } }
#endif
#ifndef SAFE_RELEASE
#define SAFE_RELEASE(p)      { if (p) { (p)->Release(); (p) = nullptr; } }
#endif

#ifndef COMMON_PROPERTY
#define COMMON_PROPERTY(__type__,__name__) \
protected: \
__type__ m_##__name__; \
public: \
const __type__ & get_##__name__() const{ \
return m_##__name__; \
} \
__type__ & get_##__name__(){ \
return m_##__name__; \
} \
void set_##__name__##(__type__ & _##__name__##_){ \
m_##__name__ = _##__name__##_; \
}
#endif

#ifndef COMMON_PROPERTY_POINTER
#define COMMON_PROPERTY_POINTER(__type__,__name__) \
protected: \
__type__ * mp_##__name__ = nullptr; \
public: \
const __type__ * get_##__name__() const{ \
return mp_##__name__; \
} \
__type__ * get_##__name__(){ \
return mp_##__name__; \
} \
void set_##__name__##(__type__ * _##__name__##_){ \
mp_##__name__ = _##__name__##_; \
}
#endif

#ifndef STATIC_PROPERTY
#define STATIC_PROPERTY(__type__,__name__) \
protected: \
static __type__ __name__; \
public: \
static __type__ & get_##__name__(){ \
return __name__; \
} \
static void set_##__name__##(const __type__ & _##__name__##_){ \
__name__ = _##__name__##_; \
}
#endif


#ifndef STATIC_PROPERTY_POINTER
#define STATIC_PROPERTY_POINTER(__type__,__name__) \
protected: \
static __type__ * __name__; \
public: \
static __type__ * get_##__name__(){ \
return __name__; \
} \
static void set_##__name__##(const __type__ * _##__name__##_){ \
__name__ = _##__name__##_; \
}
#endif

#ifndef ReportErrorString
#define ReportErrorString\
    (std::string("Error! ") + itoa(__LINE__, _errline, 20) + "of file \"" + __FILE__ + "\".\n")

#endif

#ifndef ReportError
#define ReportError(text)\
    {printf ("Error! %d of file \"%s\": %s\n",__LINE__, __FILE__, text);assert(0);}

#endif
