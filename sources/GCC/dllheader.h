#ifndef DLLHEADER_H_INCLUDED 
#define DLLHEADER_H_INCLUDED

#ifdef DLL_EXPORT
# define EXPORT __declspec (dllexport)
#else
# define EXPORT
#endif

extern EXPORT void foo ();

#endif
