
#ifndef __FDELAUN2D_
#define  __FDELAUN2D_

#ifdef __cplusplus
extern "C" {
#endif

void triangulate(void* kbd, void* ktj, void* kcm, 
	void* nex, void* ibex, void* ibno, void* nidm, void* nbk, void* ibreak, void* nbreak, void* nob, void* nib, 
	void* px, void* py, void* pd, void* dpp, void* stl, void* node, void* nelm, void* mtj, void* jac, void* idm,void* ierrcode);

#ifdef __cplusplus
}
#endif
	
#endif
