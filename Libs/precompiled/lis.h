/* Copyright (C) 2004-2007 The SSI Project. All rights reserved.
   http://ssi.is.s.u-tokyo.ac.jp/

   This software is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY.  See the LICENSE file for the precise
   terms and conditions for the use of this software.
*/

#ifndef __LIS_H__
#define __LIS_H__
/**************************************/
/* #define LIS_VERSION	"2005/07/01"  */
/* #define LIS_VERSION	"2005/07/29"  */
/* #define LIS_VERSION	"2005/09/08"  */
/* #define LIS_VERSION	"2005/09/16"  */
/* #define LIS_VERSION	"1.0.0"       */
/* #define LIS_VERSION	"1.0.1"       */
/* #define LIS_VERSION	"1.0.2"       */
/* #define LIS_VERSION	"1.1.0-beta1" */
/* #define LIS_VERSION	"1.1.0-beta2" */
#define LIS_VERSION	"1.1.0"
/**********************************/
#include <stdio.h>

#define _max(a,b) ((a) >= (b) ? (a) : (b))
#define _min(a,b) ((a) <= (b) ? (a) : (b))


#define LIS_FMT_AUTO				0
#define LIS_FMT_PLAIN				1
#define LIS_FMT_MM					2
#define LIS_FMT_LIS					3
#define LIS_FMT_LIS_ASCII			3
#define LIS_FMT_LIS_BINARY			4
#define LIS_FMT_FREE				5
#define LIS_FMT_ITBL				6
#define LIS_FMT_HB					7
#define LIS_FMT_MMB					8

#define LIS_BINARY_BIG				0
#define LIS_BINARY_LITTLE			1


#define LIS_OPTIONS_LEN					26
#define LIS_OPTIONS_SOLVER				0
#define LIS_OPTIONS_PRECON				1
#define LIS_OPTIONS_MAXITER				2
#define LIS_OPTIONS_OUTPUT				3
#define LIS_OPTIONS_RESTART				4
#define LIS_OPTIONS_ELL					5
#define LIS_OPTIONS_SCALE				6
#define LIS_OPTIONS_FILL				7
#define LIS_OPTIONS_M					8
#define LIS_OPTIONS_PSOLVER				9
#define LIS_OPTIONS_PMAXITER			10
#define LIS_OPTIONS_PRESTART			11
#define LIS_OPTIONS_PELL				12
#define LIS_OPTIONS_PPRECON				13
#define LIS_OPTIONS_ISLEVEL				14
#define LIS_OPTIONS_INITGUESS_ZEROS		15
#define LIS_OPTIONS_ADDS				16
#define LIS_OPTIONS_ADDS_ITER			17
#define LIS_OPTIONS_PRECISION			18
#define LIS_OPTIONS_USE_AT				19
#define LIS_OPTIONS_SWITCH_MAXITER		20
#define LIS_OPTIONS_SAAMG_UNSYM			21
#define LIS_OPTIONS_STORAGE				22
#define LIS_OPTIONS_STORAGE_BLOCK		23
#define LIS_OPTIONS_CONV_COND			24
#define LIS_OPTIONS_INIT_SHADOW_RESID	25

#define LIS_PARAMS_LEN				14
#define LIS_PARAMS_RESID			LIS_OPTIONS_LEN+0
#define LIS_PARAMS_W				LIS_OPTIONS_LEN+1
#define LIS_PARAMS_RELAX			LIS_OPTIONS_LEN+2
#define LIS_PARAMS_DROP				LIS_OPTIONS_LEN+3
#define LIS_PARAMS_ALPHA			LIS_OPTIONS_LEN+4
#define LIS_PARAMS_TAU				LIS_OPTIONS_LEN+5
#define LIS_PARAMS_SIGMA			LIS_OPTIONS_LEN+6
#define LIS_PARAMS_GAMMA			LIS_OPTIONS_LEN+7
#define LIS_PARAMS_SSOR_W			LIS_OPTIONS_LEN+8
#define LIS_PARAMS_PRESID			LIS_OPTIONS_LEN+9
#define LIS_PARAMS_PW				LIS_OPTIONS_LEN+10
#define LIS_PARAMS_SWITCH_RESID		LIS_OPTIONS_LEN+11
#define LIS_PARAMS_RATE				LIS_OPTIONS_LEN+12
#define LIS_PARAMS_RESID_WEIGHT		LIS_OPTIONS_LEN+13


#define LIS_OPTIONS_FILE			-1
#define LIS_OPTIONS_HELP			-2
#define LIS_OPTIONS_VER				-3


#define LIS_PRINT_NONE				0
#define LIS_PRINT_MEM				1
#define LIS_PRINT_OUT				2
#define LIS_PRINT_ALL				3

#define LIS_SCALE_NONE				0
#define LIS_SCALE_JACOBI			1
#define LIS_SCALE_SYMM_DIAG			2

#define LIS_CONV_COND_DEFAULT		0
#define LIS_CONV_COND_NRM2_R		0
#define LIS_CONV_COND_NRM2_B		1
#define LIS_CONV_COND_NRM1_B		2

#define LIS_SOLVER_LEN				20
#define LIS_SOLVER_CG				1
#define LIS_SOLVER_BICG				2
#define LIS_SOLVER_CGS				3
#define LIS_SOLVER_BICGSTAB			4
#define LIS_SOLVER_BICGSTABL		5
#define LIS_SOLVER_GPBICG			6
#define LIS_SOLVER_QMR				7
#define LIS_SOLVER_TFQMR			7
#define LIS_SOLVER_ORTHOMIN			8
#define LIS_SOLVER_GMRES			9
#define LIS_SOLVER_JACOBI			10
#define LIS_SOLVER_GS				11
#define LIS_SOLVER_SOR				12
#define LIS_SOLVER_BICGSAFE			13
#define LIS_SOLVER_CR				14
#define LIS_SOLVER_BICR				15
#define LIS_SOLVER_CRS				16
#define LIS_SOLVER_BICRSTAB			17
#define LIS_SOLVER_GPBICR			18
#define LIS_SOLVER_BICRSAFE			19
#define LIS_SOLVER_FGMRES			20

#define LIS_INS_VALUE				0
#define LIS_ADD_VALUE				1
#define LIS_SUB_VALUE				2

#define LIS_MATRIX_LOWER			0
#define LIS_MATRIX_UPPER			1
#define LIS_MATRIX_SSOR				2

#define LIS_ORIGIN_0				0
#define LIS_ORIGIN_1				1

#define LIS_RESID					0
#define LIS_RANDOM					1

#define	LIS_PRECISION_DEFAULT		0
#define LIS_PRECISION_DOUBLE		0
#define LIS_PRECISION_QUAD			1
#define LIS_PRECISION_SWITCH		2

#define LIS_LABEL_VECTOR			0
#define LIS_LABEL_MATRIX			1

#define LIS_VEC_TMP_PADD			128

#define LIS_VECTOR_NULL				-1
#define LIS_VECTOR_ASSEMBLING		0
#define LIS_VECTOR_ASSEMBLED		1

#define LIS_PRECONNAME_MAX			16
#define LIS_PRECON_REGISTER_MAX		10

#define LIS_PRECON_TYPE_LEN			12
#define LIS_PRECON_TYPE_NONE		0
#define LIS_PRECON_TYPE_JACOBI		1
#define LIS_PRECON_TYPE_ILU			2
#define LIS_PRECON_TYPE_SSOR		3
#define LIS_PRECON_TYPE_HYBRID		4
#define LIS_PRECON_TYPE_IS			5
#define LIS_PRECON_TYPE_SAI			6
#define LIS_PRECON_TYPE_SAAMG		7
#define LIS_PRECON_TYPE_ILUC		8
#define LIS_PRECON_TYPE_ILUT		9
#define LIS_PRECON_TYPE_BJACOBI		10
#define LIS_PRECON_TYPE_ADDS		11
#define LIS_PRECON_TYPE_USERDEF		LIS_PRECON_TYPE_LEN

#define LIS_MATRIX_ASSEMBLING	0
#define LIS_MATRIX_CRS			1
#define LIS_MATRIX_CSR			1
#define LIS_MATRIX_CCS			2
#define LIS_MATRIX_CSC			2
#define LIS_MATRIX_MSR			3
#define LIS_MATRIX_DIA			4
#define LIS_MATRIX_CDS			4
#define LIS_MATRIX_ELL			5
#define LIS_MATRIX_JDS			6
#define LIS_MATRIX_JAD			6
#define LIS_MATRIX_BSR			7
#define LIS_MATRIX_BSC			8
#define LIS_MATRIX_VBR			9
#define LIS_MATRIX_COO			10
#define LIS_MATRIX_DENSE		11
#define LIS_MATRIX_DNS			11
#define LIS_MATRIX_RCO			255

#define LIS_MATRIX_TJDS			12
#define LIS_MATRIX_BJDS			13
#define LIS_MATRIX_BCR			14
#define LIS_MATRIX_CJDS			15
#define LIS_MATRIX_PCRS			16
#define LIS_MATRIX_LCRS			17
#define LIS_MATRIX_LJDS			18
#define LIS_MATRIX_LBSR			19
#define LIS_MATRIX_CDIA			20
#define LIS_MATRIX_MSC			21
#define LIS_MATRIX_DECIDING_SIZE 	-(LIS_MATRIX_RCO+1)
#define LIS_MATRIX_NULL				-(LIS_MATRIX_RCO+2)

#define LIS_MATRIX_DEFAULT		LIS_MATRIX_CRS
#define LIS_MATRIX_POINT		LIS_MATRIX_CRS
#define LIS_MATRIX_BLOCK		LIS_MATRIX_BSR


#if defined(_DEBUG)
#define LIS_DEBUG_FUNC_IN	lis_debug_trace_func(1,__FUNC__)
#define LIS_DEBUG_FUNC_OUT	lis_debug_trace_func(0,__FUNC__)
#else
#define LIS_DEBUG_FUNC_IN
#define LIS_DEBUG_FUNC_OUT
#endif
/****************************************/
/****************************************/
typedef struct
{
	double	hi;
	double	lo;
} LIS_DOUBLE_DOUBLE;

typedef struct
{
	double	hi[2];
	double	lo[2];
} LIS_DOUBLE_DOUBLE_PD;

typedef struct
{
	double	*hi;
	double	*lo;
} LIS_DOUBLE_DOUBLE_PTR;

typedef double					LIS_SCALAR;
typedef double					LIS_REAL;
typedef LIS_DOUBLE_DOUBLE		LIS_QUAD;
typedef LIS_DOUBLE_DOUBLE_PD	LIS_QUAD_PD;
typedef LIS_DOUBLE_DOUBLE_PTR	LIS_QUAD_PTR;

#ifdef USE_MPI
	#include <mpi.h>
	typedef MPI_Comm			LIS_Comm;
	#define LIS_COMM_WORLD		((LIS_Comm)MPI_COMM_WORLD)
#else
	typedef int					LIS_Comm;
	#define LIS_COMM_WORLD		((LIS_Comm)0x1)
#endif

struct LIS_COMMTABLE_STRUCT
{
	LIS_Comm	comm;
	int			pad;
	int			neibpetot;
	int			imnnz;
	int			exnnz;
	int			wssize;
	int			wrsize;
	int			*neibpe;
	int			*import_ptr;
	int			*import_index;
	int			*export_ptr;
	int			*export_index;
	LIS_SCALAR	*ws;
	LIS_SCALAR	*wr;
#ifdef USE_MPI
	MPI_Request *req1,*req2;
	MPI_Status	*sta1,*sta2;
#endif
};

typedef struct LIS_COMMTABLE_STRUCT *LIS_COMMTABLE;

struct LIS_VECTOR_STRUCT
{
	int			label;
	int			status;
	int			precision;
	int			gn;
	int			n;
	int			np;
	int			pad;
	int			origin;
	int			is_copy;
	int			is_destroy;
	int			is_scaled;
	int			my_rank;
	int			nprocs;
	LIS_Comm	comm;
	int			is;
	int			ie;
	int			*ranges;
	LIS_SCALAR	*value;
	LIS_SCALAR	*value_lo;
	LIS_SCALAR	*work;
};
typedef struct LIS_VECTOR_STRUCT *LIS_VECTOR;

struct LIS_VECTOR_S_STRUCT
{
	int			label;
	int			status;
	int			precision;
	int			gn;
	int			n;
	int			np;
	int			pad;
	int			origin;
	int			is_copy;
	int			is_destroy;
	int			is_scaled;
	int			my_rank;
	int			nprocs;
	LIS_Comm	comm;
	int			is;
	int			ie;
	int			*ranges;
	LIS_SCALAR	*value;
	LIS_SCALAR	*value_lo;
	LIS_QUAD_PTR	value_hl;
	LIS_SCALAR	*work;
	int			*index;
	int			nnz;
};
typedef struct LIS_VECTOR_S_STRUCT *LIS_VECTOR_S;

#define LIS_MATRIX_OPTION_LEN		10

struct LIS_MATRIX_CORE_STRUCT
{
	int			nnz;
	int			ndz;
	int			bnr;
	int			bnc;
	int			nr;
	int			nc;
	int			bnnz;
	int			nnd;
	int			maxnzr;
	int			*ptr;
	int			*row;
	int			*col;
	int			*index;
	int			*bptr;
	int			*bindex;
	LIS_SCALAR	*value;
	LIS_SCALAR	*work;
};
typedef struct LIS_MATRIX_CORE_STRUCT *LIS_MATRIX_CORE;

struct LIS_MATRIX_DIAG_STRUCT
{
	int			label;
	int			status;
	int			precision;
	int			gn;
	int			n;
	int			np;
	int			pad;
	int			origin;
	int			is_copy;
	int			is_destroy;
	int			is_scaled;
	int			my_rank;
	int			nprocs;
	LIS_Comm	comm;
	int			is;
	int			ie;
	int			*ranges;
	LIS_SCALAR	*value;
	LIS_SCALAR	*work;

	int			bn;
	int			nr;
	int			*bns;
	int			*ptr;
	LIS_SCALAR	**v_value;
};
typedef struct LIS_MATRIX_DIAG_STRUCT *LIS_MATRIX_DIAG;

struct LIS_MATRIX_STRUCT
{
	int			label;
	int			status;
	int			precision;
	int			gn;
	int			n;
	int			np;
	int			pad;
	int			origin;
	int			is_copy;
	int			is_destroy;
	int			is_scaled;
	int			my_rank;
	int			nprocs;
	LIS_Comm	comm;
	int			is;
	int			ie;
	int			*ranges;

	int			matrix_type;
	int			nnz;			/* CRS,CCS,MSR,JDS,VBR,COO */
	int			ndz;			/* MSR */
	int			bnr;			/* BSR,BSC */
	int			bnc;			/* BSR,BSC */
	int			nr;				/* BSR,BSC,VBR */
	int			nc;				/* BSR,BSC,VBR */
	int			bnnz;			/* BSR,BSC,VBR */
	int			nnd;			/* DIA */
	int			maxnzr;			/* ELL,JDS */
	int			*ptr;			/* CRS,CCS,JDS */
	int			*row;			/* JDS,VBR,COO */
	int			*col;			/* JDS,VBR,COO */
	int			*index;			/* CRS,CCS,MSR,DIA,ELL,JDS */
	int			*bptr;			/* BSR,BSC,VBR */
	int			*bindex;		/* BSR,BSC,VBR */
	LIS_SCALAR	*value;			/* CRS,CCS,MSR,DIA,ELL,JDS,BSR,BSC,VBR,DNS,COO */
	LIS_SCALAR	*work;

	LIS_MATRIX_CORE	L;
	LIS_MATRIX_CORE	U;
	LIS_MATRIX_DIAG	D;
	LIS_MATRIX_DIAG WD;

	int			is_block;
	int			pad_comm;
	int			is_pmat;
	int			is_sorted;
	int			is_splited;
	int			is_save;
	int			is_comm;
	int			use_wd;
	int			conv_bnr;
	int			conv_bnc;
	int			*conv_row;
	int			*conv_col;
	int			options[LIS_MATRIX_OPTION_LEN];

	int			w_annz;
	int			*w_nnz;
	int			*w_row;
	int			**w_index;
	LIS_SCALAR	**w_value;
	LIS_SCALAR	***v_value;

	int			*l2g_map;
	LIS_COMMTABLE	commtable;
};
typedef struct LIS_MATRIX_STRUCT *LIS_MATRIX;


struct LIS_MATRIX_ILU_STRUCT
{
	int			n;
	int			bs;
	int			*nnz_ma;
	int			*nnz;
	int			*bsz;
	int			**index;
	LIS_SCALAR	**value;
	LIS_SCALAR	***values;
};
typedef struct LIS_MATRIX_ILU_STRUCT *LIS_MATRIX_ILU;

struct LIS_PRECON_STRUCT
{
	int			precon_type;
	LIS_MATRIX	A;							/* SSOR */
	LIS_MATRIX	At;
	LIS_MATRIX_ILU	L;						/* ilu(k),ilut,iluc,sainv */
	LIS_MATRIX_ILU	U;						/* ilu(k),ilut,iluc,sainv */
	LIS_MATRIX_DIAG WD;						/* bilu(k),bilut,biluc,bjacobi */
	LIS_VECTOR	D;							/* ilu(k),ilut,iluc,jacobi,sainv */
	LIS_VECTOR	Pb;							/* i+s */
	LIS_VECTOR	temp;						/* saamg */
	LIS_VECTOR	*work;						/* adds */
	struct LIS_SOLVER_STRUCT	*solver;	/* hybrid */
	int			worklen;					/* adds */
	int			level_num;					/* saamg */
	int			wsize;						/* saamg */
	int			solver_comm;				/* saamg */
	int			my_rank;					/* saamg */
	int			nprocs;						/* saamg */
	int			is_copy;
	LIS_COMMTABLE	commtable;				/* saamg */
};
typedef struct LIS_PRECON_STRUCT *LIS_PRECON;


struct LIS_SOLVER_STRUCT
{
	LIS_MATRIX	A,At;
	LIS_VECTOR	b,x,xx,d;
	LIS_MATRIX_DIAG WD;
	LIS_PRECON	precon;
	LIS_VECTOR	*work;
	LIS_SCALAR	*residual;
	int			worklen;
	int			options[LIS_OPTIONS_LEN];
	LIS_SCALAR	params[LIS_PARAMS_LEN];
	int			retcode;
	int			iter;
	int			iter2;
	LIS_SCALAR	resid;
	double		times;
	double		itimes;
	double		ptimes;
	double		p_c_times;
	double		p_i_times;
	int			precision;
	LIS_REAL	bnrm;
	LIS_REAL	tol;
	LIS_REAL	tol_switch;
};
typedef struct LIS_SOLVER_STRUCT *LIS_SOLVER;

struct LIS_CONV_OPTIONS_STRUCT
{
	int			bnr;
	int			bnc;
	int			*row;
	int			*col;
};
typedef struct LIS_CONV_OPTIONS_STRUCT LIS_CONV_OPTIONS;

typedef int (*LIS_PRECON_CREATE_XXX)(LIS_SOLVER solver, LIS_PRECON precon);
typedef int (*LIS_PSOLVE_XXX)(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
typedef int (*LIS_PSOLVET_XXX)(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);

typedef struct LIS_PRECON_REGISTER_STRUCT
{
	int						precon_type;
	char					name[LIS_PRECONNAME_MAX+1];
	LIS_PRECON_CREATE_XXX	pcreate;
	LIS_PSOLVE_XXX			psolve;
	LIS_PSOLVET_XXX			psolvet;
} LIS_PRECON_REGISTER;


#ifdef __cplusplus
extern "C"
{
#endif

/****************************/
/* Vector                   */
/****************************/
	extern int lis_vector_create(LIS_Comm comm, LIS_VECTOR *vec); 
	extern int lis_vector_set_size(LIS_VECTOR vec, int local_n, int global_n); 
	extern int lis_vector_destroy(LIS_VECTOR vec);
	extern int lis_vector_duplicate(void *vin, LIS_VECTOR *vout);
	extern int lis_vector_get_size(LIS_VECTOR v, int *local_n, int *global_n);
	extern int lis_vector_get_range(LIS_VECTOR v, int *is, int *ie);
	extern int lis_vector_get_value(LIS_VECTOR v, int i, LIS_SCALAR *value);
	extern int lis_vector_get_values(LIS_VECTOR v, int start, int count, LIS_SCALAR value[]);
	extern int lis_vector_set_value(int flag, int i, LIS_SCALAR value, LIS_VECTOR v);
	extern int lis_vector_set_values(int flag, int count, int index[], LIS_SCALAR value[], LIS_VECTOR v);
	extern int lis_vector_print(LIS_VECTOR x);
	extern int lis_vector_is_null(LIS_VECTOR v);
	extern int lis_vector_axpy(LIS_SCALAR alpha, LIS_VECTOR vx, LIS_VECTOR vy);
	extern int lis_vector_xpay(LIS_VECTOR vx, LIS_SCALAR alpha, LIS_VECTOR vy);
	extern int lis_vector_axpyz(LIS_SCALAR alpha, LIS_VECTOR vx, LIS_VECTOR vy, LIS_VECTOR vz);
	extern int lis_vector_set_all(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern int lis_vector_copy(LIS_VECTOR vsrc, LIS_VECTOR vdst);
	extern int lis_vector_scale(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern int lis_vector_abs(LIS_VECTOR vx);
	extern int lis_vector_reciprocal(LIS_VECTOR vx);
	extern int lis_vector_shift(LIS_SCALAR alpha, LIS_VECTOR vx);
	extern int lis_vector_pmul(LIS_VECTOR vx,LIS_VECTOR vy,LIS_VECTOR vz);
	extern int lis_vector_pdiv(LIS_VECTOR vx,LIS_VECTOR vy,LIS_VECTOR vz);
	extern int lis_vector_dot(LIS_VECTOR vx, LIS_VECTOR vy, LIS_SCALAR *val);
	extern int lis_vector_nrm2(LIS_VECTOR vx, LIS_REAL *val);
	extern int lis_vector_sum(LIS_VECTOR vx, LIS_SCALAR *val);
/****************************/
/* Matrix                   */
/****************************/
	extern int lis_matrix_create(LIS_Comm comm, LIS_MATRIX *Amat);
	extern int lis_matrix_destroy(LIS_MATRIX Amat);
	extern int lis_matrix_assemble(LIS_MATRIX A);
	extern int lis_matrix_duplicate(LIS_MATRIX Ain, LIS_MATRIX *Aout);
	extern int lis_matrix_set_size(LIS_MATRIX A, int local_n, int global_n);
	extern int lis_matrix_get_size(LIS_MATRIX A, int *local_n, int *global_n);
	extern int lis_matrix_get_range(LIS_MATRIX A, int *is, int *ie);
	extern int lis_matrix_set_type(LIS_MATRIX A, int matrix_type);
	extern int lis_matrix_get_type(LIS_MATRIX A, int *matrix_type);
	extern int lis_matrix_set_value(int flag, int i, int j, LIS_SCALAR value, LIS_MATRIX A);
	extern int lis_matrix_malloc(LIS_MATRIX A, int nnz_row, int nnz[]);
	extern int lis_matrix_get_diagonal(LIS_MATRIX A, LIS_VECTOR d);
	extern int lis_matrix_scaling(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR D, int action);
	extern int lis_matrix_convert(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_copy(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_set_blocksize(LIS_MATRIX A, int bnr, int bnc, int row[], int col[]);

	extern int lis_matrix_malloc_crs(int n, int nnz, int **ptr, int **index, LIS_SCALAR **value);
	extern int lis_matrix_set_crs(int nnz, int *row, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_ccs(int n, int nnz, int **ptr, int **index, LIS_SCALAR **value);
	extern int lis_matrix_set_ccs(int nnz, int *row, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_bsr(int n, int bnr, int bnc, int bnnz, int **bptr, int **bindex, LIS_SCALAR **value);
	extern int lis_matrix_set_bsr(int bnr, int bnc, int bnnz, int *bptr, int *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_set_msr(int nnz, int ndz, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_msr(int n, int nnz, int ndz, int **index, LIS_SCALAR **value);
	extern int lis_matrix_set_ell(int maxnzr, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_ell(int n, int maxnzr, int **index, LIS_SCALAR **value);
	extern int lis_matrix_set_jds(int nnz, int maxnzr, int *perm, int *ptr, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_jds(int n, int nnz, int maxnzr, int **perm, int **ptr, int **index, LIS_SCALAR **value);
	extern int lis_matrix_set_dia(int nnd, int *index, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_dia(int n, int nnd, int **index, LIS_SCALAR **value);
	extern int lis_matrix_malloc_bsc(int n, int bnr, int bnc, int bnnz, int **bptr, int **bindex, LIS_SCALAR **value);
	extern int lis_matrix_set_bsc(int bnr, int bnc, int bnnz, int *bptr, int *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_get_vbr_rowcol(LIS_MATRIX Ain, int *nr, int *nc, int **row, int **col);
	extern int lis_matrix_malloc_vbr(int n, int nnz, int nr, int nc, int bnnz, int **row, int **col, int **ptr, int **bptr, int **bindex, LIS_SCALAR **value);
	extern int lis_matrix_set_vbr(int nnz, int nr, int nc, int bnnz, int *row, int *col, int *ptr, int *bptr, int *bindex, LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_set_dns(LIS_SCALAR *value, LIS_MATRIX A);
	extern int lis_matrix_malloc_dns(int n, int gn, LIS_SCALAR **value);
/****************************/
/* Matrix Vector Product    */
/****************************/
	extern int lis_matvec(LIS_MATRIX A, LIS_VECTOR x, LIS_VECTOR y);
	extern int lis_matvect(LIS_MATRIX A, LIS_VECTOR x, LIS_VECTOR y);
/****************************/
/* Solvers                  */
/****************************/
	extern int lis_solver_create(LIS_SOLVER *solver);
	extern int lis_solver_destroy(LIS_SOLVER solver);
	extern int lis_solver_get_iters(LIS_SOLVER solver, int *iters);
	extern int lis_solver_get_itersex(LIS_SOLVER solver, int *iters, int *iters_double, int *iters_quad);
	extern int lis_solver_get_times(LIS_SOLVER solver, double *times);
	extern int lis_solver_get_timesex(LIS_SOLVER solver, double *times, double *itimes, double *ptimes, double *p_c_times, double *p_i_times);
	extern int lis_solver_get_residualnorm(LIS_SOLVER solver, LIS_REAL *residual);
	extern int lis_solver_get_solver(LIS_SOLVER solver, int *nsol);
	extern int lis_solver_get_status(LIS_SOLVER solver, int *status);
	extern int lis_solver_set_option(char *text, LIS_SOLVER solver);
	extern int lis_solver_set_optionC(LIS_SOLVER solver);
	extern int lis_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver);
	extern LIS_PRECON_REGISTER *precon_register_top;
	extern int precon_register_type;
	extern int lis_precon_register(char *name, LIS_PRECON_CREATE_XXX pcreate, LIS_PSOLVE_XXX psolve, LIS_PSOLVET_XXX psolvet);
	extern int lis_precon_register_free(void);
	extern int lis_get_solvername(int solver, char *name);
/****************************/
/* I/O                      */
/****************************/
	extern int lis_input(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, char *filename);
	extern int lis_input_vector(LIS_VECTOR v, char *filename);
	extern int lis_output(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_vector(LIS_VECTOR v, int format, char *filename);
	extern int lis_output_residual_history(LIS_SOLVER solver, char *filename);
/****************************/
/* Utilities                */
/****************************/
	extern int lis_initialize(int* argc, char** argv[]);
	extern int lis_finalize(void);
	extern double lis_wtime(void);
	extern void CHKERR(int err);
	extern int lis_debug_trace_func(int flag, char *func);
	extern void *lis_malloc( size_t size, char *tag );
	extern void *lis_calloc( size_t size, char *tag );
	extern void *lis_realloc( void *p, size_t size );
	extern void lis_free(void *p);
	extern void lis_free2(int n, ...);
	extern int lis_is_malloc( void *p );
	extern void lis_date(char *date);

#ifdef __cplusplus
}
#endif

#define LIS_TRUE							1
#define LIS_FALSE							0
#define LIS_FAILS							-1
#define LIS_SUCCESS							0
#define LIS_ILL_OPTION						1
#define LIS_ERR_ILL_ARG						1
#define LIS_BREAKDOWN						2
#define LIS_OUT_OF_MEMORY					3
#define LIS_ERR_OUT_OF_MEMORY				3
#define LIS_MAXITER							4
#define LIS_ERR_NOT_IMPLEMENTED				5
#define LIS_ERR_FILE_IO						6


#if 1
#define LIS_GET_ISIE(id,nprocs,n,is,ie) \
		if( (id) < (n)%(nprocs) )				\
		{								\
			(ie) = (n)/(nprocs)+1;			\
			(is) = (ie)*(id);					\
		}								\
		else							\
		{								\
			(ie) = (n)/(nprocs);				\
			(is) = (ie)*(id) + (n)%(nprocs);		\
		}								\
		(ie) = (ie)+(is);
#else
#define LIS_GET_ISIE(id,nprocs,n,is,ie) \
		(ie) = (n)/(nprocs)+1;			\
		(is) = (id)*(ie)<((n)+1)?(id)*(ie):(n)+1;	\
		(ie) = (is)+(ie)-1<(n)?(is)+(ie)-1:(n);
#endif
#endif
