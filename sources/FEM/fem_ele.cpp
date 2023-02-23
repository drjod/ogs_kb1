/*
   The members of class Element definitions.
   Designed and programmed by WW, 06/2004
 */

//#include "makros.h"
//#include <iostream>
#include "fem_ele_std.h"
#include "fem_ele.h"
#include <cfloat>
/* Objekte */
#include "rf_pcs.h" 

#include "femlib.h"
#include "mathlib.h"
#include <algorithm>
//#include "matrix_class.h"
// MSHLib
//#include "msh_elem.h"
// Will be removed when new FEM is ready
//=============================================
FiniteElement::CElement* elem_dm = NULL;
//=============================================

using FiniteElement::CElement;
using MeshLib::CElem;
using MeshLib::CEdge;
using MeshLib::CNode;
using Math_Group::vec;

extern vector<CMediumProperties*>mmp_vector;


namespace FiniteElement
{
/**************************************************************************
   FEMLib-Method:
   Task: Constructor of class CElement
   Programing:
   01/2005 WW Implementation
   01/2005 OK 1D case
   01/2006 WW Axisymmetry
   Last modified:
**************************************************************************/
CElement::CElement(int CoordFlag, const int order)//, bool _2D_mesh_with_line_elements)
	: MeshElement(NULL), Order(order), ele_dim(1), nGaussPoints(1), nGauss(1),
	  ShapeFunction(NULL), ShapeFunctionHQ(NULL),
	  GradShapeFunction(NULL), GradShapeFunctionHQ(NULL),
	  T_Flag(false), C_Flag(false), F_Flag(false), D_Flag(0), RD_Flag(false),
	  extrapo_method(ExtrapolationMethod::EXTRAPO_LINEAR)
{
	int i;
	//
	nGauss = 3;
	//flag_2D_mesh_with_line_elements = _2D_mesh_with_line_elements;
	//
	if(CoordFlag < 0)  // Axisymmetry
	{
		CoordFlag *= -1;
		axisymmetry = true;
	}
	else
		axisymmetry = false;
	//
	dim = CoordFlag / 10;
	coordinate_system = CoordFlag;
	for(i = 0; i < 4; i++)
		unit[i] = 0.0;
	switch(dim)
	{
	case 1:                               //OK
		// Memory allocated for maxium 3 nodes elements
		Jacobian = new double[1];
		invJacobian = new double[1];
		shapefct = new double[2];
		shapefctHQ = new double[3];
		dshapefct = new double[6];
		dshapefctHQ = new double[9];
		break;
	case 2:
		// Memory allocated for maxium 9 nodes elements
		Jacobian = new double[4];
		invJacobian = new double[4];
		shapefct = new double[4];
		shapefctHQ = new double[9];
		dshapefct = new double[18];
		dshapefctHQ = new double[18];
		break;
	case 3:
		// Memory allocated for maxium 20 nodes elements
		Jacobian = new double[9];
		invJacobian = new double[9];
		shapefct = new double[8];
		shapefctHQ = new double[20];
		dshapefct = new double[24];
		dshapefctHQ = new double[60];
		//
		break;
	}
	time_unit_factor = 1.0;

	if(M_Process)
		D_Flag = 4;
	if(MH_Process)
		D_Flag = 41;
	F_Flag = H_Process;
	T_Flag = T_Process;
	C_Flag = MASS_TRANSPORT_Process;
	PT_Flag = 0;                          // PCH Initialize to be no RWPT.
	RD_Flag = RD_Process;
	MCF_Flag = MULTI_COMPONENTIAL_FLOW_Process;

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	idxm = NULL;  //> global indices of local matrix rows 
	idxn = NULL;  //> global indices of local matrix columns 
	local_idx = NULL; //> local index for local assemble
	//local_matrix = NULL; //>  local matrix 
	//local_vec = NULL; //>  local vector  
#endif


}

//  Destructor of class Element
CElement::~CElement()
{
	delete  [] Jacobian;
	delete  [] invJacobian;
	delete  [] shapefct;
	delete  [] shapefctHQ;
	delete  [] dshapefct;
	delete  [] dshapefctHQ;
	Jacobian = NULL;
	shapefct = NULL;
	dshapefct = NULL;
	dshapefctHQ = NULL;
	shapefctHQ = NULL;

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	if(idxm)
	  delete [] idxm;  
	if(idxn)
	  delete [] idxn;  
	if (local_idx)
	  delete [] local_idx;
	//if (local_idx)
	//  delete [] local_matrix;
	//if (local_idx)
	//  delete [] local_vec;
#endif

}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   05/2007 WW 1D in 2D
   Last modified:
**************************************************************************/
void CElement::ConfigElement(CElem* MElement, const int nquadrature_points,
		                     bool FaceIntegration)
{
	CNode* a_node = NULL;                 //07.04.2009. WW
	MeshElement = MElement;
	Index = MeshElement->GetIndex();
	nnodes = MeshElement->nnodes;
	nnodesHQ = MeshElement->nnodesHQ;
	bool done = false;
	ConfigNumerics(MeshElement->GetElementType(), nquadrature_points);
	if (MeshElement->quadratic)
		nNodes = nnodesHQ;
	else
		nNodes = nnodes;
	// Node indices
	for(int i = 0; i < nNodes; i++)
		nodes[i] = MeshElement->nodes_index[i];
	// Put coordinates of nodes to buffer to enhance the computation
	if(!FaceIntegration)
	{
		if(dim != ele_dim)
		{
//            a_node0 = MeshElement->nodes[0];      //07.04.2007. WW
			double const* const coords_node_0 (MeshElement->nodes[0]->getData());
			for(int i = 0; i < nNodes; i++)
			{
				double const* const coords_node_i (MeshElement->nodes[i]->getData());
//               a_node = MeshElement->nodes[i];    //07.04.2007. WW
//               dy = dz = 0.;
//               dx = a_node->X()-a_node0->X();
//               dy = a_node->Y()-a_node0->Y();
//               dz = a_node->Z()-a_node0->Z();
				double dx (coords_node_i[0] - coords_node_0[0]);
				double dy (coords_node_i[1] - coords_node_0[1]);
				double dz (coords_node_i[2] - coords_node_0[2]);

				X[i] =  (*MeshElement->transform_tensor)(0,0) * dx
				       + (*MeshElement->transform_tensor)(1,0) * dy
				       + (*MeshElement->transform_tensor)(2,0) * dz;
				Y[i] =  (*MeshElement->transform_tensor)(0,1) * dx
				       + (*MeshElement->transform_tensor)(1,1) * dy
				       + (*MeshElement->transform_tensor)(2,1) * dz;
//               Z[i] =  a_node->Z();
				Z[i] = coords_node_i[2];
			}
			done = true;
		}
		else
		{
			switch(dim)
			{
			case 1:
				if(coordinate_system % 10 == 1)
				{
					for(int i = 0; i < nNodes; i++)
					{
						//07.04.2007. WW
//                        a_node = MeshElement->nodes[i];
//                        X[i] = a_node->Y();
//                        Y[i] = a_node->X();
//                        Z[i] = a_node->Z();
						double const* const coords_node_i (
						        MeshElement->nodes[i]->getData());
						X[i] = coords_node_i[1];
						Y[i] = coords_node_i[0];
						Z[i] = coords_node_i[2];
					}
					done = true;
				}
				else if(coordinate_system % 10 == 2)
				{
					for(int i = 0; i < nNodes; i++)
					{
						//07.04.2007. WW
//                        a_node = MeshElement->nodes[i];
//                        X[i] = a_node->Z();
//                        Y[i] = a_node->Y();
//                        Z[i] = a_node->X();
						double const* const coords_node_i (
						        MeshElement->nodes[i]->getData());
						X[i] = coords_node_i[2];
						Y[i] = coords_node_i[1];
						Z[i] = coords_node_i[0];
					}
					done = true;
				}
				break;
			case 2:
				if(coordinate_system % 10 == 2)
				{
					for(int i = 0; i < nNodes; i++)
					{
						//07.04.2007. WW
//                        a_node = MeshElement->nodes[i];
//                        X[i] = a_node->X();
//                        Y[i] = a_node->Z();
//                        Z[i] = a_node->Y();
						double const* const coords_node_i (
						        MeshElement->nodes[i]->getData());
						X[i] = coords_node_i[0];
						Y[i] = coords_node_i[2];
						Z[i] = coords_node_i[1];
					}
					done = true;
				}
				break;
			}
		}
	}
	//
	if(!done)
	{
		for(int i = 0; i < nNodes; i++)
		{
			a_node = MeshElement->nodes[i]; //07.04.2007. WW
			double const* const coords (a_node->getData());
			X[i] = coords[0];
			Y[i] = coords[1];
			Z[i] = coords[2];
		}
	}

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
   if(!FaceIntegration)
   {
        if(MeshElement->g_index) // ghost nodes pcs->pcs_number_of_primary_nvals
        {
            act_nodes = MeshElement->g_index[0];
            act_nodes_h = MeshElement->g_index[1];

            for(int i = 0; i < act_nodes_h; i++)
            {
                local_idx[i] = MeshElement->g_index[i+2];
            }
        }
        else
        {
            act_nodes = nnodes;
            act_nodes_h = nnodesHQ;
            for(int i = 0; i < act_nodes_h; i++)
            {
                local_idx[i] = i;
            }
        }
   }
#endif

}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   Last modified:
**************************************************************************/
void CElement::CalculateRadius()
{
	Radius = 0.0;
	ComputeShapefct(1);
	for(int i = 0; i < nnodes; i++)
		Radius += shapefct[i] * X[i];
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2005 WW Implementation
   Last modified:
**************************************************************************/
void CElement::setOrder(const int order)
{
	Order = order;
	if(Order == 1)
		nNodes = nnodes;
	else if (Order == 2)
		nNodes = nnodesHQ;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   02/2005 OK Case 1: line elements
   01/2010 NW Higher order line elements
   Last modified:
**************************************************************************/
void CElement::ConfigNumerics(MshElemType::type ele_type, const int nquadrature_points)
{
	assert(nquadrature_points>0);
	// nGauss = GetNumericsGaussPoints(ElementType);
	switch(ele_type)
	{
	case MshElemType::LINE:
		ele_dim = 1;
		nGauss = 2;
		nGaussPoints = nGauss;
		ShapeFunction = ShapeFunctionLine;
		ShapeFunctionHQ = ShapeFunctionLineHQ;
		GradShapeFunction = GradShapeFunctionLine;
		GradShapeFunctionHQ = GradShapeFunctionLineHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
		return;
	case MshElemType::QUAD:
		ele_dim = 2;
		nGauss = nquadrature_points;
		nGaussPoints = nGauss * nGauss;
		ShapeFunction = ShapeFunctionQuad;
		ShapeFunctionHQ = ShapeFunctionQuadHQ;
		GradShapeFunction = GradShapeFunctionQuad;
		GradShapeFunctionHQ = GradShapeFunctionQuadHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
		return;
	case MshElemType::HEXAHEDRON:
		ele_dim = 3;        
		nGauss = nquadrature_points;
		nGaussPoints = nGauss * nGauss * nGauss;
		ShapeFunction = ShapeFunctionHex;
		ShapeFunctionHQ = ShapeFunctionHexHQ;
		GradShapeFunction = GradShapeFunctionHex;
		GradShapeFunctionHQ = GradShapeFunctionHexHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
		return;
	case MshElemType::TRIANGLE:
		ele_dim = 2;
		nGaussPoints = nGauss = 3; // Fixed to 3
		ShapeFunction = ShapeFunctionTri;
		ShapeFunctionHQ = ShapeFunctionTriHQ;
		GradShapeFunction = GradShapeFunctionTri;
		GradShapeFunctionHQ = GradShapeFunctionTriHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
		return;
	case MshElemType::TETRAHEDRON:
		ele_dim = 3;
		//	   nGaussPoints = nGauss = 15;  // Fixed to 15
		nGaussPoints = nGauss = 5; // Fixed to 5
		ShapeFunction = ShapeFunctionTet;
		ShapeFunctionHQ = ShapeFunctionTetHQ;
		GradShapeFunction = GradShapeFunctionTet;
		GradShapeFunctionHQ = GradShapeFunctionTetHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
		return;
	case MshElemType::PRISM:
		ele_dim = 3;
		nGaussPoints = 6;         // Fixed to 6
		nGauss = 3;               // Fixed to 3
		ShapeFunction = ShapeFunctionPri;
		ShapeFunctionHQ = ShapeFunctionPriHQ;
		GradShapeFunction = GradShapeFunctionPri;
		GradShapeFunctionHQ = GradShapeFunctionPriHQ;
		extrapo_method = ExtrapolationMethod::EXTRAPO_AVERAGE;
		return;
	case MshElemType::PYRAMID:
		ele_dim = 3;
		if (Order == 1)
			nGaussPoints = nGauss = 5;
		else
			nGaussPoints = nGauss = 8;  //13;
		ShapeFunction = ShapeFunctionPyra;
		ShapeFunctionHQ = ShapeFunctionPyraHQ13;
		GradShapeFunction = GradShapeFunctionPyra;
		GradShapeFunctionHQ = GradShapeFunctionPyraHQ13;
		extrapo_method = ExtrapolationMethod::EXTRAPO_AVERAGE;
		return;
	case MshElemType::INVALID:
		std::cerr << "[CElement::ConfigNumerics] invalid element type" << "\n";
		break;
	default:
		std::cerr << "[CElement::ConfigNumerics] unknown element type" << "\n";
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   Last modified:
**************************************************************************/
double CElement::interpolate(double const * const nodalVal, const int order) const
{
	int nn = nnodes;
	double* inTerpo = shapefct;
	if(order == 2)
	{
		nn = nnodes;
		inTerpo = shapefctHQ;
	}
	double val = 0.0;
	for(int i = 0; i < nn; i++)
		val += nodalVal[i] * inTerpo[i];
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   Last modified:
**************************************************************************/
double CElement::interpolate(const int idx, CRFProcess* m_pcs, const int order)
{
	int i;
	int nn = nnodes;
	double* inTerpo = shapefct;
	double val = 0.0;
	if(order == 2)
	{
		nn = nnodes;
		inTerpo = shapefctHQ;
	}
	//
	for(i = 0; i < nn; i++)
		node_val[i] = m_pcs->GetNodeValue(nodes[i], idx);
	for(int i = 0; i < nn; i++)
		val += node_val[i] * inTerpo[i];
	return val;
} 
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   Last modified:
**************************************************************************/
//double CElement::elemnt_average (const int idx, const int order)
double CElement::elemnt_average (const int idx, CRFProcess* m_pcs, const int order )
{
	int i;
	int nn = nnodes;
	double val = 0.0;
	//WW    double* inTerpo = shapefct;
	if(order == 2)
		nn = nnodes;
	//WW       inTerpo = shapefctHQ;
	//
	for(i = 0; i < nn; i++)
		node_val[i] = m_pcs->GetNodeValue(nodes[i], idx);
	return val / (double)nn;
}


/**************************************************************************
   The generalized Jacobian caculation

   Arguments:
    const double *unit:           Local coordiantes
   return
        The determinate of Jacobian
   Programmaenderungen:
   06/2006     WW
   02/2005 OK case 1, line elements
   01/2006 WW Axisymmtry
   09/2006 WW 1D element in 3D
**************************************************************************/
double CElement::computeJacobian(const int order)
{
	int k = 0;
	int nodes_number = nnodes;
	double DetJac = 0.0;
	double* dN = dshapefct;
	//    double *sh = shapefct;
	double dx,dy,dz;
	dx = dy = dz = 0.0;

	if(order == 2)                        //OK4104
	{
		nodes_number = nnodesHQ;
		dN = dshapefctHQ;
		GradShapeFunctionHQ(dN, unit);
	}
	else
		GradShapeFunction(dN, unit);
	for(size_t i = 0; i < ele_dim*ele_dim; i++)
		Jacobian[i] = 0.0;
	//--------------------------------------------------------------------
	switch(ele_dim)
	{
	//................................................................
	case 1:
		// If Line in X or Z direction, coordinate is saved in local X
		// If Line in 3D space, a transform is applied and cast coordinate in local X
		dx = X[1] - X[0];         //+Y[1]-Y[0];
		Jacobian[0] = 0.5 * dx;
		invJacobian[0] = 2.0 / dx;
		DetJac = Jacobian[0];
		//WW
		//if(MeshElement->area>0)
		DetJac *= MeshElement->area;
		//WW          DetJac*=MeshElement->GetFluxArea();//CMCD
		if(axisymmetry && !flag_2D_mesh_with_line_elements)  // JOD
		{
			CalculateRadius();
			DetJac *= Radius; //2.0*pai*Radius;
		}
		break;
	//................................................................
	case 2:
		for(int i = 0,j = nodes_number; i < nodes_number; i++,j++)
		{
			Jacobian[0] += X[i] * dN[i];
			Jacobian[1] += Y[i] * dN[i];
			Jacobian[2] += X[i] * dN[j];
			Jacobian[3] += Y[i] * dN[j];
		}
		DetJac =  Jacobian[0] * Jacobian[3] - Jacobian[1] * Jacobian[2];
		if (fabs(DetJac) < MKleinsteZahl)
		{
			std::cout << "\n*** Jacobian: Det == 0 " << DetJac << "\n";
			abort();
		}
		invJacobian[0] = Jacobian[3];
		invJacobian[1] = -Jacobian[1];
		invJacobian[2] = -Jacobian[2];
		invJacobian[3] = Jacobian[0];
		for(size_t i = 0; i < ele_dim * ele_dim; i++)
			invJacobian[i] /= DetJac;
		//
		//By WW
		//if(MeshElement->area>0)
		DetJac *= MeshElement->area;
		//WW          DetJac*=MeshElement->GetFluxArea();//CMCD
		if(axisymmetry
			&&
				std::find_if(X, X+nodes_number, // JOD 2018-10-16
						[](double x) { return x<-1.e-10; }) == X+nodes_number) //no node has coordinate x<0.
		{
			CalculateRadius();
			DetJac *= Radius * 6.2831853072; //2.0*pai*Radius;
		}
		break;
	//................................................................
	case 3: {
		int j;
		for(int i = 0; i < nodes_number; i++)
		{
			j = i + nodes_number;
			k = i + 2 * nodes_number;

			Jacobian[0] += X[i] * dN[i];
			Jacobian[1] += Y[i] * dN[i];
			Jacobian[2] += Z[i] * dN[i];

			Jacobian[3] += X[i] * dN[j];
			Jacobian[4] += Y[i] * dN[j];
			Jacobian[5] += Z[i] * dN[j];

			Jacobian[6] += X[i] * dN[k];
			Jacobian[7] += Y[i] * dN[k];
			Jacobian[8] += Z[i] * dN[k];
		}
		DetJac = Jacobian[0] *
		         (Jacobian[4] * Jacobian[8] - Jacobian[7] * Jacobian[5])
		         + Jacobian[6] *
		         (Jacobian[1] * Jacobian[5] - Jacobian[4] * Jacobian[2])
		         + Jacobian[3] *
		         (Jacobian[2] * Jacobian[7] - Jacobian[8] * Jacobian[1]);

		if (fabs(DetJac) < MKleinsteZahl)
		{
			std::cout << "\n*** Jacobian: DetJac == 0 " << DetJac << "\n";
			abort();
		}
		invJacobian[0] =  Jacobian[4] * Jacobian[8] - Jacobian[7] * Jacobian[5];
		invJacobian[1] =  Jacobian[2] * Jacobian[7] - Jacobian[1] * Jacobian[8];
		invJacobian[2] =  Jacobian[1] * Jacobian[5] - Jacobian[2] * Jacobian[4];
		//
		invJacobian[3] =  Jacobian[5] * Jacobian[6] - Jacobian[8] * Jacobian[3];
		invJacobian[4] =  Jacobian[0] * Jacobian[8] - Jacobian[6] * Jacobian[2];
		invJacobian[5] =  Jacobian[2] * Jacobian[3] - Jacobian[5] * Jacobian[0];
		//
		invJacobian[6] =  Jacobian[3] * Jacobian[7] - Jacobian[6] * Jacobian[4];
		invJacobian[7] =  Jacobian[1] * Jacobian[6] - Jacobian[7] * Jacobian[0];
		invJacobian[8] =  Jacobian[0] * Jacobian[4] - Jacobian[3] * Jacobian[1];
		for(size_t i = 0; i < ele_dim * ele_dim; i++)
			invJacobian[i] /= DetJac;
		break;
	} // end case 3
	}
	//--------------------------------------------------------------------
	// Use absolute value (for grids by gmsh, whose orientation is clockwise)
	return fabs(DetJac);
}
/***************************************************************************
   GeoSys - Funktion: CElement::RealCoordinates

   Aufgabe:
        Mapping to real coordaintes from the local ones of quadratic traingle
   element.
   Formalparameter:
           E:
             double * x         : Array of size 3, real coordiantes
             const double *u    : Array of size 2, unit coordiantes

   Programming:
   06/2003     WW        Erste Version
   07/2005     WW        Change due to geometry element object
 **************************************************************************/
void CElement::RealCoordinates(double* realXYZ)
{
	int i;
	double* df = shapefct;
	if(Order == 2)
		df = shapefctHQ;
	for(i = 0; i < 3; i++)
		realXYZ[i] = 0.0;

	for(i = 0; i < nNodes; i++)
	{
		realXYZ[0] += df[i] * X[i];
		realXYZ[1] += df[i] * Y[i];
		realXYZ[2] += df[i] * Z[i];
	}
}
/***************************************************************************
   GeoSys - Funktion: CElement::UnitCoordinates

   Aufgabe:
        Get unit coodinates from the real ones
   element.
   Formalparameter:
           E:
             double * x         : Array of size 3, real coordiantes

   Programming:
   06/2003     WW        Erste Version
 **************************************************************************/
void CElement::UnitCoordinates(double* realXYZ)
{
	setOrder(Order);

	x1buff[0] = X[0];
	x1buff[1] = Y[0];
	x1buff[2] = Z[0];

	for(int i = 1; i < nNodes; i++)
	{
		x1buff[0] += X[i];
		x1buff[1] += Y[i];
		x1buff[2] += Z[i];
	}
	for(size_t i = 0; i < 3; i++)
		x1buff[i] /= (double)nNodes;

	for(size_t i = 0; i < ele_dim; i++)
		realXYZ[i] -= x1buff[i];

	for(size_t i = 0; i < ele_dim; i++)
	{
		unit[i] = 0.0;
		for(size_t j = 0; j < ele_dim; j++)
			unit[i] += invJacobian[j * ele_dim + i] * realXYZ[j];
	}

	for(size_t i = 0; i < ele_dim; i++)
		realXYZ[i] = unit[i];
}
/***************************************************************************

   08/2005     WW        Prism element
 **************************************************************************/
void CElement::SetGaussPoint(const int gp, int& gp_r, int& gp_s, int& gp_t)
{
	switch(MeshElement->GetElementType())
	{
	case MshElemType::LINE:               // Line
		gp_r = gp;
		unit[0] = MXPGaussPkt(nGauss, gp_r);
		return;
	case MshElemType::QUAD:               // Quadralateral
		gp_r = (int)(gp / nGauss);
		gp_s = gp % nGauss;
		unit[0] = MXPGaussPkt(nGauss, gp_r);
		unit[1] = MXPGaussPkt(nGauss, gp_s);
		return;
	case MshElemType::HEXAHEDRON:         // Hexahedra
		gp_r = (int)(gp / (nGauss * nGauss));
		gp_s = (gp % (nGauss * nGauss));
		gp_t = gp_s % nGauss;
		gp_s /= nGauss;
		unit[0] = MXPGaussPkt(nGauss, gp_r);
		unit[1] = MXPGaussPkt(nGauss, gp_s);
		unit[2] = MXPGaussPkt(nGauss, gp_t);
		return;
	case MshElemType::TRIANGLE:           // Triangle
		SamplePointTriHQ(gp, unit);
		break;
	case MshElemType::TETRAHEDRON:        // Tedrahedra
		//To be flexible          SamplePointTet15(gp, unit);
		SamplePointTet5(gp, unit);
		return;
	case MshElemType::PRISM:              // Prism
		gp_r = gp % nGauss;
		SamplePointTriHQ(gp_r, unit);
        //
		gp_s = nGaussPoints/nGauss;
		gp_t = (int)(gp / nGauss);
		unit[2] = MXPGaussPkt(gp_s,  gp_t);
		return;
	case MshElemType::PYRAMID: // Pyramid
		if (Order == 1)
			SamplePointPyramid5(gp, unit);
		else
			SamplePointPyramid8(gp, unit);  //SamplePointPyramid13(gp, unit);
		return;
	default:
		throw std::runtime_error("CElement::SetGaussPoint invalid mesh element type given");
	}
}
/***************************************************************************
   GeoSys - Funktion:
           CElement:: GetGaussData(const int gp)

   Aufgabe:
          Get Gauss points and weights, compute Jacobian
   Formalparameter:
           E:
             const int gp   : Gauss point index

   Programming:
   06/2004     WW        Erste Version
   08/2005     WW        Prism element
   02/2005 OK case 1
   02/2007 WW Abstract the calcultion of Gauss point in one function
 **************************************************************************/
double CElement::GetGaussData(int gp, int& gp_r, int& gp_s, int& gp_t)
{
	double fkt = 0.0;
	SetGaussPoint(gp, gp_r, gp_s, gp_t);
	switch(MeshElement->GetElementType())
	{
	case MshElemType::LINE:               // Line
		fkt = computeJacobian(Order) * MXPGaussFkt(nGauss, gp_r);
		break;
	case MshElemType::QUAD:               // Quadralateral
		fkt = computeJacobian(Order);
		fkt *= MXPGaussFkt(nGauss, gp_r) * MXPGaussFkt(nGauss, gp_s);
		break;
	case MshElemType::HEXAHEDRON:         // Hexahedra
		fkt = computeJacobian(Order);
		fkt *=   MXPGaussFkt(nGauss, gp_r) * MXPGaussFkt(nGauss, gp_s)
		       * MXPGaussFkt(nGauss, gp_t);
		break;
	case MshElemType::TRIANGLE:           // Triangle
		fkt = computeJacobian(Order);
		fkt *= unit[2];           // Weights
		break;
	case MshElemType::TETRAHEDRON:        // Tedrahedra
		//To be flexible          SamplePointTet15(gp, unit);
		fkt = computeJacobian(Order);
		fkt *= unit[3];           // Weights
		break;
	case MshElemType::PRISM:              // Prism
		fkt = computeJacobian(Order);
		// Weights
		fkt *= MXPGaussFktTri(nGauss, gp_r) * MXPGaussFkt(gp_s, gp_t);
		break;
	case MshElemType::PYRAMID: // Pyramid
		fkt = computeJacobian(Order);
		fkt *= unit[3];           // Weights
		break;
	default:
		throw std::runtime_error("CElement::GetGaussData invalid mesh element type given");
	}
	return fkt;
}

/***************************************************************************
   GeoSys - Funktion: FaceIntegration(const double *NodeVal)
   Task:   Used to treat Nuemann type boundary conditions (3D)
   Augument
       double *NodeVal : input, values of boundary conditions at all face node
                         Output, integration results of all face nodes
   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CElement::FaceIntegration(double* NodeVal)
{
	int i, gp, gp_r, gp_s;
	double fkt = 0.0, det, val;
	double* sf = shapefct;

	setOrder(Order);
	if(Order == 2)
	{
		sf = shapefctHQ;
		if(MeshElement->GetElementType() == MshElemType::QUAD)
			ShapeFunctionHQ = ShapeFunctionQuadHQ8;
	}

	det = MeshElement->GetVolume();
	for (i = 0; i < nNodes; i++)
		dbuff[i] = 0.0;
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		switch(MeshElement->GetElementType())
		{
		case MshElemType::LINE:   // Line
			gp_r = gp;
			unit[0] = MXPGaussPkt(nGauss, gp_r);
			fkt = 0.5* det* MXPGaussFkt(nGauss, gp_r);
			break;
		case MshElemType::TRIANGLE: // Triangle
			SamplePointTriHQ(gp, unit);
			fkt = 2.0 * det * unit[2]; // Weights
			break;
		case MshElemType::QUAD:   // Quadralateral
			gp_r = (int)(gp / nGauss);
			gp_s = gp % nGauss;
			unit[0] = MXPGaussPkt(nGauss, gp_r);
			unit[1] = MXPGaussPkt(nGauss, gp_s);
			fkt = 0.25* det* MXPGaussFkt(nGauss, gp_r) * MXPGaussFkt(nGauss, gp_s);
			break;
		default:
			throw std::runtime_error("CElement::FaceIntegration element type not handled");
		}

		ComputeShapefct(Order);
		val = 0.0;
		// Interpolation of value at Gauss point
		for (i = 0; i < nNodes; i++)
			val += NodeVal[i] * sf[i];
		// Integration
		for (i = 0; i < nNodes; i++)
			dbuff[i] += val * sf[i] * fkt;
	}
	for (i = 0; i < nNodes; i++)
		NodeVal[i] = dbuff[i];
}

/***************************************************************************
   GeoSys - Funktion:
           CElement::ComputeShapefct(const double *unit, const int order)

   Aufgabe:
         Compute values of shape function at integral point unit.
   Formalparameter:
           E:
             const double *u    : Array of size 2, unit coordiantes
             const int order    : 1, linear
                               2, quadratic

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CElement::ComputeShapefct(const int order)
{
	if(order == 1)
		ShapeFunction(shapefct, unit);
	else if(order == 2)
		ShapeFunctionHQ(shapefctHQ, unit);
}

/***************************************************************************
   GeoSys - Funktion:
           CElement::ComputeGradShapefct(const double *unit, const int order)

   Aufgabe:
         Compute values of shape function at integral point unit.
   Formalparameter:
           E:
             const double *unit    : Array of size 2, unit coordiantes
             const int order    : 1, linear
                               2, quadratic

   Programming:
   06/2004     WW        Erste Version
   10/2005     WW        2D element transform in 3D space
   06/2007     WW        1D in 2D
 **************************************************************************/
void CElement::ComputeGradShapefct(int order)
{
	int j_times_ele_dim_plus_k, j_times_nNodes_plus_i;
	static double Var[3];
	double* dN = dshapefct;

	if(order == 2)
		dN = dshapefctHQ;

	setOrder(order);
	for(int i = 0; i < nNodes; i++)
	{
		size_t j(0);
		for (j = 0, j_times_nNodes_plus_i = i; j < ele_dim; j++, j_times_nNodes_plus_i += nNodes) {
			Var[j] = dN[j_times_nNodes_plus_i];
			dN[j_times_nNodes_plus_i] = 0.0;
		}
		for (j = 0, j_times_ele_dim_plus_k = 0, j_times_nNodes_plus_i = i; j < ele_dim; j++, j_times_nNodes_plus_i
						+= nNodes) {
			for (size_t k = 0; k < ele_dim; k++, j_times_ele_dim_plus_k++) {
				dN[j_times_nNodes_plus_i] += invJacobian[j_times_ele_dim_plus_k] * Var[k];
			}
		}
	}
	// 1D element in 3D
	if((dim == 3 && ele_dim == 1) || (dim == 2 && ele_dim == 1))
		for(int i = 0; i < nNodes; i++)
		{
			for(size_t j = 1; j < dim; j++)
				dN[j * nNodes + i] = (*MeshElement->transform_tensor)(j) * dN[i];
			dN[i] *= (*MeshElement->transform_tensor)(0);
		}
	// 2D element in 3D
	if (dim == 3 && ele_dim == 2) {
		const size_t n_nodes_times_ele_dim( nNodes * ele_dim);
		for (size_t i = 0; i < n_nodes_times_ele_dim; i++)
			dShapefct[i] = dN[i];
		for (int i = 0; i < nNodes; i++)
			for (size_t j = 0; j < dim; j++) {
				dN[j * nNodes + i] = 0.0;
				for (size_t k = 0; k < ele_dim; k++)
					dN[j * nNodes + i] += (*MeshElement->transform_tensor)(j, k) * dShapefct[k
									* nNodes + i];
			}
	}
}
/***************************************************************************
   Center of reference element
   Programming:
   09/2005     WW        Erste Version
 **************************************************************************/
void CElement::SetCenterGP()
{
	// Center of the reference element
	unit[0] = unit[1] = unit[2] = 0.0;
	if(MeshElement->GetElementType() == MshElemType::TRIANGLE)
		unit[0] = unit[1] = 1.0 / 3.0;
	else if(MeshElement->GetElementType() == MshElemType::TETRAHEDRON)
		unit[0] = unit[1] = unit[2] = 0.25;
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::GetLocalIndex()
           For quadralateral and hexahedra element on the assumption that
           selected Gauss points form a quadralateral or hexahedra
   Aufgabe:
           Accumulate stress at each nodes
   Formalparameter:

   Programming:
   06/2004   WW
 **************************************************************************/
int CElement::GetLocalIndex(const int gp_r, const int gp_s, int gp_t)
{
	int LoIndex = -1;
	double r,s,t;

	//---------------------------------------------------------
	// Accumulate strains
	//---------------------------------------------------------
	switch (MeshElement->GetElementType())
	{
	case MshElemType::QUAD:               // Quadralateral
		r = MXPGaussPkt(nGauss, gp_r);
		s = MXPGaussPkt(nGauss, gp_s);
		if (r > 0.0 && s > 0.0)
			LoIndex = 0;
		else if (r < 0.0 && s > 0.0)
			LoIndex = 1;
		else if (r < 0.0 && s < 0.0)
			LoIndex = 2;
		else if (r > 0.0 && s < 0.0)
			LoIndex = 3;
		else if (fabs(r) < MKleinsteZahl && s > 0.0)
			LoIndex = 4;
		else if (r < 0.0 && fabs(s) < MKleinsteZahl)
			LoIndex = 5;
		else if (fabs(r) < MKleinsteZahl && s < 0.0)
			LoIndex = 6;
		else if (r > 0.0 && fabs(s) < MKleinsteZahl)
			LoIndex = 7;
		else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
			LoIndex = 8;
		break;
	case MshElemType::HEXAHEDRON:         // Hexahedra
		r = MXPGaussPkt(nGauss, gp_r);
		s = MXPGaussPkt(nGauss, gp_s);
		t = MXPGaussPkt(nGauss, gp_t);

		if (t > 0.0)
		{
			if (r > 0.0 && s > 0.0)
				LoIndex = 0;
			else if (r < 0.0 && s > 0.0)
				LoIndex = 1;
			else if (r < 0.0 && s < 0.0)
				LoIndex = 2;
			else if (r > 0.0 && s < 0.0)
				LoIndex = 3;
			else if (fabs(r) < MKleinsteZahl && s > 0.0)
				LoIndex = 8;
			else if (r < 0.0 && fabs(s) < MKleinsteZahl)
				LoIndex = 9;
			else if (fabs(r) < MKleinsteZahl && s < 0.0)
				LoIndex = 10;
			else if (r > 0.0 && fabs(s) < MKleinsteZahl)
				LoIndex = 11;
			else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
				return -1;
		}
		else if (fabs(t) < MKleinsteZahl)
		{
			if (fabs(r) < MKleinsteZahl || fabs(s) < MKleinsteZahl)
				return -1;
			if (r > 0.0 && s > 0.0)
				LoIndex = 16;
			else if (r < 0.0 && s > 0.0)
				LoIndex = 17;
			else if (r < 0.0 && s < 0.0)
				LoIndex = 18;
			else if (r > 0.0 && s < 0.0)
				LoIndex = 19;
		}
		if (t < 0.0)
		{
			if (r > 0.0 && s > 0.0)
				LoIndex = 4;
			else if (r < 0.0 && s > 0.0)
				LoIndex = 5;
			else if (r < 0.0 && s < 0.0)
				LoIndex = 6;
			else if (r > 0.0 && s < 0.0)
				LoIndex = 7;
			else if (fabs(r) < MKleinsteZahl && s > 0.0)
				LoIndex = 12;
			else if (r < 0.0 && fabs(s) < MKleinsteZahl)
				LoIndex = 13;
			else if (fabs(r) < MKleinsteZahl && s < 0.0)
				LoIndex = 14;
			else if (r > 0.0 && fabs(s) < MKleinsteZahl)
				LoIndex = 15;
			else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
				return -1;
		}
		break;
	default:
		throw std::runtime_error("CElement::GetLocalIndex invalid mesh element type given");
	}
	return LoIndex;
}
/***************************************************************************
   GeoSys - Funktion:
   Programming:
   02/2007   WW
 **************************************************************************/
void CElement::SetExtropoGaussPoints(const int i)
{
	int j = 0;
	MshElemType::type ElementType = MeshElement->GetElementType();
	//
	switch (ElementType)
	{
	case MshElemType::TRIANGLE:           // Triangle
		// Compute values at verteces
		// Compute values at verteces
		switch (i)
		{
		case 0:
			unit[0] = -0.1666666666667;
			unit[1] = -0.1666666666667;
			break;
		case 1:
			unit[0] = 1.6666666666667;
			unit[1] = -0.1666666666667;
			break;
		case 2:
			unit[0] = -0.1666666666667;
			unit[1] = 1.6666666666667;
			break;
		}
		break;
	case MshElemType::QUAD:               // Quadralateral element
		// Extropolation over nodes
		switch (i)
		{
		case 0:
			unit[0] = Xi_p;
			unit[1] = Xi_p;
			break;
		case 1:
			unit[0] = -Xi_p;
			unit[1] = Xi_p;
			break;
		case 2:
			unit[0] = -Xi_p;
			unit[1] = -Xi_p;
			break;
		case 3:
			unit[0] = Xi_p;
			unit[1] = -Xi_p;
			break;
		}
		break;
	case MshElemType::HEXAHEDRON:         // Hexahedra
		if (i < 4)
		{
			j = i;
			unit[2] = Xi_p;
		}
		else
		{
			j = i - 4;
			unit[2] = -Xi_p;
		}
		switch (j)
		{
		case 0:
			unit[0] = Xi_p;
			unit[1] = Xi_p;
			break;
		case 1:
			unit[0] = -Xi_p;
			unit[1] = Xi_p;
			break;
		case 2:
			unit[0] = -Xi_p;
			unit[1] = -Xi_p;
			break;
		case 3:
			unit[0] = Xi_p;
			unit[1] = -Xi_p;
			break;
		}
		break;
	case MshElemType::TETRAHEDRON:        // Tedrahedra
		// Compute values at verteces
		switch (i)
		{
		case 0:
			unit[0] = -0.166666666666667;
			unit[1] = -0.166666666666667;
			unit[2] = -0.166666666666667;
			break;
		case 1:
			unit[0] = 1.5;
			unit[1] = -0.166666666666667;
			unit[2] = -0.166666666666667;
			break;
		case 2:
			unit[0] = -0.166666666666667;
			unit[1] = 1.5;
			unit[2] = -0.166666666666667;
			break;
		case 3:
			unit[0] = -0.166666666666667;
			unit[1] = -0.166666666666667;
			unit[2] = 1.5;
		}
		break;
	case MshElemType::LINE:
		break;
	case MshElemType::PYRAMID: // WW. 09.2012. WW
		SamplePointPyramid5(i, unit);
		break;
	default:
		unit[0] = unit[1] = unit[2] = 0.; //07.01.2011. WW
		break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
   Programming:
   05/2011   NW
 **************************************************************************/
double CElement::CalcXi_p()
{
	double Xi_p = 0.0;
	MshElemType::type ElementType = MeshElement->GetElementType();
	if (ElementType == MshElemType::QUAD || ElementType == MshElemType::HEXAHEDRON)
	{
		double r = .0;
		for (gp = 0; gp < nGauss; gp++)
		{
			r = MXPGaussPkt(nGauss, gp);
			if(fabs(r) > Xi_p)
				Xi_p = fabs(r);
		}
		r = 1.0 / Xi_p;
		Xi_p = r;
	}

	return Xi_p;
}

/***************************************************************************
   GeoSys - Funktion:
   Programming:
   05/2011   NW Implementation
 **************************************************************************/
double CElement::CalcAverageGaussPointValues(double* GpValues)
{
	// average
	double avg = .0;
	for(int j = 0; j < nGauss; j++)
		avg += GpValues[j];
	avg /= nGauss;

	return avg;
}

/**************************************************************************
   ElementMatrix::AllocateMemory

   Arguments:
    const int EleIndex:  Element index,
    int type          : type
                         used to get element type.
                         type = 0, for Possion type
                         type = 1, for Possion equation with deformation coupling
                         type = 2, for Navier equation
                         type = 3, for Navier equation with pressure coupling
   type = 4, Monlithic scheme of u-p coupling
   type = 5, Mass Transport
   default = 0.

   Programmaenderungen:
   01/2005     WW

**************************************************************************/
void ElementMatrix::AllocateMemory(CElem* ele, int type)
{
	int nnodes, nnodesHQ, dim, size;
	size = 0;
	// The following two lines will be updated when new FEMGEO is ready
	nnodes = ele->GetVertexNumber();
	nnodesHQ = ele->GetNodesNumber_H();
	dim = ele->GetDimension();
	switch(type)
	{
	case 0:                               // H || T Process
		Mass = new Matrix(nnodes, nnodes);
		//        Laplace = new SymMatrix(nnodes);
		Laplace = new Matrix(nnodes, nnodes);
		RHS = new Vec(nnodes);
		break;
	case 1:                               // HM Partioned scheme, Flow
		Mass = new Matrix(nnodes, nnodes);
		//        Laplace = new SymMatrix(nnodes);
		Laplace = new Matrix(nnodes, nnodes);
		RHS = new Vec(nnodes);
		CouplingB = new Matrix(nnodes, dim * nnodesHQ);
		break;
	case 2:                               // M_Process only
		size = dim * nnodesHQ;
		Stiffness = new Matrix(size, size);
		RHS = new Vec(size);
		break;
	case 3:                               // MH Partioned scheme, M_Process
		size = dim * nnodesHQ;
		Stiffness = new Matrix(size, size);
		RHS = new Vec(size);
		CouplingA = new Matrix(dim * nnodesHQ, nnodes);
		break;
	case 4:                               // HM monothlic scheme
		Mass = new Matrix(nnodes, nnodes);
		//        Laplace = new SymMatrix(nnodes);
		Laplace = new Matrix(nnodes, nnodes);
		size = dim * nnodesHQ;
		Stiffness = new Matrix(size, size);
		RHS = new Vec(size + nnodes);
		CouplingA = new Matrix(dim * nnodesHQ, nnodes);
		CouplingB = new Matrix(nnodes, dim * nnodesHQ);
		break;
	case 5:                               // Mass Transport process
		Mass = new Matrix(nnodes, nnodes);
		Laplace = new Matrix(nnodes, nnodes);
		Advection = new Matrix(nnodes, nnodes);
		Storage = new Matrix(nnodes, nnodes);
		Content = new Matrix(nnodes, nnodes);
		RHS = new Vec(nnodes);
		break;
	}
}

/**************************************************************************
   ElementMatrix::ElementMatrix

   Arguments:
      const int EleIndex:  Element index,
                           used to get element type.
   Programmaenderungen:
   01/2005     WW

**************************************************************************/
ElementMatrix::~ElementMatrix()
{
	if(Mass)
		delete Mass;
	if(Laplace)
		delete Laplace;
	if(Advection)
		delete Advection;
	if(Storage)
		delete Storage;
	if(Content)
		delete Content;
	if(RHS)
		delete RHS;
	if(CouplingA)
		delete CouplingA;
	if(CouplingB)
		delete CouplingB;
	Mass = NULL;
	Laplace = NULL;
	Advection = NULL;
	Storage = NULL;
	Content = NULL;
	RHS = NULL;
	CouplingA = NULL;
	CouplingB = NULL;
}

/**************************************************************************
CElement::FaceNormalFluxIntegration

Used for TOTAL_FLUX calculation

Programming:
11/2014   JOD

**************************************************************************/

void CElement::FaceNormalFluxIntegration(long element_index, double *NodeVal, double *NodeVal_adv, CRFProcess* m_pcs, double* normal_vector)
{

	int gp, gp_r, gp_s;
	double fkt = 0.0, det;
	double *sf = shapefct;
	double normal_diff_flux_interpol, normal_adv_flux_interpol;
	double dbuff_adv[10], flux[3];

	setOrder(Order);
	if (Order == 2)
	{
		sf = shapefctHQ;
		if (MeshElement->GetElementType() == MshElemType::QUAD)
			ShapeFunctionHQ = ShapeFunctionQuadHQ8;
	}

	for (int i = 0; i < nNodes; i++) {
		dbuff[i] = 0.0;
		dbuff_adv[i] = 0.0;
	}

	det = MeshElement->GetVolume();
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		switch (MeshElement->GetElementType())
		{
		case MshElemType::LINE:   // Line
			gp_r = gp;
			unit[0] = MXPGaussPkt(nGauss, gp_r);
			fkt = 0.5* det* MXPGaussFkt(nGauss, gp_r);
			break;
		case MshElemType::TRIANGLE: // Triangle
			SamplePointTriHQ(gp, unit);
			fkt = 2.0 * det * unit[2]; // Weights
			break;
		case MshElemType::QUAD:   // Quadrilateral
			gp_r = (int)(gp / nGauss);
			gp_s = gp % nGauss;
			unit[0] = MXPGaussPkt(nGauss, gp_r);
			unit[1] = MXPGaussPkt(nGauss, gp_s);
			fkt = 0.25* det* MXPGaussFkt(nGauss, gp_r) * MXPGaussFkt(nGauss, gp_s);
			break;
		default:
			std::cerr << "Error in mass balance calculation: CElement::FaceIntegration element type not supported" << "\n";
		}
		//---------------------------------------------------------
		ComputeShapefct(Order);
		
		normal_diff_flux_interpol = 0.0;
		
		if ((m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW) || (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)) {

			for (int i = 0; i < nNodes; i++) 	// Darcy flux 
				normal_diff_flux_interpol += NodeVal[i] * sf[i];

			for (int i = 0; i < nNodes; i++)  // Integration
				dbuff[i] += normal_diff_flux_interpol * sf[i] * fkt;
		}
		else if ((m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT) || (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)) { 

			normal_adv_flux_interpol = 0.0;

			for (int i = 0; i < nNodes; i++) 	
				normal_adv_flux_interpol += NodeVal_adv[i] * sf[i];

			for (int i = 0; i < nNodes; i++) { // Integration
				// Fick or Fourier diffusion
				for (int l = 0; l < 3; l++)
					flux[l] = ele_gp_value[element_index]->Velocity(l, gp);// TransportFlux(l, gp);
				normal_diff_flux_interpol = PointProduction(flux, normal_vector);   //    fabs(PointProduction(flux, normal_vector));
				dbuff[i] += normal_diff_flux_interpol * sf[i] * fkt;
				// advection
				dbuff_adv[i] += normal_adv_flux_interpol * sf[i] * fkt;
			}
		} // end transport
	} // end gauss points

	for (int i = 0; i < nNodes; i++) {
		NodeVal[i] = dbuff[i];
		NodeVal_adv[i] = dbuff_adv[i];
	}
}
/***************************************************************************
GeoSys - Funktion: CalcJTC_st_ele
Task:   Used to calculate Joule-Thomson ST for elements
Augument

Output: ST for one whole element
Programming: 11/2014 WTP
**************************************************************************/
void CElement::CalcJTC_st_ele(long ele_index, double* JTCoeff, double* sourceterm_data, double k, double* kr_data, double* s_data, double n, bool multi)
{
    int i, gp, gp_r, gp_s, gp_t;
	double fkt = 0.0;
	//double det, val;
	double flux[3];
    double JTCoeff_val = 0.;
	double kr = 0.;
	double GasSat = 0.;
	double const g = 9.80665;
    double* dbuff_flux = new double[nNodes];
    double* sf = shapefct;
    //MeshLib::CElem* m_ele;
    //CFEMesh* m_msh = fem_msh_vector[0];
    //m_ele = m_msh->ele_vector[ele_index];
    CFluidProperties* m_mfp;
	if (multi)
		m_mfp = mfp_vector[1]; // 0=water, 1=gas
	else
		m_mfp = mfp_vector[0];
        
    setOrder(Order);

    //det = MeshElement->GetVolume();
    for (i = 0; i < nNodes; i++)
        dbuff_flux[i] = 0.;

    // Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		ComputeShapefct(Order);
		ComputeGradShapefct(1);   // Linear interpolation function

		//val = 0.0;
		JTCoeff_val = 0.;
		kr = 0.;
		GasSat = 0.;
		// Interpolation of values at Gauss point
		for (i = 0; i < nNodes; i++)
		{
			JTCoeff_val += (JTCoeff[i]*-1.0) * sf[i];
			GasSat += s_data[i] * sf[i];
			kr += kr_data[i] * sf[i];
		}
		
		for (int l = 0; l < 3; l++)
			if (multi)
				flux[l] = ele_gp_value[ele_index]->Velocity_g(l, gp);
			else
				flux[l] = ele_gp_value[ele_index]->Velocity(l, gp);
		double scalar_flux = MSkalarprodukt(flux, flux, 3);
		// Integration
		JTCoeff_val *= m_mfp->SpecificHeatCapacity() * m_mfp->Density() * (scalar_flux / ((kr * k) / m_mfp->Viscosity()) +  m_mfp->Density() * g * flux[2]);
		JTCoeff_val *= n * GasSat;
		
		//std::cout << "Ele: " << ele_index << " GP: " << gp << " flux[0]: " << flux[0] << " flux[1]: " << flux[1] << " flux[2]: " << flux[2] << " skalar_flux: " << scalar_flux << "\n";
		//std::cout << " C: " << m_mfp->SpecificHeatCapacity() << " dens: " << m_mfp->Density() << " Visc: " << m_mfp->Viscosity() << " k: " << k << " kr: " << kr <<  " n: " << n << "\n";
		//std::cout << "\n";

        for (i = 0; i < nNodes; i++)
		    dbuff_flux[i] += JTCoeff_val * sf[i] * fkt;
    }
    for (i = 0; i < nNodes; i++)
        sourceterm_data[i] = dbuff_flux[i];
}


/**************************************************************************
CFiniteElementStd::FaceNormalFluxIntegration

Used for TOTAL_FLUX calculation

Programming:
12/2014   JOD

**************************************************************************/

void CElement::CalculateFluxThroughFace(double factor, double *NodeVal, double *NodeVal_adv, double* normal_diff_flux, double* normal_adv_flux)
{

	int i, gp, gp_r=0, gp_s=0, gp_t=0;
	double fkt = 0.0, det;
	double Gauss_val_normal_diff_flux, Gauss_val_normal_adv_flux;

	setOrder(Order);

	det = MeshElement->GetVolume();
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		SetGaussPoint(gp, gp_r, gp_s, gp_t);
		//--------------------------------------------------------- 
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		switch (MeshElement->GetElementType())
		{
		case MshElemType::LINE:   // Line
			//gp_r = gp;
			//unit[0] = MXPGaussPkt(nGauss, gp_r);
			fkt = 0.5* det* MXPGaussFkt(nGauss, gp_r);
			break;
		case MshElemType::TRIANGLE: // Triangle
			//SamplePointTriHQ(gp, unit);
			fkt = 2.0 * det * unit[2]; // Weights
			break;
		case MshElemType::QUAD:   // Quadralateral
			/*gp_r = (int)(gp / nGauss);
			gp_s = gp % nGauss;
			unit[0] = MXPGaussPkt(nGauss, gp_r);
			unit[1] = MXPGaussPkt(nGauss, gp_s);*/
			fkt = 0.25* det* MXPGaussFkt(nGauss, gp_r) * MXPGaussFkt(nGauss, gp_s);
			break;
		default:
			throw std::runtime_error("CElement::FaceIntegration element type not handled");
		}

		fkt *= factor; // factor= 0.5 if face in domain; 1 if face at somain boundary

		ComputeShapefct(Order);

		Gauss_val_normal_diff_flux = 0;
		Gauss_val_normal_adv_flux = 0;

		for (i = 0; i < nNodes; i++) {	// Interpolation to Gauss point
			Gauss_val_normal_diff_flux += NodeVal[i] * shapefct[i];
			Gauss_val_normal_adv_flux += NodeVal_adv[i] * shapefct[i];
		}
		//////
		*normal_diff_flux += fkt * Gauss_val_normal_diff_flux;
		*normal_adv_flux += fkt * Gauss_val_normal_adv_flux;

	} // end gauss points


}

void FaceIntegration(CFEMesh* msh, std::vector<long> const &nodes_on_sfc,
		std::vector<double>&node_value_vector, Surface* m_surface, FiniteElement::DistributionType disType, int ele_gauss_points)
{
   if (!msh)
   {
      throw std::runtime_error("Warning in FaceIntegration: no MSH data, function doesn't function");
   }

   std::cout << "\tFace integration\n";
   long i, j, k, l;
   long this_number_of_nodes;
   int nfaces, nfn;
   int nodesFace[8];
   double nodesFVal[8];

   bool Const = false;
   if (disType == FiniteElement::CONSTANT
		   || disType == FiniteElement::CONSTANT_NEUMANN
		   || disType == FiniteElement::RECHARGE)	//MW
      //	if (dis_type_name.find("CONSTANT") != std::string::npos)
      Const = true;
   //----------------------------------------------------------------------
   // Interpolation of polygon values to nodes_on_sfc
   if (!Const)                                    // Get node BC by interpolation with surface
   {
      int nPointsPly = 0;
      double Area1, Area2;
      double Tol = 1.0e-9;
      bool Passed;
      const int Size = (int) nodes_on_sfc.size();
      double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];

      CGLPolyline* m_polyline = NULL;
      //Surface *m_surface = NULL;
      //m_surface = GEOGetSFCByName(geo_name);      //CC

      // list<CGLPolyline*>::const_iterator p = m_surface->polyline_of_surface_list.begin();
      std::vector<CGLPolyline*>::iterator p =
         m_surface->polyline_of_surface_vector.begin();

      for (j = 0; j < Size; j++)
      {
    	  double const*const pn (msh->nod_vector[nodes_on_sfc[j]]->getData());
//         pn[0] = msh->nod_vector[nodes_on_sfc[j]]->X();
//         pn[1] = msh->nod_vector[nodes_on_sfc[j]]->Y();
//         pn[2] = msh->nod_vector[nodes_on_sfc[j]]->Z();
         node_value_vector[j] = 0.0;
         Passed = false;
         // nodes close to first polyline
         p = m_surface->polyline_of_surface_vector.begin();
         while (p != m_surface->polyline_of_surface_vector.end())
         {
            m_polyline = *p;
            // Grativity center of this polygon
            for (i = 0; i < 3; i++)
               gC[i] = 0.0;
            vn[2] = 0.0;
            nPointsPly = (int) m_polyline->point_vector.size();
            if (m_polyline->point_vector.front() == m_polyline->point_vector.back())
               nPointsPly -= 1;
            for (i = 0; i < nPointsPly; i++)
            {
               gC[0] += m_polyline->point_vector[i]->x;
               gC[1] += m_polyline->point_vector[i]->y;
               gC[2] += m_polyline->point_vector[i]->z;

               vn[2] += m_polyline->point_vector[i]->getPropert();
            }
            for (i = 0; i < 3; i++)
               gC[i] /= (double) nPointsPly;
            // BC value at center is an average of all point values of polygon
            vn[2] /= (double) nPointsPly;

            // Area of this polygon by the grativity center
            for (i = 0; i < nPointsPly; i++)
            {
               p1[0] = m_polyline->point_vector[i]->x;
               p1[1] = m_polyline->point_vector[i]->y;
               p1[2] = m_polyline->point_vector[i]->z;
               k = i + 1;
               if (i == nPointsPly - 1)
                  k = 0;
               p2[0] = m_polyline->point_vector[k]->x;
               p2[1] = m_polyline->point_vector[k]->y;
               p2[2] = m_polyline->point_vector[k]->z;

               vn[0] = m_polyline->point_vector[i]->getPropert();
               vn[1] = m_polyline->point_vector[k]->getPropert();

               Area1 = fabs(ComputeDetTri(p1, gC, p2));

               Area2 = 0.0;
               // Check if pn is in the triangle by points (p1, gC, p2)
               Area2 = fabs(ComputeDetTri(p2, gC, pn));
               unit[0] = fabs(ComputeDetTri(gC, p1, pn));
               unit[1] = fabs(ComputeDetTri(p1, p2, pn));
               Area2 += unit[0] + unit[1];
               if (fabs(Area1 - Area2) < Tol)
               {
                  // Intopolation whin triangle (p1,p2,gC)
                  // Shape function
                  for (l = 0; l < 2; l++)
                     unit[l] /= Area1;
                  ShapeFunctionTri(NTri, unit);
                  for (l = 0; l < 3; l++)
                     node_value_vector[j] += vn[l] * NTri[l];
                  Passed = true;
                  break;
               }

            }
            //
            p++;
            if (Passed)
               break;
         }                                        // while
      }                                           //j
   }  // end !Const

   int Axisymm = 1;                               // ani-axisymmetry
   //CFEMesh* msh = m_pcs->m_msh;
   if (msh->isAxisymmetry())
      Axisymm = -1;                               // Axisymmetry is true
   CElem* elem = NULL;
   CElem* face = new CElem(1);
   CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());
   CNode* e_node = NULL;
   CElem* e_nei = NULL;
   //vec<CNode*> e_nodes(20);
   // vec<CElem*> e_neis(6);

   face->SetFace();
   this_number_of_nodes = (long) nodes_on_sfc.size();
   int nSize = (long) msh->nod_vector.size();
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

   for (i = 0; i < nSize; i++)
   {

	   msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   for (i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      k = nodes_on_sfc[i];
      G2L[k] = i;
   }

   //----------------------------------------------------------------------
   // NW 15.01.2010
   // 1) search element faces on the surface
   // 2) face integration


   std::set<long> set_nodes_on_sfc;               //unique set of node id on the surface
   for (i = 0; i < (long) nodes_on_sfc.size(); i++)
   {
      set_nodes_on_sfc.insert(nodes_on_sfc[i]);
   }

   //filtering elements: elements should have nodes on the surface
   //Notice: node-elements relation has to be constructed beforehand
   // CB THMBM
   //this->getProcess()->CheckMarkedElement(); // CB added to remove bug with deactivated Subdomains
   std::vector<long> vec_possible_elements;

   std::vector<long> elements_at_geo;
   std::vector<long> nodes_on_sfc2(nodes_on_sfc);

   msh->GetConnectedElements(nodes_on_sfc2, elements_at_geo);


   //init
   for (i = 0; i < (long) msh->ele_vector.size(); i++)
   {
      msh->ele_vector[i]->selected = 0;           //TODO can use a new variable
   }

   for (i = 0; i < this_number_of_nodes; i++)
   {
      k = nodes_on_sfc[i];
      for (j = 0; j < (long) msh->nod_vector[k]->getConnectedElementIDs().size(); j++)
      {
         l = msh->nod_vector[k]->getConnectedElementIDs()[j];
         if (msh->ele_vector[l]->selected == 0)
            vec_possible_elements.push_back(l);
         msh->ele_vector[l]->selected += 1;       // remember how many nodes of an element are on the surface
      }
   }

   //for (i = 0; i < msh->ele_vector.size(); i++)
	 //  std::cout << "----- " << msh->ele_vector[l]->selected << "\n";
   //search elements & face integration
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
      const size_t id_act_l_max = static_cast<size_t>(msh->getNumNodesLocal());
      const size_t id_act_h_min = msh-> GetNodesNumber(false);
      const size_t id_act_h_max = msh->getLargestActiveNodeID_Quadratic();
#endif

   int count;
   double fac = 1.0;
   for (i = 0; i < (long) vec_possible_elements.size(); i++)
   {
      elem = msh->ele_vector[vec_possible_elements[i]];
      if (!elem->GetMark())
         continue;
      nfaces = elem->GetFacesNumber();
      elem->SetOrder(msh->getOrder());
      for (j = 0; j < nfaces; j++)
      {
         e_nei = elem->GetNeighbor(j);
         nfn = elem->GetElementFaceNodes(j, nodesFace);
         //1st check
         if (elem->selected < nfn)
            continue;

         if(elem->GetDimension() != 3)
            continue;

         //2nd check: if all nodes of the face are on the surface
         count = 0;
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);
            if (set_nodes_on_sfc.count(e_node->GetIndex()) > 0)
            {
               count++;
            }
         }
         if (count != nfn)
            continue;
         // face integration
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);
            nodesFVal[k] = node_value_vector[G2L[e_node->GetIndex()]];
         }
         fac = 1.0;
	 // Not a surface face
         //   if (elem->GetDimension() == e_nei->GetDimension())  // removed by BW
         //     fac = 0.5;
         // BW 10.01.2020 ->if the neighbour of this face is not at the surface of a model then the source term is half of that
         // Add one more indicator for that
         if (elem->GetDimension() == e_nei->GetDimension())
             if (!m_surface->surface_at_model_surface)
             {
                 fac = 0.5;
                 std::cout
                     << "\t\t\tWarning in FaceIntegration: Surface is inside the model domain, area is multiplied with 0.5 and is calculated twice! (this line has to appear twice for each node)"
                     <<'\n';
             }

         face->SetFace(elem, j);
         face->SetOrder(msh->getOrder());
         face->ComputeVolume();
         fem->setOrder(msh->getOrder() + 1);
         fem->ConfigElement(face, ele_gauss_points, true);
         fem->FaceIntegration(nodesFVal);
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);

#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
            if(     (e_node->GetIndex() < id_act_l_max)
                ||  (    e_node->GetIndex() >= id_act_h_min
			 &&  e_node->GetIndex() < id_act_h_max)
             )
#endif
            NVal[G2L[e_node->GetIndex()]] += fac * nodesFVal[k];
         }
      }
   }

   for (i = 0; i < this_number_of_nodes; i++)
      node_value_vector[i] = NVal[i];
   for (i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);

   NVal.clear();
   G2L.clear();
   delete fem;
   delete face;
}

/**************************************************************************
 ROCKFLOW - Funktion: DomainIntegration
 Task:  Translate distributed source term within elements into nodes value
 for all kinds of element
 Programming:
 08/2005 WW Re-Implementation
 09/2010 TF re structured some things
 **************************************************************************/
void DomainIntegration(CRFProcess* m_pcs, const std::vector<long>&nodes_in_dom,
std::vector<double>&node_value_vector)
{
   CFEMesh* msh = m_pcs->m_msh;
   double nodesFVal[8];

   std::cout << "\tDomain integration\n";

   int Axisymm = 1;                               // ani-axisymmetry
   if (msh->isAxisymmetry())
      Axisymm = -1;                               // Axisymmetry is true
   CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());
   vec<CNode*> e_nodes(20);

   const size_t this_number_of_nodes (nodes_in_dom.size());
   const size_t nSize (msh->nod_vector.size());
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

   for (size_t i = 0; i < nSize; i++)
   {
      msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   for (size_t i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      G2L[nodes_in_dom[i]] = i;
   }

   size_t count = 0;
   for (size_t i = 0; i < msh->ele_vector.size(); i++)
   {
      CElem* elem (msh->ele_vector[i]);
      if (!elem->GetMark())
         continue;
      elem->GetNodes(e_nodes);
      size_t nn = elem->GetNodesNumber(msh->getOrder());
      count = 0;
      for (size_t j = 0; j < nn; j++)
      {
         for (size_t k = 0; k < this_number_of_nodes; k++)
         {
            if (*e_nodes[j] == *msh->nod_vector[nodes_in_dom[k]])
            {
               count++;
               break;
            }
         }
      }
      if (count != nn)
         continue;
      for (size_t j = 0; j < nn; j++)
         nodesFVal[j] = node_value_vector[G2L[e_nodes[j]->GetIndex()]];
      fem->ConfigElement(elem, m_pcs->m_num->ele_gauss_points, true);
      fem->setOrder(msh->getOrder() + 1);
      fem->FaceIntegration(nodesFVal);
      for (size_t j = 0; j < nn; j++)
         NVal[G2L[e_nodes[j]->GetIndex()]] += nodesFVal[j];
   }

   for (size_t i = 0; i < this_number_of_nodes; i++)
      node_value_vector[i] = NVal[i];
   for (size_t i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);

   NVal.clear();
   G2L.clear();
   e_nodes.resize(0);
   delete fem;
}

void EdgeIntegration(CFEMesh* msh, const std::vector<long>&nodes_on_ply,
std::vector<double>&node_value_vector,
FiniteElement::DistributionType dis_type, FiniteElement::PrimaryVariable prim_val,
bool flag_ignore_axisymmetry, bool flag_is_bc, int scaling_mode)
{
   size_t i, j; 
   int k, l;
   size_t this_number_of_nodes;
   size_t elemsCnode;
   int nedges, ii;
   vec<CNode*> e_nodes(3);
   vec<CEdge*> e_edges(12);

   double Jac = 0.0;
   double Weight = 0.0;
   double eta = 0.0;
   double v1, v2, radius = 0.0;
   double Shfct[3];

   double area_projection (1.0); //for projection of element areas for edges not parallel to the coordinate axes
   bool Const = false;
   if (dis_type == FiniteElement::CONSTANT || dis_type == FiniteElement::CONSTANT_NEUMANN)
      Const = true;

   CElem* elem = NULL;
   CEdge* edge = NULL;
   CNode* node = NULL;

   size_t nSize = msh->nod_vector.size();
   this_number_of_nodes = nodes_on_ply.size();
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

// CB THMBM
   // CB added to remove bug with deactivated Subdomains
   //for(i=0;i<(long)pcs_vector.size();i++){
   //  if(pcs_vector[i]->getProcessType()==this->getProcessType())
   //  {
   //    pcs_vector[i]->CheckMarkedElement();
   //    break;
   //  }
   //}

   // Unmakr edges.
   for (i = 0; i <  msh->edge_vector.size(); i++)
      msh->edge_vector[i]->SetMark(false);
   for (i = 0; i < nSize; i++)
   {
      msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   // Search edges on polyline
   for (i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      k = nodes_on_ply[i];
      G2L[k] = i;
      node = msh->nod_vector[k];
      elemsCnode = (int) node->getConnectedElementIDs().size();

      for (j = 0; j < elemsCnode; j++)
      {
         l = msh->nod_vector[k]->getConnectedElementIDs()[j];
         elem = msh->ele_vector[l];
         nedges = elem->GetEdgesNumber();
         elem->GetEdges(e_edges);
         for (ii = 0; ii < nedges; ii++)
         {
            edge = e_edges[ii];

            if (edge->GetMark())
               continue;
            edge->GetNodes(e_nodes);
            // Edge A
            if (*node == *e_nodes[0])
               e_nodes[0]->SetMark(true);
            // Edge B
            if (*node == *e_nodes[1])
               e_nodes[1]->SetMark(true);
            if (msh->getOrder())                  // Quadratic
            {
               if (*node == *e_nodes[2])
                  e_nodes[2]->SetMark(true);
            }
            if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
            {
               if (msh->getOrder())
               {
                  if (e_nodes[2]->GetMark())
                  {
                     edge->SetMark(true);
                  }
               }
               else
               {
            	   edge->SetMark(true);
               }
            }
         } // e_edges
      }  //
   }

#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
   const size_t id_act_l_max = static_cast<size_t>(msh->getNumNodesLocal());
   const size_t id_act_h_min =  msh->GetNodesNumber(false);
   const size_t id_act_h_max =  msh->getLargestActiveNodeID_Quadratic();

   struct loc_function
   {
      static bool isIDinRange( const size_t n_id, const size_t id_max0,
                                   const size_t id_min1,  const size_t id_max1 )
      {
	return (n_id < id_max0) || (n_id >= id_min1 && n_id <  id_max1);
      }
   };
#endif

   for (i = 0; i < msh->edge_vector.size(); i++)
   {
	   edge = msh->edge_vector[i];
	   if (!edge->GetMark())
		   continue;
	   edge->GetNodes(e_nodes);

	   if (flag_is_bc)
	   {
		   area_projection=AreaProjection(edge, prim_val);

	   }

      if (msh->getOrder())                        // Quadradic shape functions
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()
            && e_nodes[2]->GetMark())
         {
            Jac = 0.5 * edge->getLength()*area_projection;
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (Const && (!msh->isAxisymmetry() || flag_ignore_axisymmetry))
            {
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
              if(loc_function::isIDinRange( e_nodes[0]->GetIndex(), id_act_l_max, id_act_h_min, id_act_h_max))
#endif
               NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1 / 3.0;
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	      if(loc_function::isIDinRange( e_nodes[1]->GetIndex(), id_act_l_max, id_act_h_min, id_act_h_max))
#endif
               NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1 / 3.0;
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	      if(loc_function::isIDinRange( e_nodes[2]->GetIndex(), id_act_l_max, id_act_h_min, id_act_h_max))
#endif
               NVal[G2L[e_nodes[2]->GetIndex()]] += 4.0 * Jac * v1 / 3.0;

            }
            else
            {
               for (k = 0; k < 3; k++)            // Three nodes
               {
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
                 if( !loc_function::isIDinRange(  e_nodes[k]->GetIndex(), id_act_l_max, id_act_h_min, id_act_h_max))
                    continue;
#endif
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLineHQ(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry())
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 3; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->getData()[0];
                        Weight *= radius * 6.283185307179586;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;

                  }

               }
            }
         }
      } else                                      // Linear shape functions
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
         {
            Jac = 0.5 * edge->getLength()*area_projection;
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (!msh->isAxisymmetry()  && !flag_ignore_axisymmetry)
            {
               if (Const)
               {
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	          if(loc_function::isIDinRange( e_nodes[0]->GetIndex(), id_act_l_max,
                                            id_act_h_min, id_act_h_max))
#endif
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1;
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	          if(loc_function::isIDinRange( e_nodes[1]->GetIndex(), id_act_l_max,
                                            id_act_h_min, id_act_h_max))
#endif
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1;
               }
               else
               {
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	         if(loc_function::isIDinRange( e_nodes[0]->GetIndex(), id_act_l_max,
                                            id_act_h_min, id_act_h_max))
#endif
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * (2.0 * v1
                     + v2) / 3.0;
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	         if(loc_function::isIDinRange( e_nodes[1]->GetIndex(), id_act_l_max,
                                                id_act_h_min, id_act_h_max))
#endif
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * (v1 + 2.0
                     * v2) / 3.0;
               }
            } else                                // Axisymmetry
            {

               for (k = 0; k < 2; k++)            // Three nodes
               {
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
                  if( !loc_function::isIDinRange( e_nodes[k]->GetIndex(), id_act_l_max,
                                                  id_act_h_min, id_act_h_max) )
                    continue;
#endif
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLine(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry() && !flag_ignore_axisymmetry)
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 2; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->getData()[0];
                        Weight *= radius * 6.283185307179586;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;
                  }
               }
            }                                     // End of is (!axi)
         }
      }
   }
   for (i = 0; i < this_number_of_nodes; i++)
   {
      node_value_vector[i] = NVal[i];
      node = msh->nod_vector[nodes_on_ply[i]];
   }
   for (i = 0; i < msh->edge_vector.size(); i++)
      msh->edge_vector[i]->SetMark(true);
   for (i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);
   NVal.clear();
   G2L.clear();
   e_nodes.resize(0);
   e_edges.resize(0);

   //////////////////////////////////////////////////////////////////////////
   if(scaling_mode == 1)  // JOD 2021-08-05
   {
	   std::vector<double> scaling_vector(this_number_of_nodes, 0.);

	   for (i = 0; i < this_number_of_nodes; i++)
	   {
			 node = msh->nod_vector[nodes_on_ply[i]];
			 int nelem = 0;

			 for (j = 0; j < node->getConnectedElementIDs().size(); j++)
			 {
				elem = msh->ele_vector[node->getConnectedElementIDs()[j]];

				std::vector<size_t> node_indices;
				elem->getNodeIndices(node_indices);
				//std::cout << "ndx: " << elem->GetPatchIndex() << std::endl;
				int nnodes = 0;

				for(size_t ndx=0; ndx<node_indices.size(); ndx++)
				{
					if(std::find(nodes_on_ply.begin(), nodes_on_ply.end(), node_indices[ndx]) != nodes_on_ply.end())
						nnodes++;
				}

				if(nnodes == 2)
				{
					// take permeability first component (x-direction) for scaling
					scaling_vector[i] += mmp_vector[elem->GetPatchIndex()]->permeability_tensor[0];
					nelem++;
				}

			 }
			 scaling_vector[i] /= nelem;
	   }  // end for nodes

	   double divisor= 0.;
	   for(i=0;  i < this_number_of_nodes; i++)
	   {
		   divisor += scaling_vector[i] * node_value_vector[i];
	   }
	   divisor /= std::accumulate(node_value_vector.begin(),  node_value_vector.end(), 0.);  // divide by total st


	   for (i = 0; i < this_number_of_nodes; i++)
	   {
		   node_value_vector[i] *= scaling_vector[i] / divisor;
	   }
   }  // end if scaling mode 1 (scaling with permeability)

}

double AreaProjection(MeshLib::CEdge *edge, FiniteElement::PrimaryVariable primaryVariable)
{
	double area_projection (0);
//	const double epsilon (1.0e-8); //tolerance to decide whether projection is used
	//compute edge normal vector
	double edge_normal[3];//edge integration is 2 dimensional
	double elemNormalVector[3];

	elemNormalVector[0] = 0;
	elemNormalVector[1] = 0;
	elemNormalVector[2] = 1;

	edge->SetNormalVector(elemNormalVector, edge_normal); // Attenzione! This only works if nodes are numbered counter-clockwise.
	//TODO: distinguish coordinate systems!
	//if coordinate axes are at an angle with edges
//	if (fabs(edge_normal[0]) > epsilon && fabs(edge_normal[1]) > epsilon) {
		if (primaryVariable == FiniteElement::DISPLACEMENT_X)
			area_projection = edge_normal[0];
		else if (primaryVariable == FiniteElement::DISPLACEMENT_Y)
			area_projection = edge_normal[1];
//	}

	return area_projection;
}

}                                                 // end namespace FiniteElement
