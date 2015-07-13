
#include "FEFLOWInterface.h"

#include "StringTools.h"
#include "files0.h"
#include <QtXml>
#include <fstream>

#include "Polygon.h"
#include "msh_mesh.h"

using namespace MeshLib;
using namespace std;

/// FEFLOW Model Ascii File Format (*.fem) - Version 5.4
CFEMesh* FEFLOWInterface::readFEFLOWModelFile(const std::string &filename)
{
	std::ifstream feflow_file(filename.c_str());
	if (!feflow_file)
	{
		std::cerr << "error opening stream " << "\n";
		return NULL;
	}

	CFEMesh* m_msh = NULL;
	m_msh = new CFEMesh();

	FEFLOW_FEM_CLASS fem_class = {0,0,0,0,0,0,0,0};
	FEFLOW_FEM_DIM fem_dim = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	std::vector<GEOLIB::Point*>* points = NULL;
	std::vector<GEOLIB::Polyline*>* lines = NULL;

	std::string line_string;
	std::stringstream line_stream;
	CNode* m_nod = NULL;
	double x[12] = {};
	char a;

	bool isXZplane = false;

	while(!feflow_file.eof())
	{
		line_string = GetLineFromFile1(&feflow_file);
		//....................................................................
		// CLASS
		if(line_string.find("CLASS") != string::npos)
		{
			line_stream.str(GetLineFromFile1(&feflow_file));
			// problem class, time mode, problem orientation, dimension, nr. layers for 3D, saturation switch, precision of results, precision of coordinates
			line_stream >> fem_class.problem_class >> fem_class.time_mode >>
			fem_class.orientation >> fem_class.dimension >> fem_class.n_layers3d;
			line_stream.clear();
		}
		//....................................................................
		// DIMENS
		else if(line_string.compare("DIMENS") == 0)
		{
			// DIMENS
			line_stream.str(GetLineFromFile1(&feflow_file));
			line_stream >> fem_dim.n_nodes >> fem_dim.n_elements >>
			fem_dim.n_nodes_of_element >> ws;
			line_stream.clear();
			// Create nodes
			for(long i = 0; i < fem_dim.n_nodes; i++)
			{
				m_nod = new CNode(i);
				m_msh->nod_vector.push_back(m_nod);
			}
			m_msh->InitialNodesNumber();
		}
		//....................................................................
		// NODE
		else if(line_string.compare("NODE") == 0)
		{
			MshElemType::type t = MshElemType::INVALID;
			if (fem_dim.n_nodes_of_element==2) t = MshElemType::LINE;
			else if (fem_dim.n_nodes_of_element==3) t = MshElemType::TRIANGLE;
			else if (fem_dim.n_nodes_of_element==4 && fem_class.dimension==2) t = MshElemType::TRIANGLE;
			else if (fem_dim.n_nodes_of_element==4 && fem_class.dimension==3) t = MshElemType::TETRAHEDRON;
			else if (fem_dim.n_nodes_of_element==6 && fem_class.dimension==3) t = MshElemType::PRISM;
			else if (fem_dim.n_nodes_of_element==8 && fem_class.dimension==3) t = MshElemType::HEXAHEDRON;
			CElem* m_ele = NULL;
			for(long i = 0; i < fem_dim.n_elements; i++)
			{
				m_ele = new CElem();
				m_ele->SetElementType(t);
				m_ele->Read(feflow_file, 6);
				m_msh->ele_vector.push_back(m_ele);
			}
		}
		//....................................................................
		// COOR
		else if(line_string.compare("COOR") == 0)
		{
			const long no_nodes_per_layer =
			        (fem_class.dimension ==
			2) ? fem_dim.n_nodes : fem_dim.n_nodes / (fem_class.n_layers3d + 1);
			const long i_read = no_nodes_per_layer / 12 + 1;
			const long no_nodes_per_layer_1 = fem_class.n_layers3d * no_nodes_per_layer;
			// x, y
			for (int k = 0; k < 2; k++)
				// each point
				for(long i = 0; i < i_read; i++)
				{
					line_string = GetLineFromFile1(&feflow_file);
					line_stream.clear();
					line_stream.str(line_string);
					// each column
					for(int j = 0; j < 12; j++)
					{
						line_stream >> x[j] >> a;
						if (fem_class.dimension == 2)
						{
							long n = i * 12 + j;
							if (n >= fem_dim.n_nodes)
								break;
							m_nod = m_msh->nod_vector[n];
							if (k == 0)
								m_nod->SetX(x[j]);
							else
								m_nod->SetY(x[j]);
						}
						else
							for(int l = 0; l < fem_class.n_layers3d + 1;
							    l++)
							{
								long n = i * 12 + l *
								         no_nodes_per_layer + j;
								long nn = i * 12 +
								          no_nodes_per_layer_1 + j;
								if(nn >= fem_dim.n_nodes)
									break;
								m_nod = m_msh->nod_vector[n];
								if (k == 0)
									m_nod->SetX(x[j]);
								else
									m_nod->SetY(x[j]);
							}
					}
					line_stream.clear();
				}
		}
		//....................................................................
		// ELEV_I
		else if(line_string.compare("ELEV_I") == 0)
		{
			if (fem_class.dimension == 2)
				continue;
			const long no_nodes_per_layer = fem_dim.n_nodes / (fem_class.n_layers3d + 1);
			double z = .0;
			long n0 = 0;
			for(int l = 0; l < fem_class.n_layers3d + 1; l++)
			{
				if (l>0) {
					line_stream.clear();
					line_string = GetLineFromFile1(&feflow_file);
				}

				line_stream.clear();
				line_string = GetLineFromFile1(&feflow_file);
				line_stream.str(line_string);
				line_stream >> z >> n0;

				for(long i = 0; i < no_nodes_per_layer; i++)
				{
					//n = i+l*no_nodes_per_layer;
					long n = n0 - 1 + i + l * no_nodes_per_layer;
					m_nod = m_msh->nod_vector[n];
					m_nod->SetZ(z);
				}
			}
		}
		//....................................................................
		// EXTENTS
		else if(line_string.compare("EXTENTS") == 0)
		{
/*
      m_strInfo = "FEFLOW import: EXTENTS";
      pWin->SendMessage(WM_SETMESSAGESTRING,0,(LPARAM)(LPCSTR)m_strInfo);
      //0.00000000000000e+000,0.00000000000000e+000,7.99587000000000e+003,5.29027651306167e+003,
      line_string = GetLineFromFile1(&feflow_file);
      line_stream.clear();
      line_stream.str(line_string);
      line_stream >> x[0] >> x[1] >> x[2] >> x[3];
      //3.54101621147850e+003,3.26511822506216e+003,4.38757239488187e+003,3.81153176162252e+003,
      line_string = GetLineFromFile1(&feflow_file);
      line_stream.clear();
      line_stream.str(line_string);
      line_stream >> x[4] >> x[5] >> x[6] >> x[7];
      //
      m_pnt = new CGLPoint();
      m_pnt->id = (long)gli_points_vector.size();
      m_pnt->x = x[2];
      m_pnt->y = x[4];
      m_pnt->z = x[6];
      gli_points_vector.push_back(m_pnt);
      m_pnt = new CGLPoint();
      m_pnt->id = (long)gli_points_vector.size();
      m_pnt->x = x[3];
      m_pnt->y = x[5];
      m_pnt->z = x[7];
      gli_points_vector.push_back(m_pnt);
 */
		}
		//....................................................................
		// GRAVITY
		else if(line_string.compare("GRAVITY") == 0)
		{
			line_stream.str(GetLineFromFile1(&feflow_file));
			// vector x,y,z
			double vec[3] = {};
			line_stream >> vec[0] >> vec[1] >> vec[2];
			line_stream.clear();

			if (vec[0] == 0.0 && vec[1] == -1.0 && vec[2] == 0.0)
				// x-z plane
				isXZplane = true;
		}
		//....................................................................
		// SUPERMESH
		else if(line_string.compare("SUPERMESH") == 0)
			this->readSuperMesh(feflow_file, fem_class, &points, &lines);
		//....................................................................
	} // is line empty
	feflow_file.close();

	if (lines && lines->size() > 1)
	{
		m_msh->ConstructGrid();
		// set material ID
		this->setMaterialID(m_msh, lines);
	}

	if (isXZplane)
	{
		if (m_msh)
			for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
			{
				CNode* nod = m_msh->nod_vector[i];
				nod->SetZ(nod->getData()[1]);
				nod->SetY(0.0);
			}
		if (_geoObjects && points)
		{
			//const std::vector<GEOLIB::Point*> *points = _geoObjects->getPointVec("Feflow");
			for (size_t i = 0; i < points->size(); i++)
			{
				GEOLIB::Point* pt = (*points)[i];
				(*pt)[2] = (*pt)[1];
				(*pt)[1] = .0;
			}
			string str("Feflow");
			_geoObjects->addPointVec(points, str);
			_geoObjects->addPolylineVec(lines, str);
		}
	}

	std::cout << "done" << "\n";

	return m_msh;
}

void FEFLOWInterface::readPoints(QDomElement &nodesEle,  const std::string &tag, int dim, std::vector<GEOLIB::Point*> &points)
{
	QDomElement xmlEle = nodesEle.firstChildElement(QString::fromStdString(tag));
	if (xmlEle.isNull())
		return;
	QString str_pt_list1 = xmlEle.text().simplified();
	istringstream ss(str_pt_list1.toStdString());
	while (!ss.eof())
	{
		int pt_id = 0;
		double pt_xyz[3] = {};
		ss >> pt_id;
		for (int i = 0; i < dim; i++)
			ss >> pt_xyz[i];
		GEOLIB::Point* pnt = new GEOLIB::Point( pt_xyz[0],
		                                        pt_xyz[1],
		                                        pt_xyz[2]); //id?
		points[pt_id - 1] = pnt;
	}
}

//
void FEFLOWInterface::readSuperMesh(std::ifstream &feflow_file,
                                    FEFLOW_FEM_CLASS &fem_class,
                                    std::vector<GEOLIB::Point*>** p_points,
                                    std::vector<GEOLIB::Polyline*>** p_lines)
{
	/* Example
	   <?xml version="1.0" encoding="utf-8" standalone="no" ?>
	   <supermesh>
	   <nodes count="15">
	    <fixed>
	       1                 18   45.6363636363636
	       2   42.7878787878788   33.5151515151515
	       4   27.9393939393939   19.7575757575758
	       6   13.0909090909091   19.6363636363636
	       9   53.2727272727273   42.7878787878788
	      11   65.1515151515152   25.6363636363636
	      13   51.6969696969697    17.030303030303
	    </fixed>
	    <linear>
	       3   30.3939393939394   39.5757575757576
	       5   35.3636363636364   26.6363636363636
	       7   20.5151515151515   19.6969696969697
	       8   15.5454545454545   32.6363636363636
	      10    48.030303030303   38.1515151515152
	      12   59.2121212121212   34.2121212121212
	      14   58.4242424242424   21.3333333333333
	      15   39.8181818181818   18.3939393939394
	    </linear>
	   </nodes>
	   <polygons>
	    <polygon>
	      <nodes count="8"> 1  8  6  7  4  5  2  3 </nodes>
	    </polygon>
	    <polygon>
	      <nodes count="10"> 2  5  4 15 13 14 11 12  9 10 </nodes>
	    </polygon>
	   </polygons>
	   </supermesh>
	 */

	// extract XML strings
	ostringstream oss;
	string line_string;
	while (true)
	{
		line_string = GetLineFromFile1(&feflow_file);
		oss << line_string << "\n";
		if (line_string.find("</supermesh>") != string::npos)
			break;
	}
	QString strXML(oss.str().c_str());

	// convert string to XML
	QDomDocument doc;
	if (!doc.setContent(strXML))
	{
		// error
		std::cerr << "error illegal XML format for supermesh" << "\n";
		return;
	}

	// get geometry data from XML
	QDomElement docElem = doc.documentElement(); // #supermesh
	// #nodes
	// TODO: what is fixed/linear?
	*p_points = new std::vector<GEOLIB::Point*>();
	std::vector<GEOLIB::Point*>* points = *p_points;
	//std::vector<GEOLIB::Point*> *points = new std::vector<GEOLIB::Point*>();
	QDomElement nodesEle = docElem.firstChildElement("nodes");
	if (nodesEle.isNull())
		return;
	{
		QString str = nodesEle.attribute("count");
		const long n_points = str.toLong();
		points->resize(n_points);
		//fixed
		readPoints(nodesEle, "fixed", fem_class.dimension, *points);
		readPoints(nodesEle, "linear", fem_class.dimension, *points);
		readPoints(nodesEle, "parabolic", fem_class.dimension, *points);
	}
	//_geoObjects->addPointVec(points, string("Feflow"));

	// #polygons
	*p_lines = new std::vector<GEOLIB::Polyline*>();
	std::vector<GEOLIB::Polyline*>* lines = *p_lines;
	//std::vector<GEOLIB::Polyline*> *lines = new std::vector<GEOLIB::Polyline*>();
	QDomElement polygonsEle = docElem.firstChildElement("polygons");
	if (polygonsEle.isNull())
		return;
	{
		//GEORemoveAllPolylines(); //
		QDomNode child = polygonsEle.firstChild();
		while(!child.isNull())
		{
			if (child.nodeName() != "polygon")
			{
				child = child.nextSibling();
				continue;
			}
			QDomElement xmlEle = child.firstChildElement("nodes");
			if (xmlEle.isNull())
				continue;
			QString str = xmlEle.attribute("count");
			const long n_points = str.toLong();
			QString str_ptId_list = xmlEle.text().simplified();
			{
				GEOLIB::Polyline* line = new GEOLIB::Polyline(*points);
				lines->push_back(line);
				istringstream ss(str_ptId_list.toStdString());
				for (long i = 0; i < n_points; i++)
				{
					int pt_id = 0;
					ss >> pt_id;
					line->addPoint(pt_id - 1);
				}
				line->addPoint(line->getPointID(0));
			}
			child = child.nextSibling();
		}
		//_geoObjects->addPolylineVec(lines, "Feflow");
	}
}

void FEFLOWInterface::setMaterialID( CFEMesh* m_msh, std::vector<GEOLIB::Polyline*>* lines )
{
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		CElem* e = m_msh->ele_vector[i];
		double const* gpt = e->GetGravityCenter();
		int matId = 0;
		for (size_t j = 0; j < lines->size(); j++)
		{
			GEOLIB::Polyline* poly = (*lines)[j];
			if (!poly->isClosed())
				continue;

			GEOLIB::Polygon polygon(*poly, true);
			if (polygon.isPntInPolygon(gpt[0], gpt[1], gpt[2]))
			{
				matId = j;
				break;
			}
		}
		e->setPatchIndex(matId);
	}
}
