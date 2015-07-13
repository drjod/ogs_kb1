/**
 * \file XmlGmlInterface.cpp
 * 2011/11/23 KR as derived class from XMLInterface
 */

#include "XmlGmlInterface.h"

#include <QFile>
#include <QTextCodec>
#include <QtXml/QDomDocument>

namespace FileIO
{

XmlGmlInterface::XmlGmlInterface(ProjectData* project, const std::string &schemaFile)
: XMLInterface(project, schemaFile)
{
}

int XmlGmlInterface::readFile(const QString &fileName)
{
	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();
	std::string gliName("[NN]");

	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XmlGmlInterface::readFile() - Can't open xml-file " <<
		fileName.toStdString() << "." << "\n";
		delete file;
		return 0;
	}
	if (!checkHash(fileName))
	{
		delete file;
		return 0;
	}

	std::vector<GEOLIB::Point*>* points    = new std::vector<GEOLIB::Point*>;
	std::vector<GEOLIB::Polyline*>* polylines = new std::vector<GEOLIB::Polyline*>;
	std::vector<GEOLIB::Surface*>* surfaces  = new std::vector<GEOLIB::Surface*>;

	std::map<std::string, size_t>* pnt_names  = new std::map<std::string, size_t>;
	std::map<std::string, size_t>* ply_names  = new std::map<std::string, size_t>;
	std::map<std::string, size_t>* sfc_names  = new std::map<std::string, size_t>;

	QDomDocument doc("OGS-GLI-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //OpenGeoSysGLI
	if (docElement.nodeName().compare("OpenGeoSysGLI"))
	{
		std::cout << "XmlGmlInterface::readFile() - Unexpected XML root." << "\n";
		delete file;
		return 0;
	}

	QDomNodeList geoTypes = docElement.childNodes();

	for (int i = 0; i < geoTypes.count(); i++)
	{
		const QDomNode type_node(geoTypes.at(i));
		if (type_node.nodeName().compare("name") == 0)
			gliName = type_node.toElement().text().toStdString();
		else if (type_node.nodeName().compare("points") == 0)
		{
			readPoints(type_node, points, pnt_names);
			geoObjects->addPointVec(points, gliName, pnt_names);
		}
		else if (type_node.nodeName().compare("polylines") == 0)
			readPolylines(type_node, polylines, points,
			              geoObjects->getPointVecObj(gliName)->getIDMap(), ply_names);
		else if (type_node.nodeName().compare("surfaces") == 0)
			readSurfaces(type_node, surfaces, points,
			             geoObjects->getPointVecObj(gliName)->getIDMap(), sfc_names);
		else
			std::cout << "Unknown XML-Node found in file." << "\n";
	}
	delete file;

	if (!polylines->empty())
		geoObjects->addPolylineVec(polylines, gliName, ply_names);
	if (!surfaces->empty())
		geoObjects->addSurfaceVec(surfaces, gliName, sfc_names);
	return 1;
}

void XmlGmlInterface::readPoints( const QDomNode &pointsRoot,
                               std::vector<GEOLIB::Point*>* points,
                               std::map<std::string, size_t>* pnt_names )
{
	char* pEnd;
	QDomElement point = pointsRoot.firstChildElement();
	while (!point.isNull())
	{
		if (point.hasAttribute("id") && point.hasAttribute("x") && point.hasAttribute("y"))
		{
			_idx_map.insert (std::pair<size_t,
			                           size_t>(strtol((point.attribute("id")).
			                                          toStdString().c_str(), &pEnd,
			                                          10), points->size()));
			double zVal = (point.hasAttribute("z")) ? strtod((point.attribute(
			                                                          "z")).toStdString(
			                                                         ).c_str(),
			                                                 0) : 0.0;
			GEOLIB::Point* p =
			        new GEOLIB::Point(strtod((point.attribute("x")).toStdString().c_str(),
			                                 0),
			                          strtod((point.attribute("y")).
			                                 toStdString().c_str(), 0),
			                          zVal);
			if (point.hasAttribute("name"))
				pnt_names->insert( std::pair<std::string,
				                             size_t>(point.attribute("name").
				                                     toStdString(), points->size()) );
			points->push_back(p);
		}
		else
			std::cout <<
			"XmlGmlInterface::readPoints() - Attribute missing in <point> tag ..." <<
			"\n";

		point = point.nextSiblingElement();
	}
	if (pnt_names->empty())
		pnt_names = NULL;             // if names-map is empty, set it to NULL because it is not needed
}

void XmlGmlInterface::readPolylines( const QDomNode &polylinesRoot,
                                  std::vector<GEOLIB::Polyline*>* polylines,
                                  std::vector<GEOLIB::Point*>* points,
                                  const std::vector<size_t> &pnt_id_map,
                                  std::map<std::string, size_t>* ply_names )
{
	size_t idx(0);
	QDomElement polyline = polylinesRoot.firstChildElement();
	while (!polyline.isNull())
	{
		if (polyline.hasAttribute("id"))
		{
			idx = polylines->size();
			polylines->push_back(new GEOLIB::Polyline(*points));

			if (polyline.hasAttribute("name"))
				ply_names->insert( std::pair<std::string,
				                             size_t>(polyline.attribute("name").
				                                     toStdString(), idx) );

			QDomElement point = polyline.firstChildElement();
			while (!point.isNull())
			{
				(*polylines)[idx]->addPoint(pnt_id_map[_idx_map[atoi(point.text().
				                                                     toStdString().
				                                                     c_str())]]);
				point = point.nextSiblingElement();
			}
		}
		else
			std::cout <<
			"XmlGmlInterface::readPolylines() - Attribute missing in <polyline> tag ..."
			          <<
			"\n";

		polyline = polyline.nextSiblingElement();
	}
	if (ply_names->empty())
		ply_names = NULL;             // if names-map is empty, set it to NULL because it is not needed
}

void XmlGmlInterface::readSurfaces( const QDomNode &surfacesRoot,
                                 std::vector<GEOLIB::Surface*>* surfaces,
                                 std::vector<GEOLIB::Point*>* points,
                                 const std::vector<size_t> &pnt_id_map,
                                 std::map<std::string, size_t>* sfc_names )
{
	QDomElement surface = surfacesRoot.firstChildElement();
	while (!surface.isNull())
	{
		if (surface.hasAttribute("id"))
		{
			surfaces->push_back(new GEOLIB::Surface(*points));

			if (surface.hasAttribute("name"))
				sfc_names->insert( std::pair<std::string,
				                             size_t>(surface.attribute("name").
				                                     toStdString(),
				                                     surfaces->size() -
				                                     1) );

			QDomElement element = surface.firstChildElement();
			while (!element.isNull())
			{
				if (element.hasAttribute("p1") && element.hasAttribute("p2") &&
				    element.hasAttribute("p3"))
				{
					size_t p1 =
					        pnt_id_map[_idx_map[atoi((element.attribute("p1")).
					                                 toStdString().c_str())]];
					size_t p2 =
					        pnt_id_map[_idx_map[atoi((element.attribute("p2")).
					                                 toStdString().c_str())]];
					size_t p3 =
					        pnt_id_map[_idx_map[atoi((element.attribute("p3")).
					                                 toStdString().c_str())]];
					surfaces->back()->addTriangle(p1,p2,p3);
				}
				else
					std::cout <<
					"XmlGmlInterface::readSurfaces() - Attribute missing in <element> tag ..."
					          << "\n";
				element = element.nextSiblingElement();
			}
		}
		else
			std::cout <<
			"XmlGmlInterface::readSurfaces() - Attribute missing in <surface> tag ..." <<
			"\n";

		surface = surface.nextSiblingElement();
	}
	if (sfc_names->empty())
		sfc_names = NULL;             // if names-map is empty, set it to NULL because it is not needed
}

int XmlGmlInterface::write(std::ostream& stream)
{
	if (this->_exportName.empty())
	{
		std::cout << "Error in XmlGmlInterface::write() - No geometry specified..." << "\n";
		return 0;
	}

	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();
	size_t nPoints = 0, nPolylines = 0, nSurfaces = 0;

	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysGLI.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-GML-DOM");
	QDomElement root = doc.createElement("OpenGeoSysGLI");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.org" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" );

	doc.appendChild(root);

	QDomElement geoNameTag = doc.createElement("name");
	root.appendChild(geoNameTag);
	QDomText geoNameText = doc.createTextNode(QString::fromStdString(_exportName));
	geoNameTag.appendChild(geoNameText);

	// POINTS
	QDomElement pointsListTag = doc.createElement("points");
	root.appendChild(pointsListTag);

	const GEOLIB::PointVec* pnt_vec (geoObjects->getPointVecObj(_exportName));
	if (pnt_vec)
	{
		const std::vector<GEOLIB::Point*>* points (pnt_vec->getVector());

		if (!points->empty())
		{
			nPoints = points->size();
			for (size_t i = 0; i < nPoints; i++)
			{
				QDomElement pointTag = doc.createElement("point");
				pointTag.setAttribute("id", QString::number(i));
				pointTag.setAttribute("x", QString::number((*(*points)[i])[0], 'f'));
				pointTag.setAttribute("y", QString::number((*(*points)[i])[1], 'f'));
				pointTag.setAttribute("z", QString::number((*(*points)[i])[2], 'f'));

				std::string point_name;
				if (pnt_vec->getNameOfElementByID(i, point_name))
					pointTag.setAttribute("name", QString::fromStdString(point_name));

				pointsListTag.appendChild(pointTag);
			}
		}
		else
		{
			std::cout << "Point vector empty, abort writing geometry." << "\n";
			return 0;
		}
	}
	else
	{
		std::cout << "Point vector empty, abort writing geometry." << "\n";
		return 0;
	}

	// POLYLINES
	const GEOLIB::PolylineVec* ply_vec (geoObjects->getPolylineVecObj(_exportName));
	if (ply_vec)
	{
		const std::vector<GEOLIB::Polyline*>* polylines (ply_vec->getVector());

		if (polylines)
		{
			if (!polylines->empty())
			{
				QDomElement plyListTag = doc.createElement("polylines");
				root.appendChild(plyListTag);
				nPolylines = polylines->size();
				for (size_t i = 0; i < nPolylines; i++)
				{
					QDomElement polylineTag = doc.createElement("polyline");
					polylineTag.setAttribute("id", QString::number(i));

					std::string ply_name("");
					if (ply_vec->getNameOfElementByID(i, ply_name))
						polylineTag.setAttribute("name", QString::fromStdString(ply_name));

					plyListTag.appendChild(polylineTag);

					nPoints = (*polylines)[i]->getNumberOfPoints();
					for (size_t j = 0; j < nPoints; j++)
					{
						QDomElement plyPointTag = doc.createElement("pnt");
						polylineTag.appendChild(plyPointTag);
						QDomText plyPointText = doc.createTextNode(QString::number(((*polylines)[i])->getPointID(j)));
						plyPointTag.appendChild(plyPointText);
					}
				}
			}
			else
				std::cout << "Polyline vector empty, no polylines written to file." << "\n";
		}
	}
	else
		std::cout << "Polyline vector empty, no polylines written to file." << "\n";


	// SURFACES
	const GEOLIB::SurfaceVec* sfc_vec (geoObjects->getSurfaceVecObj(_exportName));
	if (sfc_vec)
	{
		const std::vector<GEOLIB::Surface*>* surfaces (sfc_vec->getVector());

		if (surfaces)
		{
			if (! surfaces->empty())
			{
				QDomElement sfcListTag = doc.createElement("surfaces");
				root.appendChild(sfcListTag);
				nSurfaces = surfaces->size();
				for (size_t i = 0; i < nSurfaces; i++)
				{
					QDomElement surfaceTag = doc.createElement("surface");
					surfaceTag.setAttribute("id", QString::number(i));

					std::string sfc_name("");
					if (sfc_vec->getNameOfElementByID(i, sfc_name))
						surfaceTag.setAttribute("name", QString::fromStdString(sfc_name));

					sfcListTag.appendChild(surfaceTag);

					// writing the elements compromising the surface
					size_t nElements = ((*surfaces)[i])->getNTriangles();
					for (size_t j = 0; j < nElements; j++)
					{
						QDomElement elementTag = doc.createElement("element");
						elementTag.setAttribute("p1", QString::number((*(*(*surfaces)[i])[j])[0]));
						elementTag.setAttribute("p2", QString::number((*(*(*surfaces)[i])[j])[1]));
						elementTag.setAttribute("p3", QString::number((*(*(*surfaces)[i])[j])[2]));
						surfaceTag.appendChild(elementTag);
					}
				}
			}
			else
				std::cout << "Surface vector empty, no surfaces written to file." << "\n";
		}
	}
	else
		std::cout << "Surface vector empty, no surfaces written to file." << "\n";


	//insertStyleFileDefinition(filename);
	std::string xml = doc.toString().toStdString();
	stream << xml;

	return 1;
}

}
