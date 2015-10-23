#ifndef MISC_FEMTOOLKITS_H
#define MISC_FEMTOOLKITS_H

#include <string>

namespace  MeshLib
{
   class CFEMesh;
}

void writeOGSMesh(const MeshLib::CFEMesh &mesh, const std::string file_base_name);
#endif
