#pragma once



#include "Types.h"
#include "TinyObjLoader/tiny_obj_loader.h"

/*
std::string inputfile = "cornell_box.obj";
tinyobj::attrib_t attrib;
std::vector<tinyobj::shape_t> shapes;
std::vector<tinyobj::material_t> materials;
  
std::string err;
bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());
  
if (!err.empty()) { // `err` may contain warning message.
  std::cerr << err << std::endl;
}

if (!ret) {
  exit(1);
}
*/

class CScene
{
	glm::vec3 routerPos;

	glm::vec3 FindMaxRoomCoord();
	glm::vec3 FindMinRoomCoord();
	glm::vec3 GenerateRandomVector();

public:
	float signalStrength;
	SMesh mesh;
	SGrid voxelgrid;

	bool LoadScene(std::string filename, glm::vec3 routerCoords, float sstr);
	void CreateGrid();
	void ComputeSignalStrength();
	void DeleteRoof();
	void ApplyBoxFilter();
	void ApplyMedianFilter();

};





