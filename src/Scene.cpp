#include "Scene.h"
#include <time.h> 
/*
bool CScene::LoadScene(std::string filename, glm::vec3 routerCoords, float sstr)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	std::vector<glm::vec3> glmvertices;

	routerPos = routerCoords;
	signalStrength = sstr;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str());


	if(ret) //reparse tinyobj format into glm format
	{
		//std::cout << attrib.vertices.size() << std::endl;
		for (size_t i = 0; i < attrib.vertices.size()/3; ++i)
		{
			tinyobj::real_t vx = attrib.vertices[3*i+0];
 			tinyobj::real_t vy = attrib.vertices[3*i+1];
			tinyobj::real_t vz = attrib.vertices[3*i+2];
			glm::vec3 temp(vx,vy,vz);
			glmvertices.push_back(temp);
		}
		mesh.m_vertices = glmvertices;
		return ret;
	}
	else
	{
		return ret;
	}
}
*/

bool CScene::LoadScene(std::string filename, glm::vec3 routerCoords, float sstr)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	std::vector<glm::vec3> glmvertices;

	routerPos = routerCoords;
	signalStrength = sstr;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str());
	if(ret == false)
		throw(1);

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) 
	{
	  // Loop over faces(polygon)
	  size_t index_offset = 0;
	  for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) 
	  {
	    int fv = shapes[s].mesh.num_face_vertices[f];
	    SPolygon polygon;
	    // Loop over vertices in the face.
	    for (size_t v = 0; v < fv; v++) 
	    {
	      // access to vertex
	      tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
	      tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
	      tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
	      tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
	      polygon.m_vertices.push_back(glm::vec3(vx,vy,vz));
	    }
	    //std::cout << "Number of vertices in polygon " << polygon.m_vertices.size() << std::endl;
	    index_offset += fv;
	    mesh.m_polygons.push_back(polygon);
	  }
	}
	//std::cout << "Number of polygons" << mesh.m_polygons.size() << std::endl;
	return ret;
}

glm::vec3 CScene::FindMaxRoomCoord()
{
	float xmax = std::numeric_limits<float>::min(), ymax = std::numeric_limits<float>::min(), zmax = std::numeric_limits<float>::min();
	for (size_t i = 0; i < mesh.m_polygons.size(); ++i)
	{
		for(size_t j = 0; j < mesh.m_polygons[i].m_vertices.size(); ++j)
		{
			if (mesh.m_polygons[i].m_vertices[j].x > xmax) xmax = mesh.m_polygons[i].m_vertices[j].x;
			if (mesh.m_polygons[i].m_vertices[j].y > ymax) ymax = mesh.m_polygons[i].m_vertices[j].y;
			if (mesh.m_polygons[i].m_vertices[j].z > zmax) zmax = mesh.m_polygons[i].m_vertices[j].z;
		}
	}
	//std::cout << xmax << ", " << ymax << ", " << zmax << std:: endl;
	glm::vec3 ret(xmax,ymax,zmax);
	return ret;
}

glm::vec3 CScene::FindMinRoomCoord()
{
	float xmax = std::numeric_limits<float>::max(), ymax = std::numeric_limits<float>::max(), zmax = std::numeric_limits<float>::max();
	for (size_t i = 0; i < mesh.m_polygons.size(); ++i)
	{
		for(size_t j = 0; j < mesh.m_polygons[i].m_vertices.size(); ++j)
		{
			if (mesh.m_polygons[i].m_vertices[j].x < xmax) xmax = mesh.m_polygons[i].m_vertices[j].x;
			if (mesh.m_polygons[i].m_vertices[j].y < ymax) ymax = mesh.m_polygons[i].m_vertices[j].y;
			if (mesh.m_polygons[i].m_vertices[j].z < zmax) zmax = mesh.m_polygons[i].m_vertices[j].z;
		}
	}
	//std::cout << xmax << ", " << ymax << ", " << zmax << std:: endl;
	glm::vec3 ret(xmax,ymax,zmax);
	return ret;
}

void CScene::CreateGrid()
{
	glm::vec3 min = FindMinRoomCoord();
	glm::vec3 max = FindMaxRoomCoord();
	//centering the room to start from (0,0,0)
	//(need to substract min coord from current position to get corresponding grid voxel)
	glm::vec3 len = (max-min);
	float cubeSideLen = std::min(std::min(len.x,len.y), len.z)/100.0;
	len.x = len.x / cubeSideLen;
	len.y = len.y / cubeSideLen;
	len.z = len.z / cubeSideLen;

	voxelgrid.strech = cubeSideLen;
	voxelgrid.offset = min;

	std::vector<std::vector<std::vector<float> > > temp(len.x+1, std::vector<std::vector<float> >(len.y+1, std::vector<float>(len.z+1)));

	#pragma omp parallel for collapse(3)
	for(size_t i = 0; i < temp.size(); ++i)
		for(size_t j = 0; j < temp[0].size(); ++j)
			for(size_t k = 0; k < temp[0][0].size(); ++k)
				temp[i][j][k] = 0;

	voxelgrid.grid = temp;
}


void CScene::ComputeSignalStrength() //magick
{
	const float EPS = 0.00001;

	#pragma omp parallel for//what about strength grid? rand() - reinitialise seed with thread num;
	for(int i = 0; i < 50000; i++) //MAKING 10000 RAYS, TODO - PARALLEL
	{
		//std::cout << std::endl;	
		float rayStrength = signalStrength;
		float d = 1;

		glm::vec3 curPos = routerPos;
		glm::vec3 nextPos;
		glm::vec3 intersectPos;
		glm::vec3 tempintersectPos;
		glm::vec3 direction = GenerateRandomVector();

		//std::cout << "INITIAL("<< curPos[0] << ", " << curPos[1] << ", " << curPos[2] << ")" << std::endl;
		//std::cout << "DIRECTION("<< direction[0] << ", " << direction[1] << ", " << direction[2] << ")" << std::endl;
		try
		{
			while(rayStrength > EPS)
			{
				bool isIntersect = false;
				size_t j1;
				float distance;
				float mindistance = std::numeric_limits<float>::max();
				//rayStrength = rayStrength - (EPS*d*d);
				//d+=0.05;
				for(size_t j = 0; j < mesh.m_polygons.size(); j++) //CHECKING FOR INTERSECTION
				{
					//float distance;
					//float mindistance = 1000000;
					if (glm::intersectRayTriangle(curPos, direction, mesh.m_polygons[j].m_vertices[0], mesh.m_polygons[j].m_vertices[1], mesh.m_polygons[j].m_vertices[2], tempintersectPos) == true) 
					{


						//intersectPos[2] = 1 - intersectPos[0] - intersectPos[1];
						isIntersect = true;
						tempintersectPos = curPos + direction*tempintersectPos[2];
						distance = glm::length(tempintersectPos - curPos);
						if(distance < mindistance)
						{
							j1=j;
							mindistance = distance;
							intersectPos = tempintersectPos;
						}
						//intersectPos = intersectPos[0]*mesh.m_vertices[3*j+0] + intersectPos[1]*mesh.m_vertices[3*j+1] + intersectPos[2]*mesh.m_vertices[3*j+2];
				
						//TODO Change direction;
					}
				}

				if(isIntersect) //TODO: decrease ray strength at reflection
				{
					//std::cout << "INTERSECTION FOUND at point " << "("<< intersectPos[0] << ", " << intersectPos[1] << ", " << intersectPos[2] << ")" << std::endl;
					while(rayStrength > EPS)
					{
						if(rayStrength > voxelgrid(curPos.x, curPos.y, curPos.z))
							voxelgrid(curPos.x, curPos.y, curPos.z) = rayStrength;

						rayStrength = rayStrength - (EPS*d*d);
						d+=0.05;
						if(glm::length(intersectPos - curPos) <= 1) //BUG (FIXED)- direction reflection does not work - ALL NORMALS SHOULD BE kinda (0,0,1), seems like parsing is incorrect...
						{
							rayStrength = 0.7 * rayStrength;
							curPos = intersectPos;
							glm::vec3 trNormal = glm::triangleNormal(mesh.m_polygons[j1].m_vertices[0], mesh.m_polygons[j1].m_vertices[1], mesh.m_polygons[j1].m_vertices[2]);
							//std::cout << "NORMAL at " << "("<< trNormal[0] << ", " << trNormal[1] << ", " << trNormal[2] << ")" << std::endl;
							direction = glm::reflect(direction, trNormal);

							curPos += direction;
							//std::cout << "REFLECTING at " << "("<< direction[0] << ", " << direction[1] << ", " << direction[2] << ")" << std::endl;
							break;
						}
						curPos += direction;
					}
				}
				else
				{
					while(rayStrength > EPS)
					{
						if(rayStrength > voxelgrid(curPos.x, curPos.y, curPos.z))
							voxelgrid(curPos.x, curPos.y, curPos.z) = rayStrength;
						rayStrength = rayStrength - (EPS*d*d);
						d+=0.05;
						curPos += direction;
					}
				}
			} 
		}
		catch(std::out_of_range)
		{
			//std::cout << "Ray is out_of_bounds\n";
		}
	}
	/*
	for(size_t i = 0; i < voxelgrid.grid.size(); ++i)
		for(size_t j = 0; j < voxelgrid.grid[0].size(); ++j)
			for(size_t k = 0; k < voxelgrid.grid[0][0].size(); ++k)
				if(voxelgrid.grid[i][j][k] != 0)
					pix++;
	std::cout << pix << std::endl;
	*/
}

float randMToN(float M, float N, float seed) // skommunizdil s SOF - TODO: REWORK RAND LOGIC
{
	//srand (time(NULL));
	//srand(int(time(NULL)) ^ omp_get_thread_num());
    return M + (rand() / ( RAND_MAX / (N-M) ) );  
}

glm::vec3 CScene::GenerateRandomVector()
{
	float xrand = randMToN(-1000,1000,1);
	float yrand = randMToN(-1000,1000,xrand);
	float zrand = randMToN(-1000,1000,yrand);
	return glm::normalize(glm::vec3(xrand, yrand, zrand));
}

void CScene::DeleteRoof()
{
	//1) find triangle with maximum z coord
	//2) if normal is parallel to (0,0,1) - delete this triangle from mesh
	glm::vec3 max = FindMaxRoomCoord();
	for(size_t j = 0; j < mesh.m_polygons.size(); j++)
	{
		if(mesh.m_polygons[j].m_vertices[0].z == max.z && mesh.m_polygons[j].m_vertices[1].z == max.z && mesh.m_polygons[j].m_vertices[2].z == max.z)
		{
				mesh.m_polygons.erase(mesh.m_polygons.begin() + j); 
				j--;
		}
	}
}

void CScene::ApplyBoxFilter()
{
	std::vector<std::vector<std::vector<float> > > grid = voxelgrid.grid;
	for(size_t i = 1; i < voxelgrid.grid.size()-1; ++i)
		for(size_t j = 1; j < voxelgrid.grid[0].size()-1; ++j)
			for(size_t k = 1; k < voxelgrid.grid[0][0].size()-1; ++k)			
			{
				float temp = 0;
				//#pragma omp parallel for collapse(3) reduction(+:temp);
				for(int i1 = -1; i1 <= 1; ++i1)
					for(int j1 = -1; j1 <= 1; ++j1)
						for(int k1 = -1; k1 <= 1; ++k1)
							temp += voxelgrid.grid[i+i1][j+j1][k+k1];
				temp = temp/27.0f;
				//std::cout << temp << "\n";
				grid[i][j][k] = temp;
			}

	voxelgrid.grid = grid;
}


void sort(float a[], int len)
{
	//cout << "hello" << endl;
	for(int i=0;i<=len;i++)
		for(int j=i;j<=len;j++)
			if(a[i]<a[j])
			{
				int t = a[i];
				a[i]=a[j];
				a[j]=t;
			}
}

void CScene::ApplyMedianFilter() //TODO
{
	int radius = 1;
	int radius1=100*radius;
	float med[radius1]; 
	std::vector<std::vector<std::vector<float> > > grid = voxelgrid.grid;
	for(size_t i = radius; i < voxelgrid.grid.size()-radius; ++i)
		for(size_t j = radius; j < voxelgrid.grid[0].size()-radius; ++j)
			for(size_t k = radius; k < voxelgrid.grid[0][0].size()-radius; ++k)
			{
				int h=0;
				for(int i1=-radius; i1<=radius; ++i1)
					for(int j1=-radius; j1<=radius; ++j1)
						for(int k1=-radius; k1<=radius; ++k1)
						{
							med[h] = voxelgrid.grid[i+i1][j+j1][k+k1];
							h++;
						}
				h--;
				sort(med, h);
				grid[i][j][k] = med[h/2];
			}
	voxelgrid.grid = grid;
}






