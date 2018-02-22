#pragma once

#include "glm/glm.hpp"
#include "glm/ext.hpp"
#include "vector"
#include <omp.h>
#include <iostream>

struct SRay
{
  glm::vec3 m_start;
  glm::vec3 m_dir;
};

struct SCamera
{
  glm::vec3 m_pos;          // Camera position and orientation
  glm::vec3 m_forward;      // Orthonormal basis
  glm::vec3 m_right;
  glm::vec3 m_up;

  glm::vec2 m_viewAngle;    // View angles, rad
  glm::vec2 m_resolution;  // Image resolution: w, h

  std::vector<std::vector<glm::vec3> > m_pixels;  // Pixel array

  
};

struct SPolygon
{
	std::vector<glm::vec3> m_vertices;
};

struct SMesh
{
  std::vector<glm::vec3> m_vertices;  // vertex positions
  std::vector<SPolygon> m_polygons;
 // std::vector<glm::uvec3> m_triangles;  // vetrex indices
  /*
  tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
  tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
  tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
  */
};

struct SGrid
{
	float strech;
	glm::vec3 offset;
	std::vector<std::vector<std::vector<float> > > grid;

	float& operator()(float x, float y, float z) 
	{
      return grid.at((x - offset.x)/strech).at((y - offset.y)/strech).at((z - offset.z)/strech);
    }
};

