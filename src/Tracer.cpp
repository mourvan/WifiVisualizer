#include "Tracer.h"

using namespace glm;

//glm::perspective
//glm::lookAt

SRay CTracer::MakeRay(glm::uvec2 pixelPos) //optimise
{
    SRay ray;
    ray.m_start = m_camera.m_pos;
    ray.m_dir = m_camera.m_forward + ((pixelPos.x + 0.5f) / m_camera.m_resolution.x - 0.5f)* m_camera.m_right + ((pixelPos.y + 0.5f) / m_camera.m_resolution.y - 0.5f) * m_camera.m_up;
    ray.m_dir = glm::normalize(ray.m_dir);
    return ray;
}

glm::vec3 CTracer::TraceRay(SRay ray) //if no intersection return constant color
{
  const float alpha = 0.001;
  glm::vec3 intersectPos;
  glm::vec3 tempintersectPos;
  glm::vec3 color(200, 200, 200);
  glm::vec3 curPos;


  bool isIntersect = false;
  size_t j1;
  float distance;
  float mindistance = std::numeric_limits<float>::max();

  for(size_t j = 0; j < m_pScene->mesh.m_polygons.size(); j++) //CHECKING FOR INTERSECTION
  {
        
          if (glm::intersectRayTriangle(ray.m_start, ray.m_dir, m_pScene->mesh.m_polygons[j].m_vertices[0], m_pScene->mesh.m_polygons[j].m_vertices[1], m_pScene->mesh.m_polygons[j].m_vertices[2], tempintersectPos) == true) 
          {
            isIntersect = true;
            tempintersectPos = ray.m_start + ray.m_dir*tempintersectPos[2];
            distance = glm::length(tempintersectPos - ray.m_start);
            if(distance < mindistance)
            {
              j1=j;
              mindistance = distance;
              intersectPos = tempintersectPos;
            }
          }
  }

  if(isIntersect)
  {
    //determining polygon color - camera is light source
    //return a vec3 of pixel color
    // find out own color of the last wall

    //backwards tracing
    glm::vec3 clr(1, 1, 1);
    float brightness = 100.0;
    glm::vec3 trNormal = glm::triangleNormal(m_pScene->mesh.m_polygons[j1].m_vertices[0], m_pScene->mesh.m_polygons[j1].m_vertices[1], m_pScene->mesh.m_polygons[j1].m_vertices[2]);
    float prod = glm::dot(ray.m_dir, trNormal);
    brightness = std::abs(brightness*prod);
    clr = clr * brightness;
    color = clr;
    curPos = intersectPos;

      try
      {
        while(glm::length(ray.m_start - curPos) > 1)
        {
          float red = m_pScene->voxelgrid(curPos.x, curPos.y, curPos.z) / m_pScene->signalStrength * 255.0f;
          float blue = (255.0f-red);
          float green = 0.2f * color.y;

          glm::vec3 clr(red,green,blue);
          color = (1-alpha) * color + alpha * clr;
          curPos = curPos - ray.m_dir;
        }
      }
      catch(std::out_of_range)
      {
        
      }
  }
  return color;
}

void CTracer::RenderImage(int xRes, int yRes)
{
  // Rendering
  m_camera.m_resolution = uvec2(xRes, yRes);
  std::vector<std::vector<glm::vec3> > pixels(yRes, std::vector<glm::vec3>(xRes));

  #pragma omp parallel for collapse(2)
  for(int i = 0; i < yRes; i++)
    for(int j = 0; j < xRes; j++)
    {
      SRay ray = MakeRay(uvec2(j, i));
      pixels[i][j] = TraceRay(ray);
    }
   m_camera.m_pixels = pixels;
}

void CTracer::SaveImageToFile(std::string fileName)
{
  BMP out;
    out.SetSize(m_camera.m_pixels[0].size(), m_camera.m_pixels.size());

    unsigned int r, g, b;
    RGBApixel p;
    p.Alpha = 255;
    #pragma omp parallel for collapse(2) private(r,g,b,p)
    for (size_t i = 0; i < m_camera.m_pixels.size(); ++i) 
    {
        for (size_t j = 0; j < m_camera.m_pixels[0].size(); ++j) 
        {
            r = m_camera.m_pixels[i][j].x;
            g = m_camera.m_pixels[i][j].y;
            b = m_camera.m_pixels[i][j].z;
            p.Red = r; p.Green = g; p.Blue = b;
            out.SetPixel(j, i, p);
        }
    }

    if (!out.WriteToFile(fileName.c_str()))
        std::cout << "Error writing file " << std::endl;
}


void CTracer::SetCamera(glm::vec3 mpos, glm::vec3 mforward, glm::vec3 mup, glm::vec2 mviewAngle)
{
  m_camera.m_pos = mpos;
  m_camera.m_forward = mforward;
  m_camera.m_up = mup;
  m_camera.m_viewAngle = mviewAngle;

  //glm::vec3 ViewDir = m_camera.m_forward;
  //glm::vec3 Up = m_camera.m_up;
  m_camera.m_right = glm::cross(m_camera.m_forward, m_camera.m_up);
 // std::cout << "RIGHT at point " << "("<< Right[0] << ", " << Right[1] << ", " << Right[2] << ")" << std::endl;
  m_camera.m_up = glm::cross(m_camera.m_forward, m_camera.m_right);
  //Streching vectors
  float ViewLen = glm::length(m_camera.m_forward);
  m_camera.m_up = (ViewLen * glm::tan(glm::radians(m_camera.m_viewAngle.y/2.0f)) / glm::length(m_camera.m_up)) *  m_camera.m_up;
  m_camera.m_right = (ViewLen * glm::tan(glm::radians(m_camera.m_viewAngle.x/2.0f)) / glm::length(m_camera.m_right)) * m_camera.m_right;

}

float normalise(float con)
{
  if (con > 255.0f)
    con = 255.0f;
  if (con < 0.0f)
    con = 0.0f;
  return con;
}

void CTracer::Autocontrast(float fraction) //no need to compute histogram, checking borders only
{
  //m_camera.m_pixels
   float left=255, right=0;
   float r,g,b;
    // 1) COMPUTE BORDERING VALUES
    for (size_t i = 0; i<m_camera.m_pixels.size(); ++i)
    {
        for (size_t j = 0; j<m_camera.m_pixels[0].size(); ++j)
        {
            r = m_camera.m_pixels[i][j].x;
            g = m_camera.m_pixels[i][j].y;
            b = m_camera.m_pixels[i][j].z;
            float y = 0.2125f*r + 0.7154f*g + 0.0721f*b; //find max and min
            if (y < left) 
              left = y;
            if (y > right)
              right = y;
        }

    }
    // 2) COMPUTE COEFFICIENT

    //http://spatial-analyst.net/ILWIS/htm/ilwisapp/stretch_algorithm.htm - assisted

    float interval = right - left;
    right = right - interval*fraction;
    left = left + interval*fraction;
    interval = right - left;
    float c = 255.0/interval;
    // 3) MAKE CHANGES TO IMAGE
    for (size_t i = 0; i<m_camera.m_pixels.size(); ++i)
    {
        for (size_t j = 0; j<m_camera.m_pixels[0].size(); ++j)
        {
            float rd,gd,bd; 
            r = m_camera.m_pixels[i][j].x;
            g = m_camera.m_pixels[i][j].y;
            b = m_camera.m_pixels[i][j].z;
            rd = (r - left)*c;
            gd = (g - left)*c;
            bd = (b - left)*c;
            r = normalise(rd);
            g = normalise(gd);
            b = normalise(bd);
            m_camera.m_pixels[i][j].x = r;
            m_camera.m_pixels[i][j].y = g;
            m_camera.m_pixels[i][j].z = b;
        }
    }
}
