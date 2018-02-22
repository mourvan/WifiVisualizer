#include "Tracer.h"
#include "stdio.h"
#include <fstream>

//glm::IntersectLineTriangle -  для проверки пересечения луча с поверхностью - тупо цикл по всем поверхностям
//TODO: C++ main - check if router is outside room;

int main(int argc, char** argv)
{
  CTracer tracer;
  CScene scene;

  int xRes;//Image resolution
  int yRes;
  float rx,ry,rz; //router coords
  float strength=100000; //signal strength
  std::string modelfile; //obj file path
  float cx=12500,cy=1000,cz=25000; //camera coord
  //4000 -7000 25000
  // - 
  //12500 1000 0
  //85 80 -250
  float ox=0,oy=0,oz=-1; //camera direction
  float ix=0,iy=1,iz=0; //"up" camera vector - determines incline
  float fovx = 90, fovy = 74; //fov
  bool autocontrast = false; float fraction = 0.01;
  bool medianfilter = false;
  if(argc != 2)
  {
    std::cout << "No config file specified, aborting" <<std::endl;
    return 0;
  }

  try
  {
  
    std::ifstream fin(argv[1]);
    std::string temp1, temp2;

    fin >> xRes >> yRes;
    fin >> modelfile;
    fin >> rx >> ry >> rz;
    fin >> strength;
    fin >> cx >> cy >> cz;
    fin >> ox >> oy >> oz;
    fin >> ix >> iy >> iz;
    fin >> fovx >> fovy;
    fin >> temp1;
    fin >> temp2;
    if (fin.fail())
      throw(1);

    if(temp1 == "TRUE") autocontrast = true;
    if(temp2 == "TRUE") medianfilter = true; 
  }
  catch(...)
  {
    std::cout << "Error while reading from file, aborting\n";
    return 0;
  }

  /*================================================================================================*/
  /*================================================================================================*/
  /*================================================================================================*/

  try
  {
    std::cout << "Loading scene . . . " << std::flush;
    scene.LoadScene(modelfile, glm::vec3(rx,ry,rz), strength);
    std::cout << "OK!\n" << std::flush;

    std::cout << "Creating voxel grid . . . " << std::flush;
    scene.CreateGrid();
    std::cout << "OK!\n"<< std::flush;

    std::cout << "Tracing Wi-Fi rays across the room . . . " << std::flush;
    scene.ComputeSignalStrength();
    std::cout << "OK!\n" << std::flush;

    if (medianfilter)
    {
      std::cout << "Applying median filter to voxel grid . . . " << std::flush;
      scene.ApplyMedianFilter();
      std::cout << "OK!\n" << std::flush;
    }
    else
    {
      std::cout << "Applying box filter to voxel grid . . . " << std::flush;
    	scene.ApplyBoxFilter();
      std::cout << "OK!\n" << std::flush;
    }

    scene.DeleteRoof();
   
    tracer.m_pScene = &scene;
    tracer.SetCamera(glm::vec3(cx,cy,cz) , glm::vec3(ox,oy,oz) , glm::vec3(ix,iy,iz) , glm::vec2(fovx, fovy));

    std::cout << "Rendering image from camera . . . " << std::flush;
    tracer.RenderImage(xRes, yRes);
    std::cout << "OK!\n" << std::flush;

    if (autocontrast)
    {
      std::cout << "Applying autocontrast filter to the image . . . " << std::flush;
    	tracer.Autocontrast(fraction);
      std::cout << "OK!\n" << std::flush;
    }
    std::cout << "Saving image . . . " << std::flush;
    std::string cfpath(argv[1]);
    tracer.SaveImageToFile("img/Result(" + cfpath + ").bmp");
    std::cout << "OK!\n" << std::flush;
  }
  catch(...)
  {
    std::cout << "Something went wrong, aborting\n" << std::flush;
  }
  return 0;
}

