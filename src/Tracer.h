#pragma once

#include "glm/glm.hpp"
#include "EasyBMP/EasyBMP.h"
#include "Types.h"
#include "Scene.h"

#include "string"
//#include "atlimage.h"

class CTracer
{
public:
	
  SRay MakeRay(glm::uvec2 pixelPos);  // Create ray for specified pixel - 
  glm::vec3 TraceRay(SRay ray); // Trace ray, compute its color
  void RenderImage(int xRes, int yRes);
  void SaveImageToFile(std::string fileName);
  void Autocontrast(float fraction);
  void SetCamera(glm::vec3 m_pos, glm::vec3 m_forward, glm::vec3 m_up, glm::vec2 m_viewAngle);
  //CImage* LoadImageFromFile(std::string fileName);

public:
  SCamera m_camera;
  CScene* m_pScene;
};

/*

Визуализируем модель помещения и полученную воксельную сетку с помощью обратной трассировки лучей из камеры. 
Для этого в каждом пикселе камеры генерируем луч. 

Его направление определяется ориентацией камеры, её пирамидой видимости и положением пикселя на матрице. 

Каждому лучу присваивается цвет RGB (изначально возьмём его равным 0), который затем запишется в пиксель. 
Далее трассируем луч через сцену до его пересечения с первым полигональным объектом. 
Цвет этого объекта, рассчитанный по модели Фонга, присваивается лучу. 

После этого с помощью того же метода Ray marching проходим вдоль луча от пересечения с объектом до камеры. 
В каждом очередном вокселе нужно взять значение мощности сигнала, найти соответствующий ему цвет из цветовой схемы 
и произвести альфа-смешивание между текущим цветом луча и полученным цветом. 
Так как эта операция будет выполняться много раз вдоль луча, коэффициент непрозрачности α у воксельной сетки должен быть достаточно мал, 
иначе сетка будет выглядеть абсолютно непрозрачной.


*/