/* Alicia Klinvex
April 22, 2008

Stores a set of coordinates in three dimensions
*/

#ifndef POINT_H
#define POINT_H

#include<cmath>

struct Point
{
   Point(const double _x=0, const double _y=0, const double _z=0)
   {
	  x = _x;
	  y = _y;
	  z = _z;
   }

   double distanceTo(const Point &p) const
   {
	  return sqrt(pow(x-p.x, 2) + pow(y-p.y, 2) + pow(z-p.z, 2));
   }

   double x, y, z;
};

#endif
