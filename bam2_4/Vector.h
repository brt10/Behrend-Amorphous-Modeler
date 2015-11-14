/* Alicia Klinvex
April 22, 2008

Stores vector in three dimensions
*/

#ifndef VECTOR_H
#define VECTOR_H

struct Vector
{
   Vector(const double _x = 0, const double _y = 0, const double _z = 0)
   {
	  x = _x;
	  y = _y;
	  z = _z;
   }

   Vector operator+(const Vector &v) const
   {
	  return Vector(x + v.x, y + v.y, z + v.z);
   }

   Vector operator-(const Vector &v) const
   {
	  return Vector(x - v.x, y - v.y, z - v.z);
   }

   double getLength() const
   {
	  return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
   }

   double x, y, z;
};

#endif
