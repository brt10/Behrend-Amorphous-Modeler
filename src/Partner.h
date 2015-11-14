/* Alicia Klinvex
   April 22, 2008

   Stores the data for one bond partner
*/

#ifndef PARTNER_H
#define PARTNER_H

struct Partner
{
   Partner(const Point &p = Point(0,0,0), const unsigned int _id = 0, const double _dist = 0)
   {
	  coordinates = p;
	  id = _id;
	  dist = _dist;
   }

   bool operator==(const Partner &p) const
   {
	  return (id == p.id);
   }

   Point coordinates;
   unsigned int id;
   double dist;
};

#endif
