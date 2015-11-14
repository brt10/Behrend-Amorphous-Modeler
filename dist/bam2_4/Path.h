/* Roberto Guerra
April 22, 2013

Stores four integers to identify a path suitable for bond switching:

   c-d    c d
   |     /|/
 a-b    a b

instead of storing the four atoms id we store only the id of b,
the bond indexes of b corresponding to a and c, and the
bond index of c corresponding to d.

After switching the bonds the path will change but it will be still
valid for future switches.

 
*/

#ifndef PATH_H
#define PATH_H

struct Path
{
   Path(size_t _b, size_t _bond_ba, size_t _bond_bc, size_t _bond_cd)
   {
	  b = _b;
	  bond_ba = _bond_ba;
	  bond_bc = _bond_bc;
	  bond_cd = _bond_cd;
   }

   size_t b,bond_ba,bond_bc,bond_cd;
};

#endif
