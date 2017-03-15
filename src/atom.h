#ifndef ATOM_H
#define ATOM_H

#include <math.h>

using namespace std;

class point {
	protected:
		double x, y, z;
	public:
		point();
		point(const int x, const int y, const int z);
		double distanceTo(const point &p) const;
		double distanceTo(const int x, const int y, const int z) const;
		double getX() const;
		double getY() const;
		double getZ() const;
};

class atom: public point {
	private:
		int id, type;
	public:
		atom();
		atom(const int x, const int y, const int z, const int id, const int el);
};

#endif
