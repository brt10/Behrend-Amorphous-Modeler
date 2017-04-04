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
		void setX(double x);
		void setY(double y);
		void setZ(double z);
		double getX() const;
		double getY() const;
		double getZ() const;
};

class atom: public point {
	private:
		int id, type;
		atom *nn;
	public:
		atom();
		atom(const double x, const double y, const double z, const int id, const int el);
		void setID(int id);
		void setType(int type);
		int getID();
		int getType();
		void setNN(atom* nn);
};

#endif
