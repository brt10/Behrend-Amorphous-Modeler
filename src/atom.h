#ifndef ATOM_H
#define ATOM_H

#include <math.h>
#include <string>
#include <vector>

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
		int id;
		string type;
		vector<atom> nn;
		vector<atom> bonds;
	public:
		atom();
		atom(const double x, const double y, const double z, const int id, const string type);
		void setID(int id);
		void setType(string type);
		int getID();
		string getType();
		void addNN(atom nn);
		vector<atom> getNN();
		int getBonds();
		void bond(atom b);
		bool operator==(atom b);
		string atomString();
};

#endif
