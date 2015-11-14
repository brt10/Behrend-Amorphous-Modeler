/* Alicia Klinvex
April 22, 2008

Stores the data for one atom
*/

#ifndef ATOM_H
#define ATOM_H

#include "Partner.h"
#include "Vector.h"
#include <vector>
using namespace std;

struct Atom
{
	// default and explicit value constructor
	Atom(const Point & p = Point(0,0,0), const int _id = 0, const unsigned int el = 0)
	{
		coordinates = p;
		id = _id;
		atomsInShell = 0;
		force = Vector(0,0,0);
		acceleration = Vector(0,0,0);
		radialEnergy = 0;
		angularEnergy = 0;
		repulsiveEnergy = 0;
		type = el;
		shellNumber = 0;
	}

	/*Atom(const Atom & a)
	{
		coordinates = a.coordinates;
		id = a.id;
		atomsInShell = a.atomsInShell;
		force = a.force;
		acceleration = a.acceleration;
		radialEnergy = a.radialEnergy;
		angularEnergy = a.angularEnergy;
		repulsiveEnergy = a.repulsiveEnergy;
		type = a.type;
		shellNumber = a.shellNumber;
	}*/

	// returns true if this has a lower ID than atom a
	bool hasLowerID(const Atom &a) const
	{
		return (id < a.id);
	}

	// returns true if this is in a sparser area than atom a
	bool isMoreDesparate(const Atom & a) const
	{
		return (atomsInShell < a.atomsInShell);
	}

	// returns true if two atoms are the same
	bool operator==(const Atom & a) const
	{
		return (id == a.id);
	}

	Point coordinates;
	unsigned int id;
	int atomsInShell; // describes density
	int shellNumber;  // describes proximity to bond swap
	vector<Partner> bonds;
	Vector force;
	vector<size_t> closestNeighbors; // can contribute a repulsive force
	//element type;
	unsigned int type;

	double angularEnergy, radialEnergy, repulsiveEnergy;
	Vector acceleration;
};

#endif
