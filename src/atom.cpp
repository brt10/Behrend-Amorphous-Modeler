#include <math.h>
#include "atom.h"

using namespace std;

point::point() : point(0.0, 0.0, 0.0) {}

point::point(const int x, const int y, const int z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

double point::distanceTo(const point &p) const {
	return sqrt(pow(x - p.getX(), 2) + pow(y - p.getY(), 2) + pow(z - p.getZ(), 2));
}

double point::distanceTo(const int x, const int y, const int z) const {
	return sqrt(pow(x - this->x, 2) + pow(y - this->y, 2) + pow(z - this->z, 2));
}

void point::setX(double x) {
	this->x = x;
}

void point::setY(double y) {
	this->y = y;
}

void point::setZ(double z) {
	this->z = z;
}

double point::getX() const {
	return x;
}

double point::getY() const {
	return y;
}

double point::getZ() const {
	return z;
}

atom::atom() : atom(0, 0, 0, 0, 0) {}

atom::atom(const double x = 0, const double y = 0, const double z = 0, const int id = 0, const int el = 0) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = id;
	this->type = el;
}

void atom::setID(int id) {
	this->id = id;
}

void atom::setType(int type) {
	this->type = type;
}

int atom::getID() {
	return this->id;
}

int atom::getType() {
	return this->type;
}
