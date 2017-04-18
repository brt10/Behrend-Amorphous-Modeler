#include <math.h>
#include <string>
#include <algorithm> // find
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

atom::atom() : atom(0, 0, 0, 0, "") {}

atom::atom(const double x = 0, const double y = 0, const double z = 0, const int id = 0, const string type = "") {
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = id;
	this->type = type;
}

bool atom::operator==(atom b) {
    if(this->x == b.getX() && this->y == b.getY() && this->z == b.getZ() && this->type == b.getType())
	return true;
    return false;
}

void atom::setID(int id) {
	this->id = id;
}

void atom::setType(string type) {
	this->type = type;
}

int atom::getID() {
	return this->id;
}

string atom::getType() {
	return this->type;
}

void atom::addNN(atom nn) {
    this->nn.push_back(nn);
}

vector<atom> atom::getNN() {
    return nn;
}

int atom::getBonds() {
    return (int)bonds.size();	   
}

void atom::bond(atom b) {
    vector<atom>::iterator a;
    a = find(bonds.begin(), bonds.end(), b);
    if(a != bonds.end()) { // if b isn't found in a's bonds, add it
	bonds.push_back(b);
	b.bond(*this); // recursive call to have the other atom push_back the same bond
    }
}

string atom::atomString() {
    string s = string("Atom ") + string(this->id) + string(": ") + string(this->type) + string(" (") + string(this->x) + string(", ") + string(this->y) + string(", ") + string(this->z) + string(")");
    return s;
}
