#include "molecule.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstdio>

using namespace std;


int main()
{
  Molecule mol("geom", 0);
  mol.printCoord();
  mol.calcRmatrix();
  mol.printRmatrix();
  mol.printRlist();
  mol.calcUnitMatrix();
  //mol.printUnit();
  mol.printAngles();
  //mol.printOop();
  //mol.printDihedral();
  mol.centerOfMass();

	return 0;
}