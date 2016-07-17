#include <string>
 
class Molecule
{
  public:
    int natom;
    int charge;
    int *an;
    double **coord; //coord[atomnumber][xy or z]
    double **r;
    double **ex;
    double **ey;
    double **ez;
    std::string point_group;
 
    void printCoord();
    void printRmatrix();
    void printRlist();
    void printUnit();
    void printAngles();
    void printOop();
    void printDihedral();
    void rotate(double phi);
    void translate(double x, double y, double z);
    void calcRmatrix();
    void calcUnitMatrix();
    double angle(int atom1, int atom2, int atom3);
    double oop(int atom1, int atom2, int atom3, int atom4);
    double dihedral(int atom1, int atom2, int atom3, int atom4);
    void centerOfMass();
 
    Molecule(const char*, int);
    ~Molecule();
};