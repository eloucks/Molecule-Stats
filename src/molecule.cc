#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <math.h>
#include <string>

void Molecule::printCoord() // 0 = x, 1 = y, 2 = z
{
  printf("\n");
  printf("Atomic coordinates (cartesian): \n");
  for (int i = 0; i<natom; i++)
  {
  	printf("%d %20.12f %20.12f %20.12f\n", (int) an[i], coord[i][0], coord[i][1], coord[i][2]);
  }
}

void Molecule::printRmatrix() //r values are by default in borh
{
  printf("\n");
  printf("Internuclear distances (Angstroms): \n");
  printf("  ");
  for (int i = 1; i<natom+1; i++)
  {
    printf("%10.d", i);
  }
  printf("\n");
  for (int a2 = 0; a2<natom; a2++)
  {
    printf("%2.d", a2+1); 
    for (int a1 = 0; a1<natom; a1++) //was <i
    {
      printf("%10.6f", r[a1][a2]*0.529177); // 0.529177 converts to angstroms
    }
    printf("\n");
  }
}

void Molecule::printUnit()
{
  printf("\n");
  printf("Internuclear unit vectors: \n");
  printf("X:\n");
  for (int i = 1; i<natom+1; i++)
  {
    printf("%10.d", i);
  }
  printf("\n");
  for (int a1 = 0; a1<natom; a1++)
  {
    printf("%2.d", a1+1); 
    for (int a2 = 0; a2<natom; a2++)
    {
      printf("%10.6f", ex[a1][a2]);
    }
    printf("\n");
  }
  printf("\nY:\n");
  for (int i = 1; i<natom+1; i++)
  {
    printf("%10.d", i);
  }
  printf("\n");
  for (int a1 = 0; a1<natom; a1++)
  {
    printf("%2.d", a1+1); 
    for (int a2 = 0; a2<natom; a2++)
    {
      printf("%10.6f", ey[a1][a2]);
    }
    printf("\n");
  }
  printf("\nZ:\n");
  for (int i = 1; i<natom+1; i++)
  {
    printf("%10.d", i);
  }
  printf("\n");
  for (int a1 = 0; a1<natom; a1++)
  {
    printf("%2.d", a1+1); 
    for (int a2 = 0; a2<natom; a2++)
    {
      printf("%10.6f", ez[a1][a2]);
    }
    printf("\n");
  }
}

void Molecule::printRlist() 
{
  printf("\n");
  printf("Internuclear distances (Angstroms): \n");
  for (int a1 = 0; a1<natom; a1++) 
  {
    for (int a2 = a1+1; a2<natom; a2++)
    {
      printf("%3.d %3.d %10.6f \n", a1+1, a2+1, r[a1][a2]*0.529177); //0.529177 converts to angstroms
    }
  }
}

void Molecule::printAngles()
{
  printf("\nUnique Bond angles (degrees):\n");
  for(int i=0; i < natom; i++) {
    for(int j=0; j < natom; j++) {
      for(int k=0; k < natom; k++) {
        if(i != j && i != k && j != k && i<k && r[i][j]*0.529177 > 2 && r[j][k]*0.529177 > 2)
        {
          printf("%2d - %2d - %2d %12.6f\n", i+1, j+1, k+1, angle(i,j,k)*(180.0/acos(-1.0)));
        }
      }
    }
  }
}

void Molecule::printOop() //currently excludes angles containing bonds over 2 A
{
  printf("\nUnique OOP angles (degrees):\n");
  for(int i=0; i < natom; i++) {
    for(int j=0; j < natom; j++) {
      for(int k=0; k < natom; k++) {
        for (int l =0; l < natom; l++) {
          if(i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && j<l && r[k][j]*0.529177 > 2 && r[k][l]*0.529177 > 2)
          {
            printf("%2d in plane (%2d - %2d - %2d) %12.6f\n", i+1, j+1, k+1, l+1, oop(i,j,k,l)*(180.0/acos(-1.0)));
          }
        }
      }
    }
  }
}

void Molecule::printDihedral() //currently excludes angles containing bonds over 2 A
{
  printf("\nUnique Dihedral angles (degrees):\n");
  for(int i=0; i < natom; i++) {
    for(int j=0; j < natom; j++) {
      for(int k=0; k < natom; k++) {
        for (int l =0; l < natom; l++) {
          if(i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && r[i][j]*0.529177 > 2 && r[j][k]*0.529177 > 2 && r[k][l]*0.529177 > 2 && i<l) 
          {
            printf("(%2d - %2d - %2d - %2d) %12.6f\n", i+1, j+1, k+1, l+1, dihedral(i,j,k,l)*(180.0/acos(-1.0)));
          }
        }
      }
    }
  }
}

void Molecule::translate(double x, double y, double z) //whole molecule
{
  for(int i=0; i < natom; i++) {
     coord[i][0] += x;
     coord[i][1] += y;
     coord[i][2] += z;
  }
}

void Molecule::calcRmatrix() // Number of unique combinations of 2 is = natom(natom-1)(natom-2)(natom-3)/2
{
  for (int a1 = 0; a1 < natom; a1++)
  {
    for(int a2=a1+1; a2 < natom; a2++) 
    {
      r[a1][a2] = r[a2][a1] = sqrt(
            (coord[a1][0]-coord[a2][0])*(coord[a1][0]-coord[a2][0])
          + (coord[a1][1]-coord[a2][1])*(coord[a1][1]-coord[a2][1])
          + (coord[a1][2]-coord[a2][2])*(coord[a1][2]-coord[a2][2])
                      );
    }
  }
}

void Molecule::calcUnitMatrix()
{
  for(int i=0; i < natom; i++) {
    for(int j=i+1; j < natom; j++) {
      ex[i][j] = -(coord[i][0] - coord[j][0])/r[i][j];
      ey[i][j] = -(coord[i][1] - coord[j][1])/r[i][j];
      ez[i][j] = -(coord[i][2] - coord[j][2])/r[i][j];
      ex[j][i] = -1*ex[i][j];
      ey[j][i] = -1*ey[i][j];
      ez[j][i] = -1*ez[i][j];
    }
  }
}
 

double Molecule::angle(int a1, int a2, int a3) //comes out in radians
{
  double ans = acos(
                  ex[a2][a1]*ex[a2][a3]+
                  ey[a2][a1]*ey[a2][a3]+
                  ez[a2][a1]*ez[a2][a3]);
  return ans;
}

double Molecule::oop(int i, int j, int k, int l) //comes out in radians 
{
  double exx = (ey[k][j] * ez[k][l] - ez[k][j] * ey[k][l]);// exx is e cross for x
  double exy = (ez[k][j] * ex[k][l] - ex[k][j] * ez[k][l]);
  double exz = (ex[k][j] * ey[k][l] - ey[k][j] * ex[k][l]);
  double numerator = (exx * ex[k][i] + exy * ey[k][i] + exz * ez[k][i]);
  double pre = numerator/(sin(angle(j, k, l)));
  if(pre < -1.0) pre = asin(-1.0);
  else if(pre > 1.0) pre = asin(1.0);
  else pre = asin(pre);
  return pre;
}

double Molecule::dihedral(int i, int j, int k, int l) //comes out in radians
{
  double exx1 = (ey[i][j] * ez[j][k] - ez[i][j] * ey[j][k]);// exx is e cross for x, in the first cross product
  double exy1 = (ez[i][j] * ex[j][k] - ex[i][j] * ez[j][k]);
  double exz1 = (ex[i][j] * ey[j][k] - ey[i][j] * ex[j][k]);
  double exx2 = (ey[j][k] * ez[k][l] - ez[j][k] * ey[k][l]);
  double exy2 = (ez[j][k] * ex[k][l] - ex[j][k] * ez[k][l]);
  double exz2 = (ex[j][k] * ey[k][l] - ey[j][k] * ex[k][l]);
  double numerator = (exx1 * exx2 + exy1 * exy2 + exz1 * exz2);
  double ans = numerator/(sin(angle(i,j,k))*sin(angle(j,k,l)));
  if(ans < -1.0) ans = acos(-1.0);
  else if(ans > 1.0) ans = acos(1.0);
  else ans = acos(ans);
  if (oop(l,i,j,k)*(180.0/acos(-1.0))>0) ans = ans * -1;
  return ans;
}

void Molecule::centerOfMass()
{
  double mass[] = {0, 1.00782503223, 4.00260325413, 7.01600343659, 9.01218306, 11.009305355, 12, 14.00307400443, 15.99491461957, 18.99840316273, 19.99244017617, 22.98976928196, 23.985041697, 26.981538531, 27.97692653465, 30.97376199842, 31.97207117441, 34.968852682, 39.96238312372};
  double x, y, z, m = 0;
  for(int i = 0; i<natom; i++)
  {
    x += coord[i][0]*mass[an[i]];
    y += coord[i][1]*mass[an[i]];
    z += coord[i][2]*mass[an[i]];
    m += mass[an[i]];
  }
  x = x/m;
  y = y/m;
  z = z/m;
  printf("\nMolecular center of mass (bohr): %11.8f, %11.8f, %11.8f\n", x, y, z);
  translate(-x, -y, -z);
}




Molecule::Molecule(const char *filename, int q)
{
  charge = q;

  // open filename
  std::ifstream input(filename);
  assert(input.good());

  // read the number of atoms from filename
  input >> natom;

  // allocate space
  an = new int[natom];
  coord = new double*[natom];
  for(int i=0; i < natom; i++)
  {
  	coord[i]= new double[3];
  }
  for(int i=0; i < natom; i++)
  {
   	input >> an[i] >> coord[i][0] >> coord[i][1] >> coord[i][2];
  }

  //make bond lenghs matrix
  r = new double*[natom];
  for(int i=0; i < natom; i++)
  {
    r[i]= new double[natom];
  }


  //make unit vectors matricies
  ex = new double* [natom];
  ey = new double* [natom];
  ez = new double* [natom];
  for(int i=0; i < natom; i++) {
    ex[i] = new double[natom];
    ey[i] = new double[natom];
    ez[i] = new double[natom];
  }


  input.close();
}

Molecule::~Molecule()
{
  delete[] an;
  for(int i=0; i < natom; i++)
    delete[] coord[i];
  delete[] coord;
  for(int i=0; i < natom; i++)
    delete[] r[i];
  delete[] r;
  for(int i=0; i < natom; i++) {
    delete[] ex[i]; delete[] ey[i]; delete[] ez[i];
  }
  delete[] ex; delete[] ey; delete[] ez;
}
