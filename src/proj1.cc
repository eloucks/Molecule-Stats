#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstdio>

using namespace std;

int main()
{
  ifstream input("geom");
  int natom;
  input >> natom;

  int *an = new int[natom];
  double **coord = new double*[natom];
  for(int i=0; i < natom; i++)
  {
  	coord[i]= new double[3];
  }


  for(int i=0; i < natom; i++)
  {
   	input >> an[i] >> coord[i][0] >> coord[i][1] >> coord[i][2];
  }
  input.close();

  for (int i = 0; i<natom; i++)
  {
  	printf("%d %20.12f %20.12f %20.12f\n", (int) an[i], coord[i][0], coord[i][1], coord[i][2]);
  }

	return 0;
}
