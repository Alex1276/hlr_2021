#include <stdio.h>

// Definieren Sie ein enum cardd
typedef enum{
  N=1,          //    1
  E=2,          //   10
  S=4,          //  100
  W=8 } cardd;  // 1000

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ cardd enthält
static int map[3][3] = {};

// Die Funktion set_dir soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  if ((x >=0 && x<= 2) && (y >= 0 && y<= 2) && (dir != 5) && (dir != 10) && (dir <= 12))
  {
    map[x][y] = dir;
  }
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  for(int i=0;i<3;i++) 
  {
    for(int j=0;j<3;j++)
    {
      
      switch(map[i][j])
      {
        case 1 :
          printf("N ");
          break;
        case 2 :
          printf("E ");
          break;
        case 4 :
          printf("S ");
          break;
        case 8 :
          printf("W ");
          break;
        case 3 :
          printf("NE");
          break;
        case 6 :
          printf("SE");
          break;
        case 12 :
          printf("SW");
          break;
        case 9 :
          printf("NW");
          break;
        default:
          printf("0 ");
      }
        if (j != 2) 
        {
          printf("  ");
        }
    }
    printf("\n");
  }
}

// In dieser Funktion darf nichts verändert werden!
int main (void)
{
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);
	set_dir(1, 3, N|S|E);
	set_dir(1, 1, N|S|E|W);

	show_map();

	return 0;
}
