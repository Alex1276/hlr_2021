/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

// PROBLEM: buf wird allocated wenn wir aus dem scope gehen
int *
mistake1 (void) 
{
  int buf[] = { 1, 1, 2, 3, 4, 5 };
  //int *bufpointer = (int *)&buf; 
  int * buffspace = malloc(6 * sizeof(int)); 
  
  for (int i = 0;i < 6;i++) {
    buffspace[i] = buf[i];
  }
  
  return buffspace;
}

int *
mistake2 (void) 
{
  int *buf = malloc (sizeof (int) * 2);
  buf[1] = 2;
  return buf;
}

int *
mistake3 (void)
{
  /* In dieser Funktion darf kein Speicher direkt allokiert werden. */
  //int mistake2_ = 0;
  int *buf = mistake2();
  buf[0] = 3;
  return buf;
}

int *
mistake4 (void) // DONE
{ 
  int *buf = malloc (sizeof (int));
  buf[0] = 4;
  //free (buf);
  return buf;
}

int
main (void)
{
  /* Modifizieren Sie diese Zeile nicht! */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1: %d\n", *p[0]);
  printf ("2: %d\n", *p[1]); // DONE
  printf ("3: %d\n", *p[2]);
  printf ("4: %d\n", *p[3]); // DONE

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free (p[0]-1);
  free (p[1]-1);
  free (p[2]);
  free (p[3]);

  			/* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */

  return 0;
}
