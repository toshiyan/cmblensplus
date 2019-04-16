/* chopst: copy the first [n] characters of the string [s] into
           new string and return; used to chop off last character
	   of a string passed by Fortran routines;
   M. Bartelmann, Dec. 17, 1998 */

#include <stdlib.h>
#include "chopst.h"

char *chopst (char *s, int n) {
  int i;
  char *t, *p, *q;
  
  t = (char *)malloc((n+1)*sizeof(char));

  p = &s[0];
  q = &t[0];
  for (i = 0; i < n; i++) {*q++ = *p++;}
  *q = '\0';

  return t;
}
