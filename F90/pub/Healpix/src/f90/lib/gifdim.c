/* gifdim: routine built upon GD library returning the dimensions of
           a GIF image file;
           intended for exclusive use within F90 wrapper [gifget];
   M. Bartelmann, Dec. 10, 1998 */

#include "chopst.h"
#include "gd.h"
#include <stdio.h>
#include <stddef.h>

#ifdef RS6000
void gifdim (
#else
void gifdim_ (
#endif
	      int *nx, int *ny, int *nc, char *fn, int n) {
  FILE       *in;
  gdImagePtr im;
  char       *f;

  /* chop last character off filename */

  f  = (char *)chopst(fn, n);

  /* open GIF file */

  in = fopen(f, "rb");
  im = gdImageCreateFromGif(in);
  fclose(in);

  /* get dimensions */

  *nx = gdImageSX(im);
  *ny = gdImageSY(im);
  *nc = gdImageColorsTotal(im);

  /* destroy image */

  gdImageDestroy(im);
}
