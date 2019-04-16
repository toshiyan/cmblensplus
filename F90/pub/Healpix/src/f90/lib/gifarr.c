/* gifarr: GIF input routine built upon GD library;
           intended for exclusive use within F90 wrapper [gifget];
   M. Bartelmann, Dec. 10, 1998 */

#include "chopst.h"
#include "gd.h"
#include <stdio.h>
#include <stddef.h>

#ifdef RS6000
void gifarr (
#else
void gifarr_ (
#endif
	      int *x, int *nx, int *ny,
	      int *r, int *g, int *b, int *nc, char *fn, int n) {
  FILE       *in;
  gdImagePtr im;
  int        i, j;
  char       *f;

  /* chop off last character from file name */

  f  = (char *)chopst(fn, n);

  /* read gif image from file */

  in = fopen(f, "rb");
  im = gdImageCreateFromGif(in);
  fclose(in);

  /* read color table */

  for (i = 0; i < *nc; i++) {
    r[i] = gdImageRed(im, i);
    g[i] = gdImageGreen(im, i);
    b[i] = gdImageBlue(im, i);
  }

  for (i = 0; i < *nx; i++) {
    for (j = 0; j < *ny; j++) {
      x[i + j * *nx] = gdImageGetPixel(im, i, j);
    }
  }

  /* destroy image */

  gdImageDestroy(im);
}
