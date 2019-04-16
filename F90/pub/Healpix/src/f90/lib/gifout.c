/* gifout: GIF output routine built upon GD library;
           intended for exclusive use within F90 wrapper [gifmap];
   M. Bartelmann, Nov. 25, 1998 */

#include "gd.h"
#include "chopst.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef RS6000
void gifout  (
#else
void gifout_ (
#endif
	      int *x, int *nx, int *ny,
	      int *r, int *g, int *b, int *nc, char *fn, int n) {

  gdImagePtr im;
  FILE       *out;
  int        i, j, c;
  int        *ct;
  char       *f;

  /* chop off last character from file name */

  f  = (char *)chopst(fn, n);

  /* create a new image with appropriate dimensions */

  im = gdImageCreate(*nx, *ny);

  /* allocate color table */

  ct = (int *)malloc(*nc*sizeof(int));
  for (i = 0; i < *nc; i++) {
    ct[i] = gdImageColorAllocate(im, r[i], g[i], b[i]);
  }

  /* set pixel colors in image */

  for (i = 0; i < *nx; i++) {
    for (j = 0; j < *ny; j++) {
      c = x[i + j * *nx];
      gdImageSetPixel(im, i, j, ct[c]);
    }
  }

  /* write image to gif file */

  out = fopen(f, "wb");
  gdImageGif(im, out);
  fclose(out);

  /* destroy image */

  gdImageDestroy(im);
}
