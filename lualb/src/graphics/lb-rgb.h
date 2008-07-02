#ifndef LB_RGB_H
#define LB_RGB_H

#include <stdio.h>
#include "lb-macros.h"
#include "lb-colormap.h"

LB_BEGIN_DECLS

typedef struct _LBRGB LBRGB;

LBRGB* LBRGB_new(int width, int height);
void LBRGB_destroy(LBRGB*);

int LBRGB_width(const LBRGB*);
int LBRGB_height(const LBRGB*);

void _LBRGB_fill(LBRGB*, const double rgb[3]);
void _LBRGB_get_pixel(const LBRGB*, int x, int y, double rgb[3]);
void _LBRGB_set_pixel(LBRGB*, int x, int y, const double rgb[3]);
void _LBRGB_set_pixel_rgba(LBRGB*, int, int, const double rgba[4]);
void _LBRGB_save(const LBRGB*, FILE*);

void LBRGB_map_value(LBRGB*, const LBColormap*, int, int, double);

LB_END_DECLS

#endif /* LB_RGB_H */
