#ifndef LB_RGB_H
#define LB_RGB_H

#include <stdio.h>
#include "lb-macros.h"
#include "lb-colormap.h"

LB_BEGIN_DECLS

typedef struct _LBRGB lb_rgb;

lb_rgb* lb_rgb_new(int width, int height);
void lb_rgb_destroy(lb_rgb*);
int lb_rgb_width(const lb_rgb*);
int lb_rgb_height(const lb_rgb*);
void lb_rgb_fill(lb_rgb*, const double rgb[3]);
void lb_rgb_get_pixel(const lb_rgb*, int x, int y, double rgb[3]);
void lb_rgb_set_pixel(lb_rgb*, int x, int y, const double rgb[3]);
void lb_rgb_set_pixel_rgba(lb_rgb*, int, int, const double rgba[4]);
void lb_rgb_save(const lb_rgb*, char*);
void lb_rgb_map_value(lb_rgb*, const lb_colormap*, int, int, double);

LB_END_DECLS

#endif /* LB_RGB_H */
