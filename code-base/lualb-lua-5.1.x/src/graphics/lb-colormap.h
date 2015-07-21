#ifndef LB_COLORMAP_H
#define LB_COLORMAP_H

#include "lb-macros.h"

LB_BEGIN_DECLS

typedef struct _LBColormap lb_colormap;

lb_colormap* lb_colormap_new(void);
void lb_colormap_destroy(lb_colormap*);
int lb_colormap_num_colors(const lb_colormap*);
void lb_colormap_get_color(const lb_colormap*, int, double c[3]);
void lb_colormap_set_color(lb_colormap*, int, const double c[3]);
void lb_colormap_append_color(lb_colormap*, const double c[3]);
void lb_colormap_map_value(const lb_colormap*, double, double c[3]);

LB_END_DECLS

#endif /* LB_COLORMAP_H */
