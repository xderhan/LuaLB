#ifndef LB_COLORMAP_H
#define LB_COLORMAP_H

#include "lb-macros.h"

LB_BEGIN_DECLS

typedef struct _LBColormap LBColormap;

LBColormap* LBColormap_new(void);
void LBColormap_destroy(LBColormap*);

int LBColormap_num_colors(const LBColormap*);

void _LBColormap_get_color(const LBColormap*, int, double c[3]);
void _LBColormap_set_color(LBColormap*, int, const double c[3]);
void _LBColormap_append_color(LBColormap*, const double c[3]);

void _LBColormap_map_value(const LBColormap*, double, double c[3]);

LB_END_DECLS

#endif /* LB_COLORMAP_H */
