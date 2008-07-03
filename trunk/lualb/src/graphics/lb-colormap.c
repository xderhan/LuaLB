#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <math.h>
#include <string.h>

#include "lb-memory.h"
#include "lb-messages.h"
#include "lb-colormap.h"

struct _LBColormap {
	int     ncolors;
	int	allocated;
	double* colors;
};

LBColormap*
LBColormap_new(void)
{
	LBColormap* self = lb_new(LBColormap, 1);

	self->ncolors = 0;
	self->allocated = 32;
	self->colors = lb_new(double, 3*self->allocated);

	return self;
}

void
LBColormap_destroy(LBColormap* self)
{
	lb_assert(self != NULL);

	lb_free(self->colors);
	lb_free(self);
}

int
LBColormap_num_colors(const LBColormap* self)
{
	lb_assert(self != NULL);

	return self->ncolors;
}

void
_LBColormap_get_color(const LBColormap* self, int c, double color[3])
{
	lb_assert(self != NULL);
	lb_assert(c < self->ncolors);
	lb_assert(color != NULL);

	c *= 3;

	color[0] = self->colors[c++];
	color[1] = self->colors[c++];
	color[2] = self->colors[c++];
}

void
_LBColormap_set_color(LBColormap* self, int c, const double color[3])
{
	lb_assert(self != NULL);
	lb_assert(c < self->ncolors);
	lb_assert(color != NULL);

	c *= 3;

	self->colors[c++] = color[0];
	self->colors[c++] = color[1];
	self->colors[c++] = color[2];
}

void
LBColormap_append_color(LBColormap* self, const double color[3])
{
	lb_assert(self != NULL);
	lb_assert(color != NULL);

	if (self->ncolors == self->allocated) {
		self->allocated *= 2;
		self->colors = lb_renew(double, self->colors,
					  3*self->allocated);
	}

	self->colors[3*self->ncolors]     = color[0];
	self->colors[3*self->ncolors + 1] = color[1];
	self->colors[3*self->ncolors + 2] = color[2];

	self->ncolors += 1;
}

void
_LBColormap_map_value(const LBColormap* self, double v, double color[3])
{
	int bin;

	lb_assert(self != NULL);
	lb_assert(color != NULL);
	lb_assert(self->ncolors > 1);

	v = LB_CLAMP(v, 0.0, 1.0);

        bin = floor(v*(self->ncolors - 1));
        bin = LB_CLAMP(bin, 0, self->ncolors - 1);

        {
                const double* c0 = self->colors + 3*bin;
                const double* c1 = c0 + 3;
                const double t = v*(self->ncolors - 1) - bin;
                const double t1 = 1.0 - t;

                color[0] = c0[0]*t1 + c1[0]*t;
                color[1] = c0[1]*t1 + c1[1]*t;
                color[2] = c0[2]*t1 + c1[2]*t;
        }
}
