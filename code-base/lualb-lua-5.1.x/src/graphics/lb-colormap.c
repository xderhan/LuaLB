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

lb_colormap*
lb_colormap_new(void)
{
	lb_colormap* self = lb_new(lb_colormap, 1);

	self->ncolors = 0;
	self->allocated = 32;
	self->colors = lb_new(double, 3*self->allocated);

	return self;
}

void
lb_colormap_destroy(lb_colormap* self)
{
	lb_assert(self != NULL);

	lb_free(self->colors);
	lb_free(self);
}

int
lb_colormap_num_colors(const lb_colormap* self)
{
	lb_assert(self != NULL);

	return self->ncolors;
}

void
lb_colormap_get_color(const lb_colormap* self, int c, double color[3])
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
lb_colormap_set_color(lb_colormap* self, int c, const double color[3])
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
lb_colormap_append_color(lb_colormap* self, const double color[3])
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
lb_colormap_map_value(const lb_colormap* self, double v, double color[3])
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
