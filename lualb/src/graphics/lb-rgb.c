#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <png.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "lb-rgb.h"
#include "lb-memory.h"
#include "lb-messages.h"

struct _LBRGB {
	int width;
	int height;
	int rowstride;
	double* pixels;
};

LBRGB*
LBRGB_new(int width, int height)
{
	LBRGB* self;

	lb_assert(width > 0 && height > 0);

	self = lb_new(LBRGB, 1);

	self->width = width;
	self->height = height;
	self->rowstride = 3*width;
	self->pixels = lb_new(double, self->rowstride*height);

	return self;
}

void
LBRGB_destroy(LBRGB* self)
{
	lb_assert(self != NULL);

	lb_free(self->pixels);
	lb_free(self);
}

int
LBRGB_width(const LBRGB* self)
{
	lb_assert(self != NULL);
	return self->width;
}

int
LBRGB_height(const LBRGB* self)
{
	lb_assert(self != NULL);
	return self->height;
}

void
LBRGB_fill(LBRGB* self, const double c[3])
{
	double* end;
	double* ptr;
	double color[3];

	lb_assert(self != NULL);
	lb_assert(c != NULL);

	color[0] = LB_CLAMP(c[0], 0.0, 1.0);
	color[1] = LB_CLAMP(c[1], 0.0, 1.0);
	color[2] = LB_CLAMP(c[2], 0.0, 1.0);

	end = self->pixels + self->rowstride;

	for (ptr = self->pixels; ptr != end; ptr += 3) {
		ptr[0] = color[0];
		ptr[1] = color[1];
		ptr[2] = color[2];
	}

	end = self->pixels + self->height*self->rowstride;
	for (; ptr != end; ptr += self->rowstride)
		memcpy(ptr, self->pixels, self->rowstride*sizeof(double));
}

void
LBRGB_get_pixel(const LBRGB* self, int x, int y, double pixel[3])
{
	const double* ptr;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < self->width);
	lb_assert(y >= 0 && y < self->height);
	lb_assert(pixel != NULL);

	ptr = self->pixels + y*self->rowstride + 3*x;

	pixel[0] = ptr[0];
	pixel[1] = ptr[1];
	pixel[2] = ptr[2];
}

void
LBRGB_set_pixel(LBRGB* self, int x, int y, const double pixel[3])
{
	double* ptr;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < self->width);
	lb_assert(y >= 0 && y < self->height);
	lb_assert(pixel != NULL);

	ptr = self->pixels + y*self->rowstride + 3*x;

	ptr[0] = LB_CLAMP(pixel[0], 0.0, 1.0);
	ptr[1] = LB_CLAMP(pixel[1], 0.0, 1.0);
	ptr[2] = LB_CLAMP(pixel[2], 0.0, 1.0);
}

void
LBRGB_set_pixel_rgba(LBRGB* self, int x, int y, const double rgba[4])
{
	double* ptr;

	lb_assert(self != NULL);
	lb_assert(x >= 0 && x < self->width);
	lb_assert(y >= 0 && y < self->height);
	lb_assert(rgba != NULL);

	ptr = self->pixels + y*self->rowstride + 3*x;

	{
		const double a1 = LB_CLAMP(rgba[3], 0.0, 1.0);
		const double a2 = 1.0 - a1;

		ptr[0] = a2*ptr[0] + a1*LB_CLAMP(rgba[0], 0.0, 1.0);
		ptr[1] = a2*ptr[1] + a1*LB_CLAMP(rgba[1], 0.0, 1.0);
		ptr[2] = a2*ptr[2] + a1*LB_CLAMP(rgba[2], 0.0, 1.0);
	}
}

#ifndef png_jmpbuf
#define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif /* png_jmpbuf */

void
LBRGB_save(const LBRGB* self, char* name)
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_color_8 sig_bit;

	FILE *file;

	lb_assert(self != NULL);
	lb_assert(name != NULL);
	lb_assert( (file=(FILE *)fopen(name, "w")) != NULL );

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                          NULL, NULL, NULL);
	if (png_ptr == NULL)
		lb_fatal(LB_STRLOC);

	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		lb_fatal(LB_STRLOC);
	}

	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		lb_fatal(LB_STRLOC);
	}

	png_init_io(png_ptr, file);
        png_set_IHDR(png_ptr, info_ptr, self->width, self->height, 8,
                     PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	sig_bit.red = 8;
	sig_bit.green = 8;
	sig_bit.blue = 8;
	sig_bit.alpha = 8;

	png_set_sBIT(png_ptr, info_ptr, &sig_bit);
	png_write_info(png_ptr, info_ptr);
	png_set_shift(png_ptr, &sig_bit);
	png_set_packing(png_ptr);

	{
		const double* ptr = self->pixels;
		char* bytes = lb_newa(unsigned char, self->rowstride);
		int x, y;

		for (y = 0; y < self->height; ++y) {
    			png_bytep row_ptr;

			for (x = 0; x < self->rowstride; ++x)
				bytes[x] = (unsigned char) (255.0*ptr[x]);

    			row_ptr = (png_bytep) bytes;
			png_write_rows(png_ptr, &row_ptr, 1);
            		ptr += self->rowstride;
		}
	}

	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
}

void
LBRGB_map_value(LBRGB* self, const LBColormap* cmap,
		   int x, int y, double v)
{
        double* ptr;

        lb_assert(self != NULL);
        lb_assert(cmap != NULL);
        lb_assert(x >= 0 && x < self->width);
        lb_assert(y >= 0 && y < self->height);

	ptr = self->pixels + y*self->rowstride + 3*x;
	LBColormap_map_value(cmap, v, ptr);
}
