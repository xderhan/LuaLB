#ifndef LB_MEMORY_H
#define LB_MEMORY_H

#include <sys/types.h>
#include <lb-macros.h>

LB_BEGIN_DECLS

void* lb_malloc(size_t);
void* lb_realloc(void*, size_t);
void  lb_free(void*);

#define lb_new(type, n) ((type *) lb_malloc((n)*sizeof(type)))
#define lb_renew(type, ptr, n) ((type *) lb_realloc((ptr),(n)*sizeof(type)))

#ifdef __GNUC__
# undef alloca
# define alloca(size) __builtin_alloca(size)
#else
# include <alloca.h>
#endif /* __GNUC__ */

#define lb_alloca(size) alloca(size)
#define lb_newa(type, n) ((type *) lb_alloca((n)*sizeof(type)))

void lb_scpy(double*, const double*, int, int, int);

LB_END_DECLS

#endif /* LB_MEMORY_H */
