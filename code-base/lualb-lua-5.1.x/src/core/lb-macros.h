#ifndef LB_MACROS_H
#define LB_MACROS_H

#include <stddef.h>

#ifdef __cplusplus
#  define LB_BEGIN_DECLS extern "C" {
#  define LB_END_DECLS } /* extern "C" */
#else /* not C++ */
#  define LB_BEGIN_DECLS
#  define LB_END_DECLS
#endif /* __cplusplus */

#define LB_STRINGIFY(macro_or_string)	LB_STRINGIFY_ARG(macro_or_string)
#define	LB_STRINGIFY_ARG(contents)	#contents

#define LB_STRLOC __FILE__ ":" LB_STRINGIFY(__LINE__)

#ifndef FALSE
# define FALSE (0)
#endif

#ifndef TRUE
# define TRUE (!FALSE)
#endif

#ifdef __GNUC__
# define LB_GNUC_NORETURN __attribute__((__noreturn__))
# define LB_GNUC_UNUSED __attribute__((__unused__))
# define LB_RESTRICT __restrict__
#else
# define LB_GNUC_NORETURN
# define LB_GNUC_UNUSED
# define LB_RESTRICT
#endif /* __GNUC__ */

#define LB_MAX(a, b) ((a) > (b) ? (a) : (b))
#define LB_MIN(a, b) ((a) < (b) ? (a) : (b))

#define LB_CLAMP(a, n, x) ((a) < (n) ? (n) : ((a) > (x) ? (x) : (a)))

#define LB_N_ELEMENTS(a) (sizeof(a)/sizeof(a[0]))

#endif /* LB_MACROS_H */
