#ifndef LB_MESSAGES_H
#define LB_MESSAGES_H

#include <lb-macros.h>

LB_BEGIN_DECLS

#ifndef LB_DISABLE_ASSERTS
# ifdef __GNUC__
#  define lb_assert(expr) if (!(expr))				\
     lb_fatal("file %s: line %d (%s): assertion failed: (%s)\n",	\
                __FILE__, __LINE__, __PRETTY_FUNCTION__, #expr)
#  define lb_assert_not_reached()					\
     lb_fatal("file %s: line %d (%s): should not be reached\n",	\
                __FILE__, __LINE__, __PRETTY_FUNCTION__)
# else /* !__GNUC__ */
#  define lb_assert(expr) if (!(expr))				\
     lb_fatal("file %s: line %d: assertion failed: (%s)\n",		\
                __FILE__, __LINE__, #expr)
#  define lb_assert_not_reached()					\
     lb_fatal("file %s: line %d: should not be reached\n",		\
                __FILE__, __LINE__)
# endif /* __GNUC__ */
#else /* LB_DISABLE_ASSERTS */
# define lb_assert(expr) ((void) 0)
# define lb_assert_not_reached() ((void) 0)
#endif /* LB_ASSERT */

void lb_abort(void) LB_GNUC_NORETURN;
void lb_fatal(const char*, ...) LB_GNUC_NORETURN;
void lb_warning(const char*, ...);
void lb_info(const char*, ...);

void lb_info_enable(void);
void lb_info_disable(void);
int  lb_info_enabled(void);

int lb_is_parallel(void);

LB_END_DECLS

#endif /* LB_MESSAGES_H */
