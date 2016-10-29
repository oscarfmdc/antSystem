extern void assert_msg(const char *expr, const int line,
                             const char *filename);
extern void info(const char *fmt, ...);

#ifdef NDEBUG
#  define assert(_cond_) (void)0
#else
#  define assert(_cond_) \
            if (_cond_) { /* do nothing */ } \
            else { assert_msg(#_cond_, __LINE__, __FILE__); _info_; abort(); }

#endif /* NO_DBG */
       /* --------------------------------------------------------- */

       /* ---------------------- my_assert.c ---------------------- */
       /* file to define helper functions used by customized assert.
       */
       extern void info(const char *fmt, ...)
       {
            va_list ap;

            va_start(ap, fmt);
            (void) vfprintf(stderr, fmt, ap);
            va_end(ap);
            (void)fprintf(stderr, "\n");

       }

       extern void assert_msg(const char *expr, const int line,
                             const char *filename)
       {
            /* Cleanup stuff, e.g closing file handles, releasing
             * memory etc., should be put here.
             */
            fflush(NULL);

            (void)fprintf(stderr, "\nASSERTION (%s) FAILED: %s %d\n",
                          expr, filename, line);
       }

       /* --------------------------------------------------------- */
