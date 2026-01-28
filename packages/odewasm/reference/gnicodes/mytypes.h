/* The GNU compiler predefines M_PIl,
 * whereas some other platforms do not.
 * The point is, it's a /long/ double constant.
 */
#ifndef MI_PIl
#  define M_PIl           3.141592653589793238462643383279502884L /* pi */
#endif
typedef long double FLP;
typedef int32_t INT;