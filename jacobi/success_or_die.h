#ifndef SUCCESS_OR_DIE_H
#define SUCCESS_OR_DIE_H
#include <GASPI_Ext.h>

#define SUCCESS_OR_DIE(f, args...)              \
  do                                            \
  {                                             \
    gaspi_return_t const r = f (args);          \
                                                \
    if (r != GASPI_SUCCESS)                     \
    {                                           \
      gaspi_printf ( "Error[%s:%i]: %s\n"       \
                   , __FILE__                   \
                   , __LINE__                   \
                   , gaspi_error_str (r)        \
                   );                           \
                                                \
      exit (-1);                      \
    }                                           \
  } while (0)

#endif
