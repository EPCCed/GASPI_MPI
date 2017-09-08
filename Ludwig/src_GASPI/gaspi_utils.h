#include <GASPI.h>
#include <stdlib.h>
#include <stdio.h>
void gaspicheck ( const char* file, const int line, const int ec);
#define ASSERT(x) if (!(x))						\
    {									 \
      fprintf (stderr, "Error: '%s' [%s:%i]\n", #x, __FILE__, __LINE__); \
      exit (EXIT_FAILURE);						 \
    }
#define GASPIERROR(ec) gaspicheck (__FILE__, __LINE__, ec);
void wait_or_die ( gaspi_segment_id_t, gaspi_notification_id_t,
		   gaspi_notification_t expected );
