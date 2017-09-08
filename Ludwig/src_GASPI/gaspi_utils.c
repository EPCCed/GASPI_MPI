#include "gaspi_utils.h"


void gaspicheck ( const char* file, const int line, const int ec)
{
  if (ec != GASPI_SUCCESS)
    {
      gaspi_printf ("Assertion failed in %s[%i]:%d\n", file, line, ec);
      exit (EXIT_FAILURE);
    }
}

void wait_or_die(gaspi_segment_id_t segment_id, gaspi_notification_id_t notification_id,
		 gaspi_notification_t expected ) {
  gaspi_notification_id_t id;
  GASPIERROR( gaspi_notify_waitsome (segment_id, notification_id, 1,
				     &id, GASPI_BLOCK) );
  
  ASSERT(id == notification_id);
  
  gaspi_notification_t value;
  GASPIERROR( gaspi_notify_reset (segment_id, id, &value) );

  ASSERT(value == expected);
}
