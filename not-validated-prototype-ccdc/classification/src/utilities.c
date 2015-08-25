
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h>


#include "utilities.h"


/*****************************************************************************
  NAME:  write_message

  PURPOSE:  Writes a formatted log message to the specified file handle.

  RETURN VALUE:  None

  NOTES:
      - Log Message Format:
            yyyy-mm-dd HH:mm:ss pid:module [filename]:line message
*****************************************************************************/
void write_message
(
    const char *message, /* I: message to write to the log */
    const char *module,  /* I: module the message is from */
    const char *type,    /* I: type of the error */
    char *file,          /* I: file the message was generated in */
    int line,            /* I: line number in the file where the message was
                               generated */
    FILE *fd             /* I: where to write the log message */
)
{
    time_t current_time;
    struct tm *time_info;
    int year;
    pid_t pid;

    time (&current_time);
    time_info = localtime (&current_time);
    year = time_info->tm_year + 1900;

    pid = getpid ();

    fprintf (fd, "%04d:%02d:%02d %02d:%02d:%02d %d:%s [%s]:%d [%s]:%s\n",
             year,
             time_info->tm_mon,
             time_info->tm_mday,
             time_info->tm_hour,
             time_info->tm_min,
             time_info->tm_sec,
             pid, module, basename (file), line, type, message);
}
