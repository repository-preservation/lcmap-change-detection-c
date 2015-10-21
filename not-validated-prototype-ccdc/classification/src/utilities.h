
#ifndef UTILITIES_H
#define UTILITIES_H


#include <stdio.h>


#define LOG_MESSAGE(message, module) \
            write_message((message), (module), "INFO", \
                          __FILE__, __LINE__, stdout);


#define WARNING_MESSAGE(message, module) \
            write_message((message), (module), "WARNING", \
                          __FILE__, __LINE__, stdout);


#define ERROR_MESSAGE(message, module) \
            write_message((message), (module), "ERROR", \
                          __FILE__, __LINE__, stdout);


#define RETURN_ERROR(message, module, status) \
           {write_message((message), (module), "ERROR", \
                          __FILE__, __LINE__, stdout); \
            return (status);}


void write_message
(
    const char *message, /* I: message to write to the log */
    const char *module,  /* I: module the message is from */
    const char *type,    /* I: type of the error */
    char *file,          /* I: file the message was generated in */
    int line,            /* I: line number in the file where the message was
                               generated */
    FILE * fd            /* I: where to write the log message */
);


#endif /* UTILITIES_H */
