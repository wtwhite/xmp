#include "error.h"

int p_error(char *file, int line, char *err_msg)
{
  perror("File:%s,line$d: %s\n",file,line,err_msg);
   exit(-1);
}

