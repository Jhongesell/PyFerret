#include "EF_Util.h"

void FORTRAN(ef_get_mres)(int *mres)
{
  *mres = *GLOBAL_mres_ptr;
}
