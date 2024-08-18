#include "data.h"

#if LOGQ == 24
#include "data24.c"
#elif LOGQ == 32
#include "data32.c"
#elif LOGQ == 40
#include "data40.c"
#elif LOGQ == 48
#include "data48.c"
#elif LOGQ == 56
#include "data56.c"
#elif LOGQ == 64
#include "data64.c"
#endif
