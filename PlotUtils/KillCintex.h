// this disables Cintex in ROOT 6 - you should use it as follows  HMS 8-8-2020
/*
 #if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,00)
 #include "PlotUtils/KillCintex.h"
 #else
 #include "Cintex/Cintex.h"
 #endif
 */
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,00)
namespace ROOT {
namespace Cintex {
namespace Cintex {
inline  void Enable() {}
}
}
}
#endif
