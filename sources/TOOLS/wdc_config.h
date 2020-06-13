#define LOGGING 1

#include <iostream>


#if LOGGING == 1
        #define WDC_LOG(x) std::cout << x << "\n"
#else
        #define WDC_LOG(x)
#endif

