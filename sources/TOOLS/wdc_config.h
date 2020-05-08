#define LOGGING 1


#if LOGGING == 1
        #define LOG(x) std::cout << x << "\n"
#else
        #define LOG(x)
#endif

