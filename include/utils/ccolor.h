// ccolor.hpp

#include <stdio.h>

#define CC_CONSOLE_COLOR_DEFAULT "\033[0m"
#define CC_FORECOLOR(C) "\033[" #C "m"
#define CC_BACKCOLOR(C) "\033[" #C "m"
#define CC_ATTR(A) "\033[" #A "m"

namespace zkr {
    class cc {
    public:
        class fore {
        public:
            static const char *black;
            static const char *blue;
            static const char *red;
            static const char *magenta;
            static const char *green;
            static const char *cyan;
            static const char *yellow;
            static const char *white;
            static const char *console;

            static const char *lightblack;
            static const char *lightblue;
            static const char *lightred;
            static const char *lightmagenta;
            static const char *lightgreen;
            static const char *lightcyan;
            static const char *lightyellow;
            static const char *lightwhite;
        };

        class back {
        public:
            static const char *black;
            static const char *blue;
            static const char *red;
            static const char *magenta;
            static const char *green;
            static const char *cyan;
            static const char *yellow;
            static const char *white;
            static const char *console;

            static const char *lightblack;
            static const char *lightblue;
            static const char *lightred;
            static const char *lightmagenta;
            static const char *lightgreen;
            static const char *lightcyan;
            static const char *lightyellow;
            static const char *lightwhite;
        };

        static char *color(char attr, char fg, char bg);

        static const char *console;
        static const char *underline;
        static const char *bold;
    };
}
