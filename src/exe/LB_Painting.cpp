
#include "Painter.h"

int main() {
    Painter<2, double> p(dBox<2, double>({-1, -1}, {1, 1}));
    p.save("test");
}
