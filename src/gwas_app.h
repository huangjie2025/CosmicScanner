#pragma once
#include "config.h"

namespace cosmic {
class GwasApp {
public:
    GwasApp(const GwasConfig& config);
    int run();
private:
    GwasConfig cfg;
};
}
