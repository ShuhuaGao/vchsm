#pragma once
#include "StructDefinitions.h"
#include <vector>
#include <utility>

/**
 * Implementation of automaticth2wfw in automaticth2wfw.m.
*/
FwxyStructArray automaticth2wfw(const ThxyStructArray& th, int fs2 = 8000, int fmax = 5000);