#pragma once

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "pcg_basic.h"

int* forRandOrder(int l, pcg32_random_t* rng);
bool allTrue(bool* a, int la);
