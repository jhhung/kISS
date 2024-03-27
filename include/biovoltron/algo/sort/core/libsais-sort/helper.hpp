#pragma once

#define prefetchr(address) __builtin_prefetch((const void *)(address), 0, 3)
#define prefetchw(address) __builtin_prefetch((const void *)(address), 1, 3)
