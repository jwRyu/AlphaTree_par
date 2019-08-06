#pragma once

extern size_t memuse, max_memuse;

void* Malloc(size_t size);
void* Realloc(void* ptr, size_t size);
void Free(void* ptr);
