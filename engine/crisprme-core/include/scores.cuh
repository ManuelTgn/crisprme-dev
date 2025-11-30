#pragma once

#include <common.cuh>

/// Calculate edit distance scores
void scores(const u8* query, const u8* strings, u8* result, int qlen, int slen, int n);
