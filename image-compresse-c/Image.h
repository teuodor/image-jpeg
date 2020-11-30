//
// Created by Teuodor on 22/10/2020.
//

#ifndef IMAGE_COMPRESSE_C_IMAGE_H
#define IMAGE_COMPRESSE_C_IMAGE_H

#include "Block.h"
#include <utility>
#include <iostream>
#include <cstdint>

using std::pair;
using std::make_pair;
using std::uint8_t;

class Image {
public:
    string file_name = "nt-P3.ppm", mode, comment;
    int width, height, max;
    vector<vector<float>> y, u, v;
    vector<Block> y_blocks, u_blocks, v_blocks;
    vector<int8_t> zigzag_array;
    void populate_block_lists();
    void forward_dct();
    void inverse_dct();
    void encode_zigzag();
    void decode_zigzag();
    void decode();
};



#endif //IMAGE_COMPRESSE_C_IMAGE_H
