#include <iostream>
#include "Image.h"

int main() {
    Image image = Image();
    image.populate_block_lists();
    image.forward_dct();
    image.encode_zigzag();
    image.decode_zigzag();
    image.inverse_dct();
    image.decode();
}