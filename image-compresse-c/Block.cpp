//
// Created by Teuodor on 22/10/2020.
//

#include "Block.h"

#include <utility>

Block::Block(vector<vector<float>> blockStore, string type, int *position) : block_store(std::move(blockStore)),
                                                                                         type(std::move(type)),
                                                                                         position(position) {}
