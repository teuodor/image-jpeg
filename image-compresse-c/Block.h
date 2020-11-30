//
// Created by Teuodor on 22/10/2020.
//

#ifndef IMAGE_COMPRESSE_C_BLOCK_H
#define IMAGE_COMPRESSE_C_BLOCK_H
#include<vector>
#include<string>

using std::vector;
using std::string;

class Block {
public:
    vector<vector<float>> block_store;
    string type;
    int *position;

    Block(vector<vector<float>> blockStore, string type, int *position);


};




#endif //IMAGE_COMPRESSE_C_BLOCK_H
