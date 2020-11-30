#include "Image.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <lcms2_plugin.h>

using std::ifstream;
using std::ofstream;
using std::sqrt;
using std::cos;

#define PI 3.14159265359

float clamp(float x){
    if(x < 0)
        return 0;
    if(x > 255)
        return 255;
    return x;
}

void subtract128(vector<vector<float>> &v) {
    for(int i = 0; i < 8; i++)
        for(int j = 0;j < 8; j++)
            v[i][j] -= 128;
}

void add128(vector<vector<float>> &v) {
    for(int i = 0; i < 8; i++)
        for(int j = 0;j < 8; j++)
            v[i][j] += 128;
}

void quantization(vector<vector<float>> &g) {
    float q[8][8] = {{  6,  4,  4,  6, 10, 16, 20, 24 },
                     {  5,  5,  6,  8, 10, 23, 24, 22 },
                     {  6,  5,  6, 10, 16, 23, 28, 22 },
                     {  6,  7,  9, 12, 20, 35, 32, 25 },
                     {  7,  9, 15, 22, 27, 44, 41, 31 },
                     { 10, 14, 22, 26, 32, 42, 45, 37 },
                     { 20, 26, 31, 35, 41, 48, 48, 40 },
                     { 29, 37, 38, 39, 45, 40, 31, 40 }};

    for(int i = 0; i < 8; i++)
        for(int j = 0; j < 8; j++) {
            g[i][j] /= q[i][j];
            g[i][j] = (float) ((int) g[i][j]);
        }
}

void dequantization(vector<vector<float>> &g) {
    float q[8][8] = {{  6,  4,  4,  6, 10, 16, 20, 24 },
                     {  5,  5,  6,  8, 10, 23, 24, 22 },
                     {  6,  5,  6, 10, 16, 23, 28, 22 },
                     {  6,  7,  9, 12, 20, 35, 32, 25 },
                     {  7,  9, 15, 22, 27, 44, 41, 31 },
                     { 10, 14, 22, 26, 32, 42, 45, 37 },
                     { 20, 26, 31, 35, 41, 48, 48, 40 },
                     { 29, 37, 38, 39, 45, 40, 31, 40 }};

    for(int i = 0; i < 8; i++)
        for(int j = 0; j < 8; j++) {
            g[i][j] *= q[i][j];
        }
}

float forward_dct(vector<vector<float>> const& g, int u, int v) {
    double guv = 0.25;
    if (u == 0)
        guv *= 1 / sqrt(2);
    if (v == 0)
        guv *= 1 / sqrt(2);

    double sum = 0;

    for (int x = 0; x < 8; x++)
        for (int y = 0; y < 8; y++) {
            sum += g[x][y] * cos(((2 * x + 1) * u * PI) / 16) * cos(((2 * y + 1) * v * PI) / 16);
        }
    auto to_return = (float) (guv * sum);
    return to_return;
}

float inverse_dct(vector<vector<float>> const& f, int x, int y) {
    double fxy = 0.25;
    double sum = 0;
    for(int u = 0; u < 8; u++)
        for(int v = 0; v < 8; v++) {
            double to_sum = 1;
            if(u == 0)
                to_sum *= 1 / sqrt(2);
            if(v == 0)
                to_sum *= 1 / sqrt(2);

            to_sum *= f[u][v] * cos(((2 * x + 1) * u * PI) /16) * cos(((2 * y + 1) * v * PI) /16);
            sum += to_sum;
        }
    auto to_return = (float) (fxy * sum);
    return to_return;
}

int determine_size(float n){
    if(n < 0)
        n *= -1;
    int nr = (int) n, to_return = 0;
    while(nr){
        to_return++;
        nr /= 2;
    }
    return to_return;
}

void next_position(int &i, int &j, bool &reverse, bool &margin){
    if(i == 7 && reverse && margin) {
        j++;
        reverse = !reverse;
        margin = false;
        return;
    }
    if(j == 7 && !reverse && margin) {
        i++;
        reverse = !reverse;
        margin = false;
        return;
    }

    if(i == 0 && !reverse && margin) {
        j++;
        reverse = !reverse;
        margin = false;
        return;
    }
    if(j == 0 && reverse && margin) {
        i++;
        reverse = !reverse;
        margin = false;
        return;
    }

    margin = true;
    if(reverse) {
        i++;
        j--;
    }
    else {
        i--;
        j++;
    }

}

vector<pair<int, float>> encode_zigzag(vector<vector<float>> matrix) {
    vector<pair<int, float>> to_return;
    int i = 0, j = 0;
    bool reverse = false, margin = true;
    auto pair = make_pair(determine_size(matrix[i][j]), matrix[i][j]);
    to_return.push_back(pair);
    while(i != matrix.size() - 1 || j != matrix[i].size() - 1) {
        next_position(i, j, reverse, margin);
        pair = make_pair(determine_size(matrix[i][j]), matrix[i][j]);
        to_return.push_back(pair);
    }
    return to_return;
}

vector<vector<float>> decode_zigzag(vector<pair<int, float>> array) {
    vector<vector<float>> to_return(8, vector<float>(8, 0));

    int i = 0, j = 0;
    bool reverse = false, margin = true;
    int pos = 0;
    to_return[i][j] = array[pos].second;
    while(i != 7 || j != 7) {
        pos++;
        next_position(i, j,reverse, margin);
        to_return[i][j] = array[pos].second;
    }
    return to_return;
}
int8_t trim128(float x) {
    if (x < -128)
        return -128;
    if (x > 127)
        return 127;
    return x;
}

void final_zigzag(vector<int8_t> &array, vector<pair<int, float>> initial) {
    array.push_back(initial[0].first);
    array.push_back(trim128(initial[0].second));
    int y0 = 0;
    for(int i = 1; i < initial.size(); i++) {
        if(initial[i].second == 0) {
            y0++;
        }
        else {
            array.push_back(y0);
            array.push_back(initial[i].first);
            array.push_back(trim128(initial[i].second));
            y0 = 0;
        }
    }
    if(y0 > 0) {
        array.push_back(0);
        array.push_back(0);
    }
}
vector<pair<int, float>> initial_zigzag(vector<int8_t> &array, int &start){
    vector<pair<int, float>> initial;
    initial.emplace_back(array[start], (int8_t) array[start + 1]);
    start += 2;
    while(start != 63 * 3){
        if(array[start] == 0 && array[start + 1] == 0) {
            start += 2;
            break;
        }

        for(int j = 0; j < array[start]; j++)
            initial.emplace_back(0, 0);
        initial.emplace_back((int) array[start + 1],(int8_t) array[start + 2]);
        start += 3;
    }
    int lastZero = (int)(64 - initial.size());
    while (lastZero != 0){
        initial.emplace_back(0, 0);
        lastZero--;
    }
    return initial;
}
void Image::populate_block_lists() {
    ifstream input = ifstream(file_name);
    int nr = 0;
    string str;
    int r = 0, g = 0, b = 0;
    int x = 0;
    while(std::getline(input, str)){
        if(nr == 0)
            mode = str;
        else if(nr == 1)
            comment = str;
        else if(nr == 2)
        {
            auto pos = str.find(' ');
            width = std::stoi(str.substr(0, pos));
            str.erase(0, pos + 1);
            height = std::stoi(str);
            int j = height;
            while(j --){
                vector<float> v1, v2, v3;
                y.push_back(v1);
                u.push_back(v2);
                v.push_back(v3);
            }
        } else if(nr == 3)
            max = std::stoi(str);
        else{
            int pixel = std::stoi(str);
            if(nr % 3 == 1)
                r = pixel;
            else if(nr % 3 == 2)
                g = pixel;
            else{
                b = pixel;
                if(y[x].size() == width){
                    x++;
                }

                y[x].push_back(0.299 * r + 0.587 * g + 0.114 * b);
                u[x].push_back(128 - 0.168736 * r - 0.331264 * g + 0.5 * b);
                v[x].push_back(128 + 0.5 * r - 0.418688 * g - 0.081312 * b);
            }
        }
        nr++;
    }


    for(int i = 0; i < height; i += 8)
        for(int j = 0; j < width; j += 8)
        {
            vector<vector<float>> block_store;
            vector<vector<float>> block_store_u;
            vector<vector<float>> block_store_v;
            int *position = new int[2];
            position[0] = i / 8;
            position[1] = j / 8;

            int i1 = 8;
            while(i1--){
                vector<float> vy;
                block_store.push_back(vy);
                if(i1 % 2 == 1){
                    vector<float> vb, vd;
                    block_store_u.push_back(vb);
                    block_store_v.push_back(vd);
                }
            }

            for(i1 = 0; i1 < 8; i1++)
                for(int j1 = 0; j1 < 8; j1 ++){
                    block_store[i1].push_back(y[i + i1][j + j1]);
                }

            for(i1 = 0; i1 < 8; i1 += 2){
                for(int j1 = 0; j1 < 8; j1 += 2)
                {
                    float average_b = u[i + i1][j + j1] + u[i + i1 + 1][j + j1] + u[i + i1][j + j1 + 1] + u[i + i1 + 1][j + j1 + 1];
                    average_b /= 4;
                    float average_r = v[i + i1][j + j1] + v[i + i1 + 1][j + j1] + v[i + i1][j + j1 + 1] + v[i + i1 + 1][j + j1 + 1];
                    average_r /= 4;
                    block_store_u[i1 / 2].push_back(average_b);
                    block_store_v[i1 / 2].push_back(average_r);
                }
            }

            Block block = Block(block_store, "y", position);
            Block block_u = Block(block_store_u, "u", position);
            Block block_v = Block(block_store_v, "v", position);
            y_blocks.push_back(block);
            u_blocks.push_back(block_u);
            v_blocks.push_back(block_v);
        }
    std::cout << "encode complete!!!\n";
}

void Image::decode() {
    for(int i = 0 ; i < y.size(); i++) {
        y[i].clear();
        u[i].clear();
        v[i].clear();
    }

    for(int k = 0; k <  y_blocks.size(); k++) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                y[i + y_blocks[k].position[0] * 8].push_back(y_blocks[k].block_store[i][j]);
                u[i + u_blocks[k].position[0] * 8].push_back(u_blocks[k].block_store[i][j]);
                v[i + v_blocks[k].position[0] * 8].push_back(v_blocks[k].block_store[i][j]);
            }
        }
    }

    vector<float> new_r;
    vector<float> new_g;
    vector<float> new_b;

    for(int i = 0; i < y.size(); i++)
        for(int j = 0; j < y[i].size(); j++) {
            new_r.push_back(y[i][j] + 1.402 * (v[i][j] - 128));
            new_g.push_back(y[i][j] - 0.344136 * (u[i][j] - 128) - 0.714136 * (v[i][j] - 128));
            new_b.push_back(y[i][j] + 1.772 * (u[i][j] - 128));
        }
    std::cout << "decode complete!!!\n";

    string new_file = "new_" + file_name;
    ofstream out(new_file);
    out << mode + "\n";
    out << comment + "\n";

    std::stringstream ss;
    ss << width << " " << height << "\n";
    out << ss.str();
    ss.str(string());

    ss << max << "\n";
    out << ss.str();
    ss.str(string());

    for(int i = 0; i < new_r.size(); i++){
        new_r[i] = clamp(new_r[i]);
        new_g[i] = clamp(new_g[i]);
        new_b[i] = clamp(new_b[i]);

        ss << (int) new_r[i] << "\n" << (int) new_g[i] << "\n" << (int) new_b[i] << "\n";
        out << ss.str();
        ss.str(string());
    }
    out.close();
}

void Image::forward_dct() {
    for(int k = 0; k < y_blocks.size(); k++) {
        auto y_block_store = y_blocks[k].block_store;
        vector<vector<float>> u_block_store;
        vector<vector<float>> v_block_store;

        for(int i = 0; i < 8; i++) {
            vector<float> v1, v2;
            u_block_store.push_back(v1);
            v_block_store.push_back(v2);
        }

        for(int i = 0; i < 4; i ++)
            for(int j = 0; j < 4; j++) {
                u_block_store[i * 2].push_back(u_blocks[k].block_store[i][j]);
                u_block_store[i * 2].push_back(u_blocks[k].block_store[i][j]);
                u_block_store[i * 2 + 1].push_back(u_blocks[k].block_store[i][j]);
                u_block_store[i * 2 + 1].push_back(u_blocks[k].block_store[i][j]);

                v_block_store[i * 2].push_back(v_blocks[k].block_store[i][j]);
                v_block_store[i * 2].push_back(v_blocks[k].block_store[i][j]);
                v_block_store[i * 2 + 1].push_back(v_blocks[k].block_store[i][j]);
                v_block_store[i * 2 + 1].push_back(v_blocks[k].block_store[i][j]);
            }

        subtract128(y_block_store);
        subtract128(u_block_store);
        subtract128(v_block_store);

        auto copy_y_block_store = y_block_store;
        auto copy_u_block_store = u_block_store;
        auto copy_v_block_store = v_block_store;

        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++) {
                y_block_store[i][j] = ::forward_dct(copy_y_block_store, i, j);
                u_block_store[i][j] = ::forward_dct(copy_u_block_store, i, j);
                v_block_store[i][j] = ::forward_dct(copy_v_block_store, i, j);
            }

        quantization(y_block_store);
        quantization(u_block_store);
        quantization(v_block_store);


        y_blocks[k].block_store = y_block_store;
        u_blocks[k].block_store = u_block_store;
        v_blocks[k].block_store = v_block_store;
    }
    std::cout << "forward dct complete!!!\n";
}


void Image::inverse_dct() {
    for(int k = 0; k < y_blocks.size(); k++) {
        auto y_block_store = y_blocks[k].block_store;
        auto u_block_store = u_blocks[k].block_store;
        auto v_block_store = v_blocks[k].block_store;

        dequantization(y_block_store);
        dequantization(u_block_store);
        dequantization(v_block_store);

        auto copy_y_block_store = y_block_store;
        auto copy_u_block_store = u_block_store;
        auto copy_v_block_store = v_block_store;

        for(int i = 0 ; i < 8; i++)
            for(int j = 0; j < 8; j++) {
                y_block_store[i][j] = ::inverse_dct(copy_y_block_store, i, j);
                u_block_store[i][j] = ::inverse_dct(copy_u_block_store, i, j);
                v_block_store[i][j] = ::inverse_dct(copy_v_block_store, i, j);
            }

        add128(y_block_store);
        add128(u_block_store);
        add128(v_block_store);

        y_blocks[k].block_store = y_block_store;
        u_blocks[k].block_store = u_block_store;
        v_blocks[k].block_store = v_block_store;
    }
    std::cout << "inverse dct complete!!!\n";
}

void Image::encode_zigzag() {
    std::cout << "b4 encode";
    for(int i = 0; i < y_blocks[0].block_store.size(); i++) {
        for (int j = 0; j < y_blocks[0].block_store[i].size(); j++)
            std::cout << y_blocks[0].block_store[i][j] << " ";
        std::cout << "\n";
    }
    std ::cout << "\n";
    for (auto i = 0; i < y_blocks.size(); i++) {
        auto initial_y = ::encode_zigzag(y_blocks[i].block_store);
        auto initial_u = ::encode_zigzag(u_blocks[i].block_store);
        auto initial_v = ::encode_zigzag(v_blocks[i].block_store);

        final_zigzag(zigzag_array, initial_y);
        final_zigzag(zigzag_array, initial_u);
        final_zigzag(zigzag_array, initial_v);
    }
    std::cout << "encoded zigzag!! \n";
}

void Image::decode_zigzag() {
    int start = 0;
    int pos = 0;
    while(start != zigzag_array.size()) {
        vector<pair<int, float>> initial_y = initial_zigzag(zigzag_array, start);
        auto initial_u = initial_zigzag(zigzag_array, start);
        auto initial_v = initial_zigzag(zigzag_array, start);

        y_blocks[pos].block_store = ::decode_zigzag(initial_y);
        u_blocks[pos].block_store = ::decode_zigzag(initial_u);
        v_blocks[pos].block_store = ::decode_zigzag(initial_v);
        if(pos == 0) {
            for(int i = 0; i < y_blocks[pos].block_store.size(); i++) {
                for(int j = 0; j < y_blocks[pos].block_store[i].size(); j++)
                    std::cout << y_blocks[pos].block_store[i][j] << " ";
                std::cout << "\n";
            }
        }
        pos++;
    }

    std::cout << "decoded zigzag!!! \n";
}
