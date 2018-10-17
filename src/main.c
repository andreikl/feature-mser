#include <string.h>
#include <stdio.h>

#include "khash.h"

#include "utils.h"

const char* WIDTH = "-w";
const int WIDTH_DEF = 640;
const char* HEIGHT = "-h";
const int HEIGHT_DEF = 480;
const char* DEBUG = "-d";
const int DEBUG_DEF = 0;
const char* DEBUG_X = "-dpx";
const int DEBUG_X_DEF = -1;
const char* DEBUG_Y = "-dpy";
const int DEBUG_Y_DEF = -1;
const char* ALGORITHM = "-e";
const char* ALGORITHM_MSER = "mser";
const char* HELP = "-help";

const char* INPUT = "-i";
const char* INPUT_DEF = "-";
const char* OUTPUT = "-o";
const char* OUTPUT_DEF = "-";

KHASH_MAP_INIT_STR(str, char*)
khash_t(str) *h;

void print_help() {
    printf("edge-detector [options]\n");
    printf("options:\n");
    printf("\n");
    printf("%s: help\n", HELP);
    printf("%s: width, default: %d\n", WIDTH, WIDTH_DEF);
    printf("%s: height, default: %d\n", HEIGHT, HEIGHT_DEF);
    printf("%s: show debug info, default: %d\n", DEBUG, DEBUG_DEF);
    printf("%s: show debug info for x point, default: %d\n", DEBUG_X, DEBUG_X_DEF);
    printf("%s: show debug info for y point, default: %d\n", DEBUG_Y, DEBUG_Y_DEF);
    printf("%s: input, default:  %s, -: input stream\n", INPUT, INPUT_DEF);
    printf("%s: output, default:  %s, -: output stream\n", OUTPUT, OUTPUT_DEF);
    printf("%s: algorithm's name, default: %s\n", ALGORITHM, ALGORITHM_MSER);
    printf("algorithms: %s\n", ALGORITHM_MSER);
    exit(0);
}

int main(int argc, char** argv) {
    h = kh_init(str);
    int ret;
    unsigned k;

    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            k = kh_put(str, h, argv[i], &ret);
            kh_val(h, k) = (i + 1 < argc) ? argv[i + 1] : NULL;
        }
    }

    k = kh_get(str, h, HELP);
    if (k != kh_end(h)) {
        print_help();
    }
    else {
        const char* input_file = read_input_value(INPUT, INPUT_DEF);

        int input_size;
        char* input_data = read_file(input_file, &input_size);

        printf("%s, size: %d\n", input_file, input_size);

        free(input_data);
    }

    kh_destroy(str, h);
    return 0;
}
