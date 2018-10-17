#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "khash.h"

#include "constants.h"

KHASH_MAP_INIT_STR(str, char*)
extern khash_t(str) *h;

const char* read_input_value(const char* name, const char def_value[]) {
    unsigned k = kh_get(str, h, name);
    if (k != kh_end(h)) {
        return kh_val(h, k);
    }
    return def_value;
}

char* read_file(const char* path, int* size) {
    char buffer[BUFFER_SIZE];
    FILE* fstream;
    size_t read;

    #ifdef DEBUG
    clock_t start_time = clock();
    #endif

    if (path[0] == '-') {
        fstream = stdin;
    }
    else {
        fstream = fopen(path, "r");
    }

    char* data = NULL; *size = 0;
    do {
        read = fread(buffer, sizeof(buffer[0]), BUFFER_SIZE, fstream);
        if (read > 0) {
            if (data == NULL) {
                data = malloc(read);
            }
            else {
                data = realloc(data, *size + read);
            }
            memcpy(data + *size, buffer, read);
            *size += read;
        }
    } while (read == BUFFER_SIZE);

    if (path[0] != '-') {
         fclose(fstream);
    }

    #ifdef DEBUG
    double diff = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("read_file: elapsed %f ms\n", diff);
    #endif
    
    return data;
}
