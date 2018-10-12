#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "khash.h"

#include "constants.h"

KHASH_SET_INIT_STR(str)
extern khash_t(str) *h;

const char* read_input_value(const char name[], const char def_value[]) {
    unsigned k = kh_get(str, h, name);
    if (k != kh_end(h)) {
        return kh_val(h, k);
    }
    return def_value;
}

char* read_file(const char path[], int* size) {
    char buffer[BUFFER_SIZE];
    FILE* fstream;
    size_t read;

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
    } while (read != BUFFER_SIZE);

    if (path != '-') {
         fclose(fstream);
    }
    return data;
}
