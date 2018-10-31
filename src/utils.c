#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "khash.h"

#include "constants.h"

KHASH_MAP_INIT_STR(str, char*)
extern khash_t(str) *h;

const char* read_str_value(const char* name, const char def_value[]) {
    unsigned k = kh_get(str, h, name);
    if (k != kh_end(h)) {
        return kh_val(h, k);
    }
    return def_value;
}

int read_int_value(const char name[], int def_value) {
    unsigned k = kh_get(str, h, name);
    if (k != kh_end(h)) {
        const char* value = kh_val(h, k);
        return atoi(value);
    }
    return def_value;
}

unsigned char* read_file(const char* path, int* size) {
    unsigned char buffer[BUFFER_SIZE];
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

    unsigned char* data = NULL; *size = 0;
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

void write_file(const char* path, unsigned char* data, int width, int height) {
    FILE* fstream;
    size_t written;
    unsigned char buffer[BUFFER_SIZE];
    int i = 0, j = 0, size = width * height;

#ifdef DEBUG
    clock_t start_time = clock();
#endif

    if (path[0] == '-') {
        fstream = stdout;
    }
    else {
        fstream = fopen(path, "w");
    }

    sprintf(buffer, "P6\n%d %d\n255\n\0", width, height);
    fwrite(buffer, 1, strlen(buffer), fstream);
    while (j < size) {
        if (i + 3 >= BUFFER_SIZE) {
            written += fwrite(buffer, 1, i, fstream);
            i = 0;
        }
        buffer[i] = data[j];
        buffer[i + 1] = data[j];
        buffer[i + 2] = data[j];
        i += 3; j++;
    }
    if (i > 0) {
        written += fwrite(buffer, 1, i, fstream);
    }

    if (path[0] != '-') {
        fclose(fstream);
    }

#ifdef DEBUG
    double diff = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("write_file: elapsed %f ms\n", diff);
#endif
}
