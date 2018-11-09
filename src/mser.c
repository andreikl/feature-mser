#include <math.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <kvec.h>

#include "constants.h"
#include "utils.h"

#ifdef DEBUG
extern const char* debug_file;
#endif // DEBUG

#define MSER_PIXEL_MAXVALUE 256

// MSER statistic
typedef struct _MserStatistic MserStatistic;
struct _MserStatistic {
    int num_extremal;           // number of extremal regions
    int num_unstable;           // number of unstable extremal regions
    int num_abs_unstable;       // number of regions that failed the absolute stability test
    int num_too_big;            // number of regions that failed the maximum size test
    int num_too_small;          // number of regions that failed the minimum size test
    int num_duplicates;         // number of regions that failed the duplicate test
};

// MSER: extremal region
typedef struct _MserRegion MserRegion;
struct _MserRegion {
    int level;
    int pixel;
    int area;

    MserRegion *parent; // Pointer to the parent region
    MserRegion *child;  // Pointer to the first child
    MserRegion *next;   // Pointer to the next (sister) region

    double moments[5];
};

// MSER Data
struct _MserData {
    int width;
    int height;
    int pixels;      // number of pixels
    bool* visited;   // visited pixels

    // Regions
    kvec_t(int) boundary_pixels[MSER_PIXEL_MAXVALUE];
    // Extremal regions
    kvec_t(MserRegion*) region_stack;

    // Configuration
    int delta;            // delta filter parameter
    float max_area;       // badness test parameter
    float min_area;       // badness test parameter
    float max_variation;  // badness test parameter
    float min_diversity;  // minimum diversity
    // statistic
    MserStatistic statistic;
#ifdef DEBUG
    unsigned char* debug; // debug data if it is enabled
#endif // DEBUG
};
typedef struct _MserData MserData;

// null region
#ifdef COMPILER_MSC
#define MSER_VOID_NODE ( (1ui64 << 32) - 1)
#else
#define MSER_VOID_NODE ( (1ULL << 32) - 1)
#endif

// Create a new MSER data
MserData* mser_new(int width, int height) {
    int i;
    MserData *f = (MserData *) calloc(1, sizeof(MserData));

    f->width = width;
    f->height = height;
    f->pixels = width * height;
    f->visited = (bool*)calloc(1, f->pixels * sizeof(bool));
    for (i = 0; i < MSER_PIXEL_MAXVALUE; i++) {
        kv_init(f->boundary_pixels[i]);
    }
    kv_init(f->region_stack);

    // Configuration
    f->delta = 2; // delta > 0
    f->min_area = 0.0001f; // min_area >= 0.0
    f->max_area = 0.5f; // max_area >= min_area && max_area <= 1.0
    f->max_variation = 0.5f; // max_variation > 0.0
    f->min_diversity = 0.33f; // min_diversity >= 0.0 && min_diversity < 1.0

#ifdef DEBUG
    f->debug = (unsigned char*)malloc(f->pixels);
#endif // DEBUG

    return f;
}


// Delete MSER data
void mser_delete(MserData *f) {
    int i;
    if (f) {
        if (f->visited)
            free(f->visited);

        for (i = 0; i < MSER_PIXEL_MAXVALUE; i++) {
            kv_destroy(f->boundary_pixels[i]);
        }
        kv_destroy(f->region_stack);

#ifdef DEBUG
        if (f->debug)
            free(f->debug);
#endif // DEBUG
        free(f);
    }
}

void region_accumulate(MserRegion *region, int x, int y)
{
    region->area++;
    region->moments[0] += x;
    region->moments[1] += y;
    region->moments[2] += x * x;
    region->moments[3] += x * y;
    region->moments[4] += y * y;
}

void region_merge(MserRegion *x, MserRegion *y)
{
    // Add the moments together
    x->area += y->area;
    x->moments[0] += y->moments[0];
    x->moments[1] += y->moments[1];
    x->moments[2] += y->moments[2];
    x->moments[3] += y->moments[3];
    x->moments[4] += y->moments[4];

    y->next = x->child;
    x->child = y;
    y->parent = x;
}

void process_stack(int new_level, int pixel, MserData *mser) {
    MserRegion *top, *next_top, *new_region;

    do {
        top = kv_pop(mser->region_stack);
        next_top = kv_A(mser->region_stack, mser->region_stack.n - 1);

        // 2. If new_level is smaller than the grey-level on the second component on the
        // stack, set the top of stack grey-level to new_level and return from sub-routine
        // (This occurs when the new pixel is at a grey-level for which there is not yet a component
        // instantiated, so we let the top of stack be that level by just changing its grey-level.
        if (new_level < next_top->level) {
            new_region = calloc(1, sizeof(MserRegion));
            new_region->level = new_level; new_region->pixel = pixel;
            kv_push(MserRegion*, mser->region_stack, new_region);

            region_merge(new_region, top);
            return;
        }

        // 3. Remove the top of stack and merge it into the second component on stack as follows:
        // Add the first and second moment accumulators together and/or join the pixel lists.
        // Either merge the histories of the components, or take the history from the winner. Note
        // here that the top of stack should be considered one ’time-step’ back, so its current
        // size is part of the history. Therefore the top of stack would be the winner if its
        // current size is larger than the previous size of second on stack.
        region_merge(next_top, top);
    } while (new_level > next_top->level);

}

// Process image
void mser_process(MserData *mser, unsigned char *image) {
    // local variables
    MserRegion* region;
    bool* visited = mser->visited;
    int width = mser->width, height = mser->height;
    // Make the source pixel(with its first edge) the current pixel and store the grey - level of it in the variable current level.
    int current_pixel = 0, neighbor_pixel;
    unsigned char current_edge = 0, current_level = image[0], neighbor_level, new_level;
    int priority = MSER_PIXEL_MAXVALUE, temp;

    // Push a dummy-component onto the stack, with grey-level higher than any allowed in the image.
    kv_push(MserRegion*, mser->region_stack, calloc(1, sizeof(MserRegion)));

    // mark sourse pixel as visited
    visited[current_pixel] = true;

    // 3. Push an empty component with current level onto the component stack.
step_3:
    region = calloc(1, sizeof(MserRegion));
    region->level = current_level; region->pixel = current_pixel;
    kv_push(MserRegion*, mser->region_stack, region);

    // 4. Explore the remaining edges to the neighbors of the current pixel, in order, as follows:
    // For each neighbor, check if the neighbor is already accessible. If it is not, mark it as
    // accessible and retrieve its grey-level. If the grey-level is not lower than the current one,
    // push it onto the heap of boundary pixels. If on the other hand the grey-level is lower than
    // the current one, enter the current pixel back into the queue of boundary pixels for later
    // processing (with the next edge number), consider the new pixel and its grey-level and go to step_3.
    while (true) {
        const int x = current_pixel % width;
        const int y = current_pixel / width;
        for (; current_edge < 4; current_edge++) {
            neighbor_pixel = current_pixel;

            if (current_edge == 0) {
                if (x < width - 1) neighbor_pixel = current_pixel + 1;
            }
            else if (current_edge == 1) {
                if (y < height - 1) neighbor_pixel = current_pixel + width;
            }
            else if (current_edge == 2) {
                if (x > 0) neighbor_pixel = current_pixel - 1;
            }
            else {
                if (y > 0) neighbor_pixel = current_pixel - width;
            }

            if (neighbor_pixel != current_pixel && !visited[neighbor_pixel]) {
                neighbor_level = image[neighbor_pixel];
                visited[neighbor_pixel] = true;
                if (neighbor_level >= current_level) {
                    kv_push(int, mser->boundary_pixels[neighbor_level], neighbor_pixel << 4);

                    if (neighbor_level < priority)
                        priority = neighbor_level;
                }
                else {
                    kv_push(int, mser->boundary_pixels[current_level], (current_pixel << 4) | (current_edge + 1));

                    if (current_level < priority)
                        priority = current_level;

                    current_pixel = neighbor_pixel;
                    current_edge = 0;
                    current_level = neighbor_level;

                    goto step_3;
                }
            }
        }

        // 5. Accumulate the current pixel to the component at the top of the stack (water saturates the current pixel).
        region = kv_A(mser->region_stack, mser->region_stack.n - 1);
        region_accumulate(region, x, y);

        // 6. Pop the heap of boundary pixels. If the heap is empty, we are done. If the returned
        // pixel is at the same grey-level as the previous, go to 4.
        if (priority == 256) {
            //region = mser->region_stack.a[mser->region_stack.n - 1];
            //detect (delta_, minArea_ * width * height, maxArea_ * width * height, maxVariation_, minDiversity_, regions);
            return;
        }

        temp = kv_pop(mser->boundary_pixels[priority]);
        current_pixel = temp >> 4;
        current_edge = temp & 15;

        while (mser->boundary_pixels[priority].n == 0 && (priority < 256))
            ++priority;

        new_level = image[current_pixel];
        if (new_level != current_level) {
            current_level = new_level;

            // 7. The returned pixel is at a higher grey-level, so we must now process
            // all components on the component stack until we reach the higher
            // grey-level. This is done with the processStack sub-routine, see below.
            // Then go to 4.
            process_stack(new_level, current_pixel, mser);
        }
    }
}

int mser(unsigned char *data, int width, int height) {
    int exit_code = 0;
    MserData *mser = NULL;

#ifdef DEBUG
    clock_t start_time = clock();
#endif

    mser = mser_new(width, height);
    if (!mser) {
        exit_code = 1;
        printf("\033[1;31mCould not create a MSER data.\033[0m\n");
        goto done;
    }

    mser_process(mser, (unsigned char *)data);
#ifdef DEBUG
    double diff = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("\033[1;32mmser: elapsed %f s.\033[0m\n", diff);

    write_file(debug_file, mser->debug, width, height);
#endif

done:
    // release filter
    if (mser) {
        mser_delete(mser);
    }
    return exit_code;
}
