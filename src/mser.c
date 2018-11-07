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

#ifndef M_PI
#define M_PI       3.14159265358979323846   // pi
#endif

#ifndef max
#define max(a, b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a, b)            (((a) < (b)) ? (a) : (b))
#endif

#ifdef DEBUG
extern char* debug_file;
#endif // DEBUG

#define MSER_PIXEL_MAXVALUE 256

// MSER statistic
struct _MserStatistic {
    int num_extremal;           // number of extremal regions
    int num_unstable;           // number of unstable extremal regions
    int num_abs_unstable;       // number of regions that failed the absolute stability test
    int num_too_big;            // number of regions that failed the maximum size test
    int num_too_small;          // number of regions that failed the minimum size test
    int num_duplicates;         // number of regions that failed the duplicate test
};
typedef struct _MserStatistic MserStatistic;

// MSER: extremal region
struct _MserRegion {
    int level;
    int pixel;
};
typedef struct _MserRegion MserRegion;

// MSER Pixel
struct _MserPixel {
    int index;               // index of the pixel
    unsigned char intensity; // intensity of the pixel
};
typedef struct _MserPixel MserPixel;

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

int comp(const void* x, const void* y) {
    MserPixel* vx = (MserPixel*)x;
    MserPixel* vy = (MserPixel*)y;
    if (vx->intensity > vy->intensity) return  1;
    if (vx->intensity < vy->intensity) return -1;
    return 0;
}

// Process image
void mser_process(MserData *mser, unsigned char *image) {
    // local variables
    MserRegion* region;
    bool* visited = mser->visited;
    kvec_t(int)* boundary_pixels = mser->boundary_pixels;
    int width = mser->width, height = mser->height;
    // Make the source pixel(with its first edge) the current pixel and store the grey - level of it in the variable current level.
    int current_pixel = 0, current_edge = 0, current_level = image[0];
    int neighbor_pixel, neighbor_level, new_level;
    int priority = MSER_PIXEL_MAXVALUE, temp;

    // Push a dummy-component onto the stack, with grey-level higher than any allowed in the image.
    kv_push(MserRegion, mser->region_stack, malloc(sizeof(MserRegion)));

    // mark sourse pixel as visited
    visited[current_pixel] = true;

    // 3. Push an empty component with current level onto the component stack.
step_3:
    region = malloc(sizeof(MserRegion));
    region->level = current_level; region->pixel = current_pixel;
    kv_push(MserRegion, mser->region_stack, region);

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
                    kv_push(int, boundary_pixels[neighbor_level], neighbor_pixel << 4);

                    if (neighbor_level < priority)
                        priority = neighbor_level;
                }
                else {
                    kv_push(int, boundary_pixels[current_level], (current_pixel << 4) | current_edge + 1);

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
        //region = mser->region_stack.a[mser->region_stack.n - 1];
        //accumulate(x, y);

        // 6. Pop the heap of boundary pixels. If the heap is empty, we are done. If the returned
        // pixel is at the same grey-level as the previous, go to 4.
        if (priority == 256) {
            //region = mser->region_stack.a[mser->region_stack.n - 1];
            //detect (delta_, minArea_ * width * height, maxArea_ * width * height, maxVariation_, minDiversity_, regions);
            return;
        }

        temp = kv_pop(boundary_pixels[priority]);
        current_pixel = temp >> 4;
        current_edge = temp & 15;

        while (boundary_pixels[priority].n == 0 && (priority < 256))
            ++priority;

        int new_level = image[current_pixel];
        if (new_level != current_level) {
            current_level = new_level;

            // 7. The returned pixel is at a higher grey-level, so we must now process
            // all components on the component stack until we reach the higher
            // grey-level. This is done with the processStack sub-routine, see below.
            // Then go to 4.
            //processStack(new_level, current_pixel, mser->region_stack);
        }
    }


    /* shortcuts */
    unsigned int nel = f->nel;
    int ndims = f->ndims;
    int *dims = f->dims;
    MserReg *r = f->r;
    MserPixel* image = f->image;
    MserExtrReg *er = f->er;
    unsigned int *mer = f->mer;
    int delta = f->delta;

#ifdef DEBUG
    unsigned char* debug = f->debug;
#endif

    int njoins = 0;
    // number of extremal regions
    int ner = 0;
    int nmer = 0;
    int nbig = 0;
    int nsmall = 0;
    int nbad = 0;
    int ndup = 0;

    int i, j;
    int dsubx, dsuby, subsx, subsy, idx, r_idx, n_idx, nr_idx, dx, dy;
    unsigned char val, nr_val;

    /* delete any previosuly computed ellipsoid */
    //f->nell = 0;

    // Sort pixels by intensity -----
    for (i = 0; i < nel; i++) {
        image[i].index = i;
        image[i].intensity = im[i];
    }
    qsort(image, nel, sizeof(MserPixel), comp);
    // ------------------------------

    // Compute regions and count extremal regions ---
    for (i = 0; i < (int) nel; ++i) {
        r[i].parent = MSER_VOID_NODE;
    }
    
    // process each pixel
    for (i = 0; i < nel; i++) {
        // index of the current pixel
        idx = image[i].index;
        // intensity of the current pixel
        val = im[idx];
        // index of the root of the current pixel
        r_idx = idx;

        // add the pixel to the forest as a root for now
        r[idx].parent = idx;
        r[idx].shortcut = idx;
        r[idx].area = 1;

        // neighbor index subscript
        dsubx = -1;
        dsuby = -1;
        subsx = idx % f->stride;
        subsy = idx / f->stride;
        for (dy = subsy - 1; dy <= subsy + 1; dy++) {
            for (dx = subsx - 1; dx <= subsx + 1; dx++) {
                if (dx >= 0 && dx < dims[0] && dy >= 0 && dy < dims[1]) {
                    n_idx = dy * f->stride + dx;
                    if (n_idx == idx || r[n_idx].parent == MSER_VOID_NODE)
                        continue;

//#ifdef DEBUG
//                    if (subsx == DEBUG_X && subsy == DEBUG_Y) {
//                        printf("ndx: %d\n", n_idx);
//                    }
//#endif

                    // get roots and optimise path
                    r_idx = find(r, idx);
                    nr_idx = find(r, n_idx);

                    // r_idx and nr_idx are already in the same set
                    if (r_idx == nr_idx)
                        continue;

                    nr_val = im[nr_idx];

                    //1. same intensity
                    //+(r)
                    //?(i) becames root
                    //2. different intensity
                    //+(r) ?(i) becomes root
                    //|(r)   
                    if (nr_val == val) {
                        // ROOT(IDX) becomes the child, optimize the time
                        r[r_idx].parent = nr_idx;
                        r[r_idx].shortcut = nr_idx;
                        r[nr_idx].area += r[r_idx].area;
                    }
                    else {
                        // ROOT(IDX) becomes the parent, extremal region
                        r[nr_idx].parent = r_idx;
                        r[nr_idx].shortcut = r_idx;
                        r[r_idx].area += r[nr_idx].area;

                        // count extremal region
                        ner++;
                    }
                }
            }
        }
    } // next pixel
    f->stats.num_extremal = ++ner; // 1 for root
    //----------------------------------------------------------------------------


    // Extract extremal regions --------------------------------------------------
    // Extremal regions are extracted and stored into the array ER.  The
    // structure R is also updated so that .SHORTCUT indexes the
    // corresponding extremal region if any (otherwise it is set to
    // VOID).
    // allocate memory into er for extremal regions
    if (f->rer < ner) {
        if (er)
            free(er);
        f->er = er = (MserExtrReg *) malloc(sizeof(MserExtrReg) * ner);
        f->rer = ner;
    };

    // save back
    f->nmer = ner;
    printf("mser: extremal regions %d\n", ner);

    // count again
    ner = 0;

    // fills all found extremal regions, fills shortcut in image
    if (er != NULL) {
        for (i = 0; i < (int) nel; ++i) {
            // pop next node xi
            unsigned int idx = image[i].index;
            unsigned char val = im[idx];
            unsigned int p_idx = r[idx].parent;
            unsigned char p_val = im[p_idx];

            // is extremal ?
            int is_extr = (p_val > val) || idx == p_idx;

            if (is_extr) {
                // if so, add it
                er[ner].index = idx;
                er[ner].parent = ner;
                er[ner].value = im[idx];
                er[ner].area = r[idx].area;

                // link this region to this extremal region
                r[idx].shortcut = ner;

                // increase count
                ++ner;
            } else {
                // link this region to void
                r[idx].shortcut = MSER_VOID_NODE;
            }
#ifdef DEBUG
            // draw image with regions and roots
            if (idx == p_idx) {
                debug[idx] = 255;
            }
            else {
                debug[idx] = ner % 254;
            }
#endif
        }
    }
    printf("mser: extremal regions %d\n", ner);


    // Link parent of extremal region ----------------------------------
    /*for (i = 0; i < ner; ++i) {
        unsigned int idx = er[i].index;

        do {
            idx = r[idx].parent;
        } while (r[idx].shortcut == MSER_VOID_NODE);

        er[i].parent = r[idx].shortcut;
        er[i].shortcut = i;
    }
    // -----------------------------------------------------------------

    // Compute variability of +DELTA branches --------------------------
    // For each extremal region Xi of value VAL we look for the biggest
    // parent that has value not greater than VAL+DELTA. This is dubbed
    // `top parent'.
    for (i = 0; i < ner; ++i) {
        // Xj is the current region the region and Xj are the parents
        int top_val = er[i].value + delta;
        int top = er[i].shortcut;

        // examine all parents
        while (1) {
            int next = er[top].parent;
            int next_val = er[next].value;


            // Break if:
            // - there is no node above the top or
            // - the next node is above the top value.
            if (next == top || next_val > top_val)
                break;

            // so next could be the top
            top = next;
        }

        // calculate branch variation
        {
            int area = er[i].area;
            int area_top = er[top].area;
            er[i].variation = (float) (area_top - area) / area;
            er[i].max_stable = 1;
        }


        // Optimization: since extremal regions are processed by
        // increasing intensity, all next extremal regions being processed
        // have value at least equal to the one of Xi. If any of them has
        // parent the parent of Xi (this comprises the parent itself), we
        // can safely skip most intermediate node along the branch and
        // skip directly to the top to start our search.
        {
            int parent = er[i].parent;
            int curr = er[parent].shortcut;
            er[parent].shortcut = MAX(top, curr);
        }
    }
    // ------------------------------------------------------------------------


    // maximally stable branches ----------------------------------------------
    nmer = ner;
    for (i = 0; i < ner; ++i) {
        unsigned int parent = er[i].parent;
        unsigned char val = er[i].value;
        float var = er[i].variation;
        unsigned char p_val = er[parent].value;
        float p_var = er[parent].variation;
        unsigned int loser;


        // Notice that R_parent = R_{l+1} only if p_val = val + 1. If not,
        // this and the parent region coincide and there is nothing to do.
        if (p_val > val + 1)
            continue;

        // decide which one to keep and put that in loser
        if (var < p_var)
            loser = parent;
        else loser = i;

        // make loser NON maximally stable
        if (er[loser].max_stable) {
            --nmer;
            er[loser].max_stable = 0;
        }
    }

    f->stats.num_unstable = ner - nmer;


    // Further filtering ---------------------------------------------
    // It is critical for correct duplicate detection to remove regions
    // from the bottom (smallest one first).
    {
        float max_area = (float) f->max_area * nel;
        float min_area = (float) f->min_area * nel;
        float max_var = (float) f->max_variation;
        float min_div = (float) f->min_diversity;

        // scan all extremal regions (intensity value order)
        for (i = ner - 1; i >= 0L; --i) {
            // process only maximally stable extremal regions
            if (!er[i].max_stable)
                continue;

            if (er[i].variation >= max_var) {
                ++nbad;
                goto remove;
            }
            if (er[i].area > max_area) {
                ++nbig;
                goto remove;
            }
            if (er[i].area < min_area) {
                ++nsmall;
                goto remove;
            }


            // Remove duplicates
            if (min_div < 1.0) {
                unsigned int parent = er[i].parent;
                int area, p_area;
                float div;

                // check all but the root mser
                if ((int) parent != i) {
                    // search for the maximally stable parent region
                    while (!er[parent].max_stable) {
                        unsigned int next = er[parent].parent;
                        if (next == parent)
                            break;
                        parent = next;
                    }

                    // Compare with the parent region; if the current and parent
                    // regions are too similar, keep only the parent.
                    area = er[i].area;
                    p_area = er[parent].area;
                    div = (float) (p_area - area) / (float) p_area;

                    if (div < min_div) {
                        ++ndup;
                        goto remove;
                    }
                } // remove dups end
            }
            continue;
            remove:
            er[i].max_stable = 0;
            --nmer;
        } // check next region

        f->stats.num_abs_unstable = nbad;
        f->stats.num_too_big = nbig;
        f->stats.num_too_small = nsmall;
        f->stats.num_duplicates = ndup;
    }


    // Save the result ---------------------------------------------
    // make room
    if (f->rmer < nmer) {
        if (mer)
            free(mer);
        f->mer = mer = (unsigned int *) malloc(sizeof(unsigned int) * nmer);
        f->rmer = nmer;
    }

    // save back
    f->nmer = nmer;

    j = 0;
    if (er != NULL && mer != NULL) {
        for (i = 0; i < ner; ++i) {
            if (er[i].max_stable)
                mer[j++] = er[i].index;
        }
    }*/
}

//  Fit ellipsoids
/*void mser_ell_fit(MserData *f) {
    // shortcuts
    int nel = f->nel;
    int dof = f->dof;
    int *dims = f->dims;
    int ndims = f->ndims;
    int *subs = f->subs;
    int njoins = f->njoins;
    unsigned int *joins = f->joins;
    MserReg *r = f->r;
    unsigned int *mer = f->mer;
    int nmer = f->nmer;
    mser_acc *acc = f->acc;
    mser_acc *ell = f->ell;

    int d, index, i, j;

    // already fit ?
    if (f->nell == f->nmer)
        return;

    // make room
    if (f->rell < f->nmer) {
        if (f->ell)
            free(f->ell);
        f->ell = (float *) malloc(sizeof(float) * f->nmer * f->dof);
        f->rell = f->nmer;
    }

    if (f->acc == 0) {
        f->acc = (float *) malloc(sizeof(float) * f->nel);
    }

    acc = f->acc;
    ell = f->ell;


    // Integrate moments for each dof
    for (d = 0; d < f->dof; ++d) {
        // start from the upper-left pixel (0,0,...,0)
        memset(subs, 0, sizeof(int) * ndims);

        // step 1: fill acc pretending that each region has only one pixel
        if (d < ndims) {
            // 1-order
            for (index = 0; index < nel; ++index) {
                acc[index] = (float) subs[d];
                adv(ndims, dims, subs);
            }
        } else {
            // 2-order
            // map the dof d to a second order moment E[x_i x_j]
            i = d - ndims;
            j = 0;
            while (i > j) {
                i -= j + 1;
                j++;
            }
            // initialize acc with  x_i * x_j
            for (index = 0; index < nel; ++index) {
                acc[index] = (float) (subs[i] * subs[j]);
                adv(ndims, dims, subs);
            }
        }

        // step 2: integrate
        for (i = 0; i < njoins; ++i) {
            unsigned int index = joins[i];
            unsigned int parent = r[index].parent;
            acc[parent] += acc[index];
        }

        // step 3: save back to ellpises
        for (i = 0; i < nmer; ++i) {
            unsigned int idx = mer[i];
            ell[d + dof * i] = acc[idx];
        }
    }

    // Compute central moments
    for (index = 0; index < nmer; ++index) {
        float *pt = ell + index * dof;
        unsigned int idx = mer[index];
        float area = (float) r[idx].area;

        for (d = 0; d < dof; ++d) {
            pt[d] /= area;

            if (d >= ndims) {
                // remove squared mean from moment to get variance
                i = d - ndims;
                j = 0;
                while (i > j) {
                    i -= j + 1;
                    j++;
                }
                pt[d] -= pt[i] * pt[j];
            }
        }
    }

    f->nell = nmer;
}*/

int mser(unsigned char *data, int width, int height) {
    bool err = false;
    char err_msg[1024];

    int exit_code = 0;
    MserData *filtinv = 0;

    unsigned char *datainv = NULL;
    float const *frames;
    float const *framesinv;
    int nframes = 0, nframesinv = 0;
    int i, dof;

    filtinv = mser_new(width, height);

    if (!filtinv) {
        snprintf(err_msg, sizeof(err_msg), "Could not create an MSER filter.");
        goto done;
    }

    // allocate buffer
    datainv = (unsigned char *) malloc(width * height);
    for (i = 0; i < width * height; i++) {
        datainv[i] = ~data[i]; // 255 - data[i]
    }

    if (!datainv) {
        err = false;
        snprintf(err_msg, sizeof(err_msg), "Could not allocate enough memory.");
        goto done;
    }

#ifdef DEBUG
    clock_t start_time = clock();
#endif
    mser_process(filtinv, (unsigned char *)datainv);
#ifdef DEBUG
    double diff = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("mser: elapsed %f s\n", diff);

    write_file(debug_file, filtinv->debug, width, height);
#endif

    //int  nregionsinv = mser_get_regions_num(filtinv);
    //unsigned int const *  regionsinv = mser_get_regions(filtinv);

    //for (i = 0; i < nregionsinv; ++i) {
    //    printf("%d \t ", -regionsinv[i]);
    //}

    /*mser_ell_fit(filtinv);
    nframesinv = mser_get_ell_num(filtinv);
    dof = mser_get_ell_dof(filtinv);
    const uint8_t colors[3] = {0, 0, 0};
    framesinv = mser_get_ell(filtinv);
    for (i = 0; i < nframesinv; ++i) {
        drawEllipse(framesinv, width, height, depth, data, colors);
        framesinv += dof;
    }*/

done:
    // release filter
    if (filtinv) {
        mser_delete(filtinv);
    }
    if (datainv) {
        free(datainv);
    }
    // if bad print error message
    if (err) {
        fprintf(stderr, "mser: err: %s (%d)\n", err_msg, err);
        exit_code = 1;
    }
    return exit_code;
}
