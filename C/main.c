#include<stddef.h>
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<string.h>
#include<complex.h>
#include<math.h>

/* ------------------------ */
/*          IMAGE           */
/* ------------------------ */

typedef unsigned char img_byte;

// Image is Grayscale
typedef struct {
    size_t width;
    size_t height;
    img_byte *pixels;
} Image;

static Image *make_image(size_t width, size_t height) {
    Image *img = malloc(sizeof(Image));
    img->width = width;
    img->height = height;

    img->pixels = calloc(width * height,sizeof(img_byte));
    return img;
}

static void image_set_pixel(Image *img, size_t x, size_t y, img_byte p) {
    assert(x < img->width && x >= 0);    
    assert(y < img->height && y >= 0);

    img->pixels[y * img->width + x] = p;
}

static img_byte image_get_pixel(Image *img, size_t x, size_t y) {
    assert(x < img->width && x >= 0);    
    assert(y < img->height && y >= 0);
    
    return img->pixels[y * img->width + x];
}

static void delete_image(Image *img) {
    free(img->pixels);
    free(img);
}

#define BUF_LEN 256
static Image *load_pgm(char *path) {
    // https://rosettacode.org/wiki/Bitmap/Read_a_PPM_file#C    
    char buf[BUF_LEN];

    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return NULL;

    char *err = fgets(buf, BUF_LEN, fp);
    if (err == NULL) return NULL;
    // Assert if type is not graymap, fuck color
    if (strncmp(buf, "P5\n", 3) != 0) return NULL;

    // Skip Comments
    do {
        char *err = fgets(buf, BUF_LEN, fp);
        if (err == NULL) return NULL;
    } while (strncmp(buf, "#", 1) == 0);

    // Read width and height
    int width, height;
    int c = sscanf(buf, "%u %u", &width, &height);
    if (c < 2) return NULL;

    // Read maximum value
    int max;
    c = fscanf(fp, "%u", &max);
    if (c < 2 && max != 255) return NULL;
    // Skip one spacer
    fseek(fp, 1, SEEK_CUR);

    Image *img = make_image(width, height);
    
    // Read the actual image
    size_t pc = fread(img->pixels, sizeof(img_byte), width*height, fp);
    if (pc < width*height) {
        delete_image(img);
        return NULL;
    }

    fclose(fp);

    return img;
}

static void save_pgm(char *path, Image *img) {
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return;

    fprintf(fp, "P5\n%zu %zu\n255\n", img->width, img->height);

    fwrite(img->pixels, 1, img->width*img->height, fp);

    fclose(fp);
}

/* ------------------------ */
/*           FFT            */
/* ------------------------ */

static void fft_1D(double complex buffer[], size_t n) {
    // Testing if n is a power of two
    // Explanation on StackOverflow, eventhough it should be pretty obvious
    // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
    assert((n != 0) && ((n & (n - 1)) == 0));

    // Bit-reverse the input buffer
    int j = 0;
    for (int k = 1; k < n; k++) {
        // We swap the values at k with the value at bit reversed k (j) for all k < j

        /*
          Smart reverse increment algortithm I found on StackOverflow & Github
          TLDR: We are adding one, just in reverse. And because cpus dont like to go in reverse, we do it manually
        */
        
        // Take the the zero and make the leftmost bit one, where the leftmost bit is the bit left of our new reveresed zero
        // e.g. 10000 This would be our number, lets call it n (below it's called bit, but bit gets confusing when talking)
        //       ^ This would be our new reverse LSB, which now counts to the right
        // We can directly assign n here, because n has to be a power of two, thus is our perfect single bit
        int bit = n;

        do {
            // Take n and shift it to the right by one
            bit >>= 1;
            
            // n is now basically gonna act as a carry bit and the number added at the same time, we XOR it to add it to the number j.
            // If the bit that's 1 in n and and the bit in j at the same position, are both 1, we have to carry like in normal addition, otherwise we are fine. That's what the XOR does.
            j ^= bit;

            // To test if we carried and we have to keep going or if we are finished, we just look at the conjunction of n and j
            // Because if we carried, the bit position that's one in n is now going to be zero in j (1 ^ 1 = 0), so the conjunction will be zero, otherwise n
        } while ((j & bit) == 0 && bit != 1); // (bit != 1) is the edge case

        // Only do it if k < j so we don't double swap
        if (k < j) {
            // Now swap the two numbers
            complex double temp;
            temp = buffer[j];
            buffer[j] = buffer[k];
            buffer[k] = temp;
        }
    }

    // Actual in-place radix-2 FFT
    // Or how i call it: Black magic
    // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Data_reordering,_bit_reversal,_and_in-place_algorithms
    for (int m = 2; m <= n; m <<= 1) {
        double complex w_m = cexp(-2.0 * M_PI * I / m);

        for (int k = 0; k < n; k += m) {
            double complex w = 1;
            for (int j = 0; j < m/2; j++) {
                double complex t = w * buffer[k + j + m/2];
                double complex u = buffer[k + j];
                buffer[k + j] = u + t;
                buffer[k + j + m/2] = u - t;
                w *= w_m;
            }
        }
    }
}

static void transpose(double complex buffer[], size_t width, size_t height) {
    double complex temp;

    for (int i = 0; i < height; i++) {
        for (int k = i+1; k < width; k++) {
            // buffer[i][k] <=> buffer[k][i]
            temp = buffer[i * width + k];
            buffer[i * width + k] = buffer[k * width + i];
            buffer[k * width + i] = temp;
        }
    }
}

static void fft_2D(double complex buffer[], size_t width, size_t height) {
    // Do the rows first
    for (int i = 0; i < height; i++) {
        fft_1D(buffer+i*width, width);
    }

    // Transpose the array to be able to use the same function again
    transpose(buffer, width, height);

    // Now do the columns
    for (int i = 0; i < width; i++) {
        fft_1D(buffer+i*height, height);
    }

    // Transpose back
    transpose(buffer, width, height);
}


Image *fft(Image *img) {
    // Ensure that width and height is a power of two
    assert((img->width != 0) && ((img->width & (img->width - 1)) == 0));
    assert((img->height != 0) && ((img->height & (img->height - 1)) == 0));

    double complex *buffer = calloc(img->width * img->height, sizeof(double complex));
    for (int i = 0; i < img->height * img->width; i++) buffer[i] = img->pixels[i]; 

    fft_2D(buffer, img->width, img->height);

    // Logarithmic scaling so you can actually see anything
    double max = 0;
    for (int i = 0; i < img->width * img->height; i++) {
        double amp = cabs(buffer[i]);
        if (amp > max) max = amp;
    }

    double c = 255 / log(1 + fabs(max));

    Image *transformed = make_image(img->width, img->height);

    for (int i = 0; i < img->width * img->height; i++) {
        double amp = cabs(buffer[i]);
        double processed = c * log(1 + amp);

        int m = i % img->width;
        int n = (int) (i / img->width);
        int x = m < img->width/2 ? m + img->width/2 : m - img->width/2;
        int y = n < img->height/2 ? n + img->height/2 : n - img->height/2;
        image_set_pixel(transformed, x, y, processed);
    }

    free(buffer);
    delete_image(img);

    return transformed;
}


int main() {
    Image *img = load_pgm("../cameraman.pgm");

    Image *out = fft(img);

    save_pgm("out.pgm", out);
    
    delete_image(out);
    return 0;
}
