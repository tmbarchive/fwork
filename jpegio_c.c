#include <stdio.h>
#include <stdlib.h>
#include "jpeglib.h"
#include <setjmp.h>
#include <string.h>

void c_hello(int *i) {
    printf("hello %d\n",*i);
}

void c_hello2(char *s) {
    printf("hello %s\n",s);
}

struct error_mgr {
    struct jpeg_error_mgr pub;
};

static void error_exit(j_common_ptr cinfo) {
    abort();
}

void c_read_jpeg(char *filename,JSAMPLE *image,int *width,int *height) {
    struct jpeg_decompress_struct cinfo;
    struct error_mgr jerr;
    FILE *stream;
    JSAMPARRAY buffer;
    int row_stride;

    *width = -1;
    *height = -1;
    if ((stream = fopen(filename, "rb")) == NULL) {
	fprintf(stderr, "can't open %s\n", filename);
	return;
    }

    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = error_exit;
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, stream);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);
    row_stride = cinfo.output_width * cinfo.output_components;
    buffer = (*cinfo.mem->alloc_sarray)
	((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    while (cinfo.output_scanline < cinfo.output_height) {
	jpeg_read_scanlines(&cinfo, buffer, 1);
	memcpy(image+cinfo.output_scanline*row_stride,
	       buffer[0],row_stride*sizeof(*image));
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    *width = cinfo.image_width;
    *height = cinfo.image_height;
    fclose(stream);
}

void c_write_jpeg(char *filename,int *quality,
		  unsigned char *image,int *width,int *height) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *stream;		/* target file */
    JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
    int row_stride;		/* physical row width in image buffer */
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if ((stream = fopen(filename, "wb")) == NULL) {
	fprintf(stderr, "can't open %s\n", filename);
	return;
    }
    jpeg_stdio_dest(&cinfo, stream);
    cinfo.image_width = *width; 	/* image width and height, in pixels */
    cinfo.image_height = *height;
    cinfo.input_components = 3;		/* # of color components per pixel */
    cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo,*quality, TRUE /* limit to baseline-JPEG values */);
    jpeg_start_compress(&cinfo, TRUE);
    row_stride = (*width) * 3;	/* JSAMPLEs per row in image */
    while (cinfo.next_scanline < cinfo.image_height) {
	row_pointer[0] = & image[cinfo.next_scanline * row_stride];
	jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }
    jpeg_finish_compress(&cinfo);
    fclose(stream);
    jpeg_destroy_compress(&cinfo);
}

