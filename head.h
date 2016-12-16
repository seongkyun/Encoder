#pragma warning(disable:4996)
#include <stdio.h>
#include <math.h>				  //  header file
#include <stdlib.h>
#include <string.h>

//Parameter
#define WIDTH  352                //  CIF frame size
#define HEIGHT 288

#define BLOCK_SIZE       4								//prediction block size
#define NMODE            5								//the number of intra mode
#define SR		         16								//inter mode Search Range
#define QP		         16								//Quantization Parameter
#define Transform_BLOCK  8							    //Transform Block size
#define cWIDTH		     (WIDTH>>1)						//Chroma frame size
#define cHEIGHT		     (HEIGHT>>1)					//Chroma frame size

#define cBLOCK_SIZE (BLOCK_SIZE>>1)					    //Chroma prediction block size
#define cSR			(SR>>1)							    //Chroma Search Range

#define Clip(x) ( x < 0 ? 0 : ( x > 255 ? 255 : x))     //Clipping operation

typedef unsigned char BYTE;

typedef struct MV	// motion vector structure
{
	int x, y;
}MV;

BYTE** MemAlloc_2D(int width, int height);					// 2D memory allocation
void MemFree_2D(BYTE** arr, int height);					// 2D memory free

double** MemAlloc_D_2D(int width, int height);			// 2-D memory allocation for double type
void MemFree_D_2D(double** arr, const int height);					// 2-D memory free for double type

float GetPSNR(BYTE** img_ori, BYTE** img_dist, int width, int height);	//PSNR value
void cpy_frame(BYTE** img_src, BYTE** img_dst, int width, int height);  //copy the 1 frame
int Read_Frame(FILE *fp_in, BYTE** img_in, int width, int height);											// 1 frame read from input file
void Write_Frame(FILE *fp_out, BYTE** img_in, int width, int height);										// 1 frame write on output file
void RGB_to_YUV(BYTE** img_in, BYTE** img_out, int height, int width);										// Image color conversion RGB444 to YUV444
void YUV_to_RGB(BYTE** img_in, BYTE** img_out, int width, int height);										// Image color conversion YUV444 to RGB444
void YUV444_to_420(BYTE** img_in, BYTE** img_Y, BYTE** img_U420, BYTE** img_V420, int width, int height);	// Chroma sampling   4:4:4 -> 4:2:0
void YUV420_to_444(BYTE** img_Y, BYTE** img_U420, BYTE** img_V420, BYTE** img_out, int width, int height);	// Chroma sampling   4:2:0 -> 4:4:4 

void InterPrediction(BYTE** img_ori, BYTE** img_ref, BYTE** img_pred, int** img_resi, BYTE** img_recon, int width, int height, int block_size, int search_range);

int   intra_prediction(BYTE* ori, BYTE* ref, BYTE(pred)[NMODE][BLOCK_SIZE*BLOCK_SIZE], int(resi)[NMODE][BLOCK_SIZE*BLOCK_SIZE], BYTE(recon)[NMODE][BLOCK_SIZE*BLOCK_SIZE]);
int   intra_dc(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int   intra_hor(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int   intra_ver(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int   intra_DL(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int   intra_DR(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
void Encode(BYTE** img_ori, BYTE** img_pred, int** img_resi, BYTE** img_recon, int width, int height);

// separable block-based forward DCT
void SeparableBlockFDCT_2D(const double** input, double** coeff, const int width, const int height, const int block_size, int qp);
// separable block-based inverse DCT
void SeparableBlockIDCT_2D(double** coeff, double** output, const int width, const int height, const int block_size, int qp);

void FDCT_1D(const double* input, double* coeff, const int N);			// N-point 1-D forward DCT
void IDCT_1D(const double* coeff, double* output, const int N);			// N-point 1-D inverse DCT

// matrix operation
void MatTranspose(double **mat, const int size);		// matrix transpose
