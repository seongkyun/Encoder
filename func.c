#include"head.h"


float GetPSNR(BYTE** img_ori, BYTE** img_dist, int width, int height)		//  PSNR calculation
{
	float mse = 0;
	int i, j;

	for (i = 0; i < height; i++)   // MSE calculation
	{
		for (j = 0; j < width; j++)
		{
			mse += ((img_ori[i][j] - img_dist[i][j]) * (img_ori[i][j] - img_dist[i][j])) / (float)(width*height);
		}
	}
	return 10 * (float)log10((255 * 255) / mse);		// PSNR
}

BYTE** MemAlloc_2D(int width, int height)	//  2D memory allocation
{
	BYTE** arr;
	int i;

	arr = (BYTE**)malloc(sizeof(BYTE*)* height);
	for (i = 0; i < height; i++)
	{
		arr[i] = (BYTE*)malloc(sizeof(BYTE)* width);
	}

	return arr;
}

// 2-D memory allocation for double type
double** MemAlloc_D_2D(int width, int height)
{
	double** arr;
	int i;
	arr = (double**)malloc(height * sizeof(double*));
	for (i = 0; i < height; i++)	arr[i] = (double*)malloc(width * sizeof(double));
	return arr;
}

void MemFree_2D(BYTE** arr, int height)		//  2D memory free
{
	int i;
	for (i = 0; i < height; i++)
	{
		free(arr[i]);
	}
	free(arr);
}


// 2-D memory free for double type
void MemFree_D_2D(double** arr, const int height)
{
	int i;
	for (i = 0; i < height; i++)	free(arr[i]);
	free(arr);
}

// 1 frame read from input file
int Read_Frame(FILE *fp_in, BYTE** img_in, int width, int height)
{
	int i, size = 0;

	for (i = 0; i < height; i++)
	{
		size += fread(img_in[i], sizeof(BYTE), width, fp_in);  // accumulate the reading size
	}

	return size;
}


// 1 frame write on output file
void Write_Frame(FILE* fp_out, BYTE** img_in, int width, int height)
{
	int i;

	for (i = 0; i < height; i++)
	{
		fwrite(img_in[i], sizeof(BYTE), width, fp_out);		// write on the output file
	}
}

// copy the 1 frame
void cpy_frame(BYTE** img_src, BYTE** img_dst, int width, int height)
{
	int i;

	for (i = 0; i < height; i++)
		memcpy(img_dst[i], img_src[i], sizeof(BYTE)* width);
}

void RGB_to_YUV(BYTE** img_in, BYTE** img_out, int width, int height)
{
	int i, j;
	int w[9] = { 66, 129, 25, -38, -74, 112, 112, -94, -18 };		// weight
	int temp[3] = { 0, };

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			temp[0] = w[0] * img_in[i][j] + w[1] * img_in[i + height][j] + w[2] * img_in[i + height * 2][j] + 128;
			temp[1] = w[3] * img_in[i][j] + w[4] * img_in[i + height][j] + w[5] * img_in[i + 2 * height][j] + 128;
			temp[2] = w[6] * img_in[i][j] + w[7] * img_in[i + height][j] + w[8] * img_in[i + 2 * height][j] + 128;

			img_out[i][j] = (BYTE)(temp[0] >> 8) + 16;
			img_out[i + height][j] = (BYTE)(temp[1] >> 8) + 128;
			img_out[i + 2 * height][j] = (BYTE)(temp[2] >> 8) + 128;
		}
	}
}


void YUV_to_RGB(BYTE** img_in, BYTE** img_out, int width, int height)
{
	int i, j;
	int w[5] = { 298, 409, -100, -208, 516 };	// weight
	int temp[3] = { 0, };

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			temp[0] = w[0] * (img_in[i][j] - 16) + w[1] * (img_in[i + height * 2][j] - 128) + 128;
			temp[1] = w[0] * (img_in[i][j] - 16) + w[2] * (img_in[i + height][j] - 128) + w[3] * (img_in[i + 2 * height][j] - 128) + 128;
			temp[2] = w[0] * (img_in[i][j] - 16) + w[4] * (img_in[i + height][j] - 128) + 128;

			img_out[i][j] = (BYTE)Clip((temp[0] >> 8));
			img_out[i + height][j] = (BYTE)Clip((temp[1] >> 8));
			img_out[i + 2 * height][j] = (BYTE)Clip((temp[2] >> 8));
		}
	}
}

// YUV 444 -> YUV 420
void YUV444_to_420(BYTE** img_in, BYTE** img_Y, BYTE** img_U420, BYTE** img_V420, int width, int height)
{
	int i, j;	// Loop index

	// Y component copy
	for (i = 0; i < height; i++)
	{
		memcpy(img_Y[i], img_in[i], sizeof(BYTE)* width);
	}

	//chroma sub sampling
	for (i = 0; i < height; i += 2)
	{
		for (j = 0; j < width; j += 2)
		{
			img_U420[i >> 1][j >> 1] = (BYTE)((img_in[i + height][j] + img_in[i + height + 1][j]) / 2);		// Cb calculate
			img_V420[i >> 1][j >> 1] = (BYTE)((img_in[i + height * 2][j] + img_in[i + height * 2 + 1][j]) / 2);		// Cr calculate
		}
	}
}


// YUV 420 -> YUV 444
void YUV420_to_444(BYTE** img_Y, BYTE** img_U420, BYTE** img_V420, BYTE** img_out, int width, int height)
{
	int i, j, m, n;

	// Y component copy
	for (i = 0; i < height; i++)
	{
		memcpy(img_out[i], img_Y[i], sizeof(BYTE)* width);
	}


	//chroma recon
	for (i = 0; i < height; i += 2)
	{
		for (j = 0; j < width; j += 2)
		{
			for (m = 0; m < 2; m++)
			for (n = 0; n < 2; n++)
			{
				img_out[i + m + height][j + n] = img_U420[i >> 1][j >> 1];	// Cb copy interpolation
				img_out[i + m + height * 2][j + n] = img_V420[i >> 1][j >> 1];	// Cr copy interpolation
			}
		}
	}
}

/*
Inter-prediction function
input : original image, reference image, image width & height, prediction block size, maximum search range
output: prediction image, residual image, reconstruction image
*/
void InterPrediction(BYTE** img_ori, BYTE** img_ref, BYTE** img_pred, double** img_resi, BYTE** img_recon, int width, int height, int block_size, int search_range)
{
	int i, j, m, n, x, y;										// Loop index
	int k, l;													// motion vector position
	int SR_left = 0, SR_right = 0, SR_top = 0, SR_bottom = 0;	// Search range variable

	float min_MAE;		//memory for minimum MAE value
	float temp_MAE;		//MAE temporal memory
	static MV mv[HEIGHT / BLOCK_SIZE][WIDTH / BLOCK_SIZE]; // motion vector memory

	for (i = 0; i < height; i += block_size)
	{
		for (j = 0; j < width; j += block_size)
		{
			// motion vector initialization
			k = (int)(i / block_size);
			l = (int)(j / block_size);

			mv[k][l].x = 0;
			mv[k][l].y = 0;

			//Adaptive Search Range Decision
			if (i < search_range)SR_top = -i;
			else  SR_top = -search_range;

			if (i + block_size > height - search_range)SR_bottom = height - (i + block_size);
			else  SR_bottom = search_range;

			if (j < search_range)SR_left = -j;
			else  SR_left = -search_range;

			if (j + block_size > width - search_range) SR_right = width - (j + block_size);
			else  SR_right = search_range;


			min_MAE = 255;	// Maximum MAE

			//Motion estimation
			for (m = SR_top; m <= SR_bottom; m++)
			{
				for (n = SR_left; n <= SR_right; n++)
				{
					temp_MAE = 0;

					//MAE calculation
					for (y = 0; y < block_size; y++)
					for (x = 0; x < block_size; x++)
					{
						temp_MAE += (float)abs(img_ref[i + m + y][j + n + x] - img_ori[i + y][j + x]) / (block_size*block_size);
					}

					// Best prediction block & motion vector
					if (min_MAE > temp_MAE)
					{
						min_MAE = temp_MAE;
						mv[k][l].x = n;
						mv[k][l].y = m;
					}
				}
			}

			// Best prediction & recon block copy
			for (m = 0; m < block_size; m++)
			{
				for (n = 0; n < block_size; n++)
				{
					img_pred[i + m][j + n] = img_ref[i + m + mv[k][l].y][j + n + mv[k][l].x];
					img_resi[i + m][j + n] = (double)(img_ori[i + m][j + n] - img_pred[i + m][j + n]);
				}
			}
		}
	}

	{
		double** coeff = MemAlloc_D_2D(width, height);

		//residual transform & quantization
		SeparableBlockFDCT_2D(img_resi, coeff, width, height, Transform_BLOCK, QP);

		SeparableBlockIDCT_2D(coeff, img_resi, width, height, Transform_BLOCK, QP);



		for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			if (img_resi[i][j] < 0)
				img_recon[i][j] = Clip(img_pred[i][j] + (int)-(-img_resi[i][j] + 0.5));
			else
				img_recon[i][j] = Clip(img_pred[i][j] + (int)(img_resi[i][j] + 0.5));
		}
		MemFree_D_2D(coeff, height);
	}

	{
		FILE* fp_mv = fopen("Motion_information.xls", "a");

		for (i = 0; i<height / block_size; i++)
		{
			for (j = 0; j<width / block_size; j++)
			{
				fprintf(fp_mv, "(%d,%d)\n", mv[i][j].x, mv[i][j].y);
			}
		}
	}
}



void Encode(BYTE** img_ori, BYTE** img_pred, double** img_resi, BYTE** img_recon, int width, int height)
{
	int i, j, m, n;
	int best_mode;

	static BYTE ori[BLOCK_SIZE*BLOCK_SIZE];
	static BYTE ref[BLOCK_SIZE * 3 + 1];
	static BYTE pred[NMODE][BLOCK_SIZE*BLOCK_SIZE];
	static BYTE recon[NMODE][BLOCK_SIZE*BLOCK_SIZE];
	static int  resi[NMODE][BLOCK_SIZE*BLOCK_SIZE];


	BYTE** img_padding = MemAlloc_2D(width + 1, height + 1);

	for (i = 0; i< height; i++)
	for (j = 0; j<width; j++)
		img_padding[i + 1][j + 1] = img_ori[i][j];

	for (i = 0; i<height; i++)
		img_padding[i + 1][0] = 128;

	for (i = 0; i<height + 1; i++)
		img_padding[0][i] = 128;


	//intra prediction loop
	for (i = 0; i < height; i += BLOCK_SIZE)
	{
		for (j = 0; j < width; j += BLOCK_SIZE)
		{
			// get original block 
			for (m = 0; m<BLOCK_SIZE; m++)
			for (n = 0; n<BLOCK_SIZE; n++)
			{
				ori[m * BLOCK_SIZE + n] = img_ori[i + m][j + n];
			}

			// get reference samples
			if (j != width - BLOCK_SIZE)
			{
				for (m = 0; m<2 * BLOCK_SIZE + 1; m++)
					ref[m] = img_padding[i][j + m];

			}
			else
			{
				for (m = 0; m<BLOCK_SIZE + 1; m++)
					ref[m] = img_padding[i][j + m];
				for (m = BLOCK_SIZE + 1; m<2 * BLOCK_SIZE + 1; m++)
					ref[m] = 128;
			}
			if (j != 0)
			{
				for (m = 0; m<BLOCK_SIZE; m++)
					ref[m + (2 * BLOCK_SIZE) + 1] = img_padding[(i + 1) + m][j];
			}
			else
			{
				for (m = 0; m<BLOCK_SIZE; m++)
					ref[m + (2 * BLOCK_SIZE) + 1] = 128;
			}
			// serch for best mode
			best_mode = intra_prediction(ori, ref,			// input
				pred, resi, recon);	// output

			// generate reconstructed image
			for (m = 0; m<BLOCK_SIZE; m++)
			for (n = 0; n<BLOCK_SIZE; n++)
			{
				img_pred[i + m][j + n] = pred[best_mode][m * BLOCK_SIZE + n];
				img_resi[i + m][j + n] = (double)resi[best_mode][m * BLOCK_SIZE + n];
			}

		}
	}

	{
		double** coeff = MemAlloc_D_2D(width, height);

		//residual transform & quantization
		SeparableBlockFDCT_2D(img_resi, coeff, width, height, Transform_BLOCK, QP);
		SeparableBlockIDCT_2D(coeff, img_resi, width, height, Transform_BLOCK, QP);


		//Reconstruction
		for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			if (img_resi[i][j] < 0)
				img_recon[i][j] = Clip(img_pred[i][j] + (int)-(-img_resi[i][j] + 0.5));
			else
				img_recon[i][j] = Clip(img_pred[i][j] + (int)(img_resi[i][j] + 0.5));
		}
		MemFree_D_2D(coeff, height);
	}
	MemFree_2D(img_padding, height + 1);
}

int intra_dc(BYTE* ori, BYTE* ref,							 // input
	BYTE* pred, int* resi, BYTE* recon)			 // output
{
	int n;
	int SAD = 0;
	int mean = 0;

	for (n = 0; n < BLOCK_SIZE + 1; n++)
		mean += ref[n];
	for (n = BLOCK_SIZE + 1; n < 2 * BLOCK_SIZE + 1; n++)
		mean += ref[n];
	mean /= (2 * BLOCK_SIZE + 1);

	// prediction
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		resi[n] = ori[n] - mean;

	// get predictor
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		pred[n] = ori[n] - resi[n];

	// reconstruct image
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		recon[n] = pred[n] + resi[n];

	// calculate SAD
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
	{
		SAD += abs(resi[n]);
	}

	return SAD;
}

int intra_hor(BYTE* ori, BYTE* ref, 					     // input
	BYTE* pred, int* resi, BYTE* recon)		     // output
{
	int n, m;
	int SAD = 0;

	//prediction
	for (m = 0; m<BLOCK_SIZE; m++)
	for (n = 0; n<BLOCK_SIZE; n++)
	{
		resi[m * BLOCK_SIZE + n] = ori[m * BLOCK_SIZE + n] - ref[m + 2 * BLOCK_SIZE + 1];
	}

	// get predictor
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		pred[n] = ori[n] - resi[n];

	// reconstruct image
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		recon[n] = pred[n] + resi[n];

	// calculate SAD
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
	{
		SAD += abs(resi[n]);
	}

	return SAD;

}

int intra_ver(BYTE* ori, BYTE* ref, 					  // input
	BYTE* pred, int* resi, BYTE* recon)		  // output
{
	int m, n;
	int SAD = 0;

	for (m = 0; m<BLOCK_SIZE; m++)
	for (n = 0; n<BLOCK_SIZE; n++)
	{
		resi[m * BLOCK_SIZE + n] = ori[m * BLOCK_SIZE + n] - ref[n + 1];
	}

	// get predictor
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		pred[n] = ori[n] - resi[n];

	// reconstruct image
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		recon[n] = pred[n] + resi[n];

	// calculate SAD
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
	{
		SAD += abs(resi[n]);
	}

	return SAD;
}

int intra_DL(BYTE* ori, BYTE* ref, 					      // input
	BYTE* pred, int* resi, BYTE* recon)		  // output
{
	int n, m;
	int SAD = 0;

	for (m = 0; m<BLOCK_SIZE; m++)
	for (n = 0; n<BLOCK_SIZE; n++)
	{
		resi[m * BLOCK_SIZE + n] = ori[m * BLOCK_SIZE + n] - ref[n + (m + 1) + 1];
	}

	// get predictor
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		pred[n] = ori[n] - resi[n];

	// reconstruct image
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		recon[n] = pred[n] + resi[n];

	// calculate SAD
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
	{
		SAD += abs(resi[n]);
	}

	return SAD;
}
int intra_DR(BYTE* ori, BYTE* ref, 					  // input
	BYTE* pred, int* resi, BYTE* recon)	  // output
{
	int n, m;
	int SAD = 0;

	for (m = 0; m<BLOCK_SIZE; m++)
	for (n = 0; n<BLOCK_SIZE; n++)
	{
		if (n - m >= 0)
			resi[m * BLOCK_SIZE + n] = ori[m * BLOCK_SIZE + n] - ref[n - m];
		else
		{
			resi[m * BLOCK_SIZE + n] = ori[m * BLOCK_SIZE + n] - ref[-n + m + 2 * BLOCK_SIZE + 1];
		}
	}

	// get predictor
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		pred[n] = ori[n] - resi[n];

	// reconstruct image
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
		recon[n] = pred[n] + resi[n];

	// calculate SAD
	for (n = 0; n < BLOCK_SIZE*BLOCK_SIZE; n++)
	{
		SAD += abs(resi[n]);
	}

	return SAD;
}
int intra_prediction(BYTE* ori, BYTE* ref,										 //intput
	BYTE(pred)[NMODE][BLOCK_SIZE*BLOCK_SIZE], int(resi)[NMODE][BLOCK_SIZE*BLOCK_SIZE], BYTE(recon)[NMODE][BLOCK_SIZE*BLOCK_SIZE]) // output
{

	static int SAD[NMODE];
	int min_SAD, best_mode, n;

	SAD[0] = intra_ver(ori, ref, pred[0], resi[0], recon[0]);

	SAD[1] = intra_hor(ori, ref, pred[1], resi[1], recon[1]);

	SAD[2] = intra_dc(ori, ref, pred[2], resi[2], recon[2]);

	SAD[3] = intra_DL(ori, ref, pred[3], resi[3], recon[3]);

	SAD[4] = intra_DR(ori, ref, pred[4], resi[4], recon[4]);

	best_mode = 0;
	min_SAD = SAD[0];

	for (n = 1; n < NMODE; n++)
	{
		if (min_SAD > SAD[n])
		{
			min_SAD = SAD[n];
			best_mode = n;
		}
	}
	return best_mode;
}

// separable block-based forward DCT
void SeparableBlockFDCT_2D(const double** input, double** coeff, const int width, const int height, const int block_size, int qp)
{
	int i, j;
	double **temp_hor = MemAlloc_D_2D(block_size, block_size);
	double **temp_ver = MemAlloc_D_2D(block_size, block_size);

	for (i = 0; i<height; i++)memset(coeff[i], 0, sizeof(double)*width);

	for (i = 0; i < height; i += block_size){
		for (j = 0; j < width; j += block_size){
			int k;
			for (k = 0; k < block_size; k++)	FDCT_1D(&input[i + k][j], temp_hor[k], block_size);
			MatTranspose(temp_hor, block_size);
			for (k = 0; k < block_size; k++)	FDCT_1D(temp_hor[k], temp_ver[k], block_size);
			MatTranspose(temp_ver, block_size);

			for (k = 0; k < block_size; k++)memcpy(&coeff[i + k][j], temp_ver[k], block_size*sizeof(double));

		}
	}


	// Quantization

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			coeff[i][j] = ((int)coeff[i][j] > 0 ? (int)(coeff[i][j] + qp / 2) / qp : (int)(coeff[i][j] - qp / 2) / qp);

		}
	}

	MemFree_D_2D(temp_hor, block_size);
	MemFree_D_2D(temp_ver, block_size);
}

// separable block-based inverse DCT
void SeparableBlockIDCT_2D(double** coeff, double** output, const int width, const int height, const int block_size, int qp)
{
	int i, j;
	double **temp_hor = MemAlloc_D_2D(block_size, block_size);
	double **temp_ver = MemAlloc_D_2D(block_size, block_size);

	// Inverse Quantization

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			coeff[i][j] = coeff[i][j] * qp;
		}
	}

	for (i = 0; i < height; i += block_size){
		for (j = 0; j < width; j += block_size){
			int k;
			for (k = 0; k < block_size; k++)	IDCT_1D(&coeff[i + k][j], temp_hor[k], block_size);
			MatTranspose(temp_hor, block_size);
			for (k = 0; k < block_size; k++)	IDCT_1D(temp_hor[k], temp_ver[k], block_size);
			MatTranspose(temp_ver, block_size);
			for (k = 0; k < block_size; k++)	memcpy(&output[i + k][j], temp_ver[k], block_size*sizeof(double));
		}
	}

	MemFree_D_2D(temp_hor, block_size);
	MemFree_D_2D(temp_ver, block_size);
}

// N-point 1-D forward DCT
void FDCT_1D(const double* input, double* coeff, const int N)
{
	const double PI = 3.1415926535;
	int n, k;

	for (k = 0; k < N; k++){
		double beta = k == 0 ? 1 / sqrt(2.0) : 1;
		double temp = 0;
		for (n = 0; n < N; n++){
			double basis = cos(((2 * n + 1)*PI*k) / (2.0*N));
			temp += input[n] * basis;
		}
		temp *= sqrt(2 / (double)N) * beta;

		coeff[k] = temp;
	}
}

// N-point 1-D inverse DCT
void IDCT_1D(const double* coeff, double* output, const int N)
{
	const double PI = 3.1415926535;
	int n, k;

	for (n = 0; n < N; n++){
		double temp = 0;
		for (k = 0; k < N; k++){
			double beta = k == 0 ? 1 / sqrt(2.0) : 1;
			double basis = cos(((2 * n + 1)*PI*k) / (2.0*N));
			temp += beta * coeff[k] * basis;
		}
		temp *= sqrt(2 / (double)N);

		output[n] = temp;
	}
}

// matrix transpose
void MatTranspose(double **mat, const int size)
{
	int i, j;
	for (i = 0; i < size; i++){
		for (j = i + 1; j < size; j++){
			double temp = mat[i][j];
			mat[i][j] = mat[j][i];
			mat[j][i] = temp;
		}
	}
}