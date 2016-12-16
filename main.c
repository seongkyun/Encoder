#include"head.h"


int main()
{
	FILE *fp_in0 = fopen("Suzie_CIF_150_30.rgb", "rb");			//in Frame number 1  RGB file
	FILE *fp_out0 = fopen("[predc]Suzie_CIF_13.rgb", "wb");		//Predictor RGB file
	FILE *fp_out1 = fopen("[Recon]Suzie_CIF_13.rgb", "wb");		//recon     RGB file


	BYTE **img_YUV444, **img_RGB;						//input original RGB, YUV444
	BYTE **img_ref_Y, **img_ref_U, **img_ref_V;			//input reference YUV420
	BYTE **img_ori_Y, **img_ori_U, **img_ori_V;			//input original  YUV420

	BYTE **img_recon_Y, **img_recon_U, **img_recon_V;	//recon	pointer
	BYTE **img_pred_Y, **img_pred_U, **img_pred_V;		//prediction pointer
	double **img_resi_Y, **img_resi_U, **img_resi_V;		//residual pointer

	int frame = 0, read_size = 1;

	img_YUV444 = MemAlloc_2D(WIDTH, HEIGHT * 3); // YUV 444 memory	
	img_RGB = MemAlloc_2D(WIDTH, HEIGHT * 3); // RGB memory

	img_pred_Y = MemAlloc_2D(WIDTH, HEIGHT);	 // Y component memory
	img_recon_Y = MemAlloc_2D(WIDTH, HEIGHT);	 // Y component memory
	img_resi_Y = MemAlloc_D_2D(WIDTH, HEIGHT);	 // Y component memory

	img_pred_U = MemAlloc_2D(cWIDTH, cHEIGHT);	 // U component memory
	img_recon_U = MemAlloc_2D(cWIDTH, cHEIGHT);	 // U component memory
	img_resi_U = MemAlloc_D_2D(cWIDTH, cHEIGHT);	 // U component memory

	img_pred_V = MemAlloc_2D(cWIDTH, cHEIGHT);	 // V component memory
	img_recon_V = MemAlloc_2D(cWIDTH, cHEIGHT);	 // V component memory
	img_resi_V = MemAlloc_D_2D(cWIDTH, cHEIGHT);	 // V component memory

	img_ref_Y = MemAlloc_2D(WIDTH, HEIGHT);		 // reference picture memory
	img_ref_U = MemAlloc_2D(cWIDTH, cHEIGHT);	 // reference picture memory
	img_ref_V = MemAlloc_2D(cWIDTH, cHEIGHT);	 // reference picture memory

	img_ori_Y = MemAlloc_2D(WIDTH, HEIGHT);		 // original picture memory	
	img_ori_U = MemAlloc_2D(cWIDTH, cHEIGHT);	 // original picture memory
	img_ori_V = MemAlloc_2D(cWIDTH, cHEIGHT);	 // original picture memory

	///////////////////////////////////////////////////////////////////////////////////////
	printf("		======= 2016 DSPE Simple Encoder =========\n\n");
	while (read_size = Read_Frame(fp_in0, img_RGB, WIDTH, HEIGHT * 3))		//read reference picture
	{
		if (frame == 0)// intra mode
		{
			/////////////////////////////////////////////////
			//                 Assignment 1                // ori (RGB->YUV)->(YUV444->YUV420)->(encode/intraprediction)->(predictor, recon)->(YUV420->RGB)->출력
			/////////////////////////////////////////////////

			RGB_to_YUV(img_RGB, img_YUV444, WIDTH, HEIGHT);
			YUV444_to_420(img_YUV444, img_ori_Y, img_ori_U, img_ori_V, WIDTH, HEIGHT);

			Encode(img_ori_Y, img_pred_Y, img_resi_Y, img_recon_Y, WIDTH, HEIGHT);
			Encode(img_ori_U, img_pred_U, img_resi_U, img_recon_U, cWIDTH, cHEIGHT);
			Encode(img_ori_V, img_pred_V, img_resi_V, img_recon_V, cWIDTH, cHEIGHT);
			
			YUV420_to_444(img_pred_Y, img_pred_U, img_pred_V, img_YUV444, WIDTH, HEIGHT);
			YUV_to_RGB(img_YUV444, img_RGB, WIDTH, HEIGHT);//predictor 변환
			Write_Frame(fp_out0, img_RGB, WIDTH, HEIGHT * 3);//RGB predictor 출력

			YUV420_to_444(img_recon_Y, img_recon_U, img_recon_V, img_YUV444, WIDTH, HEIGHT);
			YUV_to_RGB(img_YUV444, img_RGB, WIDTH, HEIGHT);//predictor 변환
			Write_Frame(fp_out1, img_RGB, WIDTH, HEIGHT * 3);//RGB recon 출력

			/*//Encoding Information//*/
			printf("I - Frame[POC %d] - PSNR [Y : %.2f dB, U : %.2f dB, V : %.2f dB]  \n", frame, GetPSNR(img_ori_Y, img_recon_Y, WIDTH, HEIGHT), GetPSNR(img_ori_U, img_recon_U, cWIDTH, cHEIGHT), GetPSNR(img_ori_V, img_recon_V, cWIDTH, cHEIGHT));

			cpy_frame(img_recon_Y, img_ref_Y, WIDTH, HEIGHT);
			cpy_frame(img_recon_U, img_ref_U, cWIDTH, cHEIGHT);
			cpy_frame(img_recon_V, img_ref_V, cWIDTH, cHEIGHT);

			frame++;
		}

		else // inter mode
		{
			//ori(RGB->YUV)->(YUV444->YUV420)->(inter prediction)->(predictor, recon)->(YUV420->RGB)->출력

			RGB_to_YUV(img_RGB, img_YUV444, WIDTH, HEIGHT);
			YUV444_to_420(img_YUV444, img_ori_Y, img_ori_U, img_ori_V, WIDTH, HEIGHT);

			InterPrediction(img_ori_Y, img_ref_Y, img_pred_Y, img_resi_Y, img_recon_Y, WIDTH, HEIGHT, BLOCK_SIZE, SR);
			InterPrediction(img_ori_U, img_ref_U, img_pred_U, img_resi_U, img_recon_U, cWIDTH, cHEIGHT, cBLOCK_SIZE, cSR);
			InterPrediction(img_ori_V, img_ref_V, img_pred_V, img_resi_V, img_recon_V, cWIDTH, cHEIGHT, cBLOCK_SIZE, cSR);

			YUV420_to_444(img_pred_Y, img_pred_U, img_pred_V, img_YUV444, WIDTH, HEIGHT);
			YUV_to_RGB(img_YUV444, img_RGB, WIDTH, HEIGHT);//predictor 변환
			Write_Frame(fp_out0, img_RGB, WIDTH, HEIGHT * 3);//RGB predictor 출력

			YUV420_to_444(img_recon_Y, img_recon_U, img_recon_V, img_YUV444, WIDTH, HEIGHT);
			YUV_to_RGB(img_YUV444, img_RGB, WIDTH, HEIGHT);//predictor 변환
			Write_Frame(fp_out1, img_RGB, WIDTH, HEIGHT * 3);//RGB recon 출력


			/*//Encoding Information//*/
			printf("P - Frame[POC %d] - PSNR [Y : %.2f dB, U : %.2f dB, V : %.2f dB] \n", frame, GetPSNR(img_ori_Y, img_recon_Y, WIDTH, HEIGHT), GetPSNR(img_ori_U, img_recon_U, cWIDTH, cHEIGHT), GetPSNR(img_ori_V, img_recon_V, cWIDTH, cHEIGHT));

			cpy_frame(img_recon_Y, img_ref_Y, WIDTH, HEIGHT);
			cpy_frame(img_recon_U, img_ref_U, cWIDTH, cHEIGHT);
			cpy_frame(img_recon_V, img_ref_V, cWIDTH, cHEIGHT);

			frame++;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////

	// mem free
	MemFree_2D(img_YUV444, HEIGHT * 3);
	MemFree_2D(img_RGB, HEIGHT * 3);

	MemFree_2D(img_ref_Y, HEIGHT);
	MemFree_2D(img_ref_U, cHEIGHT);
	MemFree_2D(img_ref_V, cHEIGHT);

	MemFree_2D(img_ori_Y, HEIGHT);
	MemFree_2D(img_ori_U, cHEIGHT);
	MemFree_2D(img_ori_V, cHEIGHT);

	MemFree_2D(img_pred_Y, HEIGHT);
	MemFree_2D(img_pred_U, cHEIGHT);
	MemFree_2D(img_pred_V, cHEIGHT);

	MemFree_D_2D(img_resi_Y, HEIGHT);
	MemFree_D_2D(img_resi_U, cHEIGHT);
	MemFree_D_2D(img_resi_V, cHEIGHT);

	MemFree_2D(img_recon_Y, HEIGHT);
	MemFree_2D(img_recon_U, cHEIGHT);
	MemFree_2D(img_recon_V, cHEIGHT);

	fcloseall();		//file close

	return 0;
}
