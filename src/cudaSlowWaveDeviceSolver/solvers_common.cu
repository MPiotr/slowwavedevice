extern __device__ void BiAverage(double *A, double *B, int p0, int datasize, int logsize)
{
	int stride = datasize;
	for (int q = 1; q <= logsize; q++)
	{
		stride = stride >> 1;
		if (p0 < stride)
		{
			A[p0] += A[p0 + stride];

		}
		else
		{
			if (p0 < 2 * stride)
			{

				B[p0 - stride] += B[p0];

			}
		}
		__syncthreads();
	}

}