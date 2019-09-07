#include <stdio.h>

void minmax(double *X, int N, double *min, double *max) {
	double Min = 1e10;
	double Max = -1e10;
	for (int i = 0; i < N; i++) {
		if (Min > X[i]) Min = X[i];
		if (Max < X[i]) Max = X[i];
	}
	*min = Min;
	*max = Max;
}

int getlabel(double *X, double *Y, int N, int i0, int i1, double xmin, double xmax, double ymin, double ymax) {
	int label = 3;

	if (X[i0] == xmin && X[i1] == xmin) label = 4;
	if (X[i0] == xmax && X[i1] == xmax) label = 2;
	if (Y[i0] == ymin && Y[i1] == ymin) label = 1;

	return label;
}

int main(int argc, char** argv)
{
	char filename[] = "mesh.mesh";
	char input[300];
	
	double X[10000];
	double Y[10000];
	float stepsize = -1;

	if(argc < 2) {printf("Usage: resonatorGenerator.exe [input (csv-format)] [optionally: step size]\n"); return 0;}

	sscanf(argv[1], "%s", input);
	if (argc == 3)
		sscanf(argv[2], "%g", &stepsize);

	printf("Border Inpur generation \nOpening %s file\n", filename);

	FILE* inp = fopen(input, "r");
	if (inp == nullptr) {		printf("Error opening %s\n", filename); return 1;	}

	float x, y;

	int ind = 0;
	while (fscanf(inp, "%g,%g", &x, &y) == 2) {
		X[ind] = x;
		Y[ind] = y;
		ind++;
	}
	ind--; //The shape in input is closed: i.e. first and last vertices are the same
	       //this will cause an error while mesh generation. Exclude one (last) vertex

	printf("Generating `mesh.mesh` border file...\n");
	double xmin, xmax, ymin, ymax;
	minmax(Y, ind, &ymin, &ymax);
	minmax(X, ind, &xmin, &xmax);

	if (stepsize == -1) {		
		stepsize = float((ymax - ymin) / 50.);
		printf("Mesh step is not specified, setting (Ymax - Ymin)/50 = %g\n", stepsize);
	}
	else
		printf("Mesh step is specified. %g", stepsize);

	int vert = ind;

	FILE *output = fopen(filename, "w");
	fprintf(output, "MeshVersionFormatted 0\nAngleOfCornerBound  45\nDimension  2\nVertices   %i\n", vert);
	for(int j = 0; j < vert; j++)
		fprintf(output, "%g\t%g\t1\n", X[j], Y[j]);

	fprintf(output, "Edges	%i\n", vert);

	for (int j = 0; j < vert - 1; j++) {
		fprintf(output, "%i\t%i\t%i\n", j + 1, j + 1 + 1, getlabel(X,Y,vert, j, j+1, xmin, xmax, ymin, ymax));
	}

	fprintf(output, "%i\t%i\t%i\n", vert, 1, getlabel(X, Y, vert, vert-1, 0, xmin, xmax, ymin, ymax));


	fprintf(output, "hVertices\n");
	for(int j =0; j < vert; j++)
		fprintf(output, "%g\n", stepsize);

	fclose(output);

	printf("                                     Done\n");
	return 0;	
}