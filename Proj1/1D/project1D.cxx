#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

class Triangle
{
	public:
		double         X[3];
		double         Y[3];
		double         Z[3];
		double		   colors[3][3];


	// would some methods for transforming the triangle in place be helpful?

};

class Screen
{
	public:
		unsigned char   *buffer;
		double			*zbuffer;
		int width, height;

	// would some methods for accessing and setting pixels be helpful?
	void SetPixel(int r, int c, double *col, double z){
		if(c >= 0 and r >= 0 and c < width and r < height){
			int index = 3 * (r * width + c);
			int zindex = (r * width + c);
			if(z > zbuffer[zindex]){
				buffer[index] = ceil_441(col[0] * 255);
				buffer[index + 1] = ceil_441(col[1] * 255);
				buffer[index + 2] = ceil_441(col[2] * 255);
				zbuffer[zindex] = z;
			}
		}
	}
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double findEnd(double X1, double Y1, double X2, double Y2, double row){

	//find ends of triangle

	if(X1 == X2){
		return X1;
	}

	double m = (Y2 - Y1)/(X2 - X1);
	double b = Y1 - m * X1;
	double x = (row - b)/m;

	return x;
}

double* lerpColors(double X1, double Y1, double X2, double Y2, double* F1, double* F2, double X3, double Y3){
	double t = 0;
	double* F3;

	F3 = (double*) malloc(3*sizeof(double));


	if(X1 != X2){
		t = (X3 - X1)/(X2 - X1);
	}
	else if(Y1 != Y2){
		t = (Y3 - Y1)/(Y2 - Y1);
	}

	F3[0] = F1[0] + t*(F2[0] - F1[0]);
	F3[1] = F1[1] + t*(F2[1] - F1[1]);
	F3[2] = F1[2] + t*(F2[2] - F1[2]);
	
	return F3;
}

double lerpZ(double X1, double Y1, double X2, double Y2, double Z1, double Z2, double X3, double Y3){
	double t = 0;
	double Z3 = -1;

	if(X1 != X2){
		t = (X3 - X1)/(X2 - X1);
	}
	else if(Y1 != Y2){
		t = (Y3 - Y1)/(Y2 - Y1);
	}

	Z3 = Z1 + t*(Z2 - Z1);
	
	return Z3;
}

void RasterizeGoingDownTri(Screen s, Triangle t){

	double* leftCol;
	double* rightCol;
	double* pixelCol;
	double leftZ;
	double rightZ;
	double pixelZ;

	int upperY = std::distance(t.Y, std::max_element(t.Y, t.Y+3));
	int lowerY = std::distance(t.Y, std::min_element(t.Y, t.Y+3));

	double rowMax = t.Y[upperY];
	double rowMin = t.Y[lowerY];

	for(int r = ceil_441(rowMin); r <= floor_441(rowMax); r++){

		//modding the vertices in case lowerY is not 0
		double x1 = findEnd(t.X[lowerY], t.Y[lowerY], 
					t.X[(lowerY+1)%3], t.Y[(lowerY+1)%3], r);
		double x2 = findEnd(t.X[lowerY], t.Y[lowerY], 
					t.X[(lowerY+2)%3], t.Y[(lowerY+2)%3], r);

		//credit to He He for pointing out that x1 was not always the leftEnd
		double leftEnd = std::min(x1, x2);
		double rightEnd = std::max(x1, x2);

		//color interpolation and z interpolation
		if(leftEnd == x1 and rightEnd == x2){
			//lowerY and lowerY+1 vertices for left
			leftCol = lerpColors(t.X[lowerY], rowMin, 
								 t.X[(lowerY+1)%3], rowMax, 
								 t.colors[lowerY], t.colors[(lowerY+1)%3], 
								 leftEnd, r);
			leftZ = lerpZ(t.X[lowerY], rowMin, 
						  t.X[(lowerY+1)%3], rowMax, 
						  t.Z[lowerY], t.Z[(lowerY+1)%3], 
						  leftEnd, r);

			//lowerY and lowerY+2 vertices for right
			rightCol = lerpColors(t.X[lowerY], rowMin, 
								  t.X[(lowerY+2)%3], rowMax, 
								  t.colors[lowerY], t.colors[(lowerY+2)%3], 
								  rightEnd, r);
			rightZ = lerpZ(t.X[lowerY], rowMin, 
						   t.X[(lowerY+2)%3], rowMax, 
						   t.Z[lowerY], t.Z[(lowerY+2)%3], 
						   rightEnd, r);
		}
		else{
			//lowerY and lowerY+2 vertices for left
			leftCol = lerpColors(t.X[lowerY], rowMin, 
								 t.X[(lowerY+2)%3], rowMax, 
								 t.colors[lowerY], t.colors[(lowerY+2)%3], 
								 leftEnd, r);
			leftZ = lerpZ(t.X[lowerY], rowMin, 
						  t.X[(lowerY+2)%3], rowMax, 
						  t.Z[lowerY], t.Z[(lowerY+2)%3], 
						  leftEnd, r);

			//lowerY and lowerY+1 vertices for right
			rightCol = lerpColors(t.X[lowerY], rowMin, 
								  t.X[(lowerY+1)%3], rowMax, 
								  t.colors[lowerY], t.colors[(lowerY+1)%3], 
								  rightEnd, r);
			rightZ = lerpZ(t.X[lowerY], rowMin, 
			  			  t.X[(lowerY+1)%3], rowMax, 
						  t.Z[lowerY], t.Z[(lowerY+1)%3], 
						  rightEnd, r);
		}

		for(int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++){
			pixelCol = lerpColors(leftEnd, r, rightEnd, r, leftCol, rightCol, c, r);
			pixelZ = lerpZ(leftEnd, r, rightEnd, r, leftZ, rightZ, c, r);
			s.SetPixel(r, c, pixelCol, pixelZ);
		}
	}
}



void RasterizeGoingUpTri(Screen s, Triangle t){

	double* leftCol;
	double* rightCol;
	double* pixelCol;
	double leftZ;
	double rightZ;
	double pixelZ;
	
	int upperY = std::distance(t.Y, std::max_element(t.Y, t.Y+3));
	int lowerY = std::distance(t.Y, std::min_element(t.Y, t.Y+3));

	double rowMax = t.Y[upperY];
	double rowMin = t.Y[lowerY];

	for(int r = ceil_441(rowMin); r <= floor_441(rowMax); r++){

		//modding the vertices in case upperY is not 0
		double x1 = findEnd(t.X[upperY], t.Y[upperY], 
					t.X[(upperY+1)%3], t.Y[(upperY+1)%3], r);
		double x2 = findEnd(t.X[upperY], t.Y[upperY], 
					t.X[(upperY+2)%3], t.Y[(upperY+2)%3], r);

		//credit to He He for pointing out that x1 was not always the leftEnd
		double leftEnd = std::min(x1, x2);
		double rightEnd = std::max(x1, x2);


		//color interpolation
		if(leftEnd == x1){
			//upperY and upperY+1 vertices for left
			leftCol = lerpColors(t.X[upperY], t.Y[upperY], 
								 t.X[(upperY+1)%3], t.Y[(upperY+1)%3], 
								 t.colors[upperY], t.colors[(upperY+1)%3], 
								 leftEnd, r);
			leftZ = lerpZ(t.X[upperY], t.Y[upperY], 
						  t.X[(upperY+1)%3], t.Y[(upperY+1)%3], 
						  t.Z[upperY], t.Z[(upperY+1)%3], 
						  leftEnd, r);
			//upperY and upperY+2 vertices for right
			rightCol = lerpColors(t.X[upperY], t.Y[upperY], 
								  t.X[(upperY+2)%3], t.Y[(upperY+2)%3], 
								  t.colors[upperY], t.colors[(upperY+2)%3], 
								  rightEnd, r);
			rightZ = lerpZ(t.X[upperY], t.Y[upperY], 
						  t.X[(upperY+2)%3], t.Y[(upperY+2)%3], 
						  t.Z[upperY], t.Z[(upperY+2)%3], 
						  rightEnd, r);
		}
		else{
			//upperY and upperY+2 vertices for left
			leftCol = lerpColors(t.X[upperY], t.Y[upperY], 
								 t.X[(upperY+2)%3], t.Y[(upperY+2)%3], 
								 t.colors[upperY], t.colors[(upperY+2)%3], 
								 leftEnd, r);
			leftZ = lerpZ(t.X[upperY], t.Y[upperY], 
						  t.X[(upperY+2)%3], t.Y[(upperY+2)%3], 
						  t.Z[upperY], t.Z[(upperY+2)%3], 
						  leftEnd, r);
			//upperY and upperY+1 vertices for right
			rightCol = lerpColors(t.X[upperY], t.Y[upperY], 
								  t.X[(upperY+1)%3], t.Y[(upperY+1)%3], 
								  t.colors[upperY], t.colors[(upperY+1)%3], 
								  rightEnd, r);
			rightZ = lerpZ(t.X[upperY], t.Y[upperY], 
						  t.X[(upperY+1)%3], t.Y[(upperY+1)%3], 
						  t.Z[upperY], t.Z[(upperY+1)%3], 
						  rightEnd, r);
		}

		for(int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++){
			pixelCol = lerpColors(leftEnd, r, rightEnd, r, leftCol, rightCol, c, r);
			pixelZ = lerpZ(leftEnd, r, rightEnd, r, leftZ, rightZ, c, r);
			s.SetPixel(r, c, pixelCol, pixelZ);
		}
	}
}

int main()
{
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
	double* zbuffer = (double *) malloc(1000*1000*sizeof(double));
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++)
		buffer[i] = 0;

	for (int i = 0 ; i < npixels ; i++)
		zbuffer[i] = -1;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.zbuffer = zbuffer;
	screen.width = 1000;
	screen.height = 1000;

	// YOUR CODE GOES HERE TO DEPOSIT THE COLORS FROM TRIANGLES 
	// INTO PIXELS USING THE SCANLINE ALGORITHM

	int lowerY = 0;
	int upperY = 0;
	int centerY = 0;

	double maxY = 0;
	double minY = 0;
	double midY = 0;

	for(int j = 0; j < triangles.size(); j++){

		upperY = std::distance(triangles[j].Y, std::max_element(triangles[j].Y, triangles[j].Y+3));
		lowerY = std::distance(triangles[j].Y, std::min_element(triangles[j].Y, triangles[j].Y+3));
		if(upperY == 0 and lowerY == 1)
			centerY = 2;
		else if(lowerY == 0 and upperY == 1)
			centerY = 2;
		else if(upperY == 1 and lowerY == 2)
			centerY = 0;
		else if(lowerY == 1 and upperY == 2)
			centerY = 0;
		else if(upperY == 2 and lowerY == 0)
			centerY = 1;
		else if(lowerY == 2 and upperY == 0)
			centerY = 1;

		maxY = triangles[j].Y[upperY];
		minY = triangles[j].Y[lowerY];
		midY = triangles[j].Y[centerY];

		if(midY == maxY){
			//Going down triangle

			RasterizeGoingDownTri(screen, triangles[j]);
		}
		else if(midY == minY){
			//Going up triangle

			RasterizeGoingUpTri(screen, triangles[j]);
		}
		else{
			//mid is unique, split triangle
			Triangle upT = triangles[j];
			Triangle downT = triangles[j];

			
			double splitX = findEnd(triangles[j].X[upperY], triangles[j].Y[upperY], 
							triangles[j].X[lowerY], triangles[j].Y[lowerY], midY);

			double splitZ = findEnd(triangles[j].Z[upperY], triangles[j].Y[upperY], 
							triangles[j].Z[lowerY], triangles[j].Y[lowerY], midY);

			double* splitColors = lerpColors(triangles[j].X[upperY], triangles[j].Y[upperY], 
											 triangles[j].X[lowerY], triangles[j].Y[lowerY], 
								 			 triangles[j].colors[upperY], triangles[j].colors[lowerY], 
											 splitX, midY);

			upT.Z[lowerY] = splitZ;
			upT.Y[lowerY] = midY;
			upT.X[lowerY] = splitX;
			upT.colors[lowerY][0] = splitColors[0];
			upT.colors[lowerY][1] = splitColors[1];
			upT.colors[lowerY][2] = splitColors[2];

			downT.Z[upperY] = splitZ;
			downT.Y[upperY] = midY;
			downT.X[upperY] = splitX;
			downT.colors[upperY][0] = splitColors[0];
			downT.colors[upperY][1] = splitColors[1];
			downT.colors[upperY][2] = splitColors[2];

			//Going down
			RasterizeGoingDownTri(screen, downT);

			//Going up
			RasterizeGoingUpTri(screen, upT);

		}
	}
	WriteImage(image, "allTriangles");
}
