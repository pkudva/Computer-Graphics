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
		unsigned char color[3];

	// would some methods for transforming the triangle in place be helpful?

};

class Screen
{
	public:
		unsigned char   *buffer;
		int width, height;

	// would some methods for accessing and setting pixels be helpful?
	void SetPixel(int r, int c, unsigned char* col){
		if(c >= 0 and r >= 0 and c < width and r < height){
			int index = 3 * (r * width + c);
			buffer[index] = col[0];
			buffer[index + 1] = col[1];
			buffer[index + 2] = col[2];
		}
	}
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
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
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

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

void RasterizeGoingDownTri(Screen s, Triangle t){
	int upperY = std::distance(t.Y, std::max_element(t.Y, t.Y+3));
	int lowerY = std::distance(t.Y, std::min_element(t.Y, t.Y+3));

	double topY = t.Y[upperY];
	double botY = t.Y[lowerY];

	double rowMax = floor_441(topY);
	double rowMin = ceil_441(botY);


	for(int r = rowMin; r <= rowMax; r++){

		//modding the vertices in case lowerY is not 0
		double x1 = findEnd(t.X[lowerY], t.Y[lowerY], 
					t.X[(lowerY+1)%3], t.Y[(lowerY+1)%3], r);
		double x2 = findEnd(t.X[lowerY], t.Y[lowerY], 
					t.X[(lowerY+2)%3], t.Y[(lowerY+2)%3], r);

		//credit to He He for pointing out that x1 was not always the leftEnd
		double leftEnd = std::min(x1, x2);
		double rightEnd = std::max(x1, x2);

		for(int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++){
			s.SetPixel(r, c, t.color);
		}
	}
}

void RasterizeGoingUpTri(Screen s, Triangle t){
	int upperY = std::distance(t.Y, std::max_element(t.Y, t.Y+3));
	int lowerY = std::distance(t.Y, std::min_element(t.Y, t.Y+3));

	double topY = t.Y[upperY];
	double botY = t.Y[lowerY];

	double rowMax = floor_441(topY);
	double rowMin = ceil_441(botY);

	for(int r = rowMin; r <= rowMax; r++){

		//modding the vertices in case upperY is not 0
		double x1 = findEnd(t.X[upperY], t.Y[upperY], 
					t.X[(upperY+1)%3], t.Y[(upperY+1)%3], r);
		double x2 = findEnd(t.X[upperY], t.Y[upperY], 
					t.X[(upperY+2)%3], t.Y[(upperY+2)%3], r);

		//credit to He He for pointing out that x1 was not always the leftEnd
		double leftEnd = std::min(x1, x2);
		double rightEnd = std::max(x1, x2);

		for(int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++){
			s.SetPixel(r, c, t.color);
		}
	}
}

int main()
{
	vtkImageData *image = NewImage(1786, 1344);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = 1786*1344;
	for (int i = 0 ; i < npixels*3 ; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1786;
	screen.height = 1344;

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

			upT.Y[lowerY] = midY;
			upT.X[lowerY] = splitX;

			downT.Y[upperY] = midY;
			downT.X[upperY] = splitX;

			//Going down
			RasterizeGoingDownTri(screen, downT);

			//Going up
			RasterizeGoingUpTri(screen, upT);

		}
	}
	WriteImage(image, "allTriangles");
}
