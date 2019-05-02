#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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
			int index = 3 * (r * 1000 + c);
			buffer[index] = col[0];
			buffer[index + 1] = col[1];
			buffer[index + 2] = col[2];
		}
	}
};

std::vector<Triangle>
GetTriangles(void)
{
	std::vector<Triangle> rv(100);

	unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
								{76,153,0}, {255, 204, 204}, {204, 204, 0}};
	for (int i = 0 ; i < 100 ; i++)
	{
		int idxI = i%10;
		int posI = idxI*100;
		int idxJ = i/10;
		int posJ = idxJ*100;
		int firstPt = (i%3);
		rv[i].X[firstPt] = posI;
		if (i == 50)
		   rv[i].X[firstPt] = -10;
		rv[i].Y[firstPt] = posJ+10*(idxJ+1);
		rv[i].X[(firstPt+1)%3] = posI+99;
		rv[i].Y[(firstPt+1)%3] = posJ+10*(idxJ+1);
		rv[i].X[(firstPt+2)%3] = posI+i;
		rv[i].Y[(firstPt+2)%3] = posJ;
		if (i == 5)
			rv[i].Y[(firstPt+2)%3] = -50;
		rv[i].color[0] = colors[i%6][0];
		rv[i].color[1] = colors[i%6][1];
		rv[i].color[2] = colors[i%6][2];
	}

	return rv;
}

double findEnd(double X1, double Y1, double X2, double Y2, int row){
	//find ends of triangle
	//(X1,Y1) should always be lower vertex

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
		double leftEnd = ceil_441(findEnd(t.X[lowerY], t.Y[lowerY], 
									t.X[(lowerY+1)%3], t.Y[(lowerY+1)%3], r));
		double rightEnd = floor_441(findEnd(t.X[lowerY], t.Y[lowerY], 
									t.X[(lowerY+2)%3], t.Y[(lowerY+2)%3], r));

		for(int c = leftEnd; c <= rightEnd; c++){
			//check if in screen boundaries
			s.SetPixel(r, c, t.color);
		}
	}
}

int main()
{
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1000;
	screen.height = 1000;

	// YOUR CODE GOES HERE TO DEPOSIT THE COLORS FROM TRIANGLES 
	// INTO PIXELS USING THE SCANLINE ALGORITHM

	int rowMax = 0;
	int rowMin = 0;
	int leftEnd = 0;
	int rightEnd = 0;
	int index = 0;
	double topY = 0;
	double botY = 0;

	int lowerY = 0;
	int upperY = 0;

	for(int j = 0; j < 100; j++){
		RasterizeGoingDownTri(screen, triangles[j]);
	}
	WriteImage(image, "allTriangles");
}
