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
#include <vtkDoubleArray.h>

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


class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform();
    Matrix          CameraTransform();
    Matrix          DeviceTransform(int width, int height);
};

Matrix
Camera::ViewTransform(){
	Matrix m;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			m.A[i][j] = 0;
		}
	}

	m.A[0][0] = 1/tan(angle/2);
	m.A[1][1] = 1/tan(angle/2);
	m.A[2][2] = (far + near)/(far - near);
	m.A[2][3] = -1;
	m.A[3][2] = (2*far*near)/(far - near);

	return m;
}

Matrix
Camera::CameraTransform(){
	double origin[3] = {0,0,0}; //position
	double w[3] = {0,0,0}; //origin - focus
	double u[3] = {0,0,0}; //up cross w
	double v[3] = {0,0,0}; //w cross u;

	Matrix m;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			m.A[i][j] = 0;
		}
	}

	//creating the camera frame
	origin[0] = position[0];
	origin[1] = position[1];
	origin[2] = position[2];

	double t[3] = {0 - origin[0], 0 - origin[1], 0 - origin[2]};


	double wX = origin[0] - focus[0];
	double wY = origin[1] - focus[1];
	double wZ = origin[2] - focus[2];
	double magW = sqrt(wX*wX + wY*wY + wZ*wZ);

	w[0] = wX/magW;
	w[1] = wY/magW;
	w[2] = wZ/magW;

	double uX = up[1]*w[2] - up[2]*w[1];
	double uY = up[2]*w[0] - up[0]*w[2];
	double uZ = up[0]*w[1] - up[1]*w[0];
	double magU = sqrt(uX*uX + uY*uY + uZ*uZ);

	u[0] = uX/magU;
	u[1] = uY/magU;
	u[2] = uZ/magU;

	double vX = w[1]*u[2] - w[2]*u[1];
	double vY = w[2]*u[0] - w[0]*u[2];
	double vZ = w[0]*u[1] - w[1]*u[0];
	double magV = sqrt(vX*vX + vY*vY + vZ*vZ);

	v[0] = vX/magV;
	v[1] = vY/magV;
	v[2] = vZ/magV;

	m.A[0][0] = u[0];
	m.A[0][1] = v[0];
	m.A[0][2] = w[0];
	m.A[1][0] = u[1];
	m.A[1][1] = v[1];
	m.A[1][2] = w[1];
	m.A[2][0] = u[2];
	m.A[2][1] = v[2];
	m.A[2][2] = w[2];
	m.A[3][0] = u[0]*t[0] + u[1]*t[1] + u[2]*t[2];
	m.A[3][1] = v[0]*t[0] + v[1]*t[1] + v[2]*t[2];
	m.A[3][2] = w[0]*t[0] + w[1]*t[1] + w[2]*t[2];

	m.A[3][3] = 1;

	return m;
}

Matrix
Camera::DeviceTransform(int width, int height){
	Matrix m;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			m.A[i][j] = 0;
		}
	}

	m.A[0][0] = width/2;
	m.A[1][1] = height/2;
	m.A[2][2] = 1;
	m.A[3][3] = 1;
	m.A[3][0] = width/2;
	m.A[3][1] = height/2;

	return m;
}




double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
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

Triangle
TransformTriangles(Matrix m, Triangle t){
	Triangle newT = t;
	double pointsIn[4] = {0,0,0,0};
	double pointsOut[4] = {0,0,0,0};
	for(int i = 0; i < 3; i++){
		pointsIn[0] = t.X[i];
		pointsIn[1] = t.Y[i];
		pointsIn[2] = t.Z[i];
		pointsIn[3] = 1;
		m.TransformPoint(pointsIn, pointsOut);
		newT.X[i] = pointsOut[0]/pointsOut[3];
		newT.Y[i] = pointsOut[1]/pointsOut[3];
		newT.Z[i] = pointsOut[2]/pointsOut[3];
	}

	return newT;
}

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
    rdr->SetFileName("proj1e_geometry.vtk");
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
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
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
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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

void SplitTriangles(Screen s, Triangle t){
	int lowerY = 0;
	int upperY = 0;
	int centerY = 0;

	double maxY = 0;
	double minY = 0;
	double midY = 0;

	upperY = std::distance(t.Y, std::max_element(t.Y, t.Y+3));
	lowerY = std::distance(t.Y, std::min_element(t.Y, t.Y+3));
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

	maxY = t.Y[upperY];
	minY = t.Y[lowerY];
	midY = t.Y[centerY];

	if(midY == maxY){
		//Going down triangle

		RasterizeGoingDownTri(s, t);
	}
	else if(midY == minY){
		//Going up triangle

		RasterizeGoingUpTri(s, t);
	}
	else{
		//mid is unique, split triangle
		Triangle upT = t;
		Triangle downT = t;

		
		double splitX = findEnd(t.X[upperY], t.Y[upperY], 
						t.X[lowerY], t.Y[lowerY], midY);

		double splitZ = findEnd(t.Z[upperY], t.Y[upperY], 
						t.Z[lowerY], t.Y[lowerY], midY);

		double* splitColors = lerpColors(t.X[upperY], t.Y[upperY], 
										 t.X[lowerY], t.Y[lowerY], 
							 			 t.colors[upperY], t.colors[lowerY], 
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
		RasterizeGoingDownTri(s, downT);

		//Going up
		RasterizeGoingUpTri(s, upT);
	}
}

int main()
{
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
	double* zbuffer = (double *) malloc(1000*1000*sizeof(double));
	int npixels = 1000*1000;

	Screen screen;
	screen.buffer = buffer;
	screen.zbuffer = zbuffer;
	screen.width = 1000;
	screen.height = 1000;

	// YOUR CODE GOES HERE TO DEPOSIT THE COLORS FROM TRIANGLES 
	// INTO PIXELS USING THE SCANLINE ALGORITHM

	Matrix CT;
	Matrix VT;
	Matrix DT;
	Matrix M;
	Triangle newTriangle;

	for(int k = 0; k < 4; k++){
		std::vector<Triangle> triangles = GetTriangles();

		for (int i = 0 ; i < npixels*3 ; i++)
			buffer[i] = 0;
		for (int i = 0 ; i < npixels ; i++)
			zbuffer[i] = -1;

		int f = k*250;
		Camera c = GetCamera(f, 1000);

		//Compose matrices
		CT = c.CameraTransform();
		VT = c.ViewTransform();
		DT = c.DeviceTransform(screen.width, screen.height);

		M = Matrix::ComposeMatrices(Matrix::ComposeMatrices(CT, VT), DT);

		//transform triangles to device space

		for(int j = 0; j < triangles.size(); j++){
			//printf("Triangle %d\n", j);
			newTriangle = TransformTriangles(M, triangles[j]);
			SplitTriangles(screen, newTriangle);
		}
		char temp[32];
		sprintf(temp, "frame%03d", f);
		WriteImage(image, temp);
	}
	
}
