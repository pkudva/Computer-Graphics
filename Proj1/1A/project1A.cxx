#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

vtkImageData * NewImage(int width, int height)
{
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}


int main()
{
	int width = 1024;
	int height = 1350;

	std::cerr << "In main!" << endl;
	vtkImageData *image = NewImage(width, height);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);

	int index = 0;

	unsigned char red = 0;
	unsigned char green = 0;
	unsigned char blue = 0;

	for(int x = 0; x < 27; x++){ //27 horizontal strips

		//make the color
		if(x%3 == 0){
			blue = 0;
		}
		else if (x%3 == 1){
			blue = 128;
		}
		else{
			blue = 255;
		}

		if((x/3)%3 == 0){
			green = 0;
		}
		else if((x/3)%3 == 1){
			green = 128;
		}
		else{
			green = 255;
		}

		if(x/9 == 0){
			red = 0;
		}
		else if (x/9 == 1){
			red = 128;
		}
		else if (x/9 == 2){
			red = 255;
		}

		//each strip has 50 rows of that color
		for(int r = 0; r < 50; r++){ 
			for(int c = 0; c < width; c++){
				index = 3 * ((r + x * 50) * width + c);
				buffer[index] = red;
				buffer[index+1] = green;
				buffer[index+2] = blue;
			}
		}
	}

	
	WriteImage(image, "proj1A");
}




