#include "../../core/core_headers.h"

class
NikoTestApp : public MyApp
{

	public:

	bool DoCalculation();
	void DoInteractiveUserInput();

	private:
};


IMPLEMENT_APP(NikoTestApp)

// override the DoInteractiveUserInput

void NikoTestApp::DoInteractiveUserInput()
{

}

// override the do calculation method which will be what is actually run..

bool NikoTestApp::DoCalculation()
{
	ImageFile input_stack_file;
	input_stack_file.OpenFile("Mpneumoniae_stack.mrc", false);
	Image output_image;
	Image input_image;
	//input_image.ReadSlices(&input_stack_file, 1, input_stack_file.ReturnNumberOfSlices());
	input_image.Allocate(input_stack_file.ReturnXSize(), input_stack_file.ReturnYSize(), true);
	output_image.Allocate(input_stack_file.ReturnXSize(), input_stack_file.ReturnYSize(), true);
	output_image.SetToConstant(0.0f);

	wxPrintf("number of slices = %i\n",input_stack_file.ReturnZSize());
	int i;
	for(i=1; i<= input_stack_file.ReturnZSize();i++){
		input_image.ReadSlice(&input_stack_file, i);
		output_image.AddImage(&input_image);
	}
	output_image.DivideByConstant(i);
	output_image.QuickAndDirtyWriteSlice("average_of_stack.mrc",1);

	return true;
}
