#include "../../core/core_headers.h"

class
ExtractParticlesApp : public MyApp
{

	public:

	bool DoCalculation();
	void DoInteractiveUserInput();


	private:

};

IMPLEMENT_APP(ExtractParticlesApp)



void ExtractParticlesApp::DoInteractiveUserInput()
{

	UserInput *my_input = new UserInput("ExtractParticles", 0.0);

	wxString	micrograph_filename					=	my_input->GetFilenameFromUser("Input micrograph filename","The input micrograph, in which we will look for particles","micrograph.mrc",true);
	wxString	coordinates_filename				=	my_input->GetFilenameFromUser("Coordinates (PLT) filename","The input particle coordinates, in Imagic-style PLT forlmat","coos.plt",true);
	wxString	output_stack_filename				=	my_input->GetFilenameFromUser("Filename for output stack of particles.","A stack of particles will be written to disk","particles.mrc",false);
	int			output_stack_box_size				=	my_input->GetIntFromUser("Box size for output candidate particle images (pixels)","In pixels. Give 0 to skip writing particle images to disk.","256",0);
	float       pixel_size 							=   my_input->GetFloatFromUser("pixel size","pixel size","1.0",0.1,100.0);


	delete my_input;

	my_current_job.Reset(5);
	my_current_job.ManualSetArguments("tttif",			micrograph_filename.ToStdString().c_str(),
														coordinates_filename.ToStdString().c_str(),
														output_stack_filename.ToStdString().c_str(),
														output_stack_box_size,
														pixel_size
														);


}

// override the do calculation method which will be what is actually run..

bool ExtractParticlesApp::DoCalculation()
{

	ProgressBar *my_progress_bar;
	EmpiricalDistribution my_dist;

	// Get the arguments for this job..
	wxString 	micrograph_filename 						= 	my_current_job.arguments[0].ReturnStringArgument();
	wxString 	coordinates_filename 						= 	my_current_job.arguments[1].ReturnStringArgument();
	wxString	output_stack_filename						=	my_current_job.arguments[2].ReturnStringArgument();
	int			output_stack_box_size						=	my_current_job.arguments[3].ReturnIntegerArgument();
	float		pixel_size						            =	my_current_job.arguments[4].ReturnFloatArgument();



	// Open input files so we know dimensions
	MRCFile micrograph_file(micrograph_filename.ToStdString(),false);
	MyDebugAssertTrue(micrograph_file.ReturnNumberOfSlices() == 1,"Input micrograph file should only contain one image for now");


	// Let's box particles out
	Image micrograph;
	Image micrograph_copy;
	Image box;
	box.Allocate(output_stack_box_size,output_stack_box_size,1,true);
	micrograph.ReadSlice(&micrograph_file,1);
	micrograph_copy.Allocate(micrograph.logical_x_dimension, micrograph.logical_y_dimension, 1);
	float micrograph_mean = micrograph.ReturnAverageOfRealValues();
	NumericTextFile *input_coos_file;
	input_coos_file = new NumericTextFile(coordinates_filename, OPEN_TO_READ, 3);
	int number_of_particles = input_coos_file->number_of_lines;
	MRCFile output_stack;
	output_stack.OpenFile(output_stack_filename.ToStdString(),true);
	my_progress_bar = new ProgressBar(number_of_particles);
	float plt_x, plt_y;
	int my_x, my_y;
	float temp_array[8];
	for ( int counter = 0; counter < number_of_particles; counter ++ )
	{
		input_coos_file->ReadLine(temp_array);
		plt_x = roundf(temp_array[3]/pixel_size);  // convert A to pixels in coordinate file
		//wxPrintf("pixel size = %f\n", pixel_size);
		wxPrintf("old and new plt_x = %f, %f\n", temp_array[3], plt_x);

		plt_y = roundf(temp_array[4]/pixel_size);  // convert A to pixels in coordinate file
		wxPrintf("old and new plt_y = %f, %f\n", temp_array[4], plt_y);
		/*
		 * plt_x = (my_y + phys_addr_box_center_y) + 1.0;
		 * plt_y = (logical_x_dim - (phys_addr_box_center_x + my_x)) + 1.0;
		 *
		 * my_x = -1 * ((plt_y - 1.0)  - logical_x_dim + phys_addr_box_center_x)
		 * my_y = (plt_x - 1.0)  - phys_addr_box_center_y;
		 */
		//my_x = (-1.0) * ((plt_y - 1.0) - micrograph.logical_x_dimension + micrograph.physical_address_of_box_center_x);
		//my_y = (plt_x - 1.0) - micrograph.physical_address_of_box_center_y;
		my_x = int(plt_x-micrograph.physical_address_of_box_center_x);
		my_y = int(plt_y-micrograph.physical_address_of_box_center_y);
		//wxPrintf("my_x = %i\n", my_x);
		//wxPrintf("my_y = %i\n", my_y);
		micrograph_copy.CopyFrom(&micrograph);
		micrograph_copy.RealSpaceIntegerShift(plt_x-micrograph.physical_address_of_box_center_x, plt_y-micrograph.physical_address_of_box_center_y);
		micrograph_copy.ClipInto(&box);
		//micrograph.ClipInto(&box,micrograph_mean,false,1.0,-int(plt_x-micrograph.physical_address_of_box_center_x),-int(plt_y-micrograph.physical_address_of_box_center_y),0);
		box.WriteSlice(&output_stack , counter + 1);

		//box.QuickAndDirtyWriteSlice(wxString::Format("/groups/kexin/research/scaling_mip/normality_test/real_image_147/windowed_peak_%i.mrc", counter).ToStdString(),1);

		my_progress_bar->Update(counter+1);
	}
	delete my_progress_bar;
	wxPrintf("\nExtracted %i particles\n",number_of_particles);
	delete input_coos_file;

	return true;
}
