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

	int i,j;
	int count;
	int padded_dimensions_x;
	int padded_dimensions_y;
	int pad_factor = 6;
	float sigma;
	float peak;
	float sum_of_peaks;
	ImageFile input_reconstruction_file;
	ImageFile input_search_image_file;
	Image input_image;
	Image input_reconstruction;

	AnglesAndShifts angles;
	EulerSearch global_euler_search;

	Image projection_filter;
	Image current_projection;
	Image template_reconstruction;
	Image average_density;
	Image padded_reference;

	wxString	my_symmetry = "C1";
	float	angular_step = 10.0f;
	float in_plane_angular_step = 5.0f;
	float psi_step;
	float psi_max = 360.0f;
	float psi_start = 0.0f;
	float pixel_size = 1.06f;
	float high_resolution_limit_search = 3.0f;
	float particle_radius_angstroms = 0.0f;
	float current_psi;

	float defocus1 = 10000.0f;
	float defocus2 = 10000.0f;
	float defocus_search_range = 0.0f;
	float defocus_step = 100.0f;
	float defocus_i = 0.0f;

	float variance;


	CTF input_ctf;

	float	voltage_kV = 300.0f;
	float	spherical_aberration_mm = 2.7f;
	float amplitude_contrast = 0.07f;
	float defocus_angle = 0.0f;
	float phase_shift = 0.0f;

	long pixel_counter = 0;

	ParameterMap parameter_map; // needed for euler search init
	//for (int i = 0; i < 5; i++) {parameter_map[i] = true;}
	parameter_map.SetAllTrue();
	int best_parameters_to_keep = 20;


	float mask_radius_search;
	if (particle_radius_angstroms < 1.0f) { mask_radius_search = 200.0f; } // This was the original default value.
	else mask_radius_search = particle_radius_angstroms;
	psi_step = in_plane_angular_step;

	// initialize search grid
	global_euler_search.InitGrid("C1", angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
	if (my_symmetry.StartsWith("C1")) // TODO 2x check me - w/o this O symm at least is broken
	{
		if (global_euler_search.test_mirror == true) // otherwise the theta max is set to 90.0 and test_mirror is set to true.  However, I don't want to have to test the mirrors.
		{
			global_euler_search.theta_max = 180.0f;
		}
	}
	global_euler_search.CalculateGridSearchPositions(false);
	int total_search_position = 0;
	int current_search_position = 0;
	int first_search_position = -1;
	int last_search_position = -1;
	// if running locally, search over all of them

	if (is_running_locally == true)
	{
		first_search_position = 0;
		last_search_position = global_euler_search.number_of_search_positions - 1;
	}

	input_search_image_file.OpenFile("147_Mar12_12.21.27_159_0.mrc", false);
	input_image.ReadSlice(&input_search_image_file, 1);
	input_reconstruction_file.OpenFile("6q8y_LSU_to_model5_bfactor_85_pix1.06.mrc", false);
	input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices());
	wxPrintf("template dim x = %i\n", input_reconstruction.logical_x_dimension);
	wxPrintf("template dim y = %i\n", input_reconstruction.logical_y_dimension);
	wxPrintf("template dim z = %i\n", input_reconstruction.logical_z_dimension);
	current_projection.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	projection_filter.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	padded_reference.Allocate(5832, 4096, true);
	padded_reference.SetToConstant(0.0f);
	template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
	average_density.Allocate(5832, 4096, true);
	average_density.SetToConstant(0.0f);

	input_reconstruction.ChangePixelSize(&template_reconstruction, 1.0f, 0.001f, true);
	template_reconstruction.ZeroCentralPixel();
	template_reconstruction.SwapRealSpaceQuadrants();

	input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
	input_ctf.SetDefocus((defocus1 + float(defocus_i) * defocus_step) / pixel_size, (defocus2 + float(defocus_i) * defocus_step) / pixel_size, deg_2_rad(defocus_angle));

	projection_filter.CalculateCTFImage(input_ctf);
	for (current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++)
	{
		//loop over each rotation angle
		for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
		{
			angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], current_psi, 0.0, 0.0);

			template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
			current_projection.SwapRealSpaceQuadrants();
			current_projection.MultiplyPixelWise(projection_filter);
			current_projection.BackwardFFT();
			current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges());
			variance = current_projection.ReturnSumOfSquares() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels \
					- powf(current_projection.ReturnAverageOfRealValues() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
			current_projection.DivideByConstant(sqrtf(variance));
			current_projection.ClipIntoLargerRealSpace2D(&padded_reference);
			average_density.AddImage(&padded_reference);
			total_search_position++;

			/*
			if (total_search_position == 1)
			{
				if (average_density.is_in_real_space) wxPrintf("avg in real space\n");
				average_density.DivideByConstant(2.0f);
				average_density.QuickAndDirtyWriteSlice("total.mrc",1);
				exit(0);
			}
			*/
		}
	}
	average_density.DivideByConstant(total_search_position);
	average_density.QuickAndDirtyWriteSlice("average_density.mrc",1);

	wxPrintf("total search positions = %i\n", total_search_position);


	return true;
}
