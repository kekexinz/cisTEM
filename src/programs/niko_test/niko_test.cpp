#include "../../core/core_headers.h"

class
NikoTestApp : public MyApp
{

	public:

	bool DoCalculation();
	void DoInteractiveUserInput();

	private:
};

class TemplateComparisonObject
{
public:
	Image						*input_reconstruction, *windowed_particle, *projection_filter;
	AnglesAndShifts				*angles;
	float						pixel_size_factor;
//	int							slice = 1;
};

Peak TemplateScore(void *scoring_parameters, float peak_location_x, float peak_location_y)
{
	TemplateComparisonObject *comparison_object = reinterpret_cast < TemplateComparisonObject *> (scoring_parameters);
	Image current_projection;
//	Peak box_peak;

	current_projection.Allocate(comparison_object->projection_filter->logical_x_dimension, comparison_object->projection_filter->logical_x_dimension, false);

	if (comparison_object->input_reconstruction->logical_x_dimension != current_projection.logical_x_dimension)
	{
		Image padded_projection;
		padded_projection.Allocate(comparison_object->input_reconstruction->logical_x_dimension, comparison_object->input_reconstruction->logical_x_dimension, false);
		comparison_object->input_reconstruction->ExtractSlice(padded_projection, *comparison_object->angles, 1.0f, false);
		padded_projection.SwapRealSpaceQuadrants();
		padded_projection.BackwardFFT();
		padded_projection.ChangePixelSize(&current_projection, comparison_object->pixel_size_factor, 0.001f, true);
//		padded_projection.ChangePixelSize(&padded_projection, comparison_object->pixel_size_factor, 0.001f);
//		padded_projection.ClipInto(&current_projection);
//		current_projection.ForwardFFT();
	}
	else
	{
		comparison_object->input_reconstruction->ExtractSlice(current_projection, *comparison_object->angles, 1.0f, false);
		current_projection.SwapRealSpaceQuadrants();
		current_projection.BackwardFFT();
		current_projection.ChangePixelSize(&current_projection, comparison_object->pixel_size_factor, 0.001f, true);
	}

//	current_projection.QuickAndDirtyWriteSlice("projections.mrc", comparison_object->slice);
//	comparison_object->slice++;
	current_projection.MultiplyPixelWise(*comparison_object->projection_filter);
//	current_projection.BackwardFFT();
//	current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges());
//	current_projection.Resize(comparison_object->windowed_particle->logical_x_dimension, comparison_object->windowed_particle->logical_y_dimension, 1, 0.0f);
//	current_projection.ForwardFFT();
	current_projection.ZeroCentralPixel();
	current_projection.DivideByConstant(sqrtf(current_projection.ReturnSumOfSquares()));
#ifdef MKL
	// Use the MKL
	vmcMulByConj(current_projection.real_memory_allocated/2,reinterpret_cast <MKL_Complex8 *> (comparison_object->windowed_particle->complex_values),reinterpret_cast <MKL_Complex8 *> (current_projection.complex_values),reinterpret_cast <MKL_Complex8 *> (current_projection.complex_values),VML_EP|VML_FTZDAZ_ON|VML_ERRMODE_IGNORE);
#else
	for (long pixel_counter = 0; pixel_counter < current_projection.real_memory_allocated / 2; pixel_counter ++)
	{
		current_projection.complex_values[pixel_counter] = std::conj(current_projection.complex_values[pixel_counter]) * comparison_object->windowed_particle->complex_values[pixel_counter];
	}
#endif
	current_projection.BackwardFFT();
	//current_projection.QuickAndDirtyWriteSlice("cc.mrc",1);
//	wxPrintf("ping");


	// TODO find peak at center location


	int i, j;

	float best_score;
	float sq_dist_x, sq_dist_y;

	long address_in_current_projection = 0;
	peak_location_x = peak_location_x + current_projection.physical_address_of_box_center_x;
	peak_location_y = peak_location_y + current_projection.physical_address_of_box_center_y;

	for ( j = 0; j < current_projection.logical_y_dimension; j ++ )
	{
		sq_dist_y = float(pow(j-peak_location_y, 2));
		for ( i = 0; i < current_projection.logical_x_dimension; i ++ )
		{
			sq_dist_x = float(pow(i-peak_location_x,2));
			if (sq_dist_x == 0.0f && sq_dist_y == 0.0f)
			{
				best_score = current_projection.real_values[address_in_current_projection];
				wxPrintf("at address %ld cc = %f\n", address_in_current_projection, best_score);
			}
			address_in_current_projection++;
		}
		address_in_current_projection += current_projection.padding_jump_value;
	}





	//box_peak = current_projection.FindPeakWithIntegerCoordinates();
//	wxPrintf("address = %li\n", current_projection.FindPeakWithIntegerCoordinates().physical_address_within_image);
//	box_peak.x = 0.0f;
//	box_peak.y = 0.0f;
//	box_peak.value = current_projection.real_values[33152];
//	return box_peak;
}

IMPLEMENT_APP(NikoTestApp)

// override the DoInteractiveUserInput

void NikoTestApp::DoInteractiveUserInput()
{
	wxString	input_search_images;
	wxString	input_reconstruction;

	float		pixel_size = 1.0f;
	float		voltage_kV = 300.0f;
	float		spherical_aberration_mm = 2.7f;
	float		amplitude_contrast = 0.07f;
	float 		defocus1 = 10000.0f;
	float		defocus2 = 10000.0f;;
	float		defocus_angle;
	float 		phase_shift;
	float		low_resolution_limit = 300.0f;
	float		high_resolution_limit = 8.0f;
	float		angular_range = 2.0f;
	float		angular_step = 5.0f;
	int			best_parameters_to_keep = 20;
	float 		defocus_search_range = 1000;
	float 		defocus_search_step = 10;
//	float 		defocus_step = 5;
	float		pixel_size_search_range = 0.1f;
	float		pixel_size_step = 0.001f;
//	float		pixel_size_refine_step = 0.001f;
	float		padding = 1.0;
	bool		ctf_refinement = false;
	float		mask_radius = 0.0f;
	wxString	my_symmetry = "C1";
	float 		in_plane_angular_step = 0;
	float		xy_change_threshold = 10.0f;
	bool		exclude_above_xy_threshold = false;
	//int			result_number = 1;
	wxString preexisting_particle_file_name;
	float    min_peak_radius = 0.0;

	float   input_psi = 0.0;
	float   input_theta = 0.0;
	float   input_phi = 0.0;
	float   peak_position_x = 0.0;
	float   peak_position_y = 0.0;

	int			max_threads;
	int     image_index_in_stack = 1;


	UserInput *my_input = new UserInput("Save complete cross correlation", 1.00);

	input_search_images = my_input->GetFilenameFromUser("Input images to be searched", "The input image stack, containing the images that should be searched", "image_stack.mrc", true);
	input_reconstruction = my_input->GetFilenameFromUser("Input template reconstruction", "The 3D reconstruction from which projections are calculated", "reconstruction.mrc", true);
	pixel_size = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
	voltage_kV = my_input->GetFloatFromUser("Beam energy (keV)", "The energy of the electron beam used to image the sample in kilo electron volts", "300.0", 0.0);
	spherical_aberration_mm = my_input->GetFloatFromUser("Spherical aberration (mm)", "Spherical aberration of the objective lens in millimeters", "2.7");
	amplitude_contrast = my_input->GetFloatFromUser("Amplitude contrast", "Assumed amplitude contrast", "0.07", 0.0, 1.0);
	defocus1 = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
	defocus2 = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
	defocus_angle = my_input->GetFloatFromUser("Defocus angle (degrees)", "Defocus Angle for the input image", "0.0");
	phase_shift = my_input->GetFloatFromUser("Phase shift (degrees)", "Additional phase shift in degrees", "0.0");
//	low_resolution_limit = my_input->GetFloatFromUser("Low resolution limit (A)", "Low resolution limit of the data used for alignment in Angstroms", "300.0", 0.0);
//	high_resolution_limit = my_input->GetFloatFromUser("High resolution limit (A)", "High resolution limit of the data used for alignment in Angstroms", "8.0", 0.0);
//	angular_range = my_input->GetFloatFromUser("Angular refinement range", "AAngular range to refine", "2.0", 0.1);
	angular_step = my_input->GetFloatFromUser("Out of plane angular step", "Angular step size for global grid search", "0.2", 0.01);
	in_plane_angular_step = my_input->GetFloatFromUser("In plane angular step", "Angular step size for in-plane rotations during the search", "0.1", 0.01);
//	best_parameters_to_keep = my_input->GetIntFromUser("Number of top hits to refine", "The number of best global search orientations to refine locally", "20", 1);
	defocus_search_range = my_input->GetFloatFromUser("Defocus search range (A) (0.0 = no search)", "Search range (-value ... + value) around current defocus", "200.0", 0.0);
	defocus_search_step = my_input->GetFloatFromUser("Desired defocus accuracy (A)", "Accuracy to be achieved in defocus search", "10.0", 0.0);
//	defocus_refine_step = my_input->GetFloatFromUser("Defocus refine step (A) (0.0 = no refinement)", "Step size used in the defocus refinement", "5.0", 0.0);
	pixel_size_search_range = my_input->GetFloatFromUser("Pixel size search range (A) (0.0 = no search)", "Search range (-value ... + value) around current pixel size", "0.1", 0.0);
	pixel_size_step = my_input->GetFloatFromUser("Desired pixel size accuracy (A)", "Accuracy to be achieved in pixel size search", "0.01", 0.01);
//	pixel_size_refine_step = my_input->GetFloatFromUser("Pixel size refine step (A) (0.0 = no refinement)", "Step size used in the pixel size refinement", "0.001", 0.0);
	padding = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the input volume is padded to improve projections", "2.0", 1.0);
//	ctf_refinement = my_input->GetYesNoFromUser("Refine defocus", "Should the particle defocus be refined?", "No");
	mask_radius = my_input->GetFloatFromUser("Mask radius (A) (0.0 = no mask)", "Radius of a circular mask to be applied to the input particles during refinement", "0.0", 0.0);
//	my_symmetry = my_input->GetSymmetryFromUser("Template symmetry", "The symmetry of the template reconstruction", "C1");
	xy_change_threshold = my_input->GetFloatFromUser("Moved peak warning (A)", "Threshold for displaying warning of peak location changes during refinement", "10.0", 0.0);
	exclude_above_xy_threshold = my_input->GetYesNoFromUser("Exclude moving peaks", "Should the peaks that move more than the threshold be excluded from the output MIPs?", "No");

#ifdef _OPENMP
	max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "When threading, what is the max threads to run", "1", 1);
#else
	max_threads = 1;
#endif
  image_index_in_stack = my_input->GetIntFromUser("image index in an image stack", "image index in an image stack", "1", 1);
	preexisting_particle_file_name = my_input->GetFilenameFromUser("Frealign parameter file name", "an input parameter file to match reconstruction", "myparams.star", true );

	input_psi = my_input->GetFloatFromUser("Input angle psi", "input value for angle psi", "0.0", 0.0);
	input_theta = my_input->GetFloatFromUser("Input angle theta", "input value for angle theta", "0.0", 0.0);
	input_phi = my_input->GetFloatFromUser("Input angle phi", "input value for angle phi", "0.0", 0.0);
	peak_position_x = my_input->GetFloatFromUser("Input shift x", "input value for translation in x", "0.0", 0.0);
	peak_position_y = my_input->GetFloatFromUser("Input shift y", "input value for translation in y", "0.0", 0.0);


	int first_search_position = -1;
	int last_search_position = -1;
	//int image_number_for_gui = 0;
	//int number_of_jobs_per_image_in_gui = 0;
	float threshold_for_result_plotting = 0.0f;


	delete my_input;

//	my_current_job.Reset(42);

	my_current_job.ManualSetArguments("ttfffffffffffifffffbffffbtfiiiiftfffff",
															input_search_images.ToUTF8().data(),
															input_reconstruction.ToUTF8().data(),
															pixel_size,
															voltage_kV,
															spherical_aberration_mm,
															amplitude_contrast,
															defocus1,
															defocus2,
															defocus_angle,
															low_resolution_limit,
															high_resolution_limit,
															angular_range,
															angular_step,
															best_parameters_to_keep,
															defocus_search_range,
															defocus_search_step,
//															defocus_refine_step,
															pixel_size_search_range,
															pixel_size_step,
//															pixel_size_refine_step,
															padding,
															ctf_refinement,
															mask_radius,
															phase_shift,
															min_peak_radius,
															xy_change_threshold,
															exclude_above_xy_threshold,
															my_symmetry.ToUTF8().data(),
															in_plane_angular_step,
															first_search_position,
															last_search_position,
															max_threads,
															image_index_in_stack,
															threshold_for_result_plotting,
															preexisting_particle_file_name.ToUTF8().data(),
															input_psi,
															input_theta,
															input_phi,
															peak_position_x,
															peak_position_y
															);

		wxPrintf("Input finished...\n");
}

// override the do calculation method which will be what is actually run..


bool NikoTestApp::DoCalculation()
{

	// read in inputs
	wxDateTime start_time = wxDateTime::Now();

	wxString	input_search_images_filename = my_current_job.arguments[0].ReturnStringArgument();
	wxString	input_reconstruction_filename = my_current_job.arguments[1].ReturnStringArgument();
	float		pixel_size = my_current_job.arguments[2].ReturnFloatArgument();
	float		voltage_kV = my_current_job.arguments[3].ReturnFloatArgument();
	float		spherical_aberration_mm = my_current_job.arguments[4].ReturnFloatArgument();
	float		amplitude_contrast = my_current_job.arguments[5].ReturnFloatArgument();
	float 		defocus1 = my_current_job.arguments[6].ReturnFloatArgument();
	float		defocus2 = my_current_job.arguments[7].ReturnFloatArgument();
	float		defocus_angle = my_current_job.arguments[8].ReturnFloatArgument();;
	float		low_resolution_limit = my_current_job.arguments[9].ReturnFloatArgument();
	float		high_resolution_limit_search = my_current_job.arguments[10].ReturnFloatArgument();
	float		angular_range = my_current_job.arguments[11].ReturnFloatArgument();
	float		angular_step = my_current_job.arguments[12].ReturnFloatArgument();
	int			best_parameters_to_keep = my_current_job.arguments[13].ReturnIntegerArgument();
	float 		defocus_search_range = my_current_job.arguments[14].ReturnFloatArgument();
	float 		defocus_search_step = my_current_job.arguments[15].ReturnFloatArgument();
//	float 		defocus_refine_step = my_current_job.arguments[15].ReturnFloatArgument();
	float 		defocus_refine_step = 2.0f * defocus_search_step;
	float 		pixel_size_search_range = my_current_job.arguments[16].ReturnFloatArgument();
	float 		pixel_size_search_step = my_current_job.arguments[17].ReturnFloatArgument();
//	float 		pixel_size_refine_step = my_current_job.arguments[17].ReturnFloatArgument();
	float 		pixel_size_refine_step = 2.0f * pixel_size_search_step;
	float		padding = my_current_job.arguments[18].ReturnFloatArgument();
	bool		ctf_refinement = my_current_job.arguments[19].ReturnBoolArgument();
	float		mask_radius = my_current_job.arguments[20].ReturnFloatArgument();
	float 		phase_shift = my_current_job.arguments[21].ReturnFloatArgument();
	float		min_peak_radius = my_current_job.arguments[22].ReturnFloatArgument();
	float		xy_change_threshold = my_current_job.arguments[23].ReturnFloatArgument();
	bool		exclude_above_xy_threshold = my_current_job.arguments[24].ReturnBoolArgument();
	wxString 	my_symmetry = my_current_job.arguments[25].ReturnStringArgument();
	float		in_plane_angular_step = my_current_job.arguments[26].ReturnFloatArgument();
	int 		first_search_position = my_current_job.arguments[27].ReturnIntegerArgument();
	int 		last_search_position = my_current_job.arguments[28].ReturnIntegerArgument();
	//int 		image_number_for_gui = my_current_job.arguments[29].ReturnIntegerArgument();
	//int 		number_of_jobs_per_image_in_gui = my_current_job.arguments[30].ReturnIntegerArgument();
	//int 		result_number = my_current_job.arguments[31].ReturnIntegerArgument();
	int 		max_threads = my_current_job.arguments[29].ReturnIntegerArgument();
	int     image_index_in_stack = my_current_job.arguments[30].ReturnIntegerArgument();
	float		threshold_for_result_plotting = my_current_job.arguments[31].ReturnFloatArgument();
	wxString 	preexisting_particle_file_name = my_current_job.arguments[32].ReturnStringArgument();
	float		input_psi = my_current_job.arguments[33].ReturnFloatArgument();
	float		input_theta = my_current_job.arguments[34].ReturnFloatArgument();
	float		input_phi = my_current_job.arguments[35].ReturnFloatArgument();
	float		peak_position_x = my_current_job.arguments[36].ReturnFloatArgument();
	float		peak_position_y = my_current_job.arguments[37].ReturnFloatArgument();



	// create variables and allocate memory
	//if (is_running_locally == false) max_threads = number_of_threads_requested_on_command_line;

	int i,j;
	ParameterMap parameter_map;
	parameter_map.SetAllTrue();

	float outer_mask_radius;

	float temp_float;
	double temp_double_array[5];

	int number_of_rotations;
	long total_correlation_positions;
	long current_correlation_position;
	long pixel_counter;
	float sq_dist_x, sq_dist_y;
	long address;
	long best_address;

	int current_x;
	int current_y;

	int phi_i;
	int theta_i;
	int psi_i;
	int defocus_i;
	int defocus_is = 0;
	int size_i;
	int size_is = 0;

	float current_psi;
	float psi_max;
	float psi_start;
	float psi_step;

	int current_search_position;

	AnglesAndShifts angles;
	TemplateComparisonObject template_object;

	ImageFile input_search_image_file;
	ImageFile input_reconstruction_file;

	Curve whitening_filter;
	Curve number_of_terms;

	input_search_image_file.OpenFile(input_search_images_filename.ToStdString(), false);
	input_reconstruction_file.OpenFile(input_reconstruction_filename.ToStdString(), false);

	Image input_image;
	Image windowed_particle;
	Image padded_reference;
	Image input_reconstruction;

	Image projection_filter;

	Peak template_peak;
	Peak best_peak;
	long current_address;
	long address_offset;

	float current_defocus;
	float current_pixel_size;

	float current_shift_x;
	float current_shift_y;
	float current_shift_z;

	int ii, jj, kk, ll;
	float mult_i;
	float mult_i_start;
	//float defocus_step;
	float score_adjustment;
	float offset_distance;
//	float offset_warning_threshold = 10.0f;

	int number_of_peaks_found = 0;
	int peak_number;

	int number_preexisting_particles = 1;

	EulerSearch	global_euler_search;


	input_image.ReadSlice(&input_search_image_file, image_index_in_stack);
	padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);

	// read in template and pre-process
	input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices());
	if (padding != 1.0f)
	{
		input_reconstruction.Resize(input_reconstruction.logical_x_dimension * padding, input_reconstruction.logical_y_dimension * padding, input_reconstruction.logical_z_dimension * padding, input_reconstruction.ReturnAverageOfRealValuesOnEdges());
	}
	input_reconstruction.ForwardFFT();
	//input_reconstruction.CosineMask(0.1, 0.01, true);
	//input_reconstruction.Whiten();
	//if (first_search_position == 0) input_reconstruction.QuickAndDirtyWriteSlices("/tmp/filter.mrc", 1, input_reconstruction.logical_z_dimension);
	input_reconstruction.ZeroCentralPixel();
	input_reconstruction.SwapRealSpaceQuadrants();

	CTF input_ctf;


	/*  we will not use 0 or negative angle parameters
	if (angular_step <= 0)
	{
		angular_step = CalculateAngularStep(high_resolution_limit_search, mask_radius_search);
	}

	if (in_plane_angular_step <= 0)
	{
		psi_step = rad_2_deg(pixel_size / mask_radius_search);
		psi_step = 360.0 / int(360.0 / psi_step + 0.5);
	}
	else
	*/


	psi_step = in_plane_angular_step;


	// calculate radius of template
	temp_float = (float(input_reconstruction_file.ReturnXSize()) / 2.0f - 1.0f) * pixel_size;
	if (mask_radius > temp_float) mask_radius = temp_float;

	// for now, I am assuming the MTF has been applied already.
	// work out the filter to just whiten the image..

	whitening_filter.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
	number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

	wxDateTime my_time_out;
	wxDateTime my_time_in;

	psi_start = 0.0f;
	psi_max = 360.0f;
	global_euler_search.InitGrid(my_symmetry, angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
//	wxPrintf("%s",my_symmetry);
	if (my_symmetry.StartsWith("C1")) // TODO 2x check me - w/o this O symm at least is broken
	{
		if (global_euler_search.test_mirror == true) // otherwise the theta max is set to 90.0 and test_mirror is set to true.  However, I don't want to have to test the mirrors.
		{
			global_euler_search.theta_max = 180.0f;
		}
	}

	global_euler_search.CalculateGridSearchPositions(false);


	first_search_position = 0;
	last_search_position = global_euler_search.number_of_search_positions - 1;

	number_of_rotations = 0;

	for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
	{
		number_of_rotations++;
	}
	total_correlation_positions = 0.0;
	for (current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++)
	{
		//loop over each rotation angle

		for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
		{
			total_correlation_positions++;
		}
	}

	//Loop over ever search position

	wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
	wxPrintf("Searching %i rotations per position (psi).\n", number_of_rotations);
	wxPrintf("Total number of correlation positions: %ld\n", total_correlation_positions);

	wxPrintf("Performing Search...\n\n");
	//float 		*cc_output = new float[total_correlation_positions];


//	wxPrintf("Searching %i - %i of %i total positions\n", first_search_position, last_search_position, global_euler_search.number_of_search_positions);
//	wxPrintf("psi_start = %f, psi_max = %f, psi_step = %f\n", psi_start, psi_max, psi_step);

	// remove outliers

	input_image.ReplaceOutliersWithMean(5.0f);
	input_image.ForwardFFT();

	input_image.ZeroCentralPixel();
	input_image.Compute1DPowerSpectrumCurve(&whitening_filter, &number_of_terms);
	whitening_filter.SquareRoot();
	whitening_filter.Reciprocal();
	whitening_filter.MultiplyByConstant(1.0f / whitening_filter.ReturnMaximumValue());

	input_image.ApplyCurveFilter(&whitening_filter);
	input_image.ZeroCentralPixel();
	input_image.DivideByConstant(sqrt(input_image.ReturnSumOfSquares()));
	input_image.BackwardFFT();
	wxPrintf("input image preprocessing done...\n");

	/*
	// read in star file for peak location & pose
	if (! DoesFileExist(preexisting_particle_file_name))
	{
		SendError(wxString::Format("Error: Input parameter file %s not found\n", preexisting_particle_file_name));
		exit(-1);
	}


	FrealignParameterFile  parameter_file;
	parameter_file.Open(preexisting_particle_file_name,0,29);  // CHECKME if it is 17 values
	parameter_file.ReadFile(true, number_preexisting_particles);

	wxPrintf("\nReading in translational and rotional information for %i particles from the supplied parameter file\n", number_preexisting_particles);
	// Read the first line so that all of the values are initialized in parameter_vect  TODO if there are multiple images in an stack or multiple particles in an image
	float parameter_vect[29] = {0.0f};
	parameter_file.ReadLine(parameter_vect);
	current_psi   = parameter_vect[1];   //TODO for clarification, use arguments instead of files for input angles and x y.
	current_theta = parameter_vect[2];
	current_phi   = parameter_vect[3];

	current_shift_x = parameter_vect[4];
	current_shift_y = parameter_vect[5];
	current_shift_z = 0;
	parameter_file.Close();
	*/

	wxPrintf("input psi = %f\n", input_psi);
	wxPrintf("input theta = %f\n", input_theta);
	wxPrintf("input phi = %f\n", input_phi);


	// TODO: physical_address_of_box_center_x vs. shift_x

	wxPrintf("physical address of box center x = %i\n", input_image.physical_address_of_box_center_x);
	wxPrintf("physical address of box center y = %i\n", input_image.physical_address_of_box_center_y);


	// extract particle


	// FIXME: rearrange output_file
	// CHECKME: if total_correlation_positions == output numbers?
//	std::ofstream my_file ("new.txt",std::ofstream::out);

//	long records_to_read = total_correlation_positions;
//	float *temp_cc_array = new float [records_to_read];
	bool is_finished = false;
	long counter = 0;
	int image_pixel_counter = 0;

	int tid;
	#ifdef _OPENMP
		wxPrintf("max threads = %i\n",max_threads);
	#pragma omp parallel num_threads(max_threads) default(none) shared(input_image, mask_radius, pixel_size, \
		defocus_search_range, defocus_refine_step, pixel_size_search_range, \
		pixel_size_refine_step, defocus1, defocus2, defocus_angle, angular_step, in_plane_angular_step, whitening_filter, input_reconstruction, min_peak_radius, \
		input_reconstruction_file, voltage_kV, spherical_aberration_mm, amplitude_contrast, \
		phase_shift, max_threads, xy_change_threshold, exclude_above_xy_threshold, peak_position_x, peak_position_y, last_search_position, first_search_position,\
		global_euler_search, psi_max, psi_start, psi_step, counter) \
	private(padded_reference, windowed_particle, sq_dist_x, sq_dist_y, address, template_object, input_ctf, \
		angles, temp_float, projection_filter, template_peak, i, j, best_address, current_search_position, current_psi, tid)
{

	input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));

	// create template object to store windowed particle and template projection
	windowed_particle.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), true);
	projection_filter.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	template_object.input_reconstruction = &input_reconstruction;
	template_object.windowed_particle = &windowed_particle;
	template_object.projection_filter = &projection_filter;
	template_object.angles = &angles;

	padded_reference.CopyFrom(&input_image);
	padded_reference.RealSpaceIntegerShift(peak_position_x, peak_position_y); //TODO
	padded_reference.ClipInto(&windowed_particle);
	windowed_particle.ForwardFFT();
	windowed_particle.SwapRealSpaceQuadrants();
	template_object.pixel_size_factor = 1.0f;

//	windowed_particle.QuickAndDirtyWriteSlice("windowed_particle.mrc",1);


	//peak_position_x = peak_position_x + input_image.physical_address_of_box_center_x;
	//peak_position_y = peak_position_y + input_image.physical_address_of_box_center_y;



	// CHECK OMP
	//int tid = omp_get_thread_num();
	//printf("Hello world from omp thread %d\n", tid);


	// we are now located at the peak center

	// CHECK
	//input_image.real_values[image_pixel_counter] = 0;

	input_ctf.SetDefocus(defocus1 / pixel_size, defocus2 / pixel_size, deg_2_rad(defocus_angle));
	projection_filter.CalculateCTFImage(input_ctf);
	projection_filter.ApplyCurveFilter(&whitening_filter);


	#pragma omp for schedule(dynamic,1)
	//int current_correlation_position = 0;
	for (current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++)
	{
		//loop over each rotation angle
		angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], current_psi, 0.0, 0.0);
		for (current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step)
		{
			template_peak = TemplateScore(&template_object, peak_position_x, peak_position_y);
	//		tid = ReturnThreadNumberOfCurrentThread();

		//current_correlation_position++;
		//wxPrintf("template cc = %f calculated by thread %i\n", template_peak.value, tid);

	//	temp_cc_array[counter] = float(template_peak.value);
	//		counter ++;
		//cc_output[current_correlation_position] = template_peak.value;

		}
	//image_pixel_counter++;
	}
	//	image_pixel_counter+=input_image.padding_jump_value;
	//is_finished = true;
	//break;


	//wxPrintf("total records = %ld\n", counter);
	/*
	for (int i = 0; i < 10; i++)
	{
		wxPrintf("head of temp_cc_array %f\n",temp_cc_array[i]);
	}

	if (!my_file.is_open())
	{
		wxPrintf("failed to open\n");
	}
	else
	{
		wxPrintf("File opened successfully\n");
		my_file.write((char*)temp_cc_array, records_to_read);
	}

	delete [] temp_cc_array;
	*/
	}
	#else
		wxPrintf("not using openmp\n");
	#endif

	return true;

}
/*
bool NikoTestApp::DoCalculation()
{



	ImageFile input_psi_file, input_theta_file, input_phi_file;

	input_psi_file.OpenFile("LSU_psi.mrc", false);
	input_theta_file.OpenFile("LSU_theta.mrc", false);
	input_phi_file.OpenFile("LSU_phi.mrc", false);

	Image input_psi, input_theta, input_phi;
	input_psi.ReadSlice(&input_psi_file, 1);
	input_theta.ReadSlice(&input_theta_file, 1);
	input_phi.ReadSlice(&input_phi_file, 1);

	Image angle_distance;
	angle_distance.Allocate(input_psi.logical_x_dimension, input_psi.logical_y_dimension, true);
	angle_distance.SetToConstant(0.0f);
	long address_offset,current_address;
	long address = 0;
	int shift_x, shift_y;

	for (int j = 0; j < input_psi.logical_y_dimension; j++)
	{
		for (int i = 0; i < input_psi.logical_x_dimension; i++)
		{
			for (shift_y = -1; shift_y <= 1; shift_y ++)
			{
				for (shift_x = -1; shift_x <= 1; shift_x ++)
				{
					address_offset = (input_psi.logical_x_dimension + input_psi.padding_jump_value) * shift_y + shift_x;
					current_address = address + address_offset;
					angle_distance.real_values[address] += input_psi.real_values[current_address];

				}
			}
			angle_distance.real_values[address] -= input_psi.real_values[address];
			angle_distance.real_values[address] /= 8;

			address++;
		}
		address+=input_psi.padding_jump_value;
	}
	angle_distance.QuickAndDirtyWriteSlice("angle_distance.mrc",1);
	return true;
}
*/
