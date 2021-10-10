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
	MRCFile input_img_patch("6sct_1024_1A_3_no_shift.mrc", false);
	//MRCFile input_ref_patch("ref.mrc", false);
	Image img;
	Image larger_img;
	Image mask,local_mean,local_std;
	//Image ref;

	img.Allocate(input_img_patch.ReturnXSize(), input_img_patch.ReturnYSize(), true);
	img.ReadSlice(&input_img_patch, 1);
	//ref.Allocate(input_ref_patch.ReturnXSize(), input_ref_patch.ReturnYSize(), true);
	larger_img.Allocate(2048,2048,true);
	img.ClipIntoLargerRealSpace2D(&larger_img);
	//ref.ForwardFFT();

	larger_img.ReplaceOutliersWithMean(5.0f);

	mask.Allocate(larger_img.logical_x_dimension, larger_img.logical_y_dimension, 1);
	local_mean.Allocate(larger_img.logical_x_dimension, larger_img.logical_y_dimension, 1);
	local_std.Allocate(larger_img.logical_x_dimension, larger_img.logical_y_dimension, 1);

	mask.SetToConstant(1.0f);
	float P = mask.CosineMask(100 / 2, 1, false, true, 0.0);
	wxPrintf("P=%f\n",P);
	wxPrintf("sum=%f\n", mask.ReturnSumOfRealValues());
	mask.QuickAndDirtyWriteSlice("circular_mask.mrc",1);
	 // end of circular mask definition

	larger_img.ComputeLocalMeanAndVarianceMaps(&local_mean, &local_std, &mask, long(P));

	local_mean.QuickAndDirtyWriteSlice("circular_mask/local_mean.mrc",1);
	local_std.SquareRootRealValues();
	local_std.QuickAndDirtyWriteSlice("circular_mask/local_std.mrc",1);
	exit(0);

	/*
	// override original
	// calculates average density from projections
	AnglesAndShifts angles;
	EulerSearch global_euler_search;

	Image projection_filter;
	Image current_projection;
	Image template_reconstruction;
	Image average_density;
	Image padded_reference;
	//Image input_image;
	//ImageFile input_reconstruction_file;

	Image input_reconstruction;
	wxString	my_symmetry = "C1";
	float	angular_step = 20.0f;
	float in_plane_angular_step = 5.0f;
	float psi_step;
	float psi_max = 360.0f;
	float psi_start = 0.0f;
	float pixel_size = 1.0f; // requires user input
	float high_resolution_limit_search = 3.0f;
	float particle_radius_angstroms = 0.0f;
	float current_psi;

	float defocus1 = 4293.14f;
	float defocus2 = 4024.42f;

	float variance;


	CTF input_ctf;

	float	voltage_kV = 300.0f;
	float	spherical_aberration_mm = 2.7f;
	float amplitude_contrast = 0.07f;
	float defocus_angle = 20.48f;
	float phase_shift = 0.0f;

	long pixel_counter = 0;

	ParameterMap parameter_map; // needed for euler search init
	//for (int i = 0; i < 5; i++) {parameter_map[i] = true;}
	parameter_map.SetAllTrue();
	int best_parameters_to_keep = 20;

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

	//input_search_image_file.OpenFile("6sct_chainA_D500_P40_5000eA_1A.mrc", false);
	//input_image.ReadSlice(&input_search_image_file, 1);
	MRCFile input_reconstruction_file("6q8y_LSU_to_model5_bfactor_85_pix1.06.mrc", false);
	input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnZSize());
	//input_reconstruction_file.OpenFile("6q8y_LSU_to_model5_bfactor_85_pix1.06.mrc", false);
	//input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices());


	wxPrintf("template dim x = %i\n", input_reconstruction.logical_x_dimension);
	wxPrintf("template dim y = %i\n", input_reconstruction.logical_y_dimension);
	wxPrintf("template dim z = %i\n", input_reconstruction.logical_z_dimension);
	current_projection.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);
	projection_filter.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), false);

	padded_reference.Allocate(5832, 4096, true);
	padded_reference.SetToConstant(0.0f);
	template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
	average_density.Allocate(5832, 4096, true);
//	average_density.Allocate(input_reconstruction_file.ReturnXSize(), input_reconstruction_file.ReturnXSize(), true);
	average_density.SetToConstant(0.0f);

	input_reconstruction.ChangePixelSize(&template_reconstruction, 1.0f, 0.001f, true); // true is for FFT
	template_reconstruction.ZeroCentralPixel();
	template_reconstruction.SwapRealSpaceQuadrants();

	input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
	input_ctf.SetDefocus(defocus1 / pixel_size, defocus2 / pixel_size, deg_2_rad(defocus_angle));
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
			//variance = current_projection.ReturnSumOfSquares() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels \
					- powf(current_projection.ReturnAverageOfRealValues() * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
			//current_projection.DivideByConstant(sqrtf(variance));
			current_projection.ClipIntoLargerRealSpace2D(&padded_reference);
			average_density.AddImage(&padded_reference);
			total_search_position++;

		}
	}
	average_density.DivideByConstant(total_search_position);
	average_density.QuickAndDirtyWriteSlice("average_density_LSU.mrc",1);

	wxPrintf("total search positions = %i\n", total_search_position);
	*/
	/*
	int i,j;
	int count;
	int padded_dimensions_x;
	int padded_dimensions_y;
	int pad_factor = 6;
	float sigma;
	float peak;
	float sum_of_peaks;

	MRCFile input_file_3d("input3d.mrc", false);
	MRCFile input_file_2d("input2d.mrc", false);
	MRCFile output_file("output.mrc", true);
	Image input_volume;
	Image input_image;
	Image padded_image;
	Image output_image;
	Image temp_image;
	Image temp_image2;

	padded_dimensions_x = ReturnClosestFactorizedUpper(pad_factor * input_file_2d.ReturnXSize(), 3);
	padded_dimensions_y = ReturnClosestFactorizedUpper(pad_factor * input_file_2d.ReturnYSize(), 3);
	input_volume.Allocate(input_file_3d.ReturnXSize(), input_file_3d.ReturnYSize(), input_file_3d.ReturnZSize(), true);
	input_image.Allocate(input_file_2d.ReturnXSize(), input_file_2d.ReturnYSize(), true);
	padded_image.Allocate(padded_dimensions_x, padded_dimensions_y, true);

	input_volume.ReadSlices(&input_file_3d, 1, input_file_3d.ReturnZSize());
	input_image.ReadSlice(&input_file_2d, 1);
	sigma = sqrtf(input_image.ReturnVarianceOfRealValues());
	input_image.ClipIntoLargerRealSpace2D(&padded_image);
	padded_image.AddGaussianNoise(10.0 * sigma);
	padded_image.WriteSlice(&output_file, 1);
	padded_image.ForwardFFT();

	input_image.AddSlices(input_volume);
//	input_image.QuickAndDirtyWriteSlice("junk.mrc", 1);
	count = 1;
	output_image.CopyFrom(&padded_image);
	temp_image.CopyFrom(&input_image);
	temp_image.Resize(padded_dimensions_x, padded_dimensions_y, 1);
//	temp_image2.CopyFrom(&input_image);
//	temp_image2.Resize(padded_dimensions_x, padded_dimensions_y, 1);
//	temp_image2.RealSpaceIntegerShift(input_image.logical_x_dimension, input_image.logical_y_dimension, 0);
//	temp_image.AddImage(&temp_image2);
//	temp_image.QuickAndDirtyWriteSlice("junk.mrc", 1);
	temp_image.ForwardFFT();
	output_image.ConjugateMultiplyPixelWise(temp_image);
	output_image.SwapRealSpaceQuadrants();
	output_image.BackwardFFT();
	peak = output_image.real_values[output_image.logical_x_dimension / 2 + (output_image.logical_x_dimension + output_image.padding_jump_value) * output_image.logical_y_dimension / 2];
	wxPrintf("\nPeak with whole projection = %g background = %g\n\n", peak, output_image.ReturnVarianceOfRealValues(float(2 * input_image.logical_x_dimension), 0.0, 0.0, 0.0, true));
//	wxPrintf("\nPeak with whole projection = %g\n\n", output_image.ReturnMaximumValue());
	output_image.WriteSlice(&output_file, 2);
	sum_of_peaks = 0.0;
	for (i = 1; i <= 3; i += 2)
	{
		for (j = 1; j <= 3; j += 2)
		{
			output_image.CopyFrom(&padded_image);
			temp_image.CopyFrom(&input_image);
			temp_image.SquareMaskWithValue(float(input_image.logical_x_dimension) / 2.0, 0.0, false, i * input_image.logical_x_dimension / 4, j * input_image.logical_x_dimension / 4);
			temp_image.Resize(padded_dimensions_x, padded_dimensions_y, 1);
			temp_image.ForwardFFT();
			output_image.ConjugateMultiplyPixelWise(temp_image);
			output_image.SwapRealSpaceQuadrants();
			output_image.BackwardFFT();
			peak = output_image.real_values[output_image.logical_x_dimension / 2 + (output_image.logical_x_dimension + output_image.padding_jump_value) * output_image.logical_y_dimension / 2];
//			peak = output_image.ReturnMaximumValue();
			wxPrintf("Quarter peak = %i %i %g\n", i, j, peak);
			sum_of_peaks += peak;
			count++;
			output_image.WriteSlice(&output_file, 1 + count);
		}
	}
	wxPrintf("\nSum of quarter peaks = %g\n\n", sum_of_peaks);

	sum_of_peaks = 0.0;
	for (i = 0; i < 4; i ++)
	{
		output_image.CopyFrom(&padded_image);
		input_image.AddSlices(input_volume, i * input_volume.logical_z_dimension / 4 + 1, (i + 1) * input_volume.logical_z_dimension / 4);
//		input_image.QuickAndDirtyWriteSlice("junk.mrc", i + 1);
		temp_image.CopyFrom(&input_image);
		temp_image.Resize(padded_dimensions_x, padded_dimensions_y, 1);
		temp_image.ForwardFFT();
		output_image.ConjugateMultiplyPixelWise(temp_image);
		output_image.SwapRealSpaceQuadrants();
		output_image.BackwardFFT();
		peak = output_image.real_values[output_image.logical_x_dimension / 2 + (output_image.logical_x_dimension + output_image.padding_jump_value) * output_image.logical_y_dimension / 2];
//		peak = output_image.ReturnMaximumValue();
		wxPrintf("Slice peak = %i %g\n", i + 1, peak);
		sum_of_peaks += peak;
		count++;
		output_image.WriteSlice(&output_file, 1 + count);
	}
	wxPrintf("\nSum of slice peaks = %g\n", sum_of_peaks);
 */
/*	wxPrintf("\nDoing 1000 FFTs %i x %i\n", output_image.logical_x_dimension, output_image.logical_y_dimension);
	for (i = 0; i < 1000; i++)
	{
		output_image.is_in_real_space = false;
		output_image.SetToConstant(1.0);
		output_image.BackwardFFT();
	}
	wxPrintf("\nFinished\n");
*/

/*	int i, j;
	int slice_thickness;
	int first_slice, last_slice, middle_slice;
	long offset;
	long pixel_counter;
	float bfactor = 20.0;
	float mask_radius = 75.0;
	float pixel_size = 1.237;
//	float pixel_size = 0.97;
	float bfactor_res_limit = 8.0;
	float resolution_limit = 3.8;
//	float resolution_limit = 3.0;
	float cosine_edge = 5.0;
	float bfactor_pixels;

	MRCFile input_file("input.mrc", false);
	MRCFile output_file_2D("output2D.mrc", true);
	MRCFile output_file_3D("output3D.mrc", true);
	Image input_image;
	Image output_image;
	Image output_image_3D;

	Curve power_spectrum;
	Curve number_of_terms;

	UserInput my_input("NikoTest", 1.00);
	pixel_size = my_input.GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
	mask_radius = my_input.GetFloatFromUser("Mask radius (A)", "Radius of mask to be applied to input 3D map, in Angstroms", "100.0", 0.0);
	bfactor = my_input.GetFloatFromUser("B-Factor (A^2)", "B-factor to be applied to dampen the 3D map after spectral flattening, in Angstroms squared", "20.0");
	bfactor_res_limit = my_input.GetFloatFromUser("Low resolution limit for spectral flattening (A)", "The resolution at which spectral flattening starts being applied, in Angstroms", "8.0", 0.0);
	resolution_limit = my_input.GetFloatFromUser("High resolution limit (A)", "Resolution of low-pass filter applied to final output maps, in Angstroms", "3.0", 0.0);

	slice_thickness = myroundint(resolution_limit / pixel_size);
	input_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), slice_thickness, true);
	output_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), true);
	output_image_3D.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), input_file.ReturnZSize(), true);

	wxPrintf("\nCalculating 3D spectrum...\n");

	power_spectrum.SetupXAxis(0.0, 0.5 * sqrtf(3.0), int((output_image_3D.logical_x_dimension / 2.0 + 1.0) * sqrtf(3.0) + 1.0));
	number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(3.0), int((output_image_3D.logical_x_dimension / 2.0 + 1.0) * sqrtf(3.0) + 1.0));
	output_image_3D.ReadSlices(&input_file, 1, input_file.ReturnZSize());
	output_image_3D.CosineMask(mask_radius / pixel_size, 10.0 / pixel_size);

	first_slice = int((input_file.ReturnZSize() - slice_thickness + 1) / 2.0);
	last_slice = first_slice + slice_thickness;
	pixel_counter = 0;
	for (j = first_slice; j < last_slice; j++)
	{
		offset = j * (output_image_3D.logical_x_dimension + output_image_3D.padding_jump_value) * output_image_3D.logical_y_dimension;
		for (i = 0; i < (output_image_3D.logical_x_dimension + output_image_3D.padding_jump_value) * output_image_3D.logical_y_dimension; i++) {input_image.real_values[pixel_counter] = output_image_3D.real_values[i + offset]; pixel_counter++;}
	}

	output_image_3D.ForwardFFT();
	output_image_3D.Compute1DPowerSpectrumCurve(&power_spectrum, &number_of_terms);
	power_spectrum.SquareRoot();
	wxPrintf("Done with 3D spectrum. Starting slice estimation...\n");

//	input_image.ReadSlices(&input_file, first_slice, last_slice);
	bfactor_pixels = bfactor / pixel_size / pixel_size;
	input_image.ForwardFFT();
	input_image.ApplyBFactorAndWhiten(power_spectrum, bfactor_pixels, bfactor_pixels, pixel_size / bfactor_res_limit);
//	input_image.ApplyBFactor(bfactor_pixels);
//	input_image.CosineMask(pixel_size / resolution_limit, cosine_edge / input_file.ReturnXSize());
	input_image.CosineMask(pixel_size / resolution_limit - pixel_size / 40.0, pixel_size / 20.0);
	input_image.BackwardFFT();

	middle_slice = int(slice_thickness / 2.0);
	offset = middle_slice * (input_file.ReturnXSize() + input_image.padding_jump_value) * input_file.ReturnYSize();
	pixel_counter = 0;
	for (i = 0; i < (input_file.ReturnXSize() + input_image.padding_jump_value) * input_file.ReturnYSize(); i++) {output_image.real_values[pixel_counter] = input_image.real_values[i + offset]; pixel_counter++;}
//	output_image.ForwardFFT();
//	output_image.CosineMask(pixel_size / resolution_limit - pixel_size / 40.0, pixel_size / 20.0);
//	output_image.BackwardFFT();
	output_image.WriteSlice(&output_file_2D, 1);
	wxPrintf("Done with slices. Starting 3D B-factor application...\n");

	output_image_3D.ApplyBFactorAndWhiten(power_spectrum, bfactor_pixels, bfactor_pixels, pixel_size / bfactor_res_limit);
	output_image_3D.CosineMask(pixel_size / resolution_limit - pixel_size / 40.0, pixel_size / 20.0);
	output_image_3D.BackwardFFT();
	output_image_3D.WriteSlices(&output_file_3D, 1, input_file.ReturnZSize());
	wxPrintf("Done with 3D B-factor application.\n");
/*
	int i;
	int min_class;
	int max_class;
	int count;
	float temp_float;
	float input_parameters[17];

	MRCFile input_file("input.mrc", false);
	MRCFile output_file("output.mrc", true);
	Image input_image;
	Image padded_image;
	Image ctf_image;
	Image sum_image;
	CTF ctf;
	AnglesAndShifts rotation_angle;

	input_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), true);
	padded_image.Allocate(4 * input_file.ReturnXSize(), 4 * input_file.ReturnYSize(), true);
	ctf_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), false);
	sum_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), false);
	sum_image.SetToConstant(0.0);

	FrealignParameterFile input_par_file("input.par", OPEN_TO_READ);
	input_par_file.ReadFile();

//	count = 0;
//	for (i = 1; i <= input_par_file.number_of_lines; i++)
//	{
//		if (i % 100 == 1) wxPrintf("Working on line %i\n", i);
//		input_par_file.ReadLine(input_parameters);
//		input_image.ReadSlice(&input_file, int(input_parameters[0] + 0.5));
//		count++;
//		input_image.WriteSlice(&output_file, count);
//	}

	for (i = 1; i <= input_par_file.number_of_lines; i++)
	{
		if (i % 100 == 1) wxPrintf("Rotating image %i\n", i);
		input_par_file.ReadLine(input_parameters);
		input_image.ReadSlice(&input_file, i);
		input_image.RealSpaceIntegerShift(-input_parameters[4], -input_parameters[5]);
		input_image.ForwardFFT();
		input_image.ClipInto(&padded_image);
		padded_image.BackwardFFT();
		rotation_angle.GenerateRotationMatrix2D(-input_parameters[1]);
		padded_image.Rotate2DSample(input_image, rotation_angle);
		input_image.WriteSlice(&output_file, i);
		if (input_parameters[7] == 2)
		{
			ctf.Init(300.0, 0.0, 0.1, input_parameters[8], input_parameters[9], input_parameters[10], 0.0, 0.0, 0.0, 1.0, input_parameters[11]);
			ctf_image.CalculateCTFImage(ctf);
			input_image.ForwardFFT();
			input_image.PhaseFlipPixelWise(ctf_image);
			sum_image.AddImage(&input_image);
		}
//		if (i == 1001) break;
	}
	sum_image.QuickAndDirtyWriteSlice("sum.mrc", 1);
*/

/*	FrealignParameterFile input_par_file("input.par", OPEN_TO_READ);
	FrealignParameterFile output_par_file("output.par", OPEN_TO_WRITE);
	input_par_file.ReadFile(true);
	input_par_file.ReduceAngles();
	min_class = myroundint(input_par_file.ReturnMin(7));
	max_class = myroundint(input_par_file.ReturnMax(7));
	for (i = min_class; i <= max_class; i++)
	{
		temp_float = input_par_file.ReturnDistributionMax(2, i);
		if (temp_float != 0.0) wxPrintf("theta max, sigma = %i %g %g\n", i, temp_float, input_par_file.ReturnDistributionSigma(2, temp_float, i));
//		input_par_file.SetParameters(2, temp_float, i);
		temp_float = input_par_file.ReturnDistributionMax(3, i);
		if (temp_float != 0.0) wxPrintf("phi max, sigma = %i %g %g\n", i, temp_float, input_par_file.ReturnDistributionSigma(3, temp_float, i));
//		input_par_file.SetParameters(3, temp_float, i);
	} */
//	for (i = 1; i <= input_par_file.number_of_lines; i++)
//	{
//		input_par_file.ReadLine(input_parameters);
//		output_par_file.WriteLine(input_parameters);
//	}

//	MRCFile input_file("input.mrc", false);
//	MRCFile output_file("output.mrc", true);
//	Image input_image;
//	Image filtered_image;
//	Image kernel;

//	input_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), input_file.ReturnZSize(), true);
//	filtered_image.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), input_file.ReturnZSize(), true);
//	kernel.Allocate(input_file.ReturnXSize(), input_file.ReturnYSize(), input_file.ReturnZSize(), true);
//	input_image.ReadSlices(&input_file,1,input_image.logical_z_dimension);

/*	kernel.SetToConstant(1.0);
	kernel.CosineMask(8.0, 8.0, false, true, 0.0);
//	kernel.real_values[0] = 1.0;
	temp_float = kernel.ReturnAverageOfRealValues() * kernel.number_of_real_space_pixels;
//	wxPrintf("average = %g\n", temp_float);
//	kernel.WriteSlices(&output_file,1,input_image.logical_z_dimension);
	kernel.ForwardFFT();
	kernel.SwapRealSpaceQuadrants();
	kernel.MultiplyByConstant(float(kernel.number_of_real_space_pixels) / temp_float);
//	kernel.CosineMask(0.03, 0.03, true);

	input_image.SetMinimumValue(0.0);
	filtered_image.CopyFrom(&input_image);
	filtered_image.ForwardFFT();
	filtered_image.MultiplyPixelWise(kernel);
//	filtered_image.CosineMask(0.01, 0.02);
	filtered_image.BackwardFFT();
//	filtered_image.MultiplyByConstant(0.3);
	input_image.SubtractImage(&filtered_image);
*/
//	input_image.SetToConstant(1.0);
//	input_image.CorrectSinc(45.0, 1.0, true, 0.0);
//	for (i = 0; i < input_image.real_memory_allocated; i++) if (input_image.real_values[i] < 0.0) input_image.real_values[i] = -log(-input_image.real_values[i] + 1.0);
//	input_image.WriteSlices(&output_file,1,input_image.logical_z_dimension);
//	temp_float = -420.5; wxPrintf("%g\n", fmodf(temp_float, 360.0));
//	filtered_image.WriteSlices(&output_file,1,img.logical_z_dimension);

	return true;
}
