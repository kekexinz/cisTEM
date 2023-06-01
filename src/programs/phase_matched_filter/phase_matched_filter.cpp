#include "../../core/core_headers.h"

class
        PhaseMatchedFilterApp : public MyApp {

  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );

  private:
};

IMPLEMENT_APP(PhaseMatchedFilterApp)

// override the DoInteractiveUserInput

void PhaseMatchedFilterApp::DoInteractiveUserInput( ) {
    wxString input_search_images;
    wxString input_reconstruction;

    wxString output_histogram_file;
    wxString correlation_avg_output_file;
    wxString correlation_std_output_file;
    wxString output_mip_file;

    float pixel_size              = 1.0f;
    float voltage_kV              = 300.0f;
    float spherical_aberration_mm = 2.7f;
    float amplitude_contrast      = 0.07f;
    float defocus1                = 10000.0f;
    float defocus2                = 10000.0f;

    float    defocus_angle;
    float    phase_shift;
    float    low_resolution_limit      = 300.0;
    float    high_resolution_limit     = 8.0;
    float    angular_step              = 5.0;
    int      best_parameters_to_keep   = 20;
    float    defocus_search_range      = 500;
    float    defocus_step              = 50;
    float    pixel_size_search_range   = 0.1f;
    float    pixel_size_step           = 0.02f;
    float    padding                   = 1.0;
    float    particle_radius_angstroms = 0.0f;
    wxString my_symmetry               = "C1";
    float    in_plane_angular_step     = 0;
    int      max_threads               = 1; // Only used for the GPU code

    bool  do_constrained_search   = false;
    float phi_start               = 0.0;
    float phi_max                 = 360.0;
    float theta_start             = 0.0;
    float theta_max               = 180.0;
    float psi_start               = 0.0;
    float psi_max                 = 360.0;
    bool  do_phase_matched_filter = false;
    int   whitening_option        = 1;
    bool  do_whiten_image         = true;

    UserInput* my_input = new UserInput("PhaseMatchedFilter", 1.00);

    input_search_images   = my_input->GetFilenameFromUser("Input images to be searched", "The input image stack, containing the images that should be searched", "image_stack.mrc", true);
    input_reconstruction  = my_input->GetFilenameFromUser("Input template reconstruction", "The 3D reconstruction from which projections are calculated", "reconstruction.mrc", true);
    output_histogram_file = my_input->GetFilenameFromUser("Output histogram of correlation values", "histogram of all correlation values", "histogram.txt", false);
    output_mip_file       = my_input->GetFilenameFromUser("Output mip", "mip", "mip.mrc", false);
    //correlation_avg_output_file = my_input->GetFilenameFromUser("Correlation average value", "The file for saving the average value of all correlation images", "corr_average.mrc", false);
    //correlation_std_output_file = my_input->GetFilenameFromUser("Correlation standard deviation value", "The file for saving the std value of all correlation images", "corr_variance.mrc", false);
    pixel_size              = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    voltage_kV              = my_input->GetFloatFromUser("Beam energy (keV)", "The energy of the electron beam used to image the sample in kilo electron volts", "300.0", 0.0);
    spherical_aberration_mm = my_input->GetFloatFromUser("Spherical aberration (mm)", "Spherical aberration of the objective lens in millimeters", "2.7");
    amplitude_contrast      = my_input->GetFloatFromUser("Amplitude contrast", "Assumed amplitude contrast", "0.07", 0.0, 1.0);
    defocus1                = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    defocus2                = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    defocus_angle           = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    phase_shift             = my_input->GetFloatFromUser("Phase Shift (degrees)", "Additional phase shift in degrees", "0.0");
    //    low_resolution_limit = my_input->GetFloatFromUser("Low resolution limit (A)", "Low resolution limit of the data used for alignment in Angstroms", "300.0", 0.0);
    high_resolution_limit = my_input->GetFloatFromUser("High resolution limit (A)", "High resolution limit of the data used for alignment in Angstroms", "8.0", 0.0);
    angular_step          = my_input->GetFloatFromUser("Out of plane angular step (0.0 = set automatically)", "Angular step size for global grid search", "0.0", 0.0);
    in_plane_angular_step = my_input->GetFloatFromUser("In plane angular step (0.0 = set automatically)", "Angular step size for in-plane rotations during the search", "0.0", 0.0);
    do_constrained_search = my_input->GetYesNoFromUser("Perform constraied angular search", "yes no", "no");
    if ( do_constrained_search ) {
        phi_max     = my_input->GetFloatFromUser("Constrained search Phi max", "max phi for constrained search", "360.0", 0.0, 360.0);
        phi_start   = my_input->GetFloatFromUser("Constrained search Phi min", "min phi for constrained search", "0.0", 0.0, 360.0);
        theta_max   = my_input->GetFloatFromUser("Constrained search Theta max", "max theta for constrained search", "180.0", 0.0, 180.0);
        theta_start = my_input->GetFloatFromUser("Constrained search Theta min", "max theta for constrained search", "0.0", 0.0, 180.0);
        psi_max     = my_input->GetFloatFromUser("Constrained search Psi max", "max psi for constrained search", "360.0", 0.0, 360.0);
        psi_start   = my_input->GetFloatFromUser("Constrained search Psi min", "max psi for constrained search", "0.0.0", 0.0, 360.0);
    }

    float fixed_phi   = my_input->GetFloatFromUser("Constrained search Phi", "max phi for constrained search", "360.0", 0.0, 360.0);
    float fixed_theta = my_input->GetFloatFromUser("Constrained search Theta", "max theta for constrained search", "360.0", 0.0, 360.0);
    float fixed_psi   = my_input->GetFloatFromUser("Constrained search Psi", "max psi for constrained search", "360.0", 0.0, 360.0);

    //    best_parameters_to_keep = my_input->GetIntFromUser("Number of top hits to refine", "The number of best global search orientations to refine locally", "20", 1);
    padding                   = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the input volume is padded to improve projections", "1.0", 1.0, 2.0);
    particle_radius_angstroms = my_input->GetFloatFromUser("Mask radius for global search (A) (0.0 = max)", "Radius of a circular mask to be applied to the input images during global search", "0.0", 0.0);
    my_symmetry               = my_input->GetSymmetryFromUser("Template symmetry", "The symmetry of the template reconstruction", "C1");
    do_phase_matched_filter   = my_input->GetYesNoFromUser("Perform phase-only matched filter in cc calculation", "yes no", "no");
    whitening_option          = my_input->GetIntFromUser("How to whiten ref? (1: by img  2: by self  3: no whitening)", "1: by img  2: by self  3: no whitening", "1", 1);
    do_whiten_image           = my_input->GetYesNoFromUser("Whiten image?", "yes no", "yes");

    int   first_search_position = -1;
    int   last_search_position  = -1;
    float min_peak_radius       = 10.0f;

    delete my_input;

    my_current_job.ManualSetArguments("ttffffffffffiffffffftftiifbffffffffftbib", input_search_images.ToUTF8( ).data( ),
                                      input_reconstruction.ToUTF8( ).data( ),
                                      pixel_size,
                                      voltage_kV,
                                      spherical_aberration_mm,
                                      amplitude_contrast,
                                      defocus1,
                                      defocus2,
                                      defocus_angle,
                                      low_resolution_limit,
                                      high_resolution_limit,
                                      angular_step,
                                      best_parameters_to_keep,
                                      defocus_search_range,
                                      defocus_step,
                                      pixel_size_search_range,
                                      pixel_size_step,
                                      padding,
                                      particle_radius_angstroms,
                                      phase_shift,
                                      my_symmetry.ToUTF8( ).data( ),
                                      in_plane_angular_step,
                                      output_histogram_file.ToUTF8( ).data( ),
                                      first_search_position,
                                      last_search_position,
                                      min_peak_radius,
                                      do_constrained_search,
                                      phi_start,
                                      phi_max,
                                      theta_start,
                                      theta_max,
                                      psi_start,
                                      psi_max,
                                      fixed_phi,
                                      fixed_theta,
                                      fixed_psi,
                                      output_mip_file.ToUTF8( ).data( ),
                                      do_phase_matched_filter,
                                      whitening_option,
                                      do_whiten_image
                                      //correlation_avg_output_file.ToUTF8( ).data( ),
                                      //correlation_std_output_file.ToUTF8( ).data( )

    );
}

// override the do calculation method which will be what is actually run..

bool PhaseMatchedFilterApp::DoCalculation( ) {
    wxString input_search_images_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_reconstruction_filename = my_current_job.arguments[1].ReturnStringArgument( );
    float    pixel_size                    = my_current_job.arguments[2].ReturnFloatArgument( );
    float    voltage_kV                    = my_current_job.arguments[3].ReturnFloatArgument( );
    float    spherical_aberration_mm       = my_current_job.arguments[4].ReturnFloatArgument( );
    float    amplitude_contrast            = my_current_job.arguments[5].ReturnFloatArgument( );
    float    defocus1                      = my_current_job.arguments[6].ReturnFloatArgument( );
    float    defocus2                      = my_current_job.arguments[7].ReturnFloatArgument( );
    float    defocus_angle                 = my_current_job.arguments[8].ReturnFloatArgument( );
    float    low_resolution_limit          = my_current_job.arguments[9].ReturnFloatArgument( );
    float    high_resolution_limit_search  = my_current_job.arguments[10].ReturnFloatArgument( );
    float    angular_step                  = my_current_job.arguments[11].ReturnFloatArgument( );
    int      best_parameters_to_keep       = my_current_job.arguments[12].ReturnIntegerArgument( );
    float    defocus_search_range          = my_current_job.arguments[13].ReturnFloatArgument( );
    float    defocus_step                  = my_current_job.arguments[14].ReturnFloatArgument( );
    float    pixel_size_search_range       = my_current_job.arguments[15].ReturnFloatArgument( );
    float    pixel_size_step               = my_current_job.arguments[16].ReturnFloatArgument( );
    float    padding                       = my_current_job.arguments[17].ReturnFloatArgument( );
    float    particle_radius_angstroms     = my_current_job.arguments[18].ReturnFloatArgument( );
    float    phase_shift                   = my_current_job.arguments[19].ReturnFloatArgument( );
    wxString my_symmetry                   = my_current_job.arguments[20].ReturnStringArgument( );
    float    in_plane_angular_step         = my_current_job.arguments[21].ReturnFloatArgument( );
    wxString output_histogram_file         = my_current_job.arguments[22].ReturnStringArgument( );
    int      first_search_position         = my_current_job.arguments[23].ReturnIntegerArgument( );
    int      last_search_position          = my_current_job.arguments[24].ReturnIntegerArgument( );
    float    min_peak_radius               = my_current_job.arguments[25].ReturnFloatArgument( );
    bool     do_constrained_search         = my_current_job.arguments[26].ReturnBoolArgument( );
    float    phi_start                     = my_current_job.arguments[27].ReturnFloatArgument( );
    float    phi_max                       = my_current_job.arguments[28].ReturnFloatArgument( );
    float    theta_start                   = my_current_job.arguments[29].ReturnFloatArgument( );
    float    theta_max                     = my_current_job.arguments[30].ReturnFloatArgument( );
    float    psi_start                     = my_current_job.arguments[31].ReturnFloatArgument( );
    float    psi_max                       = my_current_job.arguments[32].ReturnFloatArgument( );
    float    fixed_phi                     = my_current_job.arguments[33].ReturnFloatArgument( );
    float    fixed_theta                   = my_current_job.arguments[34].ReturnFloatArgument( );
    float    fixed_psi                     = my_current_job.arguments[35].ReturnFloatArgument( );
    wxString output_mip_file               = my_current_job.arguments[36].ReturnStringArgument( );
    bool     do_phase_matched_filter       = my_current_job.arguments[37].ReturnBoolArgument( );
    int      whitening_option              = my_current_job.arguments[38].ReturnIntegerArgument( );
    bool     do_whiten_image               = my_current_job.arguments[39].ReturnBoolArgument( );

    //wxString correlation_avg_output_file   = my_current_job.arguments[29].ReturnStringArgument( );
    //wxString correlation_std_output_file   = my_current_job.arguments[30].ReturnStringArgument( );

    // read in image file (with particle in the center, no randomness in x, y, z, defocus1, defocus2, defocus angle when simulating particles)
    // read in 3d template
    ImageFile input_search_image_file;
    ImageFile input_reconstruction_file;
    input_search_image_file.OpenFile(input_search_images_filename.ToStdString( ), false);
    input_reconstruction_file.OpenFile(input_reconstruction_filename.ToStdString( ), false);

    Image input_image;
    Image input_reconstruction;
    Image template_reconstruction;
    Image padded_projection;
    Image current_projection;
    Image projection_filter;
    Image padded_reference;

    Image max_intensity_projection;
    Image correlation_pixel_sum_image;
    Image correlation_pixel_sum_of_squares_image;

    input_image.ReadSlice(&input_search_image_file, 1);
    input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices( ));
    template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
    current_projection.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
    projection_filter.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
    padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    padded_reference.SetToConstant(0.0f);
    max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);

    correlation_pixel_sum_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    correlation_pixel_sum_of_squares_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    double* correlation_pixel_sum            = new double[input_image.real_memory_allocated];
    double* correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];

    padded_reference.SetToConstant(0.0f);
    max_intensity_projection.SetToConstant(-FLT_MAX);
    //best_psi.SetToConstant(0.0f);
    //best_theta.SetToConstant(0.0f);
    //best_phi.SetToConstant(0.0f);
    //best_defocus.SetToConstant(0.0f);

    ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
    ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);

    if ( padding != 1.0f ) {
        input_reconstruction.Resize(input_reconstruction.logical_x_dimension * padding, input_reconstruction.logical_y_dimension * padding, input_reconstruction.logical_z_dimension * padding, input_reconstruction.ReturnAverageOfRealValuesOnEdges( ));

        padded_projection.Allocate(input_reconstruction_file.ReturnXSize( ) * padding, input_reconstruction_file.ReturnXSize( ) * padding, false);
    }

    input_reconstruction.ChangePixelSize(&template_reconstruction, 1, 0.001f, true);
    template_reconstruction.ZeroCentralPixel( );
    template_reconstruction.SwapRealSpaceQuadrants( );

    // setup curve
    double sqrt_input_pixels = sqrt((double)(input_image.logical_x_dimension * input_image.logical_y_dimension));

    int   histogram_number_of_points = 512;
    float histogram_min              = -12.5f;
    float histogram_max              = 30.0f;

    float  histogram_step        = (histogram_max - histogram_min) / float(histogram_number_of_points);
    double histogram_min_scaled  = histogram_min / sqrt_input_pixels;
    double histogram_step_scaled = histogram_step / sqrt_input_pixels;

    long* histogram_data;
    int   current_bin;

    histogram_data = new long[histogram_number_of_points];

    for ( int counter = 0; counter < histogram_number_of_points; counter++ ) {
        histogram_data[counter] = 0;
    }

    // preprocess (whitening, normalization)
    Curve whitening_filter_ref, whitening_filter_img;
    Curve number_of_terms_ref, number_of_terms_img;
    whitening_filter_img.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    number_of_terms_img.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    whitening_filter_ref.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((template_reconstruction.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    number_of_terms_ref.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((template_reconstruction.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

    input_image.ReplaceOutliersWithMean(5.0f);
    input_image.ForwardFFT( );
    input_image.SwapRealSpaceQuadrants( );

    input_image.ZeroCentralPixel( );
    input_image.Compute1DPowerSpectrumCurve(&whitening_filter_img, &number_of_terms_img);
    whitening_filter_img.SquareRoot( );
    whitening_filter_img.Reciprocal( );
    whitening_filter_img.MultiplyByConstant(1.0f / whitening_filter_img.ReturnMaximumValue( ));
    //whitening_filter_img.WriteToFile("w_filter_img.txt");
    if ( do_whiten_image ) {
        input_image.ApplyCurveFilter(&whitening_filter_img);
        input_image.ZeroCentralPixel( );
        input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares( )));
    }
    else
        input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares( )));
    // set up CTF
    CTF input_ctf;
    input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
    input_ctf.SetDefocus(defocus1 / pixel_size, defocus2 / pixel_size, deg_2_rad(defocus_angle));
    projection_filter.CalculateCTFImage(input_ctf);

    if ( whitening_option == 1 )
        projection_filter.ApplyCurveFilter(&whitening_filter_img);

    AnglesAndShifts angles;
    EulerSearch     global_euler_search;
    ParameterMap    parameter_map; // needed for euler search init
    //for (int i = 0; i < 5; i++) {parameter_map[i] = true;}
    parameter_map.SetAllTrue( );
    float variance;
    float current_psi;
    float psi_step;
    int   number_of_rotations;
    long  total_correlation_positions;
    long  current_correlation_position;
    long  total_correlation_positions_per_thread;
    long  pixel_counter;

    int current_search_position;
    int current_x;
    int current_y;

    // set up search grid
    // in-plane PSI
    psi_step = in_plane_angular_step;

    //psi_start = 0.0f;
    //psi_max   = 360.0f;
    // out-of-plane PHI+THETA
    //global_euler_search.InitGrid(my_symmetry, angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
    if ( do_constrained_search )
        global_euler_search.InitConstrainedGrid(my_symmetry, angular_step, phi_start, phi_max, theta_start, theta_max, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);
    else
        global_euler_search.InitGrid(my_symmetry, angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);

    if ( my_symmetry.StartsWith("C") ) // TODO 2x check me - w/o this O symm at least is broken
    {
        if ( global_euler_search.test_mirror == true ) // otherwise the theta max is set to 90.0 and test_mirror is set to true.  However, I don't want to have to test the mirrors.
        {
            if ( ! do_constrained_search )
                global_euler_search.theta_max = 180.0f;
        }
    }

    if ( ! do_constrained_search )
        global_euler_search.CalculateGridSearchPositions(false);
    else
        global_euler_search.CalculateConstrainedGridSearchPositions(false);

    total_correlation_positions  = 0;
    current_correlation_position = 0;
    if ( is_running_locally == true ) {
        wxPrintf("running locally...\n");
        first_search_position = 0;
        last_search_position  = global_euler_search.number_of_search_positions - 1;
    }

    for ( current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
        //loop over each rotation angle

        for ( current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step ) {
            total_correlation_positions++;
        }
    }

    number_of_rotations = 0;
    for ( current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step ) {
        number_of_rotations++;
    }

    ProgressBar* my_progress;
    if ( is_running_locally == true ) {
        my_progress = new ProgressBar(total_correlation_positions);
    }

    wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
    wxPrintf("Searching %i rotations per position.\n", number_of_rotations);
    wxPrintf("There are %li correlation positions total.\n\n", total_correlation_positions);

    wxPrintf("Performing Search...\n\n");

    // start rotational search
    for ( current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
        for ( current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step ) {
            angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], current_psi, 0.0, 0.0);
            //angles.Init(fixed_phi, fixed_theta, fixed_psi, 0.0, 0.0); // 91 clathrin
            //angles.Init(-132.69, 134.20, 181.39, 0, 0); // 1 mature60S

            if ( padding != 1.0f ) {
                template_reconstruction.ExtractSlice(padded_projection, angles, 1.0f, false);
                padded_projection.SwapRealSpaceQuadrants( );
                padded_projection.BackwardFFT( );
                padded_projection.ClipInto(&current_projection);
                current_projection.ForwardFFT( );
            }
            else {
                template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
                current_projection.SwapRealSpaceQuadrants( );
            }

            current_projection.MultiplyPixelWise(projection_filter); // only the ctf

            if ( whitening_option == 2 ) {
                ///*
                // apply ref's own whitening filter
                current_projection.Compute1DPowerSpectrumCurve(&whitening_filter_ref, &number_of_terms_ref);
                whitening_filter_ref.SquareRoot( );
                whitening_filter_ref.Reciprocal( );
                whitening_filter_ref.MultiplyByConstant(1.0f / whitening_filter_ref.ReturnMaximumValue( ));
                //whitening_filter_ref.WriteToFile("w_filter_ref.txt");

                current_projection.ApplyCurveFilter(&whitening_filter_ref);
                current_projection.ZeroCentralPixel( );
                current_projection.DivideByConstant(sqrtf(current_projection.ReturnSumOfSquares( )));
                //*/
            }

            current_projection.BackwardFFT( );

            current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges( ));
            variance = current_projection.ReturnSumOfSquares( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels - powf(current_projection.ReturnAverageOfRealValues( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
            current_projection.DivideByConstant(sqrtf(variance));
            current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

            padded_reference.ForwardFFT( );
            padded_reference.ZeroCentralPixel( );

            //#ifdef MKL
            // Use the MKL
            //vmcMulByConj(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(input_image.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);
            //#else
            //NumericTextFile amp_file("amplitude.txt", OPEN_TO_WRITE, 1);

            //amp_file.WriteCommentLine("amplitude");
            //double temp_double[1];

            for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
            }

            if ( do_phase_matched_filter ) {
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                    float amplitude = abs(padded_reference.complex_values[pixel_counter]);

                    if ( amplitude == 0.0f )
                        amplitude = 1.0; //0.000001f;
                    //amp_file.WriteLine(temp_double);

                    padded_reference.complex_values[pixel_counter] /= amplitude;
                }
            }

            //amp_file.Close( );
            //#endif

            padded_reference.BackwardFFT( );
            //padded_reference.QuickAndDirtyWriteSlice("cc.mrc", 1);
            //exit(0);

            // update mip, and histogram..
            pixel_counter = 0;

            for ( current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++ ) {
                for ( current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++ ) {
                    // first mip

                    if ( padded_reference.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter] ) {
                        max_intensity_projection.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
                        //best_psi.real_values[pixel_counter]                 = current_psi;
                        //best_theta.real_values[pixel_counter]               = global_euler_search.list_of_search_parameters[current_search_position][1];
                        //best_phi.real_values[pixel_counter]                 = global_euler_search.list_of_search_parameters[current_search_position][0];
                        //best_defocus.real_values[pixel_counter]             = float(defocus_i) * defocus_step;
                        //best_pixel_size.real_values[pixel_counter]          = float(size_i) * pixel_size_step;
                    }

                    // histogram

                    current_bin = int(double((padded_reference.real_values[pixel_counter]) - histogram_min_scaled) / histogram_step_scaled);

                    if ( current_bin >= 0 && current_bin <= histogram_number_of_points ) {
                        histogram_data[current_bin] += 1;
                    }

                    pixel_counter++;
                }

                pixel_counter += padded_reference.padding_jump_value;
            }

            /*
            for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                correlation_pixel_sum[pixel_counter] += padded_reference.real_values[pixel_counter];
            }
            padded_reference.SquareRealValues( );
            for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                correlation_pixel_sum_of_squares[pixel_counter] += padded_reference.real_values[pixel_counter];
            }
			*/
            current_projection.is_in_real_space = false;
            padded_reference.is_in_real_space   = true;

            current_correlation_position++;
            if ( is_running_locally == true )
                my_progress->Update(current_correlation_position);

        } // end of IP search
    } // end of OOP search

    // manually add in correct orientation
    current_projection.is_in_real_space = false;
    padded_reference.is_in_real_space   = true;
    angles.Init(fixed_phi, fixed_theta, fixed_psi, 0.0, 0.0);
    if ( padding != 1.0f ) {
        template_reconstruction.ExtractSlice(padded_projection, angles, 1.0f, false);
        padded_projection.SwapRealSpaceQuadrants( );
        padded_projection.BackwardFFT( );
        padded_projection.ClipInto(&current_projection);
        current_projection.ForwardFFT( );
    }
    else {
        template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
        current_projection.SwapRealSpaceQuadrants( );
    }

    current_projection.MultiplyPixelWise(projection_filter);

    if ( whitening_option == 2 ) {
        ///*
        // apply ref's own whitening filter
        current_projection.Compute1DPowerSpectrumCurve(&whitening_filter_ref, &number_of_terms_ref);
        whitening_filter_ref.SquareRoot( );
        whitening_filter_ref.Reciprocal( );
        whitening_filter_ref.MultiplyByConstant(1.0f / whitening_filter_ref.ReturnMaximumValue( ));
        //whitening_filter_ref.WriteToFile("w_filter_ref.txt");

        current_projection.ApplyCurveFilter(&whitening_filter_ref);
        current_projection.ZeroCentralPixel( );
        current_projection.DivideByConstant(sqrtf(current_projection.ReturnSumOfSquares( )));
        //*/
    }

    current_projection.BackwardFFT( );

    current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges( ));
    variance = current_projection.ReturnSumOfSquares( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels - powf(current_projection.ReturnAverageOfRealValues( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
    current_projection.DivideByConstant(sqrtf(variance));
    current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

    padded_reference.ForwardFFT( );
    padded_reference.ZeroCentralPixel( );

    for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
        padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
    }

    if ( do_phase_matched_filter ) {
        for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
            float amplitude = abs(padded_reference.complex_values[pixel_counter]);

            if ( amplitude == 0.0f )
                amplitude = 1.0; //0.000001f;
            //amp_file.WriteLine(temp_double);

            padded_reference.complex_values[pixel_counter] /= amplitude;
        }
    }

    ////// update mip, and histogram..
    padded_reference.BackwardFFT( );
    padded_reference.QuickAndDirtyWriteSlice("correct_pose_cc.mrc", 1);

    pixel_counter = 0;

    for ( current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++ ) {
        for ( current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++ ) {
            // first mip

            if ( padded_reference.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter] ) {
                max_intensity_projection.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
                //best_psi.real_values[pixel_counter]                 = current_psi;
                //best_theta.real_values[pixel_counter]               = global_euler_search.list_of_search_parameters[current_search_position][1];
                //best_phi.real_values[pixel_counter]                 = global_euler_search.list_of_search_parameters[current_search_position][0];
                //best_defocus.real_values[pixel_counter]             = float(defocus_i) * defocus_step;
                //best_pixel_size.real_values[pixel_counter]          = float(size_i) * pixel_size_step;
            }

            pixel_counter++;
        }

        pixel_counter += padded_reference.padding_jump_value;
    }
    //////

    long   number_of_result_floats = 0;
    float* pointer_to_histogram_data;
    pointer_to_histogram_data = (float*)histogram_data;

    number_of_result_floats = histogram_number_of_points * sizeof(long) / sizeof(float);
    float* result           = new float[number_of_result_floats];

    for ( pixel_counter = 0; pixel_counter < histogram_number_of_points * 2; pixel_counter++ ) {
        result[pixel_counter] = pointer_to_histogram_data[pixel_counter];
    }

    // write out histogram
    float temp_float;
    temp_float = histogram_min + (histogram_step / 2.0f); // start position
    NumericTextFile histogram_file(wxString::Format("%s.txt", output_histogram_file).ToStdString( ), OPEN_TO_WRITE, 2);

    histogram_file.WriteCommentLine("SNR, histogram");

    double temp_double_array[2];

    for ( int line_counter = 0; line_counter < histogram_number_of_points; line_counter++ ) {
        temp_double_array[0] = temp_float + histogram_step * float(line_counter);
        temp_double_array[1] = histogram_data[line_counter];
        histogram_file.WriteLine(temp_double_array);
    }
    histogram_file.Close( );

    if ( is_running_locally == true ) {
        delete my_progress;

        // scale images..

        for ( pixel_counter = 0; pixel_counter < input_image.real_memory_allocated; pixel_counter++ ) {
            correlation_pixel_sum[pixel_counter] /= float(total_correlation_positions);
            correlation_pixel_sum_of_squares[pixel_counter] = correlation_pixel_sum_of_squares[pixel_counter] / float(total_correlation_positions) - powf(correlation_pixel_sum[pixel_counter], 2);
            if ( correlation_pixel_sum_of_squares[pixel_counter] > 0.0f ) {
                correlation_pixel_sum_of_squares[pixel_counter] = sqrtf(correlation_pixel_sum_of_squares[pixel_counter]) * (float)sqrt_input_pixels;
            }
            else
                correlation_pixel_sum_of_squares[pixel_counter] = 0.0f;
            correlation_pixel_sum[pixel_counter] *= (float)sqrt_input_pixels;
        }

        max_intensity_projection.MultiplyByConstant((float)sqrt_input_pixels);
        max_intensity_projection.QuickAndDirtyWriteSlice(output_mip_file.ToStdString( ), true);
    }
    return true;
}
