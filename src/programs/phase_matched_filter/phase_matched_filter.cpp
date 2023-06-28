#include "../../core/core_headers.h"
#include "../../constants/constants.h"

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

    wxString mip_output_file;
    wxString best_psi_output_file;
    wxString best_theta_output_file;
    wxString best_phi_output_file;
    wxString scaled_mip_output_file;

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

    bool use_gpu_input = false;
    int  max_threads   = 1;

    bool  do_constrained_search = false;
    float phi_start             = 0.0;
    float phi_max               = 360.0;
    float theta_start           = 0.0;
    float theta_max             = 180.0;
    float psi_start             = 0.0;
    float psi_max               = 360.0;
    float fixed_phi             = 0.0;
    float fixed_theta           = 0.0;
    float fixed_psi             = 0.0;
    bool  do_whiten_ref         = 1;
    bool  do_whiten_image       = true;
    bool  do_phase_only         = true;

    UserInput* my_input = new UserInput("PhaseMatchedFilter", 1.00);

    input_search_images         = my_input->GetFilenameFromUser("Input images to be searched", "The input image stack, containing the images that should be searched", "image_stack.mrc", true);
    input_reconstruction        = my_input->GetFilenameFromUser("Input template reconstruction", "The 3D reconstruction from which projections are calculated", "reconstruction.mrc", true);
    output_histogram_file       = my_input->GetFilenameFromUser("Output histogram of correlation values", "histogram of all correlation values", "histogram.txt", false);
    mip_output_file             = my_input->GetFilenameFromUser("Output MIP file", "The file for saving the maximum intensity projection image", "mip.mrc", false);
    scaled_mip_output_file      = my_input->GetFilenameFromUser("Output Scaled MIP file", "The file for saving the maximum intensity projection image divided by correlation variance", "mip_scaled.mrc", false);
    best_psi_output_file        = my_input->GetFilenameFromUser("Output psi file", "The file for saving the best psi image", "psi.mrc", false);
    best_theta_output_file      = my_input->GetFilenameFromUser("Output theta file", "The file for saving the best psi image", "theta.mrc", false);
    best_phi_output_file        = my_input->GetFilenameFromUser("Output phi file", "The file for saving the best psi image", "phi.mrc", false);
    correlation_avg_output_file = my_input->GetFilenameFromUser("Correlation average value", "The file for saving the average value of all correlation images", "corr_average.mrc", false);
    correlation_std_output_file = my_input->GetFilenameFromUser("Correlation SD output file", "The file for saving the standard deviation of all correlation images", "corr_standard_deviation.mrc", false);
    pixel_size                  = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    voltage_kV                  = my_input->GetFloatFromUser("Beam energy (keV)", "The energy of the electron beam used to image the sample in kilo electron volts", "300.0", 0.0);
    spherical_aberration_mm     = my_input->GetFloatFromUser("Spherical aberration (mm)", "Spherical aberration of the objective lens in millimeters", "2.7");
    amplitude_contrast          = my_input->GetFloatFromUser("Amplitude contrast", "Assumed amplitude contrast", "0.07", 0.0, 1.0);
    defocus1                    = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    defocus2                    = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    defocus_angle               = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    phase_shift                 = my_input->GetFloatFromUser("Phase Shift (degrees)", "Additional phase shift in degrees", "0.0");
    high_resolution_limit       = my_input->GetFloatFromUser("High resolution limit (A)", "High resolution limit of the data used for alignment in Angstroms", "8.0", 0.0);
    angular_step                = my_input->GetFloatFromUser("Out of plane angular step (0.0 = set automatically)", "Angular step size for global grid search", "0.0", 0.0);
    in_plane_angular_step       = my_input->GetFloatFromUser("In plane angular step (0.0 = set automatically)", "Angular step size for in-plane rotations during the search", "0.0", 0.0);
    do_constrained_search       = my_input->GetYesNoFromUser("Perform constraied angular search", "yes no", "no");
    if ( do_constrained_search ) {
        phi_max     = my_input->GetFloatFromUser("Constrained search Phi max", "max phi for constrained search", "360.0", 0.0, 360.0);
        phi_start   = my_input->GetFloatFromUser("Constrained search Phi min", "min phi for constrained search", "0.0", 0.0, 360.0);
        theta_max   = my_input->GetFloatFromUser("Constrained search Theta max", "max theta for constrained search", "180.0", 0.0, 180.0);
        theta_start = my_input->GetFloatFromUser("Constrained search Theta min", "max theta for constrained search", "0.0", 0.0, 180.0);
        psi_max     = my_input->GetFloatFromUser("Constrained search Psi max", "max psi for constrained search", "360.0", 0.0, 360.0);
        psi_start   = my_input->GetFloatFromUser("Constrained search Psi min", "max psi for constrained search", "0.0.0", 0.0, 360.0);
        fixed_phi   = my_input->GetFloatFromUser("Constrained search Phi", "max phi for constrained search", "360.0", 0.0, 360.0);
        fixed_theta = my_input->GetFloatFromUser("Constrained search Theta", "max theta for constrained search", "360.0", 0.0, 360.0);
        fixed_psi   = my_input->GetFloatFromUser("Constrained search Psi", "max psi for constrained search", "360.0", 0.0, 360.0);
    }
    //    best_parameters_to_keep = my_input->GetIntFromUser("Number of top hits to refine", "The number of best global search orientations to refine locally", "20", 1);
    padding                   = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the input volume is padded to improve projections", "1.0", 1.0, 2.0);
    particle_radius_angstroms = my_input->GetFloatFromUser("Mask radius for global search (A) (0.0 = max)", "Radius of a circular mask to be applied to the input images during global search", "0.0", 0.0);
    my_symmetry               = my_input->GetSymmetryFromUser("Template symmetry", "The symmetry of the template reconstruction", "C1");
    do_whiten_ref             = my_input->GetYesNoFromUser("Whiten reference?", "yes no", "yes");
    do_whiten_image           = my_input->GetYesNoFromUser("Whiten image?", "yes no", "yes");
    do_phase_only             = my_input->GetYesNoFromUser("do_phase_only?", "yes no", "yes");
#ifdef ENABLEGPU
    use_gpu_input = my_input->GetYesNoFromUser("Use GPU", "Offload expensive calcs to GPU", "No");
#endif
    max_threads = my_input->GetIntFromUser("Max. threads to use for calculation", "when threading, what is the max threads to run", "1", 1);

    int   first_search_position = -1;
    int   last_search_position  = -1;
    float min_peak_radius       = 10.0f;

    delete my_input;

    my_current_job.ManualSetArguments("ttffffffffffiffffffftftiifbfffffffffbbbtttttttbi", input_search_images.ToUTF8( ).data( ),
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
                                      do_whiten_ref,
                                      do_whiten_image,
                                      do_phase_only,
                                      mip_output_file.ToUTF8( ).data( ),
                                      scaled_mip_output_file.ToUTF8( ).data( ),
                                      best_psi_output_file.ToUTF8( ).data( ),
                                      best_theta_output_file.ToUTF8( ).data( ),
                                      best_phi_output_file.ToUTF8( ).data( ),
                                      correlation_avg_output_file.ToUTF8( ).data( ),
                                      correlation_std_output_file.ToUTF8( ).data( ),
                                      use_gpu_input,
                                      max_threads);
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
    int      do_whiten_ref                 = my_current_job.arguments[36].ReturnBoolArgument( );
    bool     do_whiten_image               = my_current_job.arguments[37].ReturnBoolArgument( );
    bool     do_phase_only                 = my_current_job.arguments[38].ReturnBoolArgument( );
    wxString mip_output_file               = my_current_job.arguments[39].ReturnStringArgument( );
    wxString scaled_mip_output_file        = my_current_job.arguments[40].ReturnStringArgument( );
    wxString best_psi_output_file          = my_current_job.arguments[41].ReturnStringArgument( );
    wxString best_theta_output_file        = my_current_job.arguments[42].ReturnStringArgument( );
    wxString best_phi_output_file          = my_current_job.arguments[43].ReturnStringArgument( );
    wxString correlation_avg_output_file   = my_current_job.arguments[44].ReturnStringArgument( );
    wxString correlation_std_output_file   = my_current_job.arguments[45].ReturnStringArgument( );
    bool     use_gpu                       = my_current_job.arguments[46].ReturnBoolArgument( );
    int      max_threads                   = my_current_job.arguments[47].ReturnIntegerArgument( );

    // This condition applies to GUI and CLI - it is just a recommendation to the user.
    if ( use_gpu && max_threads <= 1 ) {
        SendInfo("Warning, you are only using one thread on the GPU. Suggested minimum is 2. Check compute saturation using nvidia-smi -l 1\n");
    }
    if ( ! use_gpu ) {
        SendInfo("GPU disabled\nCan use up to 44 threads on frankfurt\n.");
    }

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

    Image best_psi;
    Image best_theta;
    Image best_phi;

    input_image.ReadSlice(&input_search_image_file, 1);
    input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices( ));
    template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
    projection_filter.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
    if ( use_gpu ) {
        current_projection.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
        padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
        padded_reference.SetToConstant(0.0f);
    }
    max_intensity_projection.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    best_psi.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    best_theta.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    best_phi.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);

    correlation_pixel_sum_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    correlation_pixel_sum_of_squares_image.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
    double* correlation_pixel_sum            = new double[input_image.real_memory_allocated];
    double* correlation_pixel_sum_of_squares = new double[input_image.real_memory_allocated];

    max_intensity_projection.SetToConstant(-FLT_MAX);
    best_psi.SetToConstant(0.0f);
    best_theta.SetToConstant(0.0f);
    best_phi.SetToConstant(0.0f);

    ZeroDoubleArray(correlation_pixel_sum, input_image.real_memory_allocated);
    ZeroDoubleArray(correlation_pixel_sum_of_squares, input_image.real_memory_allocated);

    if ( padding != 1.0f ) {
        input_reconstruction.Resize(input_reconstruction.logical_x_dimension * padding, input_reconstruction.logical_y_dimension * padding, input_reconstruction.logical_z_dimension * padding, input_reconstruction.ReturnAverageOfRealValuesOnEdges( ));
        if ( use_gpu )
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
    float amplitude;
    long  pixel_counter;
    Image temp_image;

    Curve whitening_filter;
    Curve number_of_terms;
    whitening_filter.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

    input_image.ReplaceOutliersWithMean(5.0f);
    input_image.ForwardFFT( );
    input_image.SwapRealSpaceQuadrants( );

    input_image.ZeroCentralPixel( );
    input_image.Compute1DPowerSpectrumCurve(&whitening_filter, &number_of_terms);
    whitening_filter.SquareRoot( );
    whitening_filter.Reciprocal( );
    whitening_filter.MultiplyByConstant(1.0f / whitening_filter.ReturnMaximumValue( ));
    //whitening_filter.WriteToFile("w_filter.txt");
    if ( do_whiten_image ) {
        input_image.ApplyCurveFilter(&whitening_filter);
    }

    input_image.ZeroCentralPixel( );
    input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares( )));

    // set up CTF
    CTF input_ctf;
    input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));
    input_ctf.SetDefocus(defocus1 / pixel_size, defocus2 / pixel_size, deg_2_rad(defocus_angle));
    projection_filter.CalculateCTFImage(input_ctf);

    if ( do_whiten_ref == 1 )
        projection_filter.ApplyCurveFilter(&whitening_filter);

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

    wxDateTime overall_start;
    wxDateTime overall_finish;
    overall_start = wxDateTime::Now( );

    // GPU code
    float actual_number_of_ccs_calculated = 0.0;

    bool first_gpu_loop = true;
    int  nThreads       = 2;
    int  nGPUs          = 1;
    int  nJobs          = last_search_position - first_search_position + 1;
    if ( use_gpu && max_threads > nJobs ) {
        wxPrintf("\n\tWarning, you request more threads (%d) than there are search positions (%d)\n", max_threads, nJobs);
        max_threads = nJobs;
    }

    int minPos = first_search_position;
    int maxPos = last_search_position;
    int incPos = (nJobs) / (max_threads);

#ifdef ENABLEGPU
    TemplateMatchingCore* GPU;
    DeviceManager         gpuDev;
#endif

    if ( use_gpu ) {
        total_correlation_positions_per_thread = total_correlation_positions / max_threads;

#ifdef ENABLEGPU
        //    checkCudaErrors(cudaGetDeviceCount(&nGPUs));
        GPU = new TemplateMatchingCore[max_threads];
        gpuDev.Init(nGPUs);

//    wxPrintf("Host: %s is running\nnThreads: %d\nnGPUs: %d\n:nSearchPos %d \n",hostNameBuffer,nThreads, nGPUs, maxPos);

//    TemplateMatchingCore GPU(number_of_jobs_per_image_in_gui);
#endif
    }

    if ( use_gpu ) {
#ifdef ENABLEGPU

#pragma omp parallel num_threads(max_threads)
        {
            int tIDX = ReturnThreadNumberOfCurrentThread( );
            gpuDev.SetGpu(tIDX);

            if ( first_gpu_loop ) {

                int t_first_search_position = first_search_position + (tIDX * incPos);
                int t_last_search_position  = first_search_position + (incPos - 1) + (tIDX * incPos);

                if ( tIDX == (max_threads - 1) )
                    t_last_search_position = maxPos;

                GPU[tIDX].Init(this, template_reconstruction, input_image, current_projection,
                               pixel_size_search_range, pixel_size_step, pixel_size,
                               defocus_search_range, defocus_step, defocus1, defocus2,
                               psi_max, psi_start, psi_step,
                               angles, global_euler_search,
                               float(histogram_min_scaled), float(histogram_step_scaled), histogram_number_of_points,
                               0, t_first_search_position, t_last_search_position,
                               my_progress, total_correlation_positions_per_thread, is_running_locally);

                wxPrintf("%d\n", tIDX);
                wxPrintf("%d\n", t_first_search_position);
                wxPrintf("%d\n", t_last_search_position);
                wxPrintf("Staring TemplateMatchingCore object %d to work on position range %d-%d\n", tIDX, t_first_search_position, t_last_search_position);

                first_gpu_loop = false;
            }
            else {
                GPU[tIDX].template_reconstruction.CopyFrom(&template_reconstruction);
            }
        } // end of omp block
#endif
    }

    if ( use_gpu ) {
#ifdef ENABLEGPU
        //            wxPrintf("\n\n\t\tsizeI defI %d %d\n\n\n", size_i, defocus_i);

#pragma omp parallel num_threads(max_threads)
        {
            int tIDX = ReturnThreadNumberOfCurrentThread( );
            gpuDev.SetGpu(tIDX);

            GPU[tIDX].RunInnerLoop(projection_filter, 0.0, 0.0, tIDX, current_correlation_position, true); // no pixel size or defocus search

#pragma omp critical
            {

                Image mip_buffer;
                mip_buffer.CopyFrom(&max_intensity_projection);
                Image psi_buffer;
                psi_buffer.CopyFrom(&max_intensity_projection);
                Image phi_buffer;
                phi_buffer.CopyFrom(&max_intensity_projection);
                Image theta_buffer;
                theta_buffer.CopyFrom(&max_intensity_projection);

                GPU[tIDX].d_max_intensity_projection.CopyDeviceToHost(mip_buffer, true, false);
                GPU[tIDX].d_best_psi.CopyDeviceToHost(psi_buffer, true, false);
                GPU[tIDX].d_best_phi.CopyDeviceToHost(phi_buffer, true, false);
                GPU[tIDX].d_best_theta.CopyDeviceToHost(theta_buffer, true, false);

                // mip_buffer.QuickAndDirtyWriteSlice("tmpMipBuffer.mrc", 1, 1);
                // TODO should prob aggregate these across all workers
                // TODO add a copySum method that allocates a pinned buffer, copies there then sumes into the wanted image.
                Image sum;
                Image sumSq;

                sum.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
                sumSq.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);

                sum.SetToConstant(0.0f);
                sumSq.SetToConstant(0.0f);

                GPU[tIDX].d_sum3.CopyDeviceToHost(sum, true, false);
                GPU[tIDX].d_sumSq3.CopyDeviceToHost(sumSq, true, false);

                GPU[tIDX].d_max_intensity_projection.Wait( );

                // TODO swap max_padding for explicit padding in x/y and limit calcs to that region.
                pixel_counter = 0;
                for ( current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++ ) {
                    for ( current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++ ) {
                        // first mip

                        if ( mip_buffer.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter] ) {
                            max_intensity_projection.real_values[pixel_counter] = mip_buffer.real_values[pixel_counter];
                            best_psi.real_values[pixel_counter]                 = psi_buffer.real_values[pixel_counter];
                            best_theta.real_values[pixel_counter]               = theta_buffer.real_values[pixel_counter];
                            best_phi.real_values[pixel_counter]                 = phi_buffer.real_values[pixel_counter];
                            // best_defocus.real_values[pixel_counter]             = float(defocus_i) * defocus_step;
                            // best_pixel_size.real_values[pixel_counter]          = float(size_i) * pixel_size_step;
                        }

                        correlation_pixel_sum[pixel_counter] += (double)sum.real_values[pixel_counter];
                        correlation_pixel_sum_of_squares[pixel_counter] += (double)sumSq.real_values[pixel_counter];

                        pixel_counter++;
                    }

                    pixel_counter += max_intensity_projection.padding_jump_value;
                }

                GPU[tIDX].histogram.CopyToHostAndAdd(histogram_data);

                //                    current_correlation_position += GPU[tIDX].total_number_of_cccs_calculated;
                actual_number_of_ccs_calculated += GPU[tIDX].total_number_of_cccs_calculated;

            } // end of omp critical block
        } // end of parallel block

        //continue;

#endif
    }
    else {
// start rotational search
// openmp
#pragma omp parallel for num_threads(max_threads) default(shared) private(current_search_position, angles, current_psi, variance, pixel_counter, current_x, current_y, current_projection, padded_projection, amplitude, temp_image, padded_reference)
        for ( current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
            // wxPrintf("thread id = %i position %i\n", ReturnThreadNumberOfCurrentThread( ), current_search_position);
            // private variables for each thread to update sum and sumSq of ccs
            long    current_correlation_position_private = 0;
            double* correlation_pixel_sum_private        = new double[input_image.real_memory_allocated];
            ZeroDoubleArray(correlation_pixel_sum_private, input_image.real_memory_allocated);
            double* correlation_pixel_sum_of_squares_private = new double[input_image.real_memory_allocated];
            ZeroDoubleArray(correlation_pixel_sum_of_squares_private, input_image.real_memory_allocated);
            // test mkl
            float* amplitude_array = new float[input_image.real_memory_allocated];
            ZeroFloatArray(amplitude_array, input_image.real_memory_allocated);

            for ( current_psi = psi_start; current_psi <= psi_max; current_psi += psi_step ) {
                current_projection.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
                padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
                padded_reference.SetToConstant(0.0f);
                padded_projection.Allocate(input_reconstruction_file.ReturnXSize( ) * padding, input_reconstruction_file.ReturnXSize( ) * padding, false);

                angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], current_psi, 0.0, 0.0);

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

                current_projection.BackwardFFT( );

                current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges( ));
                variance = current_projection.ReturnSumOfSquares( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels - powf(current_projection.ReturnAverageOfRealValues( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
                current_projection.DivideByConstant(sqrtf(variance));
                current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

                padded_reference.ForwardFFT( );
                padded_reference.ZeroCentralPixel( );

#ifdef MKL
                // Use the MKL
                vmcMulByConj(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(input_image.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);
                vmcAbs(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), amplitude_array, VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                    if ( amplitude_array[pixel_counter] == 0.0 )
                        amplitude_array[pixel_counter] = 0.000001f;
                }
                // vmcDiv(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), amplitude_array, reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);

                // amp_file.WriteCommentLine("amplitude");
                // float temp_double[1];
                // for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                //     temp_double[0] = amplitude_array[pixel_counter];
                //     amp_file.WriteLine(temp_double);
                // }
                // amp_file.Close( );
                //exit(0);

#else

                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                    padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
                }

                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                    amplitude = abs(padded_reference.complex_values[pixel_counter]);
                    // wxPrintf("amplitude of cc = %f\n", amplitude);

                    if ( amplitude == 0.0f )
                        amplitude = 0.000001f;
                }

#endif
                // TODO mkl for element wise complex / float
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                    padded_reference.complex_values[pixel_counter] /= amplitude_array[pixel_counter];
                }
                padded_reference.BackwardFFT( );
                //padded_reference.QuickAndDirtyWriteSlice("cc.mrc", 1);
                //exit(0);
                current_correlation_position_private++;
// update mip, and histogram..
#pragma omp critical
                {
                    pixel_counter = 0;

                    for ( current_y = 0; current_y < max_intensity_projection.logical_y_dimension; current_y++ ) {
                        for ( current_x = 0; current_x < max_intensity_projection.logical_x_dimension; current_x++ ) {
                            // first mip

                            if ( padded_reference.real_values[pixel_counter] > max_intensity_projection.real_values[pixel_counter] ) {
                                max_intensity_projection.real_values[pixel_counter] = padded_reference.real_values[pixel_counter];
                                best_psi.real_values[pixel_counter]                 = current_psi;
                                best_theta.real_values[pixel_counter]               = global_euler_search.list_of_search_parameters[current_search_position][1];
                                best_phi.real_values[pixel_counter]                 = global_euler_search.list_of_search_parameters[current_search_position][0];
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

                } // end of critical block

                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                    correlation_pixel_sum_private[pixel_counter] += padded_reference.real_values[pixel_counter];
                }
                padded_reference.SquareRealValues( );
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                    correlation_pixel_sum_of_squares_private[pixel_counter] += padded_reference.real_values[pixel_counter];
                }

                current_projection.is_in_real_space = false;
                padded_reference.is_in_real_space   = true;

                //current_correlation_position = current_correlation_position + 1;

            } // end of IP search
#pragma omp critical
            {
                current_correlation_position += current_correlation_position_private;
                for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated; pixel_counter++ ) {
                    correlation_pixel_sum[pixel_counter] += correlation_pixel_sum_private[pixel_counter];
                    correlation_pixel_sum_of_squares[pixel_counter] += correlation_pixel_sum_of_squares_private[pixel_counter];
                }
                if ( is_running_locally == true )
                    my_progress->Update(current_correlation_position);
            }
        } // end of OOP search
        // wxPrintf("end of OOP\n");
        // end of omp block

    } // end of else (not gpu)

    wxPrintf("\n\n\tTimings: Overall: %s\n", (wxDateTime::Now( ) - overall_start).Format( ));

    // calculate avg and std
    for ( pixel_counter = 0; pixel_counter < input_image.real_memory_allocated; pixel_counter++ ) {
        correlation_pixel_sum_image.real_values[pixel_counter]            = (float)correlation_pixel_sum[pixel_counter];
        correlation_pixel_sum_of_squares_image.real_values[pixel_counter] = (float)correlation_pixel_sum_of_squares[pixel_counter];
    }
    if ( is_running_locally == true ) {
        // scaling
        delete my_progress;
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
    }
    max_intensity_projection.MultiplyByConstant((float)sqrt_input_pixels);
    max_intensity_projection.QuickAndDirtyWriteSlice(mip_output_file.ToStdString( ), true);

    for ( pixel_counter = 0; pixel_counter < input_image.real_memory_allocated; pixel_counter++ ) {
        max_intensity_projection.real_values[pixel_counter] -= correlation_pixel_sum[pixel_counter];
        if ( correlation_pixel_sum_of_squares[pixel_counter] > 0.0f ) {
            max_intensity_projection.real_values[pixel_counter] /= correlation_pixel_sum_of_squares[pixel_counter];
        }
        else
            max_intensity_projection.real_values[pixel_counter] = 0.0f;
        correlation_pixel_sum_image.real_values[pixel_counter]            = correlation_pixel_sum[pixel_counter];
        correlation_pixel_sum_of_squares_image.real_values[pixel_counter] = correlation_pixel_sum_of_squares[pixel_counter];
    }

    max_intensity_projection.QuickAndDirtyWriteSlice(scaled_mip_output_file.ToStdString( ), 1, pixel_size);
    correlation_pixel_sum_image.QuickAndDirtyWriteSlice(correlation_avg_output_file.ToStdString( ), 1, pixel_size);
    correlation_pixel_sum_of_squares_image.QuickAndDirtyWriteSlice(correlation_std_output_file.ToStdString( ), 1, pixel_size);
    best_psi.QuickAndDirtyWriteSlice(best_psi_output_file.ToStdString( ), 1, pixel_size);
    best_theta.QuickAndDirtyWriteSlice(best_theta_output_file.ToStdString( ), 1, pixel_size);
    best_phi.QuickAndDirtyWriteSlice(best_phi_output_file.ToStdString( ), 1, pixel_size);

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
    NumericTextFile histogram_file(output_histogram_file.ToStdString( ), OPEN_TO_WRITE, 2);

    histogram_file.WriteCommentLine("SNR, histogram");

    double temp_double_array[2];

    for ( int line_counter = 0; line_counter < histogram_number_of_points; line_counter++ ) {
        temp_double_array[0] = temp_float + histogram_step * float(line_counter);
        temp_double_array[1] = histogram_data[line_counter];
        histogram_file.WriteLine(temp_double_array);
    }
    histogram_file.Close( );

    return true;
}
