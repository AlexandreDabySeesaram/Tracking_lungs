
import numpy as np
import os

from create_vti_from_scans import get_metada_PGM, pgm2array, array2vti, get_z_metadata_flatten_image


N_patients      = 9
Patients_Ids    = list(range(1,  N_patients + 1))
# Patients_Ids.remove(1)

# Patients_Ids = [1]

initial_path = '/Users/daby/LargeFiles/Geometries_alexandre/Patient_data/Images/'

destination_path = "./"

RAW_vti                             = False
Vti_Binary_PGM                      = False         # check +- 1 on pixels
Vti_RAW_PGM                         = False
Blur_binary_vti                     = False
Clean_BINARY_PGM                    = False
Copy_initial_image                  = False
Downscale_vti                       = False
Threshold_blurred_images            = True
original_dir = os.getcwd()                          # remember original path



if RAW_vti:
    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        folder                      = prefix+"_images/"
        file_init                   = Folder+prefix+"_M0_RAW_1.vti"
        new_folder                  = destination_path+prefix
        new_file                    = new_folder+ "/Image_01.vti"
        new_file_ref                = new_folder+ "/Image_00.vti"
        command_mkdir               = "mkdir "+new_folder
        command_cp                  = "cp "+initial_path+file_init+" "+new_file
        os.system(command_cp)
        command_cp_ref              =  "cp E2/Image_P_00.vti "+new_file_ref
        os.system(command_cp_ref)

if Clean_BINARY_PGM:
    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        # Gif
        post_gif                    = prefix+"/LUNG/"+"Pat"+str(i)+"_inspi?_????_colorcls.gif"
        command_del                 = "rm "+post_gif
        os.system(command_del)

        # textured binary
        post_textured               = prefix+"/LUNG/"+"Pat"+str(i)+"_inspi?_????_predcls.pgm"
        command_del                 = "rm "+post_textured
        os.system(command_del)

        os.chdir(prefix+"/LUNG")
        for filename in os.listdir():
            
            if filename.endswith("_lung.pgm"):
                # Construct new filename by removing "_lung"
                new_filename        = filename.replace("_lung", "")
                # Rename the file
                os.rename(filename, new_filename)
                print(f"Renamed: {filename} -> {new_filename}")
        os.chdir(original_dir)



# /!\ Need to read metadata from raw images

if Vti_Binary_PGM:
    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        input_files                 = prefix+"/LUNG/"+"Pat"+str(i)+"_inspi1"
        input_files_raw             = prefix+"/PGM/"+"Pat"+str(i)+"_inspi1"             # For metadata only here
        output = prefix+"/"+prefix+"_Binary"
        pgm_files, image_array      = pgm2array(input_files)
        pgm_files_raw, _            = pgm2array(input_files_raw)
        metadata_fields             = ["Slice_Location", "Pixel_Size"]
        metadata                    = get_metada_PGM(
                                                        input_file        = pgm_files_raw[0],
                                                        metadata_fields   = metadata_fields
                                                    )
        image_shape, pixel_size, image_pos, flatten_image_array = get_z_metadata_flatten_image(pgm_files_raw, metadata, image_array)
        array2vti(
                image_shape         = image_shape,
                pixel_size          = pixel_size,
                image_pos           = image_pos,
                input_array         = flatten_image_array,
                field_name          = 'pixel intensity',
                output              = output)
        print("Patient "+prefix+"done")



if Blur_binary_vti:
    import vtk 

    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        input_file                  = prefix+"/"+prefix+"_Binary.vti"
        # Read the .vti file
        reader                      = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_file)
        reader.Update()

        # Apply Gaussian smoothing
        gaussian                    = vtk.vtkImageGaussianSmooth()
        gaussian.SetInputConnection(reader.GetOutputPort())
        gaussian.SetStandardDeviations(5.0, 5.0, 5.0)                               # Standard deviations for the Gaussian in X, Y, Z
        gaussian.SetRadiusFactors(10.0, 10.0, 10.0)                                 # Radius factors 
        gaussian.Update()

        # Write the output to a new .vti file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(prefix+"/Image_Binary_blurred_01.vti")
        writer.SetInputConnection(gaussian.GetOutputPort())
        writer.Write()



def get_values_segmented_images(patient, blurred = False):

    from vtk.util.numpy_support import vtk_to_numpy
    import numpy as np


    prefix = "PA"+str(patient)
    input_image = prefix+"/"+prefix+ "_Binary.vti"
    # Load the binary segmented VTI file
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(input_image)
    reader.Update()

    # Get the image data
    image_data = reader.GetOutput()

    # Extract scalar data
    scalars = image_data.GetPointData().GetScalars()

    # Convert VTK scalars to a Numpy array
    scalar_array = vtk_to_numpy(scalars)

    # Get all unique values
    unique_values = np.unique(scalar_array)

    print(f"Unique values in the image: {unique_values}")


if Threshold_blurred_images:
    import vtk
    for i in Patients_Ids:
        prefix = "PA"+str(i)
        input_image = prefix+"/Image_Binary_blurred_01.vti"
        output_image = prefix+"/"+prefix+"_Image_Binary_thrshd_01.vti"
        # Load the VTI file
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_image)
        reader.Update()

        # Get the image data
        image_data = reader.GetOutput()

        # Apply a threshold to isolate the blurred zone
        threshold_filter = vtk.vtkImageThreshold()
        threshold_filter.SetInputData(image_data)

        # Set the threshold values based on your blurred range
        # These values depend on the range of values introduced by the Gaussian filter
        # Adjust the ranges as necessary (example: 0.1 to 0.9 for blurred values)
        threshold_filter.ThresholdBetween(1, 70)

        # Assign a high constant value to the blurred zone
        new_value = 500.0  # Replace with your desired high value
        threshold_filter.SetInValue(new_value)  # Value for the blurred zone
        threshold_filter.SetOutValue(0)         # Value for other zones (optional)
        threshold_filter.ReplaceInOn()
        threshold_filter.ReplaceOutOn()
        threshold_filter.Update()

        # Get the modified output
        output_data = threshold_filter.GetOutput()

        # Save the modified data back to a .vti file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(output_image)
        writer.SetInputData(output_data)
        writer.Write()


if Copy_initial_image:

    for i in Patients_Ids:
        prefix = "PA"+str(i)
        # cp_command_images = "cp PA1/Image_Binary_blurred_01.vti "+prefix+"/Image_Binary_blurred_00.vti"
        cp_command_images = "cp PA1/PA1_Image_Binary_thrshd_01.vti "+prefix+"/"+prefix+"_Image_Binary_thrshd_00.vti"
        os.system(cp_command_images)



if Downscale_vti:
    import numpy as np
    import dolfin_warp     as dwarp
    # start by looping all patients vti to get smaller image resolution
    metadatas = []
    for i in Patients_Ids:
        prefix = "PA"+str(i)
        input_files_raw             = prefix+"/PGM/"+"Pat"+str(i)+"_inspi1"             # For metadata only here
        pgm_files_raw, _            = pgm2array(input_files_raw)
        metadata                    = get_metada_PGM(
                                                        input_file        = pgm_files_raw[0],
                                                        metadata_fields   = ["Rows", "Columns"]
                                                    )
        metadata.append(len(pgm_files_raw))
        metadata = [int(meta) for meta in metadata]
        metadatas.append(metadata)


    # Compute scaling factor for each image
    meta_array = np.array(metadatas)
    min_xy = min(meta_array[:,0])
    min_z = min(meta_array[:,2])

    target = np.array([min_xy,min_xy,min_z]) 

    scaling_factors = (meta_array/target).tolist()


    # rescale all image
    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        # Create copy of original sized image
        initial_image1 = prefix+"/Image_Binary_blurred_01.vti"
        command_cp_clone1 = "cp "+initial_image1+" "+prefix+"/Image_Binary_blurred_scaled_01.vti"
        initial_image2 = prefix+"/Image_Binary_blurred_00.vti"
        command_cp_clone2 = "cp "+initial_image2+" "+prefix+"/Image_Binary_blurred_scaled_00.vti"
        os.system(command_cp_clone1)
        os.system(command_cp_clone2)
        print("Patient "+prefix+" cloned")



        dwarp.compute_downsampled_images(
                                        images_folder           = "./"+prefix+"/", 
                                        images_basename         = "Image_Binary_blurred_scaled",
                                        downsampling_factors    = scaling_factors[i-1]
                                        )