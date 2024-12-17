
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
Threshold_blurred_images            = False
Left_right_lungs                    = False
Threshold_blurred_images_LL_RR      = False
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


if Left_right_lungs:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    for i in Patients_Ids:

        prefix                      = "PA"+str(i)
        input_file                  = prefix+"/"+prefix+"_Binary.vti"
        output_file                 = prefix+"/"+prefix+"_Binary_LL_RL.vti"

        # Load the VTI file
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_file)
        reader.Update()

        # Get the image data
        image_data = reader.GetOutput()

        # Extract scalar data
        scalars = image_data.GetPointData().GetScalars()

        # Convert VTK scalars to a Numpy array
        scalar_array = vtk_to_numpy(scalars)

        # Convert the numpy array to a signed type, e.g., int32 or float32
        signed_scalar_array = scalar_array.astype(np.int32)  # Use np.float32 for floating-point

        # Replace values
        signed_scalar_array[signed_scalar_array == 200] = -100

        # Convert back to VTK array and set it
        modified_scalars = numpy_to_vtk(signed_scalar_array)
        image_data.GetPointData().SetScalars(modified_scalars)

        # Save the modified VTI file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(image_data)
        writer.Write()

        print("Sucess. "+prefix)




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



# if Blur_binary_vti:
#     import vtk 

#     for i in Patients_Ids:
#         prefix                      = "PA"+str(i)
#         input_file                  = prefix+"/"+prefix+"_Binary.vti"
#         # Read the .vti file
#         reader                      = vtk.vtkXMLImageDataReader()
#         reader.SetFileName(input_file)
#         reader.Update()

#         # Apply Gaussian smoothing
#         gaussian                    = vtk.vtkImageGaussianSmooth()
#         gaussian.SetInputConnection(reader.GetOutputPort())
#         gaussian.SetStandardDeviations(5.0, 5.0, 5.0)                               # Standard deviations for the Gaussian in X, Y, Z
#         gaussian.SetRadiusFactors(10.0, 10.0, 10.0)                                 # Radius factors 
#         gaussian.Update()

#         # Write the output to a new .vti file
#         writer = vtk.vtkXMLImageDataWriter()
#         writer.SetFileName(prefix+"/Image_Binary_blurred_01.vti")
#         writer.SetInputConnection(gaussian.GetOutputPort())
#         writer.Write()


def create_LL_RL_binary(Patients_Ids, input_name, output_name):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    for i in Patients_Ids:

        prefix                      = "PA"+str(i)
        input_file                  = prefix+"/"+prefix+input_name+".vti"
        output_file                 = prefix+"/"+prefix+output_name+".vti"

        # Load the VTI file
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_file)
        reader.Update()

        # Get the image data
        image_data = reader.GetOutput()

        # Extract scalar data
        scalars = image_data.GetPointData().GetScalars()

        # Convert VTK scalars to a Numpy array
        scalar_array = vtk_to_numpy(scalars)

        # Convert the numpy array to a signed type, e.g., int32 or float32
        signed_scalar_array = scalar_array.astype(np.int32)  # Use np.float32 for floating-point

        # Replace values
        signed_scalar_array[signed_scalar_array == 200] = -100

        # Convert back to VTK array and set it
        modified_scalars = numpy_to_vtk(signed_scalar_array)
        image_data.GetPointData().SetScalars(modified_scalars)

        # Save the modified VTI file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(image_data)
        writer.Write()

        print("Sucess. "+prefix)


def blur_vti(Patients_Ids, input_name, output_name, radius = 10.0, std = 5.0):
    import vtk 

    for i in Patients_Ids:
        prefix                      = "PA"+str(i)
        input_file                  = prefix+"/"+prefix+input_name+".vti"
        # Read the .vti file
        reader                      = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_file)
        reader.Update()

        # Apply Gaussian smoothing
        gaussian                    = vtk.vtkImageGaussianSmooth()
        gaussian.SetInputConnection(reader.GetOutputPort())
        gaussian.SetStandardDeviations(std, std, std)                               # Standard deviations for the Gaussian in X, Y, Z
        gaussian.SetRadiusFactors(radius, radius, radius)                                 # Radius factors 
        gaussian.Update()

        # Write the output to a new .vti file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(prefix+"/"+prefix+output_name+".vti")
        writer.SetInputConnection(gaussian.GetOutputPort())
        writer.Write()
        print("Done blurring. "+prefix)


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



def threshold_blurred_LL_RL(Patients_Ids, input_name, output_name, field_name = "pixel intensity"):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    for i in Patients_Ids:

        prefix                      = "PA"+str(i)
        input_file                  = prefix+"/"+prefix+input_name+".vti"

        # Load the VTI file
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(input_file)
        reader.Update()

        # Get the image data
        image_data = reader.GetOutput()

        # Extract scalar data and convert to Numpy array
        scalars = image_data.GetPointData().GetScalars()
        scalar_array = vtk_to_numpy(scalars)

        # Apply the value adjustments
        scalar_array[(scalar_array >= -60) & (scalar_array <= -45)] = -150
        scalar_array[(scalar_array >= 45) & (scalar_array <= 60)] = 150

        # Convert the modified Numpy array back to VTK format
        modified_scalars = numpy_to_vtk(scalar_array)
        modified_scalars.SetName(field_name)  # Set the name of the scalar field

        # Assign the modified scalars back to the image data
        image_data.GetPointData().SetScalars(modified_scalars)

        # Save the modified image to a new VTI file
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(prefix+"/"+prefix+output_name+".vti")
        writer.SetInputData(image_data)
        writer.Write()

        print("Done thresholding. "+prefix)




def vti_unsigned_char(patient, input_name):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    prefix = "PA"+str(patient)
    input_file                  = prefix+"/"+prefix+input_name+".vti"
    output_image                = prefix+"/"+prefix+"_INT"+input_name+".vti"


    # Load the VTI file
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Get the image data
    image_data = reader.GetOutput()

    # Extract the scalar data as a Numpy array
    scalars = image_data.GetPointData().GetScalars()
    scalar_array = vtk_to_numpy(scalars)

    # Rescale the data to fit within [0, 255]
    min_val = scalar_array.min()
    max_val = scalar_array.max()

    # Avoid division by zero if all values are the same
    if max_val > min_val:
        rescaled_array = 255 * (scalar_array - min_val) / (max_val - min_val)
    else:
        rescaled_array = np.zeros_like(scalar_array)

    # Convert the data type to unsigned char (uint8)
    rescaled_array = rescaled_array.astype(np.uint8)

    # Convert the modified Numpy array back to VTK format
    modified_scalars = numpy_to_vtk(rescaled_array, deep=True)
    # modified_scalars.SetName("RescaledScalars")  # Optional: Set a new name for the scalar field

    # Assign the new scalars back to the image data
    image_data.GetPointData().SetScalars(modified_scalars)

    # Save the modified image to a new VTI file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_image)
    writer.SetInputData(image_data)
    writer.Write()

    print("Reformating to unsigned char done."+prefix)


def vti_int8(patient, input_name):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    prefix = "PA"+str(patient)
    input_file                  = prefix+"/"+prefix+input_name+".vti"
    output_image                = prefix+"/"+prefix+"_INT"+input_name+".vti"


    # Load the VTI file
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Get the image data
    image_data = reader.GetOutput()

    # Extract the scalar data as a Numpy array
    scalars = image_data.GetPointData().GetScalars()
    scalar_array = vtk_to_numpy(scalars)

    # Rescale the data to fit within [0, 255]
    min_val = scalar_array.min()
    max_val = scalar_array.max()

    # Avoid division by zero if all values are the same
    if max_val > min_val:
        rescaled_array = -128 + 255 * (scalar_array - min_val) / (max_val - min_val)
    else:
        rescaled_array = np.zeros_like(scalar_array)

    # Convert the data type to unsigned char (uint8)
    rescaled_array = rescaled_array.astype(np.int8)

    # Convert the modified Numpy array back to VTK format
    modified_scalars = numpy_to_vtk(rescaled_array, deep=True)
    # modified_scalars.SetName("RescaledScalars")  # Optional: Set a new name for the scalar field

    # Assign the new scalars back to the image data
    image_data.GetPointData().SetScalars(modified_scalars)

    # Save the modified image to a new VTI file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_image)
    writer.SetInputData(image_data)
    writer.Write()

    print("Reformating to unsigned char done."+prefix)



if Copy_initial_image:

    for i in Patients_Ids:
        prefix = "PA"+str(i)
        # cp_command_images = "cp PA1/Image_Binary_blurred_01.vti "+prefix+"/Image_Binary_blurred_00.vti"
        cp_command_images = "cp PA1/PA1_Image_Binary_thrshd_01.vti "+prefix+"/"+prefix+"_Image_Binary_thrshd_00.vti"
        os.system(cp_command_images)

def prepare_pairs_images(Patients_Ids, input_name):

    for i in Patients_Ids:
        prefix = "PA"+str(i)
        cp_command_images_01 = "cp "+prefix+"/"+prefix+input_name+".vti "+prefix+"/"+prefix+input_name+"_01.vti" 
        os.system(cp_command_images_01)


        # cp_command_images = "cp PA1/Image_Binary_blurred_01.vti "+prefix+"/Image_Binary_blurred_00.vti"
        cp_command_images_00 = "cp PA1/PA1"+input_name+"_01.vti "+prefix+"/"+prefix+input_name+"_00.vti" 
        os.system(cp_command_images_00)

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



def from_threshold_to_gradually_threshold(patient, input_name, output_name):
    blur_vti(
            Patients_Ids    = [patient],
            input_name      = input_name,
            output_name     = input_name+"very_blurred", 
            # radius          = 150, 
            std             = 50)

    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import numpy as np
    prefix = "PA"+str(patient)
    input_file_1 = prefix+"/"+prefix+input_name+".vti"

    input_file_3 = prefix+"/"+prefix+"_Binary_LL_RL"+".vti"


    input_file_2 = prefix+"/"+prefix+input_name+"very_blurred"+".vti"

    output_file = prefix+"/"+prefix+output_name+".vti"

    # Load the VTI files (1 = threshold, 2 = very blurred)
    reader_1 = vtk.vtkXMLImageDataReader()
    reader_1.SetFileName(input_file_1)
    reader_1.Update()

    reader_2 = vtk.vtkXMLImageDataReader()
    reader_2.SetFileName(input_file_2)
    reader_2.Update()

    reader_3 = vtk.vtkXMLImageDataReader()
    reader_3.SetFileName(input_file_3)
    reader_3.Update()

    # Get the image data
    image_data_1 = reader_1.GetOutput()
    image_data_2 = reader_2.GetOutput()
    image_data_3 = reader_3.GetOutput()


    # Extract the scalar data as a Numpy array
    scalars_1 = image_data_1.GetPointData().GetScalars()
    scalar_array_1 = vtk_to_numpy(scalars_1)

    scalars_2 = image_data_2.GetPointData().GetScalars()
    scalar_array_2 = vtk_to_numpy(scalars_2)


    scalars_3 = image_data_3.GetPointData().GetScalars()
    scalar_array_3 = vtk_to_numpy(scalars_3)

    # Rescale the data of very blurred
    min_val_1 = scalar_array_1.min()
    max_val_1 = scalar_array_1.max()

    min_val_2 = scalar_array_2.min()
    max_val_2 = scalar_array_2.max()

    # Avoid division by zero if all values are the same
    # if max_val_2 > min_val_2:
    #     rescaled_array_blurred = min_val_1 + (max_val_1 - min_val_1) * (scalar_array_2 - min_val_2) / (max_val_2 - min_val_2)
    # else:
    #     rescaled_array_blurred = np.zeros_like(scalar_array)



    ## Debug
    rescaled_array_blurred = scalar_array_2
    ##
    rescaled_array_blurred[(np.abs(scalar_array_3) >= 90) ] = scalar_array_3[(np.abs(scalar_array_3) >= 90) ]

    rescaled_array_blurred[ (np.abs(scalar_array_1) >= 120)] = scalar_array_1[ (np.abs(scalar_array_1) >= 120)]


    # Convert the modified Numpy array back to VTK format
    new_values = numpy_to_vtk(rescaled_array_blurred, deep=True)
    # modified_scalars.SetName("RescaledScalars")  # Optional: Set a new name for the scalar field

    # Assign the new scalars back to the image data
    image_data_1.GetPointData().SetScalars(new_values)

    # Save the modified image to a new VTI file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(image_data_1)
    writer.Write()

    print("Done. written at"+output_name+ ".vti")















#%% Pipeline

# Put one lung to -100 and the other to +100
# create_LL_RL_binary(
#         Patients_Ids    = Patients_Ids,
#         input_name      = "_Binary",
#         output_name     = "_Binary_LL_RL")


# Blur signed lung binaries
# blur_vti(
#         Patients_Ids    = Patients_Ids,
#         input_name      = "_Binary_LL_RL",
#         output_name     = "_Binary_LL_RL_blurred")


# # thresholds blurred signed lung binaries
# threshold_blurred_LL_RL(
#         Patients_Ids    = Patients_Ids,
#         input_name      = "_Binary_LL_RL_blurred",
#         output_name     = "_Binary_LL_RL_blurred_thrshld")

#### STRAT1 -
match strategy:
    case "blurred_threshold":
        #  blur thresholds blurred signed lung binaries
        blur_vti(
                Patients_Ids    = Patients_Ids,
                input_name      = "_Binary_LL_RL_blurred_thrshld",
                output_name     = "_Binary_LL_RL_blurred_thrshld_blurred")

        # Copy initial image and target images to right names and folders
        prepare_pairs_images(Patients_Ids, "_Binary_LL_RL_blurred_thrshld_blurred")

        for patient in Patients_Ids:
        # Convert to unsigned char the vti (effectively rescaling the oxel values), otherwise the kernel crashes
            vti_unsigned_char(patient,"_Binary_LL_RL_blurred_thrshld_blurred_00")
            vti_unsigned_char(patient,"_Binary_LL_RL_blurred_thrshld_blurred_01")

    case "external_progressive_gradient":
        for patient in Patients_Ids:
            from_threshold_to_gradually_threshold(patient, "_Binary_LL_RL_blurred_thrshld", "_thrshld_external_gradient")


        # # Copy initial image and target images to right names and folders
        prepare_pairs_images(Patients_Ids, "_thrshld_external_gradient")

        for patient in Patients_Ids:
        # Convert to unsigned char the vti (effectively rescaling the voxel values), otherwise the kernel crashes
            vti_unsigned_char(patient,"_thrshld_external_gradient_00")
            vti_unsigned_char(patient,"_thrshld_external_gradient_01")

            # vti_int8(patient,"_thrshld_external_gradient_00")
            # vti_int8(patient,"_thrshld_external_gradient_01")