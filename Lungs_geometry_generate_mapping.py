
import numpy as np
import dolfin
import meshio
import vtk
import fenics as fenics
import dolfin_warp     as dwarp
import dolfin_mech as dmech
import myVTKPythonLibrary as myvtk

#%% Define initial mesh (that hold all mappings)

mesh_LL                                                     =fenics.Mesh("Initial_Config/mesh_LL.xml")
mesh_RL                                                     =fenics.Mesh("Initial_Config/mesh_RL.xml")

regul_level                                                 = 0.0
regul_type                                                  = "discrete-equilibrated-tractions-normal"

if any([_ in regul_type for _ in ["linear", "simple"]]):
    regul_model = "hooke"
else:
    regul_model = "ciarletgeymonatneohookean"




#%% Create initial scaling 
# The objective is to create an initial scaled mapping that lies within all possible rib cages to help with tracking
def initial_scaling(mesh, lung, coef =-0.4, reduced_kinematics_model = "translation+rotation+scaling+shear"):
    saving_name_initial_scalaing = "initial_scaling_"+lung+".dat"
    alpha = coef                                                                                # Scaling factor

    dim = 3
    x = dolfin.SpatialCoordinate(mesh)
    dV = dolfin.Measure(
            "dx",
            domain=mesh)

    center_gravity = [0,0,0]                                                                    # Center of mass coordinates
    for k_dim in range(dim):
        center_gravity[k_dim] = dolfin.assemble(x[k_dim]*dV)/dolfin.assemble(1*dV)

    x_mid = center_gravity[0]
    y_mid = center_gravity[1]
    z_mid = center_gravity[2]

    ## Write down reduced_kinematics initialisation file
    reduced_disp_initial_scaling_list = [-x_mid,-y_mid,-z_mid]                            # List of initial reduced displacements

    if ("rotation" in reduced_kinematics_model):
        reduced_disp_initial_scaling_list+=[0,0,0]

    reduced_disp_initial_scaling_list+=[1,1,1]

    if ("shear" in reduced_kinematics_model):
        reduced_disp_initial_scaling_list+=[0,0,0]

    reduced_disp_initial_scaling = alpha*np.array([reduced_disp_initial_scaling_list])               # Reduced displacement for the 6 modes (3 translation and 3 rotations) reduced-kinematics

    print(f"shape of initial reduced displacement is {reduced_disp_initial_scaling.shape}")
    np.savetxt(saving_name_initial_scalaing, reduced_disp_initial_scaling)                      # Save the reduced displacements

#%% Define tracking functions

destination_path = "./"                                                                         #Root path for patients solution folders

def reduced_kiematics(image_base_name, patient,lung, mesh, tol=1e-6, images_quadrature = 6, reduced_kinematics_model = "translation+rotation+scaling+shear"):

   
    prefix = "PA"+str(patient)
    dwarp.warp(
            images_char_func                            = False,
            working_folder                              = destination_path+prefix,
            # working_basename                            = "thrshd_mapping_reduced_kinematics"+'_'+prefix+'_'+lung,
            working_basename                            = "TEST_mapping_reduced_kinematics"+'_'+prefix+'_'+lung,
            images_folder                               = destination_path+prefix,
            images_basename                             = prefix+image_base_name,
            images_ext                                  = "vti",
            mesh                                        = mesh,
            kinematics_type                             = "reduced",
            reduced_kinematics_model                    = reduced_kinematics_model,
            images_quadrature                           = 6,
            n_iter_max                                  = 10,
            regul_poisson                               = 0.3,
            regul_type                                  = regul_type,
            regul_model                                 = regul_model,
            regul_level                                 = regul_level,
            relax_type                                  ="backtracking",
            tol_dU                                      = tol,
            continue_after_fail                         = 1,
            write_VTU_files                             = True,
            write_VTU_files_with_preserved_connectivity = True,
            initialize_reduced_U_from_file              = True,
            initialize_reduced_U_filename               = "initial_scaling_"+lung+".dat",
            print_iterations                            = 1,
            save_reduced_disp                           = True) 


def tracking(patient,lung, mesh, tol=1e-3, regul = 0.5, images_quadrature = 6):
    prefix = "PA"+str(patient)
    dwarp.warp(
            working_folder                              = destination_path+prefix,
            working_basename                            = "thrshd_mapping"+'_'+prefix+'_'+lung,
            images_folder                               = destination_path+prefix,
            images_basename                             = prefix+"_INT_thrshld_external_gradient_blurred",
            images_ext                                  = "vti",
            mesh                                        = mesh,
            n_iter_max                                  = 1000,
            nonlinearsolver                             = "newton",
            regul_types                                 = [ "continuous-hyperelastic"],         #"discrete-mesh",  "continuous-equilibrated", "discrete-tractions-normal", "continuous-equilibrated", "discrete-tractions-tangential", "continuous-hyperelastic" "continuous-hyperelastic",
            regul_model                                 = "ogdenciarletgeymonatneohookeanmooneyrivlin",
            regul_levels                                = [regul],
            regul_poisson                               = 0.0,
            images_quadrature                           = 8,
            relax_type                                  = "backtracking",
            tol_dU                                      = tol,
            write_VTU_files                             = True,
            write_VTU_files_with_preserved_connectivity = True,
            images_char_func                            = True,
            initialize_U_from_file                      = 1,
            initialize_U_folder                         = destination_path+prefix, 
            initialize_U_basename                       = "thrshd_mapping_reduced_kinematics"+'_'+prefix+'_'+lung,
            initialize_U_ext                            = "vtu",
            initialize_U_array_name                     = "displacement",
            initialize_U_method                         = "dofs_transfer",                      # dofs_transfer, interpolation, projection
            print_iterations                            = 1)                                    



def warp_and_blur(patient, attenuation_factors, lung, mesh):
    import Create_data
    prefix = "PA"+str(patient)
    common_basename = "_INT_thrshld_external_gradient"
    images_quadrature_progressive = np.linspace(1, 6, len(attenuation_factors))  # Generate m evenly spaced values
    images_quadrature_progressive = np.round(images_quadrature_progressive).astype(int)

    for i in range(len(attenuation_factors)):
        Create_data.gaussian_windowing(
                working_folder                              = destination_path+prefix,
                image_name                                  = prefix+common_basename,
                attenuation_factor                          = attenuation_factors[i],   
                verbose                                     = True
                )
        image_base_name = common_basename + "_downsampled="+str(attenuation_factor)
        reduced_kiematics(image_base_name, patient,lung, mesh, tol=1e-6, images_quadrature = images_quadrature_progressive[i])

        ### Update with new reduced kinematics: Check that results of reduced kinematics is saved in .dat

        if i >=1 :
            initialize_reduced_U_filename = self.working_basename+"_reduced_kinematics_inter.dat"

                ###### Check that .dat is saved in reduced kinematics for latter initialisation check name for initial .dat files
        else:
            initialize_reduced_U_filename = "initial_scaling_"+lung+".dat"

        dwarp.warp(
            images_char_func                = False,
            working_folder                  = working_folder,
            working_basename                = working_basename+"_downsampled="+str(attenuation_factor),
            images_folder                   = images_folder,
            images_basename                 = images_basename_blur_factor,
            images_quadrature               = images_quadrature_progressive[i],
            mesh                            = mesh,
            kinematics_type                 = kinematics_type,
            reduced_kinematics_model        = reduced_kinematics_model,
            normalize_energies              = normalize_energies,
            relax_type                      = relax_type,
            tol_dU                          = tol_dU,
            write_qois_limited_precision    = write_qois_limited_precision, 
            initialize_reduced_U_from_file  = initialize_reduced_U_from_file,
            initialize_reduced_U_filename   = initialize_reduced_U_filename,
            )





N_patients = 9
# Patients_Ids = list(range(1,  N_patients + 1))
# Patients_Ids.remove(1)
# # Patients_Ids.remove(3)
# # Patients_Ids.remove(4)
# # Patients_Ids.remove(6)



# Patients_Ids.remove(2)



Lungs = ['RL','LL']
Lungs = ['RL']
# Lungs = ['LL']


reduced_kinematics_model = "translation+scaling+rotation"
Patients_Ids = [5]

# for lung in Lungs:
#     match lung:
#         case 'LL':
#             mesh = mesh_LL
#         case 'RL':
#             mesh = mesh_RL
#     initial_scaling(mesh, lung, coef =0.9, reduced_kinematics_model = reduced_kinematics_model)
#     for patient in Patients_Ids:
#         reduced_kiematics(
#                 image_base_name                         = "_INT_thrshld_external_gradient_blurred",
#                 patient                                 = patient,
#                 lung                                    = lung,
#                 mesh                                    = mesh,
#                 tol                                     = 1e-6,
#                 images_quadrature                       = 3,
#                 reduced_kinematics_model                = reduced_kinematics_model
#                 )

        # tracking(
        #         patient                                 = patient,
        #         lung                                    = lung,
        #         mesh                                    = mesh,
        #         tol                                     = 1e-4,
        #         regul                                   = 0.01)



#%% Test of blur and warp from dwarp directly

image_base_name = "_INT_thrshld_external_gradient_blurred"
reduced_kinematics_model = "translation+scaling+rotation"


for lung in Lungs:
    match lung:
        case 'LL':
            mesh = mesh_LL
        case 'RL':
            mesh = mesh_RL
    initial_scaling(mesh, lung, coef =0.9, reduced_kinematics_model = reduced_kinematics_model)
    for patient in Patients_Ids:
        prefix = "PA"+str(patient)
        
        dwarp.blur_and_warp(
                    attenuation_factors                         = [8,16,32],
                    images_char_func                            = False,
                    working_folder                              = destination_path+prefix,
                    # working_basename                            = "thrshd_mapping_reduced_kinematics"+'_'+prefix+'_'+lung,
                    working_basename                            = "MB_mapping_reduced_kinematics"+'_'+prefix+'_'+lung,
                    images_folder                               = destination_path+prefix,
                    images_basename                             = prefix+image_base_name,
                    images_ext                                  = "vti",
                    mesh                                        = mesh,
                    kinematics_type                             = "reduced",
                    reduced_kinematics_model                    = reduced_kinematics_model,
                    images_quadrature                           = 6,
                    n_iter_max                                  = 3,
                    regul_poisson                               = 0.3,
                    regul_type                                  = regul_type,
                    regul_model                                 = regul_model,
                    regul_level                                 = regul_level,
                    relax_type                                  ="backtracking",
                    tol_dU                                      = 1e-2,
                    continue_after_fail                         = 1,
                    write_VTU_files                             = True,
                    write_VTU_files_with_preserved_connectivity = True,
                    initialize_reduced_U_from_file              = True,
                    initialize_reduced_U_filename               = "initial_scaling_"+lung+".dat",
                    print_iterations                            = 1,
                    save_reduced_disp                           = True) 