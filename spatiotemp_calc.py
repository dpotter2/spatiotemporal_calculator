#
#   spatiotemporal calculator
#
#       add some pithy notes here
#
#   Initial - 2019.02.13 - David Potter
#   version 3 started - 2019.02.26

import numpy as np  # math functions
import pint         # units
 
u = pint.UnitRegistry()
u.define('micron = 1 * micrometer')

#####
##### Coordinate system = X positive to the right, Y positive is into the monitor, Z postive is up
#####

### There are two coordinate systems here: 1) the sample XYZ 2) the light path XYZ
### The sample and light path Y-axes are coincedent; the X & Z may not be depending on system design

### A sample sphere moves along each compass rose axis at some velocity. It starts at the 0, 0, 0 XYZ of the sample coord

### The light sheet propogates along the X axis in the positive direction; Z is the detection axis, camera at positive end

### The light path XYZ can be rotated about the sample XYZ. Angle is measured in degrees CCW, about the Y axis, from the X axis

### Detection objective scanning systems are tretaed as if the detection is immobile and the stage scans the sample
### through the sheet. 

### The XYZ coordinate of the sphere in the sample coord is transformed to the light path coord if required by system design
###        or will be when I finish it

### The sphere and the cutting plane are both transformed to light path XYZ 0, 0, 0 to make the math simpler

### The resulting small circle center points are then transformed back into the sample coord system and mapped to the sample
### coord compass rose

### Each sphere along a compass rose axis is calculated starting from 0,0,0 but the results will be shifted some amount along
### each axis so that the spheres don't overlap each other. Duplicate (mirrored) axes will not show both for clarity.

### Take the imaging time interval and apply it to the moving sphere as follows: At each interval calculate the location
### of the sphere center in the sample coord (place into an array?). Transform the sphere location to the coordinate in the
### light path coord system.

### Consider the sheet thickness when cutting the sphere. Essentially two small circles convolved. If the sphere is hollow,
### then there will be four small circles convolved!!  (convolove = Z projection?)

### Also calculate the "smear" of the sphere during the exposure interval by cutting the sphere at the beginning and the end
### of the exposure and convolving. Each cut is as above thus there could be a whole lot of comvolved small circles if the 
### sphere is hollow

### Calculate light sheet location for each time interval. Determine "line" endpoints. Calculate intersection with the
### trasformed sphere coordinate. Create an image tiff with the cut and convolved spheres. Repeat until we have the image
### volume for each axis of the compass rose. Good lord...

#############################################################################################################################

# fluorophore data
        # my plan is to use this data to analyze resolution based on contrast changes related to emmision intensity
        # "shell_thickness" above is part of this math as well
image_background = 120                      #
fluorophore_emission = 0.525 * u.micron
fluorophore_cross_section = 0.00            # not sure how I will do this
fluorophore_relative_energy = 0.00          # or this either

# excitation lines
laser_445 = 0.445 * u.micron
laser_488 = 0.488 * u.micron
laser_560 = 0.560 * u.micron
laser_642 = 0.642 * u.micron

#laser_power = input("Laser power: ") * u.milliwatt
laser_power = 350 * u.milliwatt
#print("laser power", laser_power)

# detection objective
detection_NA = 1.1
medium_refractive_index = 1.333             # water for the moment
obj_theta_rad = np.arcsin(detection_NA/medium_refractive_index) * u.radian
obj_theta_deg = np.rad2deg(obj_theta_rad)
print("Detection theta:", round(obj_theta_rad, 4), round(obj_theta_deg, 4))
det_obj_resolution = (0.61 * fluorophore_emission) / detection_NA
print("Detection resolution", round(det_obj_resolution, 4))
det_obj_mag = 25                            # as labeled on obj body
det_man_tube_lens = 200 * u.millimeter      # manufactures required tube lens focal length
det_obj_focal_length = det_man_tube_lens / det_obj_mag      # geometric optics approximation
back_pupil_dia = det_obj_focal_length * np.tan(obj_theta_rad)
system_tube_lens = 500 * u.millimeter       # as used tube lens focal length
system_mag = det_obj_mag * (system_tube_lens / det_man_tube_lens)       # image mag at the camera chip
print("System mag:", round(system_mag, 2), "Pupil", round(back_pupil_dia, 2), "Focal length:", round(det_obj_focal_length, 2))

# camera
camera_manufacture = "Hamamatsu"
camera_model = "Orca Flash 4 v2"
camera_pix_xy = np.array([6.5, 6.5]) * u.micron       # camera pixel spacing
camera_micron_pixel = camera_pix_xy / system_mag

#######################################################################

# sample coord versus detection coord systems
detect_to_sample_Z_ang_deg = 31.8      # CCW angle relative to sample X axis, rotated about the Y axis; from the point of view of Y+
                            # 0 angle equals an orthogonal excitation to detection axis i.e. a standard light sheet
detect_to_sample_Z_ang_rad = np.deg2rad(detect_to_sample_Z_ang_deg)

# scan information
z_voxel = 0.21078 * u.micron
volume_image_number = 101
if detect_to_sample_Z_ang_deg > 0:
    scan_interval = z_voxel / np.abs(np.sin(detect_to_sample_Z_ang_rad))     # volume is skewed
else:
    scan_interval = z_voxel

scan_dist = (volume_image_number - 1) * scan_interval
print("scan distance:", round(scan_dist, 2), "scan interval",round(scan_interval, 4))

##########################################################################################

# experiment parameters
num_channels = 4
cam_exposure = 20 * u.millisecond
cam_memory_speed = (800 * u.megabytes) / (1 * u.second)                    # megabytes per second
print("Max camera memory speed:", cam_memory_speed)
system_house_keeping = 1.15                                    # a reasonable estimate for SLM flipping etc.
    # I should put a variable in here to deal with camera chip scan method as the options consume different times
cam_cycle_time = cam_exposure * system_house_keeping           # cycle time for one image of one channel
cycle_per_stage_interval = num_channels * cam_cycle_time
#print("cycle_per_stage_interval -", cycle_per_stage_interval)
cycle_per_stack = cycle_per_stage_interval * volume_image_number
#print("cycle_per_stack -", cycle_per_stack)
time_point_pause = 0.050 * u.second                             # time pause between each time point stack
#print("time_point_pause - ", time_point_pause)
stacks_taken = 5                            # I intend on adding code to deal with apparent discontinous particle movement related to the time phase
total_imaging_time_ms = ((cycle_per_stack * stacks_taken) +
                      (time_point_pause * (stacks_taken - 1)))       # in seconds I hope
total_imaging_time = total_imaging_time_ms.to(u.second)
#print("total_imaging_time -", total_imaging_time)                      

# image info
image_pix = np.array([800, 600])            # the ROI of the captured image
image_fov = image_pix * camera_micron_pixel
image_bits = 16                             # dynamic range
bit_per_byte = 8
mega_binary_factor = (2^20)                 # bytes per "mega"
image_memory = (image_pix[0] * image_pix[1] * image_bits / bit_per_byte) / mega_binary_factor    # memory size of each image in the stack in megabytes
volume_memory_size_mega = image_memory * volume_image_number #* u.megabytes               # in megabytes
volume_memory_size_giga = (volume_memory_size_mega / 1024) * u.gigabytes                    # in gigabytes
memory_per_second = volume_memory_size_giga / total_imaging_time            # in megabytes per second
#print("volume_memory_size_giga:", round(volume_memory_size_giga, 2), "Memory_per_second:", round(memory_per_second, 2))

# light sheet parameters
sheet_length = 30 * u.micron        # calculate based on system optics or bessels parameters
sheet_width = 50 * u.micron         # along the Y axis
sheet_thick = 0.4* u.micron         # calculate this properly later


#################    CHANGE THIS     ##########################

# sample sphere represented as a 2D circle (the small circle as defined by the XZ cutting plane)
# Y location doesn't change the cutting plane, so it is carried through the processing and only used when
# drawing the sphere section in the final multi tiff tiff


circle_rad = 2.50 * u.micron                # outer radius of the sphere in microns
circle_center = np.array([0, 0, 0]) * u.micron  # X, Y, Z coords of sphere location in an array
shell_thickness = 40 * u.nanometer          # sphere fluorophore is "membrance bound" or "intravesicular"
        # will need to treat as two great circles cut by the plane
if shell_thickness >= circle_rad:
    shell_thickness = 0                     # i.e. the sphere is solid i.e. flurophore is intravesicular
    print("woo hoo")
else:  # shell_thickness < circle_rad
    print("by jimminy") # add code soon

sphere_velocity = (10 * u.micron) / (1 * u.second)        # microns per second
sphere_motion_vector = np.array([1, 0, 0])              # 1, 0, 1 -> shere is moving along the +X and +Z axes i.e. 45 deg in the XZ plane

shift = circle_center
print("shift:", shift)                     # set the shift to center distances
circle_offset = circle_center - shift       # shift the great circle to XYZ 0,0,0


#detect_to_stage_transform = np.sin(detect_to_sample_Z_ang_rad)

# current stage location
scan_loc = np.array([0, 0, 0.250]) * u.micron  # current  coordinate of stage (at X = 0 and Y = 0)

transform = np.array([np.sin(detect_to_sample_Z_ang_rad), 0, np.cos(detect_to_sample_Z_ang_rad)])
pt1 = 0.5 * sheet_length * transform * (-1)  # upper end point
pt1_shift = pt1 - shift
pt2 = 0.5 * sheet_length * transform         # lower end point
pt2_shift = pt2 - shift
print("pt1_shift:", pt1_shift, "pt2_shift:", pt2_shift)

# great circle diameter and signum
line = pt2_shift - pt1_shift
print("line:", line)
if line[2] < 0.0:
    sgn_z = -1.0
else:
    sgn_z = 1.0

#line_dr = np.hypot(*line) # argument unpacking - pulls the two elements out of the array and pythags them to get the hypotnuse
line_dr = np.hypot(line[0], line[2])
print("line_dr -", line_dr)
line_D = (pt1_shift[0] * pt2_shift[2]) - (pt2_shift[0] * pt1_shift[2])
print("line D", line_D)

# discriminant -> checks for line / great circle intersection
does_intersect = (np.square(circle_rad) * np.square(line_dr)) - np.square(line_D)

if does_intersect > 0.0:            # intersect
    discriminant = 2
elif does_intersect == 0.0:         # tangent
    discriminant = 1
elif does_intersect < 0.0:          # pass like ships in the night
    discriminant = 0

# calculate both intersection points, cap diameter and center -> this is the circle in the image plane
#   note: this is the shifted cutting plane
if discriminant == 2:
    int_pt1_x = (((line_D * line[2]) + (sgn_z * line[0] * np.sqrt(np.square(circle_rad) 
                * np.square(line_dr) - np.square(line_D)))) / np.square(line_dr))
    int_pt1_x_shift = int_pt1_x + shift[0]
    int_pt1_z = (((-1 * line_D * line[0]) + (np.abs(line[2]) * np.sqrt(np.square(circle_rad)
                * np.square(line_dr) - np.square(line_D)))) / np.square(line_dr))
    int_pt1_z_shift = int_pt1_z + shift[2]
    print("intersect points 1:", round(int_pt1_x_shift, 4), round(int_pt1_z_shift, 4))
    #print("intersect_pt1 -", intersect_pt1)

    int_pt2_x = (((line_D * line[2]) - (sgn_z * line[0] * np.sqrt(np.square(circle_rad)
                * np.square(line_dr) - np.square(line_D)))) / np.square(line_dr))
    int_pt2_x_shift = int_pt2_x + shift[0]
    int_pt2_z = (((-1 * line_D * line[0]) - (np.abs(line[2]) * np.sqrt(np.square(circle_rad)
                * np.square(line_dr) - np.square(line_D)))) / np.square(line_dr))
    int_pt2_z_shift = int_pt2_z + shift[2]
    print("intersect points 2:", round(int_pt2_x_shift, 4), round(int_pt2_z_shift, 4))

    # calculate cap diameter and center coordinates
    line_x_side = int_pt2_x_shift - int_pt1_x_shift
    line_z_side = int_pt2_z_shift - int_pt1_z_shift
    cap_diameter = np.sqrt(np.square(line_x_side) + np.square(line_z_side))
    cap_center_x = int_pt1_x_shift + 0.5 * line_x_side
    cap_center_z = int_pt1_z_shift + 0.5 * line_z_side
    print("cap dia & center:", round((cap_diameter / 2), 4), round(cap_center_x, 4), round(cap_center_z, 4))
  


print("That's all folks")