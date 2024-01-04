# # [Loading IERS EOP data](@id tutorial_03_eop)
# _This example was generated on DATEOFTODAY._

# When working with frames associated with the Earth, it is imperative to incorporate 
# Earth Orientation Parameters (EOP) data. The EOP data are required for the accurate 
# construction of  various Earth associated frames within the [`FrameSystem`](@ref).

# ### Creating compatible EOP file

# To minimize dependencies on external sources, `FrameTransformations` defines a 
# standardized format for EOP data. The expected format consists of a file with 
# the '.eop.dat' extension. 
# This file should contain columns representing different Earth orientation parameters. 
# The columns include:

# 1. **J2000 UTC**: Julian Date (UTC) for J2000 epoch.
# 2. **J2000 TT**: Julian Date (Terrestrial Time) for J2000 epoch.
# 3. **X-Pole**: X-coordinate of the Celestial Intermediate Pole (CIP) in microarcseconds.
# 4. **Y-Pole**: Y-coordinate of the Celestial Intermediate Pole (CIP) in microarcseconds.
# 5. **UT1-UTC**: The difference between Coordinated Universal Time (UTC) and Universal Time 1 (UT1), in seconds.
# 6. **UT1-TT**: The difference between UT1 and TT, in seconds.
# 7. **LOD**: Length of Day, the excess length of a day in milliseconds.
# 8. **dX**: Nutation correction in the X-direction in microarcseconds.
# 9. **dY**: Nutation correction in the Y-direction in microarcseconds.

# Ensure that the provided EOP file adheres to this format. 

# If needed, use the [`Orient.prepare_eop`](@ref) utility function to generate a 
# file in the correct format before attempting to load the data.

# ### Loading EOP data 

# Once a EOP file compatible with the reader is avaliable, the data could be 
# loaded in the environment for later use by the [`FrameSystem`](@ref). 

# To initialize the EOP just run the [`Orient.init_eop`](@ref) function and 
# provide as input the previosly formatted file.