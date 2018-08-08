# skarn_model

# 
author: Evan J. Ramos

#
publication: Ramos, E.J., Hesse, M.A., Barnes, J.D., Jordan, J.S., and Lackey, J.S., 2018. Re-evaluating fluid
             sources during skarn formation: an assessment of the Empire Mountain skarn, Sierra Nevada, USA.
             Geochemistry, Geophysics, Geosystems.

#
In this repository, you will find all scripts and functions that you can use to replicate model results found in the above
publication. The main script is found in the file skarn_main.m. Additionally, we provide a more generic model where
the user can model hydrothermal fluid flow and oxygen isotope transport under the conditions they specificy (e.g.,
different host rock permeability, different intrusion size, etc.). This script is found in the file skarn_generic.m.
Lastly, to generate figures and process model outputs, you can call the script skarn_post_process.m. This file gives example
calls but can be tailored to your specification.

# NOTE 1
All other functions will be called in the scripts. Feel free to tinker with the code in the functions, but
some unbeknownst errors may arise in part because of this.

# NOTE 2
All code (in scripts and functions) include line-by-line comments on what the line (or section of code) accomplishes.

# NOTE 3 
If you have any comments, queries, or suggestions, either contact me directly through GitHub or send me an email at ejramos@utexas.edu with the subject line "Ramos et al 2018 G3 model question".
