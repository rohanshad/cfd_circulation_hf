# Accessory Scripts

Code repository for accessory scripts for CFD and Aortic Incompetence in LVAD patients paper. 
Published in Circulation Heart Failure: (add link here) 

- __234_face_are.txt__ : Tuned parameters `Rp, C, Rd` for a single case
- __234_test.flow__: Patient specific flow waveform
- __compute_rcr_values.py__: Script to tune `Rp, C, Rd` for a single case
- __position_vertices.py__: Helper script to place aortic valve geometry into aortic root model
- __solver.inp__: Simulation settings used for main flow solver setup
- __particle_track.txt__: Simulation settings for advection-diffusion solver
- __flowmapper.R__: R script for generating patient specific flow curves and calculate periodicity of flow
