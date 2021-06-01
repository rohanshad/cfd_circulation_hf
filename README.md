# Accessory Scripts

Code repository for accessory scripts for the paper "Patient-Specific Computational Fluid Dynamics Reveal Localized Flow Patterns Predictive of Post–LVAD Aortic Incompetence". 
Coming soon in Circulation Heart Failure: 2021;14:e008034. DOI: 10.1161/CIRCHEARTFAILURE.120.008034

<img width="2039" alt="Screen Shot 2021-05-05 at 12 42 56 PM" src="https://user-images.githubusercontent.com/12033026/117199756-6f153c00-ad9f-11eb-8d0d-7d55fc6b7cfc.png">


- __234_face_are.txt__ : Tuned parameters `Rp, C, Rd` for a single case
- __234_test.flow__: Patient specific flow waveform
- __compute_rcr_values.py__: Script to tune `Rp, C, Rd` for a single case
- __position_vertices.py__: Helper script to place aortic valve geometry into aortic root model
- __solver.inp__: Simulation settings used for main flow solver setup
- __particle_track.txt__: Simulation settings for advection-diffusion solver
- __flowmapper.R__: R script for generating patient specific flow curves and calculate periodicity of flow
