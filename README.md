# Topology-Optimization

The aircraft pylon is a typical cantilever structure hanging the turbine engine or other payloads to the wing of aircraft. The main objective of this study is to learn topology optimization and use the concept for aircraft pylon structure to reduce the mass maintaining its stiffness. The design is carried out in SolidWorks CAD and later imported to Altair Inspire for topology optimization. 

A simplified geometry is designed in SolidWorks as shown below with a length of the pylon around 3 m and mass around 595.11 kg. 

The geometry is simplified for the study with the fixed boundary condition at the top four nodes connected to the wing of the aircraft. The inertial load, thrust and weight of the turbine were applied at the bottom four nodes connected to the turbine engine. Five load cases were assessed including - 1) Steady flight, 2) Accent, 3) Decent, 4) Steady flight and right turn, 5) Steady flight and left turn. 

The von-mises stress results with the extreme load cases without topology optimization is shown. 

With the topology optimization, the mass of the pylon is reduced to 159.19 kg which implies 73% reduction. The von-mises stress results with the topology optimization is shown for the extreme load case. 

The topology optimization code to design a cantilever beam. The code minimizes the compliance of the cantilever.
The SIMP method with gradient filtering is employed to reduce mesh dependence and checkerboarding. 
The Sigmund's 99 - line topology optimization code is used for this case. 
