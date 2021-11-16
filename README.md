# Simulation_verify
This is used to verify the simulation repository is right.  
First, I used the differential cross section to generate sufficient events, 
and then we can obtain a list like (Energy_7Li, theta_c.m., phi_c.m.), 
where the angle phi_c.m. can be randomly filled among all the angles(0-2PI).
further saved in a Root Tree.  
Second, Judging the generated events whether enters in detector surface,
if meet the conditions, the event can be as effective, and then fill in histogram.
Third, We should consider the angular resolution, But now, I have no idea about that.

# Usage
You shoud do 5 steps as follows:  
1. mkdir build  
2. cd build   
3. cmake ..  
4. make -jN  
5. a.out  
