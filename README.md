# Microrobot-Leg-Control

<p align="center">
  <img src="https://github.com/Ali-7800/Microrobot-Leg-Control/blob/main/img/model.png" />
</p>

## Main Functions
### Get Parameters Function ```get_params()```
This is used by all the other functions and scripts to retrieve parameters like electrical resistivity or conduction coefficient, lengths such as beam width or V-shape actuator span or coordinates for the contact points.

### Forward Kinematic Function ```fcn_FK(dV1,dV2)```
This function takes the voltage differences on each V-shape actuator and calculates the leg tip displacements dx_c and dy_c.

### Inverse Kinematic Function ```fcn_IK(dx_c,dy_c)```
This function takes the voltage differences on each V-shape actuator and calculates the leg tip displacements dx_c and dy_c.

## Main Scripts
### Workspace genration ```workspace_gen.m```
This scripts generates the three workspace (No Contact, Unamplified, and Amplified) Using Either IK (Slow) or FK (Fast).
**To generate the workspaces make sure to have the correct maximum voltage value before running the script**

#### Steps to generate the workspaces
1. Run the first section to generate the No Contact workspace
2. Run either the FK section or IK section and rename the output to with contact (deltas_with_contact = deltas for FK and space_with_contact = space)
3. Disable the contact model by commenting the block highlighted in fcn_FK or fcn_IK
4. Run the same section you ran in part 2 again and rename the output to without contact (deltas_without_contact = deltas for FK and space_without_contact = space)

### Workspace plotting ```workspace.m```
This script plots the workspace and shows it just like the image below
<p align="center">
  <img src="https://github.com/Ali-7800/Microrobot-Leg-Control/blob/main/img/workspace.png" />
</p>

Run the section corresponding to section you ran in ```workspace_gen.m``` to plot the three workspaces.

In the third section of this script you can plot a trajectory onto the workspaces to see if it is within them by changing the parametric equations x and y to the desired trajectory.

### Trajectory generation ```trajectory_gen.m```
This calculates the voltages needed to make the leg tip move in a trajectory.
To calculate the voltages just run the first section to create matrix to store the voltages and then add a section following the format given in the script using the parametric equations x and y for your trajectory.

## Side Functions

### Angle Range Between Contacts Function ```fcn_contacts(dx_c,dy_c)```
Outputs the range of possible alpha values [alpha_min,alpha_max] between the contacts for a given (dx_c,dy_c), used by ```fcn_FK``` and ```fcn_IK``` to determine when to use the contact model.

### Contact Check ```fcn_contactCheck(dx_c,dy_c)```
Checks if the point (dx_c,dy_c) makes the leg touch any of the contact points.

### La(N2)
Calculates $&Lambda_{a}$ given $N^{2}$





