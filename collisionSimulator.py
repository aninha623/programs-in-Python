#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 16:05:56 2023


"""
import math

from simVis1 import visualize




def main():
    
    # initialize parameters
    
    radius =0.15
    dt =.2
    
    # simulation steps
    simSteps = 800
    

    # Default is Ax1xmin= 0,Ax1xmax = 1, Ax1ymin = 0, Ax1ymax = 1
    vis = visualize() 
    

    x1 = 0.2
    y1 = 0.25
    
    v1x = -0.05
    v1y = -0.05
    
    x2 = 0.7
    y2 = 0.7
    
    
    v2x = 0.05
    v2y = 0.05
    
    
    # run simulation
    simulate(simSteps,vis,x1,y1,x2,y2,v1x,v1y,v2x,v2y,dt,radius)
    


# DO NOT CHANGE THIS FUNCTION
def boundary_locations(vis,radius):
    
    """
    Calculate the boundary locations within the visualization window.
    This function calculates the boundary coordinates for a containing box 
    within the visualization window to ensure that particles 
    do not exceed these boundaries.

    Parameters
    ----------
    vis : visualization object
        The visualization object used for determining window boundaries.
    radius : float
        The radius of the particles.

    Returns
    -------
    xLow : float
        The lower bound of the containing box's threshold to hit the vertical left wall.
    yLow : float
        The lower bound of the containing box's threshold to hit the horizontal bottom wall.
    xHigh : float
        The upper bound of the containing box's threshold to hit the vertical right wall.
    yHigh : float
        The upper bound of the containing box's threshold to hit the horizontal top wall.

    """
    
    # Adjust the boundary according to the 
    # circle radius and display window size
    # Assume circles have the same radius.
    xLow = vis.Ax1xmin + radius 
    xHigh = vis.Ax1xmax - radius 
    
    yLow = vis.Ax1ymin + radius 
    yHigh = vis.Ax1ymax - radius 
    
    return xLow, yLow, xHigh, yHigh
    
    
def get_distance(x1,y1,x2,y2):
    
    """
    Calculate the Euclidean distance between two points based on their coordinates (x,y)
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first point
    y1: float
        y-coordinate of the first point
    x2 : float
        x-coordinate of the second point
    y2: float
        y-coordinate of the second point
    
    Returns
    ------
    distance_btwn_pts: float
        The distance between the two points rounded to 3 decimal places
    
    """
    
    delta_x_sq = (x2 - x1)**2
    delta_y_sq = (y2 - y1)**2
    distance_btwn_pts = math.sqrt(delta_x_sq + delta_y_sq)
    distance_btwn_pts = round(distance_btwn_pts, 3)
    return distance_btwn_pts

def updateX(x,vx,dt):
    
    """
    Update x coordinate of an object's position based on its velocity and time interval
    
    Parameters
    ------
    x : float
        x-coordinate of the object's position
    vx: float
        velocity of the object in the x direction
    dt: float
        time step for the update
    
    Returns
    ------
    x_updated : float
        new x-coordinate of the object's position determined by how fast it has been moving over a designated interval
    """
    
    x_updated= x + (vx * dt)
    x_updated= round(x_updated, 3)
    return x_updated

def updateY(y,vy,dt):
    
    """
    Update y coordinate of an object's position based on its velocity and time interval
    
    Parameters
    ------
    y : float
        y-coordinate of the object's position
    vy: float
        velocity of the object in the y direction
    dt: float
        time step for the update
    
    Returns
    ------
    y_updated : float
        new y-coordinate of the object's position determined by how fast it has been moving over a designated interval
    """
    
    y_updated= y + (vy * dt)
    y_updated= round(y_updated, 3)
    return y_updated

def boxCollision(x,y,vx,vy,xLow,xHigh,yLow,yHigh):
    
    """
    Check if a circle with a center at point (x,y) that is moving with velocity (vx,vy) hits the walls of a containing box
    with dimensions defined by thresholds. 
    If a collision occurs, the position and velocity of the object is recalculated based on wall reflections
    
    Parameters
    ------
    x : float
        current x coordinate of the circle
    y: float
         current y coordinate of the circle
    vx: float
        current x component of the circle's velocity
    vy: float
        current y component of the circle's velocity
    xLow: float
        The lower bound of the containing box's threshold to hit the vertical left wall.
    xHigh: float
        The upper bound of the containing box's threshold to hit the vertical right wall.
    yLow: float
        The lower bound of the containing box's threshold to hit the horizontal bottom wall.
    yHigh: float
        The upper bound of the containing box's threshold to hit the horizontal top wall.
    
    Returns
    ------
    x: float
        updated x coordinate of the circle following collision
    y: float
        updated y coordinate of the circle following collision
    vx: float
        updated x component of the circle's velocity following collision
    vy: float
        updated y component of the circle's velocity following collision
    
    """
    #if x surpasses the value of the threshold, it hits the
    #box, and velocity and position is recalculated accordingly
    if x < xLow:
        x = xLow
        vx = -vx
    elif x > xHigh:
        x = xHigh
        vx = -vx
    elif y < yLow:
        y = yLow
        vy = -vy
    elif y > yHigh:
        y = yHigh 
        vy = -vy
    return x,y,vx,vy

def Overlap(x1,y1,radius1,x2,y2,radius2):
    
    """
    Determine whether 2 circles overlap/collide or not
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first circle's center 
    y1: float
        y-coordinate of the first circle's center 
    radius1: float
        radius of the first circle
    x2 : float
        x-coordinate of the second circle's center 
    y2: float
        y-coordinate of the second circle's center 
    radius2: float
        radius of the second circle
    
    Returns
    ------
    True (Boolean) if the two circles overlap/collide
    False (Boolean) if the two circles do not
    
    """
    
    sum_radii= radius1 + radius2
    return get_distance(x1,y1,x2,y2) < sum_radii

def get_unit_direction(x1,y1,x2,y2):
   
    """
    Calcuate the unit direction vector between two points given their x and y coordinates
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first point
    y1: float
        y-coordinate of the first point
    x2 : float
        x-coordinate of the second point
    y2: float
        y-coordinate of the second point
    
    Returns
    ------
    tuple x_unit, y_unit (floats) representing the unit direction vector rounded to 3 decimal points
    
    """
    
    distance_btwn_pts = get_distance(x1, y1, x2, y2)
    if distance_btwn_pts != 0:
        x_unit= round((x1 - x2)/distance_btwn_pts, 3)
        y_unit= round((y1 - y2)/distance_btwn_pts, 3)
    else:
        #if points overlap/distance is zero unit direction 
        #is undefined = math.nan
        x_unit= math.nan
        y_unit= math.nan
    return (x_unit, y_unit)

def update_collision_velocity(x1,y1,v1x,v1y, x2,y2,v2x,v2y):
    
    """
    Calculate updated velocities of two circles following a collision
    using their initial positions and initial velocities
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first circle's center 
    y1: float
        y-coordinate of the first circle's center 
    v1x: float
        x component of the current velocity vector of the first circle
    v1y: float
        y component of the current velocity vector of the first circle
    x2 : float
        x-coordinate of the second circle's center 
    y2: float
        y-coordinate of the second circle's center 
    v2x: float
        x component of the current velocity vector of the second circle
    v2y: float
        y component of the current velocity vector of the second circle
    
    Returns
    ------
    All rounded to 3 decimal places: 
        
    V1x: float
        updated x component of the velocity of the first circle
    V1y: float
        updated y component of the velocity of the first circle
    V2x: float
        updated x component of the velocity of the second circle
    V2y: float
        updated y component of the velocity of the second circle
    
    """
    
    numerator1= ((v1x - v2x)*(x1 - x2)) + ((v1y - v2y)*(y1 - y2))
    numerator2= ((v2x - v1x)*(x2 - x1)) + ((v2y - v1y)*(y2 - y1))
    denominator= get_distance(x1, y1, x2, y2)**2
    #update velocities
    if denominator == 0:
        #updated velocities will be negative but have the same magnitude
        V1x = -v1x
        V1y = -v1y
        V2x = -v2x
        V2y = -v2y
    else: 
        V1x = v1x - ((numerator1 / denominator)*(x1 - x2))
        V1y = v1y - ((numerator1 / denominator)*(y1 - y2))
        V2x = v2x - ((numerator2 / denominator)*(x2 - x1))
        V2y = v2y - ((numerator2 / denominator)*(y2 - y1))
        V1x = round(V1x, 3)
        V1y = round(V1y, 3)
        V2x = round(V2x, 3)
        V2y = round(V2y, 3)
    return V1x, V1y, V2x, V2y
    
def dot_product(x1,y1,x2,y2):
    
    """
    Calcuate the dot product of two 2D vectors based on their x and y components
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first point
    y1: float
        y-coordinate of the first point
    x2 : float
        x-coordinate of the second point
    y2: float
        y-coordinate of the second point
    
    Returns
    ------
    final_dot_product: float
        a single number representing the dot product of the two vectors
    
    """
    
    x_product= x1 * x2
    y_product= y1 * y2
    final_dot_product= x_product + y_product
    return final_dot_product

    
def circleCollision(x1,y1,radius1,v1x, v1y, x2,y2,radius2,v2x,v2y):
     
    """
    Simulate the collision of the two circles and their subsequent reactions. 
    Calculates the new positions and velocities of the circles after they collide.
    Once they collide, the first circle is shifted over proportionally to its position 
    and velocity such that the two circles are no longer overlapping. 
    
    Parameters
    ------
    x1 : float
        x-coordinate of the first circle's center 
    y1: float
        y-coordinate of the first circle's center 
    radius1: float
        radius of the first circle
    v1x: float
        x component of the current velocity vector of the first circle
    v1y: float
        y component of the current velocity vector of the first circle
    x2 : float
        x-coordinate of the second circle's center 
    y2: float
        y-coordinate of the second circle's center 
    radius2: float
        radius of the second circle
    v2x: float
        x component of the current velocity vector of the second circle
    v2y: float
        y component of the current velocity vector of the second circle
    
    Returns
    ------
    All rounded to 3 decimal places: 
        
    updated_x1: float
        updated x coordinate of the center of the first circle
    updated_y1: float
        updated y coordinate of the center of the first circle
    uv1x: float
        updated x component of the velocity of the first circle
    uv1y: float
        updated y component of the velocity of the first circle
    x2: float
        x coordinate of the center of the second circle, does not change
    y2: float
        y coordinate of the center of the second circle, does not change
    uv2x: float
        updated x component of the velocity of the second circle
    uv2y: float
        updated y component of the velocity of the second circle
    
    """
    
    if Overlap(x1, y1, radius1, x2, y2, radius2) == True :
        sum_of_radii = radius1 + radius2
        #update velocity
        uv1x, uv1y, uv2x, uv2y = update_collision_velocity(x1,y1,v1x,v1y, x2,y2,v2x,v2y)
        #update direction
        udirection_x, udirection_y = get_unit_direction(x1, y1, x2, y2)
        #update position
        #if the circles overlap completely, shift over according to fixed value:
        if x1 == x2 and y1 == y2:
            direction_x = 0.707
            direction_y = 0.707
            updated_x1 = x1 + (math.fabs(get_distance(x1, y1, x2, y2) - sum_of_radii) * direction_x) 
            updated_x1 = round(updated_x1, 3)
            updated_y1 = y1 + (math.fabs(get_distance(x1, y1, x2, y2) - sum_of_radii) * direction_y) 
            updated_y1 = round(updated_y1, 3)
        #if the circles don't overlap completely, shift over proportionally to overlap:
        else:
            updated_x1 = x1 + (math.fabs(get_distance(x1, y1, x2, y2) - sum_of_radii) * udirection_x) 
            updated_x1 = round(updated_x1, 3)
            updated_y1 = y1 + (math.fabs(get_distance(x1, y1, x2, y2) - sum_of_radii) * udirection_y) 
            updated_y1 = round(updated_y1, 3)
    #if they do not overlap, return initial conditions
    else:
        updated_x1=x1
        updated_y1=y1
        uv1x=v1x
        uv1y=v1y
        x2=x2
        y2=y2
        uv2x=v2x
        uv2y=v2y
    return updated_x1, updated_y1, uv1x, uv1y, x2, y2, uv2x, uv2y
        

# DON NOT CHANGE THIS FUNCTION
def simulate(simSteps,vis,x1,y1,x2,y2,v1x,v1y,v2x,v2y,dt,radius):
    
    """
    Run a particle simulation for a specified number of time steps.

    Parameters
    ----------
    simSteps : int
        Number of time steps to run the simulation.
    vis : visualize
        An instance of the 'visualize' class for visualization.
    x1 : float
        Initial x-coordinate of the center of circle 1.
    y1 : float
        Initial y-coordinate of the center of circle 1.
    x2 : float
        Initial x-coordinate of the center of circle 2.
    y2 : float
        Initial y-coordinate of the center of circle 2.
    v1x : float
        Initial x-component of the velocity of circle 1.
    v1y : float
        Initial y-component of the velocity of circle 1.
    v2x : float
        Initial x-component of the velocity of circle 2.
    v2y : float
        Initial y-component of the velocity of circle 2.
    dt : float
        Time step size for the simulation.
    radius : float
        Radius of both particles.

    Returns
    -------
    None

    """
    
    # This function runs a two-particle simulation for a specified number of 
    # time steps. It updates the positions of two circles and checks for 
    # their collisions with both the container boundary and each other. 
    # The simulation is visualized using the'visualize' class instance 
    # provided as 'vis'.
    
    
    # Identify the boundary location for the containing window
    # Assumption: both particles have the same radius.
    xLow,yLow,xHigh, yHigh = boundary_locations(vis,radius) 
    
    for i in range(simSteps):
        
            
    
    
        # update x, y position of particles
        x1 = updateX(x1,v1x,dt)
        y1 = updateY(y1,v1y,dt)
        
        x2 = updateX(x2,v2x,dt)
        y2 = updateY(y2,v2y,dt)
            
        
                    
        # Check for box and circle collisions
        
        x1,y1,v1x,v1y = boxCollision(x1,y1,v1x,v1y,xLow,xHigh,yLow,yHigh)
        
        
        
        x2,y2,v2x,v2y = boxCollision(x2,y2,v2x,v2y,xLow,xHigh,yLow,yHigh)
        
        
        x1,y1,v1x,v1y, x2,y2,v2x,v2y = circleCollision\
                                      (x1,y1,radius,v1x, \
                                        v1y, x2,y2,radius,v2x,v2y)
                    
        
                    
        # draw circles
        vis.circle(x1, y1, radius, 0)
        vis.circle(x2, y2, radius, 1)
        
                    
        # pause plots and clear window axis 1
        vis.plotPause()
        
        
        vis.axis1Clear()
    
    # redraw circles after last iteration
    vis.circle(x1, y1, radius, 0)
    vis.circle(x2, y2, radius, 1)
        
if __name__ == '__main__': 
    
    # Call the main function to excute the simulation
    # For testing purpose comment main() and call another 
    # function.
    main()
    
    
    