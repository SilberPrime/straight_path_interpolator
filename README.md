# straight_path_interpolator

This function provides a straight path interpolation between two geographical points (inluding altitude) in 3Dspace in a standard coordinate system and the UTM coordinate system. 
The function was originally developed for the purpose of interpolating a meteoroid flight trajectory, but it can be used for other purposes too. 

Given the begin and end points in 3D space (latitude, longitude, altitude), the remaining points along the straight path between these two points are calculated. All results are written in a .txt file, and plotted in a 3D figure. The output also includes the horizontal range along the great circle, vertical distance, total distance, azimuth (as seen from the begin point and end point), and radiant.
