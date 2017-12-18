#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

//////////////////////
//                  //
// Constants/Macros //
//                  //
//////////////////////
#define INIT_BACK_PROJECTION_DELTATIME_S  (0.02)
#define FORWARD_PROJECTION_DELTATIME_S    (0.02)
#define SPLINE_PROJECTION                    (3)
#define FORWARD_PROJECTION_POINTS           (10)
#define MILESPERHOUR_TO_METERSPERSEC   (0.44704)
#define PATH_PROJECTION_POINTS              (25)
#define DIFF_SPEED                           (1)
#define REFERENCE_SPEED                   (48.0)
#define ALLOWABLE_DISTANCE                (40.0)
#define SPEED_BUFFER                       (2.5)
#define SPEED_REDUCTION_FACTOR             (0.8)
#define SPEED_INCREASE_FACTOR              (1.1)
#define LANE_WIDTH                         (4.0)
#define MIDDLE_LANE                          (1)
#define MAX_LANES                            (3)
#define MAX_FRENET_S                  (6945.554)

#define DEBUG_ENABLED                        (0)         


// Debug print macro
// Enabled/disabled with DEBUG_ENABLED #define
#define DEBUG(x) do {              \
  if (DEBUG_ENABLED)               \
  {                                \
    std::cerr << x << std::endl;   \
  }                                \
} while (0);



//////////////////////
//                  //
//      ENUMS       //
//                  //
////////////////////// 
typedef enum SensorFusionTypes
{
  ID = 0,
  GLOBAL_X = 1,
  GLOBAL_Y = 2,
  VELOCITY_X = 3,
  VELOCITY_Y = 4,
  FRENET_S = 5,
  FRENET_D = 6
} eSensorFusionTypes;


//////////////////////
//                  //
//    FUNCTIONS     //
//                  //
////////////////////// 

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) 
{
  string Res = "";
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    Res = "";
  }
  else if ((b1 != string::npos) && (b2 != string::npos)) 
  {
    Res = s.substr(b1, b2 - b1 + 2);
  }
  
  return Res;
}

// Calculate the distance between two cartesian points
double distance(double x1, double y1, double x2, double y2)
{
  double x_diff = x2 - x1;
  double y_diff = y2 - y1;
	return sqrt(x_diff*x_diff + y_diff*y_diff);
}

// Find the closest waypoint to an x/y cartesian point
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{
  // start with a large number for the closest waypoint distance
	double closestLen = 1e8;
  // Set the closest waypoint to invalid (i.e. not set)
	int closestWaypoint = -1;
  int i;

  // iterate through each waypoint and determine how close it is
	for(i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);

    // Compare with previously known closest and update if necessary
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}

	return closestWaypoint;
}

// Given an x/y cartesian coordinate, determine the next waypoint
// Note this is not necessarily the closest (as the closest could be behind)
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  double map_x, map_y, heading, angle;
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  if (-1 == closestWaypoint)
  {
    DEBUG("Closest waypoint not found!");
    goto Error;
  }

	map_x = maps_x[closestWaypoint];
	map_y = maps_y[closestWaypoint];

	heading = atan2((map_y - y), (map_x - x));

	angle = abs(theta - heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

Error:

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int i;
  int prev_wp;
  double n_x, n_y, x_x, x_y = 0.0;
  double proj_norm, proj_x, proj_y = 0.0;
  double frenet_d, frenet_s = 0.0;
  double center_x, center_y, centerToPos, centerToRef = 0.0;
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  if (-1 == next_wp)
  {
    DEBUG("Next waypoint not found!");
    goto Error;
  }

	prev_wp = next_wp - 1;
	if(0 == next_wp)
	{
		prev_wp  = maps_x.size() - 1;
	}

  n_x = maps_x[next_wp] - maps_x[prev_wp];
	n_y = maps_y[next_wp] - maps_y[prev_wp];
	x_x = x - maps_x[prev_wp];
	x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	proj_norm = ((x_x * n_x) + (x_y * n_y))/((n_x * n_x) + (n_y * n_y));
	proj_x = proj_norm * n_x;
	proj_y = proj_norm * n_y;

	frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point
	center_x = 1000 - maps_x[prev_wp];
	center_y = 2000 - maps_y[prev_wp];
	centerToPos = distance(center_x, center_y, x_x, x_y);
	centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	for( i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

Error:

	return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while((s > maps_s[prev_wp + 1]) && 
        (prev_wp < (int)(maps_s.size() - 1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));

	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s*cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s*sin(heading);

	double perp_heading = heading - pi()/2;

	double x = seg_x + (d * cos(perp_heading));
	double y = seg_y + (d * sin(perp_heading));

	return {x,y};
}


/* @brief Support wrapping of frenet_s coordinates
 *
 * @param input the frenet s value to wrap awhen crossing/finishing each lap.
 * @return normalised (wrapped) frenet value
 */
double normalise_frenet_s(double input)
{
  if (input > MAX_FRENET_S)
  {
    input -= MAX_FRENET_S;
  }
  return input;
}

// Determine if the car location is within a lane's tolerance
bool SameLane(double frenet_d, int lane)
{ 
  bool same_lane = false;

  double inside_edge_lane  = lane * LANE_WIDTH;
  double outside_edge_lane = (lane + 1) * LANE_WIDTH;

  DEBUG("Testing same lane...");
  DEBUG("  Inside Edge  : " << inside_edge_lane);
  DEBUG("  Outside Edge : " << outside_edge_lane);
  DEBUG("  Test Frenet  : " << frenet_d);
  // Assumes that the frenet_d value increases moving away from the center of the road
  if ((frenet_d > inside_edge_lane) && 
      (frenet_d < outside_edge_lane))
  {
    same_lane = true;
  }

  return same_lane;
}

// Determine if the focus car (FC) in front is a danger to our car (OC)
//   - Check the position of the car and the lane of the car
//   - Check the speed values of each car
bool IsCarFrontDanger(double OC_frenet_s, int OC_lane, double OC_speed, vector<double> sensor_data)
{
  // in danger?
  bool in_danger = false;
  bool forward_car = false;

  // Window for worrying about whether the other car is too close in front.
  double front_window = normalise_frenet_s(OC_frenet_s + ALLOWABLE_DISTANCE);

  double FC_frenet_s = sensor_data[FRENET_S];
  double FC_frenet_d = sensor_data[FRENET_D];

  DEBUG("Testing if car is in the same lane...");
  DEBUG("  OC Frenet S    : " << OC_frenet_s);
  DEBUG("  Forward window : " << front_window);
  DEBUG("  FC Frenet S    : " << FC_frenet_s);

  // Is the FC within the forward distance window?
  if (front_window < OC_frenet_s)
  {
    // wrapping for another lap of track
    if ((OC_frenet_s < FC_frenet_s) || 
        (front_window > FC_frenet_s))
    {
      forward_car = true;
    }
  }
  else if ((OC_frenet_s < FC_frenet_s) && 
           (front_window > FC_frenet_s))
  {
    forward_car = true;
  }

  if (forward_car)
  {
    DEBUG("## FC in Forward Window ##");
    // Are we in the same lane?
    if (SameLane(FC_frenet_d, OC_lane))
    {
      DEBUG("## FC in Same Lane ##");

      // Set the allowable speed to the vehicles speed minus a little
      double FC_speed_x = sensor_data[VELOCITY_X];
      double FC_speed_y = sensor_data[VELOCITY_Y];
      double FC_speed = sqrt((FC_speed_x * FC_speed_x) + (FC_speed_y * FC_speed_y));

      // Calculate a speed that is less than the car in front (to prevent a collision)
      double OC_speed_allow = FC_speed - DIFF_SPEED;
      
      DEBUG("Testing if car needs to slow down...")
      DEBUG("  OC Speed : " << OC_speed);
      DEBUG("  FC Speed : " << FC_speed);
      DEBUG("  OC Speed allow : " << OC_speed_allow);

      if (OC_speed_allow < OC_speed)
      {
        DEBUG("## FC going too slow ##");
        in_danger = true;
      }
      else
      {
        DEBUG("## FC going fast enough ##");
      }
    }
    else
    {
      DEBUG("## FC not in Same Lane ##")
    }
  }
  else
  {
    DEBUG("## FC not in Forward Window ##")
  }

  return in_danger;
}

// Determine if the focus car (FC) behind is a danger to our car (OC)
//   - Check the position of the car and the lane of the car
//   - Check the speed values of each car
bool IsCarBehindDanger(double OC_frenet_s, int OC_lane, double OC_speed, vector<double> sensor_data)
{
  // in danger?
  bool in_danger = false;

  // Window for worrying about whether the other car is too close behind.
  double backward_window  = (OC_frenet_s - ALLOWABLE_DISTANCE/2);

  double FC_frenet_s = sensor_data[FRENET_S];
  double FC_frenet_d = sensor_data[FRENET_D];

  DEBUG("Testing if car is in the same lane...");
  DEBUG("  OC Frenet S     : " << OC_frenet_s);
  DEBUG("  Backward window : " << backward_window);
  DEBUG("  FC Frenet S     : " << FC_frenet_s);

  // Is the FC within the backward distance window?
  if ((OC_frenet_s > FC_frenet_s) && 
      (backward_window < FC_frenet_s))
  {
    DEBUG("## FC in Backward Window ##");
    // Are we in the same lane?
    if (SameLane(FC_frenet_d, OC_lane))
    {
      DEBUG("## FC in Same Lane ##");

      in_danger = true;
    }
    else
    {
      DEBUG("## FC not in Same Lane ##")
    }
  }
  else
  {
    DEBUG("## FC not in Backward Window ##")
  }

  return in_danger;
}


// Convert to the car reference coordinate system
vector<double> convertToCarRef(double pos_x, double pos_y, double ref_x, double ref_y, double yaw)
{
  double output_x, output_y = 0.0;

  double car_ref_x = pos_x - ref_x;
  double car_ref_y = pos_y - ref_y;

  output_x = (car_ref_x * cos(-yaw)) - (car_ref_y * sin(-yaw));
  output_y = (car_ref_x * sin(-yaw)) + (car_ref_y * cos(-yaw));

  return {output_x, output_y};  
}

// Convert from the car reference coordinate system to the global coordinates
vector<double> convertFromCarRef(double pos_x, double pos_y, double ref_x, double ref_y, double yaw)
{
  double output_x = ref_x + (pos_x * cos(yaw)) - (pos_y * sin(yaw));
  double output_y = ref_y + (pos_x * sin(yaw)) + (pos_y * cos(yaw));

  return {output_x, output_y};
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0


  // Lane the vehicle is in - start in MIDDLE LANE
  int OC_lane = MIDDLE_LANE;
  double OC_reference_speed = (REFERENCE_SPEED - SPEED_BUFFER);

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);

  }
  double max_frenet_s = MAX_FRENET_S;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&OC_lane,&OC_reference_speed,&max_frenet_s](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") 
        {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = deg2rad(j[1]["yaw"]);
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          // Previous and Next vectors
          vector<double> spline_x_vals;
          vector<double> spline_y_vals;
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int i = 0;
          bool front_car_danger = false;
          bool try_lane_change = false;

          // Print some debug information
          DEBUG("OC STATE INFORMATION...");
          DEBUG("  OC Current Lane : " << OC_lane);
          DEBUG("  OC Pos X        : " << car_x);
          DEBUG("  OC Pos Y        : " << car_y);
          DEBUG("  OC Frenet S     : " << car_s);
          DEBUG("  OC Car Frenet D : " << car_d);
          DEBUG("  OC Car Yaw      : " << car_yaw);
          DEBUG("  OC Car Speed    : " << car_speed);

          // Is there a collision likely
          for (i = 0; i < sensor_fusion.size(); i++)
          {
            // check if the vehicle is in the same lane, in front and travelling too slow
            front_car_danger |= IsCarFrontDanger(car_s, OC_lane, car_speed, sensor_fusion[i]);
          }
          
          if (front_car_danger)
          {
            // Need to slow down
            DEBUG("## Slowing Down ##");
            car_speed -= 1;

            // Since we're obstructed, perhaps try a lane change
            try_lane_change = true;
          }
          else
          {
            // Can speed up if not at max speed
            if (car_speed < OC_reference_speed)
            {
              DEBUG("## Speeding up ##")
              car_speed += 1;
            }
            else
            {
              DEBUG("## Crusing ##")
            }
          }
                                    
          // Now we've altered speed, should we try a lane change?
          if (try_lane_change)
          {
            DEBUG("Try Lane Change...")
            bool left_lane_clear  = true;
            bool right_lane_clear = true;

            int OC_left_lane  = OC_lane - 1;
            int OC_right_lane = OC_lane + 1;

            DEBUG("  OC lane       : " << OC_lane);
            DEBUG("  OC left lane  : " << OC_left_lane);
            DEBUG("  OC right lane : " << OC_right_lane);

            for (i = 0; i < sensor_fusion.size(); i++)
            {
              // Check if left lane is clear
              // but only if there is a left lane
              if (OC_left_lane >= 0)
              {
                left_lane_clear &= !IsCarFrontDanger(car_s, OC_left_lane, OC_reference_speed, sensor_fusion[i]);
                left_lane_clear &= !IsCarBehindDanger(car_s, OC_left_lane, OC_reference_speed, sensor_fusion[i]);
              }
              else
              {
                left_lane_clear = false;
              }

              // Check if right lane is clean
              // but only if there is a right lane
              if (OC_right_lane < MAX_LANES)
              {
                right_lane_clear &= !IsCarFrontDanger(car_s, OC_right_lane, OC_reference_speed, sensor_fusion[i]);
                right_lane_clear &= !IsCarBehindDanger(car_s, OC_right_lane, OC_reference_speed, sensor_fusion[i]);
              } 
              else
              {
                right_lane_clear = false;
              }
            }

            
            // Can we change lanes?
            // Preferencing left lane
            if (left_lane_clear)
            {
              DEBUG("## Left lane clear ##")
              OC_lane = OC_left_lane;
            }
            else if (right_lane_clear)
            {
              DEBUG("## Right lane clear ##")
              OC_lane = OC_right_lane;
            }
          }

          // Now we've figured out the state of play, what lane we should be in and the speed
          // determine the waypoints that get us there.
          double OC_ref_x = car_x;
          double OC_ref_y = car_y;
          double OC_x_prev, OC_y_prev, OC_x_prev_prev, OC_y_prev_prev = 0.0;

          // Do we have some left over waypoints?
          if (previous_path_x.size() < 2)
          {
            // Bootstrap the "previous path" with some constructed data
            // Back-pedal the car from the initial state based on current speed and yaw
            // Projecting back assuming no acceleration
            double OC_backward_projection = OC_reference_speed * INIT_BACK_PROJECTION_DELTATIME_S * MILESPERHOUR_TO_METERSPERSEC;
            
            // move the car backward in the direction of the current steering yaw
            double OC_backward_x = OC_backward_projection * cos(car_yaw);
            double OC_backward_y = OC_backward_projection * sin(car_yaw);

            DEBUG("Re-Initialise spline...");
            DEBUG("  Dummy back projection : " << OC_backward_projection);
            DEBUG("  Dummy back position   : (" << OC_backward_x << ", " << OC_backward_y << ")");

            OC_x_prev = OC_ref_x - OC_backward_x;
            OC_y_prev = OC_ref_y - OC_backward_y;

            OC_x_prev_prev = OC_x_prev - OC_backward_x;
            OC_y_prev_prev = OC_y_prev - OC_backward_y;

            end_path_s = car_s;

          }
          else
          { 
            // The state size has to be greater than 2 unless the system is in the initial state.
            OC_ref_x = previous_path_x[previous_path_x.size() - 1];
            OC_ref_y = previous_path_y[previous_path_y.size() - 1];

            // Use the last 2 previous state values for the 
            OC_x_prev = previous_path_x[previous_path_x.size() - 1];
            OC_y_prev = previous_path_y[previous_path_y.size() - 1];

            OC_x_prev_prev = previous_path_x[previous_path_x.size() - 2];
            OC_y_prev_prev = previous_path_y[previous_path_y.size() - 2];
          }

          // prime the spline with the previous state vector
          vector<double> OC_prev      = convertToCarRef(OC_x_prev, OC_y_prev, OC_ref_x, OC_ref_y, car_yaw);
          vector<double> OC_prev_prev = convertToCarRef(OC_x_prev_prev, OC_y_prev_prev, OC_ref_x, OC_ref_y, car_yaw);
          spline_x_vals.push_back(OC_prev_prev[0]);
          spline_y_vals.push_back(OC_prev_prev[1]);

          spline_x_vals.push_back(OC_prev[0]);
          spline_y_vals.push_back(OC_prev[1]);

          DEBUG("Creating waypoint spline...");
          DEBUG("  Add spline prev 2 : (" << OC_x_prev << ", " << OC_y_prev << ")");
          DEBUG("  Add spline prev 1 : (" << OC_x_prev_prev << ", " << OC_y_prev_prev << ")");

          // Calculate point offsets to generate a spline
          int OC_frenet_d = LANE_WIDTH/2 + OC_lane * LANE_WIDTH;
          for (i = 1; i <  SPLINE_PROJECTION; i++)
          {
            vector<double> OC_next_point = getXY(end_path_s + i * 3 * FORWARD_PROJECTION_POINTS, OC_frenet_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
  
            DEBUG("  Add spline next " << i << " : (" << OC_next_point[0] << ", " << OC_next_point[1] << ")");

            vector<double> OC_next_point_cr = convertToCarRef(OC_next_point[0], OC_next_point[1], OC_ref_x, OC_ref_y, car_yaw);
            spline_x_vals.push_back(OC_next_point_cr[0]);
            spline_y_vals.push_back(OC_next_point_cr[1]);
          }


          // Generate a spline to fit the path points
          tk::spline OC_path_spline;
          OC_path_spline.set_points(spline_x_vals, spline_y_vals);
          double last_x_pos = 0;


          double targetXDistance      = FORWARD_PROJECTION_POINTS * SPLINE_PROJECTION;
          double targetYDistance      = OC_path_spline(targetXDistance);
          double targetDistance       = sqrt((targetXDistance * targetXDistance) + (targetYDistance * targetYDistance));
          double distanceRatio        = targetXDistance/targetDistance;
          double distanceInTimeStep   = distanceRatio * car_speed * FORWARD_PROJECTION_DELTATIME_S * MILESPERHOUR_TO_METERSPERSEC;

          DEBUG("Project waypoints...");
          DEBUG("  OC target distance X      : " << targetXDistance);
          DEBUG("  OC target distance Y      : " << targetYDistance);
          DEBUG("  OC target distance        : " << targetDistance);
          DEBUG("  OC distance X/total ratio : " << distanceRatio);
          DEBUG("  OC distance in time-step  : " << distanceInTimeStep);

          DEBUG("Creating waypoints...");
          // Add the previous path not yet traversed
          for (i = 0; i < previous_path_x.size(); i++)
          {
            DEBUG("  Add waypoints prev " << i << " : (" << previous_path_x[i] << ", " << previous_path_y[i] << ")");
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // fill up the projected path with the newly generated spline
          for (i = 0; i < (PATH_PROJECTION_POINTS - previous_path_x.size()); i++)
          {
            double x_offset = (i + 1) * distanceInTimeStep;
            double y_offset = OC_path_spline(x_offset);
            vector<double> OC_next_waypoint = convertFromCarRef(x_offset, y_offset, OC_ref_x, OC_ref_y, car_yaw);

            DEBUG("  Add waypoints next " << i << " : (" << previous_path_x[i] << ", " << previous_path_y[i] << ")");

            next_x_vals.push_back(OC_next_waypoint[0]);
            next_y_vals.push_back(OC_next_waypoint[1]);
          }


        	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
        	msgJson["next_x"] = next_x_vals;
        	msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";
          //std::cout << msg << std::endl;

        	//this_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } 
      else 
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
