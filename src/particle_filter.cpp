/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 20;
	
	double x_std = std[0];
	double y_std = std[1];
	double theta_std = std[2];
	for(int i = 0; i<num_particles; i++)
	{
		weights.push_back(1);
	}
	normal_distribution<double> dist_x(x,x_std);
	normal_distribution<double> dist_y(y,y_std);
	normal_distribution<double> dist_theta(theta,theta_std);
	default_random_engine gen;
	//init particles
	for(int i =0 ; i<num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;
		weights.push_back(p.weight);
		particles.push_back(p);	
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	double x_std = std_pos[0];
	double y_std = std_pos[1];
	double theta_std = std_pos[2];
	normal_distribution<double> dist_x(0,x_std);
	normal_distribution<double> dist_y(0,y_std);
	normal_distribution<double> dist_theta(0,theta_std);
	default_random_engine gen;
	std::vector<Particle>::iterator itr = particles.begin();	
	for(; itr != particles.end(); itr++)
	{
		/*bicycle mode has two mode for yaw_rate*/
		if(fabs(yaw_rate) >  0.0001)
		{
				//predict x position in particles
				itr->x = itr->x + velocity/yaw_rate*(sin(itr->theta + delta_t*yaw_rate) - sin(itr->theta));
				itr->x = itr->x + dist_x(gen);

				//predict y position in particles
				itr->y = itr->y + velocity/yaw_rate*(cos(itr->theta) - cos(itr->theta + delta_t*yaw_rate));
				itr->y = itr->y + dist_y(gen);

				//predict theta position in particles
				itr->theta = itr->theta + delta_t*yaw_rate;
				itr->theta = itr->theta + dist_theta(gen);
		}
		else
		{
				//yaw rate ==0 
				itr->theta = itr->theta + dist_theta(gen);
				itr->x = itr->x + velocity*cos(itr->theta)*delta_t;
				itr->x = itr->x + dist_x(gen);
				itr->y = itr->y + velocity*sin(itr->theta)*delta_t;
				itr->y = itr->y + dist_y(gen);

		}
	}	
}

void ParticleFilter::dataAssociation(std::vector<Map::single_landmark_s>& vaild_landmark, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i = 0; i < observations.size(); i++)
	{
		double min = 9999999;
		for(int j = 0; j < vaild_landmark.size(); j++)
		{
			double shortness = dist(observations[i].x,observations[i].y,vaild_landmark[j].x_f,vaild_landmark[j].y_f); 
			if(shortness < min)
			{
				observations[i].id = vaild_landmark[j].id_i;
				min = shortness; 
			}	
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	const std::vector<Map::single_landmark_s>& landmarks = map_landmarks.landmark_list;
	if(landmarks.empty())
		return;
	weights.clear();
    std::vector<Map::single_landmark_s> vaild_landmark;
	for(int i =0 ; i < num_particles; i++)
	{
		for(int j =0; j < landmarks.size(); j++)
		{
			if(dist(particles[i].x,particles[i].y,landmarks[j].x_f,landmarks[j].y_f) < sensor_range)
			{
				vaild_landmark.push_back(landmarks[j]);
			}
		}
		vector<LandmarkObs> obs_map;
		for(int j = 0; j < observations.size(); j++)
		{
			//transform veichel coridinate to map coridinate
			LandmarkObs o;
			double cos_theta = cos(particles[i].theta);
			double sin_theta = sin(particles[i].theta);
			o.x = particles[i].x + cos_theta*observations[j].x - sin_theta*observations[j].y;
			o.y = particles[i].y + sin_theta*observations[j].x + cos_theta*observations[j].y;
			obs_map.push_back(o);
		}	
		//associate the closet landmark and observation
		dataAssociation(vaild_landmark,obs_map);

	cout<<"---------"<<obs_map.size()<<"--"<<__LINE__<<endl;
		SetAssociations(particles[i],obs_map);
	cout<<"---------"<<__LINE__<<endl;
		double x_cft = 1.0 / (2.0 * std_landmark[0] * std_landmark[0]);
		double y_cft = 1.0 / (2.0 * std_landmark[1] * std_landmark[1]);
		double xy_cft = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);

		double prob = 1.0;
		for(int j=0;j<obs_map.size();++j)
		{       
				double tx=x_cft*pow((obs_map[j].x-landmarks[obs_map[j].id-1].x_f),2);
				double ty=y_cft*pow((obs_map[j].y-landmarks[obs_map[j].id-1].y_f),2);
				double p = xy_cft * exp(-1.0*(tx+ty));
				prob *= p;
		}
		particles[i].weight = prob;
		weights.push_back(prob);
	}	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
		std::default_random_engine gen;
		std::discrete_distribution<int> dist(weights.begin(), weights.end());

		std::vector<Particle> resample_p;

		for(int i = 0; i < num_particles; i++) {
				int index = dist(gen);
				resample_p.push_back(particles.at(index));
		}

		particles = resample_p;

}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<LandmarkObs>& obs_map)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		particle.associations.clear();
		particle.sense_x.clear();
		particle.sense_y.clear();
		for(int j=0;j<obs_map.size();++j) {
				particle.associations.push_back(obs_map[j].id);
				particle.sense_x.push_back(obs_map[j].x);
				particle.sense_y.push_back(obs_map[j].y);
		}
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
