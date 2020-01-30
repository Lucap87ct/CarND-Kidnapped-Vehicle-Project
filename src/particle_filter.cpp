/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

  // set number of particles
  num_particles_ = 1000;

  std::default_random_engine noise_gen;
  std::normal_distribution<double> distrib_x(x, std[0]);
  std::normal_distribution<double> distrib_y(y, std[1]);
  std::normal_distribution<double> distrib_theta(theta, std[2]);
  Particle init_particle;

  // init all particles with random distribution with input noise around input
  // state
  for (int i = 0; i < num_particles_; i++) {
    init_particle.x = distrib_x(noise_gen);
    init_particle.y = distrib_y(noise_gen);
    init_particle.theta = distrib_theta(noise_gen);
    init_particle.weight = 1.0;
    particles.push_back(init_particle);
  }

  std::cout
      << "--------------------- Particles initialized ---------------------"
      << std::endl;

  is_initialized_ = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {

  std::default_random_engine noise_gen;

  // predict particles with input movement using bicycle model
  for (auto &particle : particles) {
    auto x_0 = particle.x;
    auto y_0 = particle.x;
    auto theta_0 = particle.x;

    if (abs(yaw_rate) > 0.001) {
      particle.x =
          x_0 + velocity *
                    ((sin(theta_0 + yaw_rate * delta_t) - sin(theta_0))) /
                    yaw_rate;
      particle.y =
          y_0 + velocity *
                    ((-cos(theta_0 + yaw_rate * delta_t) + cos(theta_0))) /
                    yaw_rate;
      particle.theta = theta_0 + yaw_rate * delta_t;

    } else {
      particle.x = x_0 - velocity * sin(theta_0);
      particle.y = y_0 + velocity * cos(theta_0);
    }

    // add random noise to particles using input noise
    std::normal_distribution<double> distrib_x(particle.x, std_pos[0]);
    std::normal_distribution<double> distrib_y(particle.y, std_pos[1]);
    std::normal_distribution<double> distrib_theta(particle.theta, std_pos[2]);
    particle.x = distrib_x(noise_gen);
    particle.y = distrib_y(noise_gen);
    particle.theta = distrib_theta(noise_gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no
   * scaling). The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}
