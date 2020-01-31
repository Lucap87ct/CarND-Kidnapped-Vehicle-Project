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

    if (fabs(yaw_rate) > 0.001) {
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
      particle.x = x_0 + velocity * cos(theta_0);
      particle.y = y_0 + velocity * sin(theta_0);
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

void ParticleFilter::dataAssociation(const vector<LandmarkObs> landmarks,
                                     vector<LandmarkObs> &observations) {
  double dist_temp;
  for (auto &observation : observations) {
    double min_dist = 1e3;
    for (const auto &landmark : landmarks) {
      dist_temp = std::sqrt(pow(observation.x - landmark.x, 2) +
                            pow(observation.y - landmark.y, 2));
      if (dist_temp < min_dist) {
        min_dist = dist_temp;
        observation.id = landmark.id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {

  double obs_x_temp, obs_y_temp;
  double landmark_x_temp, landmark_y_temp;
  double exponent, update_prob_temp;

  // loop on all particles
  for (auto &particle : particles) {

    // transform observations in map frame at particle state
    vector<LandmarkObs> observations_in_map_frame;
    for (const auto &observation : observations) {
      obs_x_temp = particle.x + cos(particle.theta) * observation.x -
                   sin(particle.theta) * observation.y;
      obs_y_temp = particle.y + sin(particle.theta) * observation.x +
                   cos(particle.theta * observation.y);
      observations_in_map_frame.push_back(
          LandmarkObs{observation.id, obs_x_temp, obs_y_temp});
    }

    // find landmarks in range of the particle
    vector<LandmarkObs> landmarks_in_range;
    double dist_landmark_particle;
    for (const auto &landmark : map_landmarks.landmark_list) {
      dist_landmark_particle = std::sqrt(pow(landmark.x_f - particle.x, 2) +
                                         pow(landmark.y_f - particle.y, 2));
      if (dist_landmark_particle <= sensor_range) {
        landmarks_in_range.push_back(
            LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
      }
    }

    // TODO data association mapped observations and landmarks
    dataAssociation(landmarks_in_range, observations_in_map_frame);

    // TODO find closest landmark to observation
    landmark_x_temp = 0.0;
    landmark_y_temp = 0.0;

    // reset weight
    particle.weight = 1.0;

    // calculate weight
    for (const auto &observation : observations_in_map_frame) {
      exponent = -(pow(observation.x - landmark_x_temp, 2) /
                       (2 * M_PI * std_landmark[0] * std_landmark[1]) +
                   pow(observation.y - landmark_y_temp, 2) /
                       (2 * M_PI * std_landmark[0] * std_landmark[1]));
      update_prob_temp =
          exp(exponent) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      update_prob_temp = std::max(update_prob_temp, 0.01);
      particle.weight = particle.weight * update_prob_temp;
    }
  }
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
  // associations: The landmark id that goes along with each listed
  // association sense_x: the associations x mapping already converted to
  // world coordinates sense_y: the associations y mapping already converted
  // to world coordinates
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
