#pragma once
#include <lcm/lcm-cpp.hpp>
#include "messages/trajectory.hpp"
#include "messages/control.hpp"
#include "messages/plantState.hpp"

#include <atomic>
#include <mutex>
#include <vector>

#include <array>       // for std::array
#include <algorithm>   // for std::clamp

class Agent {
    public:
        Agent(int id, int mpc_iter, int alg_iter, int num_neighbors); // constructor with agent ID and algorithm iterations; Agent ID starts from 0

        // start main loop
        void run_async();
        void run_sync();
        void setParameters(const double* params);

        private:
        // Communication handlers and publishers
        void publishLoop();
        void handleTrajectoryMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const messages::trajectory* msg);
        void handlePlantMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const messages::plantState* msg);
        void publish_trajectory(int current_alg_iter, int current_mpc_iter);
        void publish_control_action();

        // Algorithmic functions
        void compute_local_update();
        double F_i(double state, double lambda, double control, const double* neighbor_state);
        double G_i(double state, double lambda, const double* neighbor_state, const double* neighbor_lambda);
        double h_i(double state, double lambda);

        lcm::LCM lcm_; 
        int id_;                                                    // ID of agent
        int alg_iter_;                                              // Maximum Algorithm iterations
        int mpc_iter_;                                              // Number of MPC steps
        int current_mpc_iter_;                                      // Current MPC iteration
        double predicted_state_[31];                                // Current predicted state of the agent
        double predicted_lambda_[31];                               // Current lambda values
        double predicted_control_[31];                              // Current predicted control actions
        int num_agents_;                                            // Number of neighboring agents
        std::vector<messages::trajectory> neighbor_trajectories_;   // Store neighbor trajectories
        std::mutex neighbor_mutex_;                                 // Mutex for neighbor trajectories
        bool latest_plant_state_received_ = false;                  // Flag to check if latest plant state is received
        int num_discretization_points_ = 31;                        // Number of discretization points in prediction horizon
        double Delta_t_ = 0.1;                                      // Time step size   
        std::array<double, 8> parameters_;                           // Parameters for dynamics and cost functions  
        // Synchronization variables for asynchronous version
        std::atomic<bool> recieved_all_trajectories_{false};        // Flag to indicate if trajectories from all neighbors have been received for current iteration
        std::atomic<int> traj_counter_{0};                          // Counter for number of trajectories received for current iteration
};  
