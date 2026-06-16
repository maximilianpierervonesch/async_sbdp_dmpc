
#include "simulator.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>        // for std::sin and std::cos
#include <fstream>      // for file output

Simulator::Simulator( int mpc_iter, double init_state[3])
        : mpc_iter_(mpc_iter), state_{init_state[0], init_state[1], init_state[2]},
	lcm_("udpm://239.255.76.67:7667?ttl=1") // define adress for multicast
          {
            if (!lcm_.good()) {
                throw std::runtime_error("LCM initialization failed");
            }
            // subscribe to plant messages 
            lcm_.subscribe("CONTROL_CHANNEL", &Simulator::handleControlMessage, this);

            // initialize state trajectory storage
            state_trajectory_.resize(mpc_iter_ + 1);
            for (int i = 0; i < mpc_iter_ + 1; ++i) {
                state_trajectory_[i].resize(3);
            }
            control_trajectory_.resize(mpc_iter_);
            for (int i = 0; i < mpc_iter_; ++i) {
                control_trajectory_[i].resize(3);
            }
        }

        void Simulator::handleControlMessage(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const messages::control* msg)
    {
        // update current control input
        control_[msg->id] = msg->control;  // assuming single control input per Agent
        std::cout << "Simulator received control message from Agent " << msg->id << " with control " << msg->control << std::endl;
        num_controls_received_++;
    }

    void Simulator::run() {
        // Save state and control trajectories
        state_trajectory_[0] = {state_[0], state_[1], state_[2]};
        for (int step = 0; step < mpc_iter_; ++step)
        {
            // 1. Publish updated state
            messages::plantState msg;
            msg.timestamp = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
            for (int i = 0; i < 3; ++i) {
                msg.states[i] = state_[i];
            }
            lcm_.publish("PLANT_CHANNEL", &msg);
            std::cout << "Simulator published plant state for step " << step << std::endl;
            
            received_control_ = false; // reset control received flag for this iteration
            num_controls_received_ = 0; // reset control count for this iteration
            while (!received_control_) {
                std::cout << "Waiting for control inputs for step " << step << std::endl;
                lcm_.handle(); 
                if(num_controls_received_ == num_agents_) {
                 received_control_ = true;  
            }     
            }; 

            // 2. Update state based on received control inputs and centralized dynamics
            std::cout << "Received controls for step " << step << ": [" << control_[0] << ", " << control_[1] << ", " << control_[2] << "]" << std::endl;
            std::cout << "Current State at step " << step << ": [" << state_[0] << ", " << state_[1] << ", " << state_[2] << "]" << std::endl;

            double new_state[3];

            new_state[0] = state_[0] +  dt_*(control_[0] + std::sin(state_[1]) + std::sin(state_[2])); 
            new_state[1] = state_[1] +  dt_*(control_[1] + std::sin(state_[0]) + std::sin(state_[2])); 
            new_state[2] = state_[2] +  dt_*(control_[2] + std::sin(state_[0]) + std::sin(state_[1])); 

            state_[0] = new_state[0];
            state_[1] = new_state[1];
            state_[2] = new_state[2];

            std::cout << "Simulator updated state for step " << step << ": [" << state_[0] << ", " << state_[1] << ", " << state_[2] << "]" << std::endl;

            // Reset for next iteration
            received_control_ = false;

            // Save state and control trajectories
            state_trajectory_[step + 1] = {state_[0], state_[1], state_[2]};
            control_trajectory_[step] = {control_[0], control_[1], control_[2]};
        }   
        std::cout << "Simulator completed " << mpc_iter_ << " iterations.\n";   
        
        // write state and control trajectories to files
        std::ofstream state_file("state_trajectory.csv");
        std::ofstream control_file("control_trajectory.csv"); 
        state_file << "Step,Agent1,Agent2,Agent3\n";
        control_file << "Step,Agent1,Agent2,Agent3\n";
        for (int i = 0; i < mpc_iter_; ++i) {
            state_file << i << "," << state_trajectory_[i][0] << "," << state_trajectory_[i][1] << "," << state_trajectory_[i][2] << "\n";
            control_file << i << "," << control_trajectory_[i][0] << "," << control_trajectory_[i][1] << "," << control_trajectory_[i][2] << "\n";
        }
        state_file.close();
        control_file.close();
    } 
        
    int main(int argc, char** argv)
    {
   
    int mpc_iter = std::stoi(argv[1]);
    double init_state[3] = {std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4])};

    std::cout << "Starting simulator for " << mpc_iter << " iterations\n";

    Simulator simulator(mpc_iter, init_state);
    simulator.run();

    return 0;
    }
