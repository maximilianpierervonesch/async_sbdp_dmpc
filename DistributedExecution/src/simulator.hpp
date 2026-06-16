#pragma once
#include <lcm/lcm-cpp.hpp>
#include "messages/control.hpp"
#include "messages/plantState.hpp"


class Simulator {
    public:
        Simulator(int mpc_iter, double init_states[3]); // constructor with number of agents
        // start listening for messages
        void run();

        private:   
        void handleControlMessage(const lcm::ReceiveBuffer* rbuf,
                                const std::string& chan,
                                const messages::control* msg);
        
        lcm::LCM lcm_;
        int mpc_iter_;                          // Number of MPC steps
        double state_[3];                       // Current state of the system
        double control_[3];                     // Current control inputs for all agents
        double dt_ = 0.1;                       // Sampling time 
        int num_controls_received_ = 0;         // Number of control inputs recieved in current iteration
        int num_agents_ = 3;                    // Number of agents in the system
        bool received_control_ = false;         // Flag to indicate if control input is received
        std::vector<std::vector<double>> state_trajectory_;                // closed-loop state trajectory of the system
        std::vector<std::vector<double>> control_trajectory_;                // closed-loop control trajectory of the system
        
};
