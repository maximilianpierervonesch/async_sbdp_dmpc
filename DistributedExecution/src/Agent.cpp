#include "Agent.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>        // for std::sin and std::cos
#include <fstream>      // for file output


Agent::Agent(int id, int mpc_iter, int alg_iter, int num_agents)
        : id_(id),
          mpc_iter_(mpc_iter),
          alg_iter_(alg_iter),
          num_agents_(num_agents),
	 lcm_("udpm://239.255.76.67:7667?ttl=1") // define adress for multicast
          {
        std::cout << "lcm_.good() = " << lcm_.good() << std::endl;

    if (!lcm_.good()) {
        throw std::runtime_error("LCM instance not initialized");
    }

            // Initialize arrays
            for (int i = 0; i < num_discretization_points_; ++i) predicted_state_[i] = 0.0;
            for (int i = 0; i < num_discretization_points_; ++i) predicted_lambda_[i] = 0.0;
            for (int i = 0; i < num_discretization_points_; ++i) predicted_control_[i] = 0.0;

            // Subscribe to trajectory messages
            lcm_.subscribe("TRAJECTORY_CHANNEL", &Agent::handleTrajectoryMessage, this);
            // subscribe to plant messages 
            lcm_.subscribe("PLANT_CHANNEL", &Agent::handlePlantMessage, this);

            // resize neighbor trajectories storage
            neighbor_trajectories_.resize(num_agents_); 
        }

    void Agent::setParameters(const double* params)
        {
            for (size_t i = 0; i < parameters_.size(); ++i) {
                parameters_[i] = params[i];
            }
        }

    void Agent::run_async() 
    {
        for (int step = 0; step < mpc_iter_; ++step)
        {
        current_mpc_iter_ = step;
        // ------------------------
        // 1. Receive latest plant state
        // ------------------------
        bool state_received = false;
        while (!state_received)
        {
            //std::cout << "Agent " << id_ << " waiting for plant state for step " << step << std::endl;

            // Non-blocking poll
            lcm_.handle();
            // Check if a new PlantState message has arrived
            if (latest_plant_state_received_)
            {   
                // save latest plant state
                latest_plant_state_received_ = false;
                state_received = true;
                //std::cout << "Agent " << id_ << " received latest plant state for step " << step << std::endl;
            }
        }
       auto t_start = std::chrono::steady_clock::now();
        // ------------------------
        // 2. Run fixed SBDP iterations
        // ------------------------
        // send initial trajectory
        publish_trajectory(0, current_mpc_iter_);      // send to neighbors
        for (int k = 0; k < alg_iter_; ++k)
        {               
            // Drain all pending messages
            while (lcm_.handleTimeout(0) == 0) {};
            compute_local_update();                          // local trajectory update
            publish_trajectory(k+1, current_mpc_iter_);      // send to neighbors
           // Drain all pending messages
            while (lcm_.handleTimeout(0) == 0) {};                
        }
        auto t_end = std::chrono::steady_clock::now();
        auto latency_us = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
        // std::cout << "Agent " << id_ << " completed " << alg_iter_ << " Algorithm iterations in "
        //     << latency_us << " us." << std::endl; 
        // write computation times to file 
            std::ofstream latency_file("latency"+std::to_string(id_)+".csv", std::ios::app);
            latency_file << latency_us << "\n";
            latency_file.close();    
        // ------------------------
        // 3. Send computed control to simulator
        // ------------------------
        publish_control_action();
    }
    }

     void Agent::run_sync() 
    {
        for (int step = 0; step < mpc_iter_; ++step)
        {
        current_mpc_iter_ = step;
        // ------------------------
        // 1. Receive latest plant state
        // ------------------------
        bool state_received = false;
        while (!state_received)
        {
            //std::cout << "Agent " << id_ << " waiting for plant state for step " << step << std::endl;

            // Non-blocking poll
            lcm_.handle();
            // Check if a new PlantState message has arrived
            if (latest_plant_state_received_)
            {   
                // save latest plant state
                latest_plant_state_received_ = false;
                state_received = true;
                //std::cout << "Agent " << id_ << " received latest plant state for step " << step << std::endl;
            }
        }
       auto t_start = std::chrono::steady_clock::now();
        // ------------------------
        // 2. Run fixed SBDP iterations
        // ------------------------
        // send initial trajectory
        publish_trajectory(0, current_mpc_iter_);      // send to neighbors
        for (int k = 0; k < alg_iter_; ++k)
        {               
            // wait for trajectories from all neighbors for current iteration
            bool all_neighbors_updated = false;
            while (!all_neighbors_updated) {
                lcm_.handle();
                if (recieved_all_trajectories_) {
                    recieved_all_trajectories_ = false; // reset for next iteration
                    all_neighbors_updated = true;
                }
            }
            compute_local_update();                          // local trajectory update
            publish_trajectory(k+1, current_mpc_iter_);      // send to neighbors
        }
        auto t_end = std::chrono::steady_clock::now();
        auto latency_us = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
        // std::cout << "Agent " << id_ << " completed " << alg_iter_ << " Algorithm iterations in "
        //     << latency_us << " us." << std::endl; 
        // write computation times to file 
            std::ofstream latency_file("latency"+std::to_string(id_)+".csv", std::ios::app);
            latency_file << latency_us << "\n";
            latency_file.close();
        // ------------------------
        // 3. Send computed control to simulator
        // ------------------------
        publish_control_action();
        }
    }

    void Agent::compute_local_update()
    {
        // make copy of neighbor trajectories since they may be updated asynchronously during computation
        double neighbor_states[num_agents_][num_discretization_points_];
        double neighbor_adjoints[num_agents_][num_discretization_points_];
        {
            std::lock_guard<std::mutex> lock(neighbor_mutex_);
            for (int j = 0; j < num_agents_; ++j) {
                if (j == id_) continue;
                for (int i = 0; i < num_discretization_points_; ++i) {
                    neighbor_states[j][i]   = neighbor_trajectories_[j].states[i];
                    neighbor_adjoints[j][i] = neighbor_trajectories_[j].adjoint[i];
                }
            }
        }
        // auxiliary variables
         double x_mid;
         double lambda_mid;
         double F_current;
         double F_mid;
         double G_current;
         double G_mid;
         double neighbor_state[num_agents_-1]; // excluding self
         double neighbor_lambda[num_agents_-1]; // excluding self
         double alpha = 0.05;                   // damping for forward-backward sweep
        // inner FB loop
        for ( int q = 0; q < 50; ++q)
            {
            // x integration forward in time    
            for ( int i = 0; i < num_discretization_points_ - 1; ++i)
                {
                predicted_control_[i] = alpha * h_i(predicted_state_[i], predicted_lambda_[i]) + (1 - alpha) * predicted_control_[i]; // control update with relaxation

                    for ( int j = 0; j < num_agents_; ++j)
                        {
                        if(j != id_) 
                            {
                             neighbor_state[j] = neighbor_states[j][i];
                            }
                        }
                F_current = F_i(predicted_state_[i], predicted_lambda_[i], predicted_control_[i], neighbor_state);
                x_mid = predicted_state_[i] + Delta_t_ * F_current;
                F_mid = F_i(x_mid, predicted_lambda_[i], predicted_control_[i], neighbor_state);
                predicted_state_[i + 1] = predicted_state_[i] + 0.5 * Delta_t_ * (F_current + F_mid);
                }

            // lambda integration backward in time
            // set terminal condition
            predicted_lambda_[num_discretization_points_ - 1] = parameters_[7] * predicted_state_[num_discretization_points_ - 1];
            for ( int i = num_discretization_points_ - 1; i > 0; --i)
                {
                     for ( int j = 0; j < num_agents_; ++j)
                        {
                        if(j != id_) 
                            {
                             neighbor_state[j] = neighbor_states[j][i];
                             neighbor_lambda[j] = neighbor_adjoints[j][i];
                            }
                        }
                    G_current = G_i(predicted_state_[i], predicted_lambda_[i], neighbor_state, neighbor_lambda);
                    lambda_mid = predicted_lambda_[i] - Delta_t_ * G_current;
                    G_mid = G_i(predicted_state_[i], lambda_mid, neighbor_state, neighbor_lambda);
                    predicted_lambda_[i - 1] = predicted_lambda_[i] - 0.5 * Delta_t_ * (G_current + G_mid);
                } 
            }
    }

    double Agent::F_i(double state, double lambda, const double control, const double* neighbor_state)
    {
         // Placeholder for actual dynamics
        double ai = parameters_[0]; // example parameter
        double ci = parameters_[2]; // example parameter
        double out = -ai * state + control; // local control term
        for ( int j = 0; j < num_agents_; ++j)
            {
                if(j != id_) 
                    {
                        out += ci * std::sin(neighbor_state[j]); // coupling term
                    }
            }   
        return out; 
    }

    double Agent::G_i(double state, double lambda, const double* neighbor_state, const double* neighbor_lambda)
    {
       double ai = parameters_[0]; // example parameter
       double ci = parameters_[2]; // example parameter
       double qi = parameters_[6]; // state cost weight
       double out = -2 * qi * state - ai * lambda;
       for ( int j = 0; j < num_agents_; ++j)
            {
                if(j != id_) 
                    {
                     out += -ci * std::cos(neighbor_state[j]) * neighbor_lambda[j] ; // coupling term
                    }
            }
       return out;
    }

    double Agent::h_i(double state, double lambda)
    {   
        double bi = parameters_[1]; // parameter 
        double umin = parameters_[3]; // lower bound 
        double umax = parameters_[4]; // upper bound
        double ri = parameters_[5]; // control cost weight


        double u_opt = - (lambda * bi) / (2 * ri) ; // unconstrained optimal control
        return std::clamp(u_opt, umin, umax);
    }

    void Agent::handleTrajectoryMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const messages::trajectory* msg)
    {
        if (msg->id == id_) return;
        // Validate message ID is within bounds
        if (msg->id < 0 || msg->id >= num_agents_) {
            std::cerr << "Error: Received trajectory from Agent " << msg->id << " but number of agents is " << num_agents_ << std::endl;
            return;
        }
        if(msg->mpc_step < current_mpc_iter_) {
            return;
        }
        {
        // save neighbor trajectory at global ID index in
        std::lock_guard<std::mutex> lock(neighbor_mutex_);
        neighbor_trajectories_[msg->id] = *msg;
        // increment count of received trajectories for current iteration
        traj_counter_++;
        if (traj_counter_ == num_agents_ - 1) {
            recieved_all_trajectories_ = true;
            traj_counter_ = 0; // reset for next iteration
        }
        }
    }

    void Agent::handlePlantMessage(const lcm::ReceiveBuffer* rbuf, const std::string& chan, const messages::plantState* msg)
    {
        // update current state
        for (int i = 0; i < 1; ++i)
        {
            predicted_state_[0] = msg->states[id_ + i];       // assuming states are packed per agent and id starts from 0
        }
        latest_plant_state_received_ = true;
    }

    void Agent::publish_trajectory(int current_alg_iter, int current_mpc_iter) 
    {
        messages::trajectory traj_msg;
        traj_msg.id = id_;
        traj_msg.step = current_alg_iter;
        traj_msg.mpc_step = current_mpc_iter; ;
        // Copy predicted_state_ to traj_msg.states
        for (int i = 0; i < num_discretization_points_; ++i) 
        {
            traj_msg.states[i] = predicted_state_[i];
            traj_msg.adjoint[i] = predicted_lambda_[i];
        }
        lcm_.publish("TRAJECTORY_CHANNEL", &traj_msg);
    }

    void Agent::publish_control_action() 
    {
        messages::control ctrl_msg;
        ctrl_msg.id = id_;
        // Classical MPC: use first control input
        ctrl_msg.control = predicted_control_[0]; // 
        //std::cout << "Agent " << id_ << " published control: " << ctrl_msg.control << std::endl;
        ctrl_msg.timestamp = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        lcm_.publish("CONTROL_CHANNEL", &ctrl_msg);
    }

    int main(int argc, char** argv)
    {
        if (argc < 5)
        {
            std::cerr << "Error: Usage is agent <id> <mpc_iterations> <alg_iterations> <num_agents>\n";
            return 1;
        }

        int id = std::stoi(argv[1]);
        int mpc_iter = std::stoi(argv[2]);
        int alg_iter = std::stoi(argv[3]);
        int num_agents = std::stoi(argv[4]);


        std::cout << "Starting agent " << id
                << " for " << mpc_iter << " MPC iterations\n"
                << " with " << alg_iter << " algorithm iterations each.\n";

        Agent agent(id, mpc_iter, alg_iter, num_agents);
        double params[8] = {0.0, 1.0, 1.0, -2.0, 2.0, 0.1, 1.0, 0};  // ai bi ci umin umax ri qi pi
        agent.setParameters(params); // example parameters
        agent.run_async();  // toggle between run_async() and run_sync()

        return 0;
    }
