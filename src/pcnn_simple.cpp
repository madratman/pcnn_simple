#include "pcnn_simple.h"	 
#include <unordered_set>
#include <iostream>
#include <algorithm>

using namespace std;

pcnn::pcnn(const unsigned int size, const conn_type connection_type, const pcnn_parameters & parameters) :
	m_oscillators(size, pcnn_oscillator()), 
	network(size, connection_type)
{
	m_params = parameters;
}

pcnn::pcnn(const unsigned int size, const conn_type connection_type, const size_t height, const size_t width, const pcnn_parameters & parameters) :
m_oscillators(size, pcnn_oscillator()),
network(size, connection_type, height, width)
{
    m_params = parameters;
}


pcnn::~pcnn() { }

void pcnn::simulate(const unsigned int steps, const pcnn_stimulus & stimulus, pcnn_dynamic & output_dynamic) {
	output_dynamic.resize(steps, size());

	for (unsigned int current_time_step = 0; current_time_step < steps; current_time_step++) {
		std::cout<<"currently simulating at step no "<< current_time_step << std::endl;;
		calculate_states(stimulus, current_time_step);
		// for (unsigned int index = 0; index < size(); index++) 
		// 	cout << "step " << current_time_step <<", test m_osc in pcnn::simulate " << m_oscillators[index].output<<endl;
		store_dynamic(current_time_step, output_dynamic);
	}
}

void pcnn::calculate_states(const pcnn_stimulus & stimulus, const unsigned int current_step) {
	// TODO initial_threshold should be a pcnn_parameter

	std::vector<double> feeding(size(), 0.0);
	std::vector<double> linking(size(), 0.0);
	std::vector<double> outputs(size(), 0.0);
	auto stimulus_max_ptr = std::max_element(stimulus.begin(), stimulus.end());
	double initial_threshold = *stimulus_max_ptr + 0.01;
	std::vector<double> internal_activity_vector;

	for (unsigned int index = 0; index < size(); index++) {
		pcnn_oscillator & current_oscillator = m_oscillators[index];
		std::vector<size_t> neighbors;
		get_neighbors(index, neighbors);

		double feeding_influence = stimulus[index];
		double linking_influence = 0.0;

		for (std::vector<size_t>::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++) {
			// output_neighbour is the value of output of the current neighbouring oscillator
			const double output_neighbor = m_oscillators[(*iter)].output;
			// std::cout<< index << ", " << *iter << endl; 

			int current = (*iter);
			const int upper_index = index - m_width;
	        const int upper_left_index = upper_index - 1;
       		const int upper_right_index = upper_index + 1;

			const int lower_index = index + m_width;
			const int lower_left_index = lower_index - 1;
			const int lower_right_index = lower_index + 1;
 
			const int left_index = index - 1;
			const int right_index = index + 1;

			const int node_row_index = std::floor(index / m_width);
        	const int upper_row_index = node_row_index - 1;
       		const int lower_row_index = node_row_index + 1;
       		// cout<< "index " << index << ", neighbor "<<  current  << ", stimulus " << stimulus[index] <<endl;

     	    if ((current == upper_left_index) && (upper_left_index >= 0) && (std::floor(upper_left_index / m_width) == upper_row_index)){
					linking_influence += output_neighbor * m_params.W[0];
					// if( output_neighbor == 1)
					// 	cout << linking_influence<<endl;
					// cout << "index " << index <<". upper_left_index. " << output_neighbor<<endl;
					// cout << "index " << index << "  " << " step " << current_step<< endl; 
				}				
			
			if ((current == upper_index) && upper_index >= 0){
					linking_influence += output_neighbor * m_params.W[1];
				}
       		
       		if ((current == upper_right_index) && (upper_right_index >= 0) && (std::floor(upper_right_index / m_width) == upper_row_index)){
					linking_influence += output_neighbor * m_params.W[2];
				}

	        if ((current == left_index) && (left_index >= 0) && (std::floor(left_index / m_width) == node_row_index)){
					linking_influence += output_neighbor * m_params.W[3];
				}
 			
 			if (index == current){
					linking_influence += output_neighbor * m_params.W[4];
				}
	        
	        if ((current == right_index) && (right_index < size()) && (std::floor(right_index / m_width) == node_row_index)){
					linking_influence += output_neighbor * m_params.W[5];
				}

     	  
     	    if ((current == lower_left_index) && (lower_right_index < size()) && (std::floor(lower_left_index / m_width) == upper_row_index)){
					linking_influence += output_neighbor * m_params.W[6];
					}				
			
			if ((current == lower_index) && (lower_index < size())){
					linking_influence += output_neighbor * m_params.W[7];
				}
       		
       		if ((current == lower_right_index) && (lower_right_index < size()) && (std::floor(lower_right_index / m_width) == upper_row_index)){
					linking_influence += output_neighbor * m_params.W[8];
				}
			
		}

   		// cout<< " linking_influence "<< linking_influence << ". feeding_influence "<< feeding_influence <<endl;
		
		// heavily simplified!
		feeding[index] = feeding_influence;
		linking[index] = linking_influence;

		/* calculate internal activity */
		double internal_activity = feeding[index] + (m_params.B * linking[index]);
		internal_activity_vector.push_back(internal_activity);

		if(current_step == 0)
		{
			current_oscillator.threshold = initial_threshold;
		}

		/* calculate output of the oscillator */
		if (internal_activity > current_oscillator.threshold) {
			outputs[index] = OUTPUT_ACTIVE_STATE;
		}
		else {
			outputs[index] = OUTPUT_INACTIVE_STATE;
		}
		// cout << "internal_activity " << internal_activity << ", current_oscillator.threshold "<<current_oscillator.threshold<< ", output "<< outputs[index]<<endl;
	}

	/* find minimum internal energy so as to set the threshold for next time step */
	auto internal_activity_min_ptr = std::min(internal_activity_vector.begin(), internal_activity_vector.end());
	auto internal_activity_max_ptr = std::max(internal_activity_vector.begin(), internal_activity_vector.end());
	double internal_activity_min = *internal_activity_min_ptr;
	double internal_activity_max = *internal_activity_max_ptr;

	/* update states of oscillators */
	for (unsigned int index = 0; index < size(); index++) {
		pcnn_oscillator & oscillator = m_oscillators[index];

		oscillator.feeding = feeding[index];
		oscillator.linking = linking[index];
		oscillator.output = outputs[index];
		// cout << "test outputs in calculate_states " << outputs[index]<<endl;
		// cout << "test m_osc in calculate_states " << oscillator.output<<endl;
		if(oscillator.output == OUTPUT_INACTIVE_STATE){
			oscillator.threshold -= m_params.step_value; // decrease threshold by step value
		}
		else{
			oscillator.threshold *= m_params.VT; // penalize if it has pulsed
			if(oscillator.threshold > 255)
				oscillator.threshold = 255;
			// oscillator.threshold = internal_activity_max; // penalize if it has pulsed
		}
	}
	// std::cout<<"calculate_states ";
}

void pcnn::store_dynamic(const unsigned int step, pcnn_dynamic & dynamic) {
	pcnn_network_state & current_state = (pcnn_network_state &) dynamic[step];
	current_state.m_output.resize(size());

	current_state.m_time = step;
	for (size_t i = 0; i < m_oscillators.size(); i++) {
		// cout <<  "inside store_dynamic, output "<< m_oscillators[i].output<<endl;
		current_state.m_output[i] = m_oscillators[i].output;
	}
	// std::cout<<"store_dynamic"<<std::endl;
}

pcnn_dynamic::pcnn_dynamic() {}

pcnn_dynamic::~pcnn_dynamic() {}

/* TODO: implementation */
pcnn_dynamic::pcnn_dynamic(const unsigned int number_oscillators, const unsigned int simulation_steps) { }

void pcnn_dynamic::allocate_sync_ensembles(ensemble_data<pcnn_ensemble> & ensembles) const {
	std::unordered_set<unsigned int> traverse_oscillators;
	traverse_oscillators.reserve(number_oscillators());

	for (const_reverse_iterator iter_state = crbegin(); iter_state != crend(); iter_state++) {
		pcnn_ensemble ensemble;
		const pcnn_network_state & state_network = (*iter_state);

		for (unsigned int i = 0; i < number_oscillators(); i++) {
			if (state_network.m_output[i] == OUTPUT_ACTIVE_STATE) {
				if (traverse_oscillators.find(i) == traverse_oscillators.end()) {
					ensemble.push_back(i);
					traverse_oscillators.insert(i);
				}
			}
		}

		if (!ensemble.empty()) {
			ensembles.push_back(ensemble);
		}
	}
}

void pcnn_dynamic::allocate_spike_ensembles(ensemble_data<pcnn_ensemble> & ensembles) const {
	for (const_iterator iter_state = cbegin(); iter_state != cend(); iter_state++) {
		pcnn_ensemble ensemble;
		const pcnn_network_state & state_network = (*iter_state);

		for (unsigned int i = 0; i < number_oscillators(); i++) {
			if (state_network.m_output[i] == OUTPUT_ACTIVE_STATE) {
				ensemble.push_back(i);
			}
		}

		if (!ensemble.empty()) {
			ensembles.push_back(ensemble);
		}
	}
}

void pcnn_dynamic::allocate_time_signal(pcnn_time_signal & time_signal) const {
	time_signal.resize(size());

	for (size_t t = 0; t < size(); t++) {
		const pcnn_network_state & state_network = (*this)[t];

		for (unsigned int i = 0; i < number_oscillators(); i++) {
			if (state_network.m_output[i] == OUTPUT_ACTIVE_STATE) {
				time_signal[t]++;
			}
		}
	}
}