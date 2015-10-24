#include "pcnn_simple.h"	 
#include <unordered_set>
#include <iostream>
using namespace std;

pcnn::pcnn(const unsigned int size, 
		   const conn_type connection_type, 
		   const pcnn_parameters & parameters,
		   const unsigned int width_oscillators,
		   const unsigned int height_oscillators) :
	m_oscillators(size, pcnn_oscillator()), 
	network(size, connection_type, width_oscillators, height_oscillators){

	cout<<"pcnn_simple constructor"<<endl;
	m_params = parameters;
}

pcnn::~pcnn() { }

void pcnn::simulate(const unsigned int steps, const pcnn_stimulus & stimulus, pcnn_dynamic & output_dynamic) {
	output_dynamic.resize(steps, size());

	for (unsigned int i = 0; i < steps; i++) {
		std::cout<<"currently simulating at step no "<< i << std::endl;;
		calculate_states(stimulus);
		store_dynamic(i, output_dynamic);
	}
}

void pcnn::calculate_states(const pcnn_stimulus & stimulus) {
	std::vector<double> feeding(size(), 0.0);
	std::vector<double> linking(size(), 0.0);
	std::vector<double> outputs(size(), 0.0);

	for (unsigned int index = 0; index < size(); index++) {
		pcnn_oscillator & current_oscillator = m_oscillators[index];
		std::vector<unsigned int> neighbors;
		get_neighbors(index, neighbors);

		double feeding_influence = stimulus[index];
		double linking_influence = 0.0;

		for (std::vector<unsigned int>::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++) {
			const double output_neighbor = m_oscillators[(*iter)].output;

			linking_influence += output_neighbor * m_params.W[index];
		}

		// heavily simplified!
		feeding[index] = feeding_influence;
		linking[index] = linking_influence;

		/* calculate internal activity */
		double internal_activity = feeding[index] + (m_params.B * linking[index]);

		/* calculate output of the oscillator */
		if (internal_activity > current_oscillator.threshold) {
			outputs[index] = OUTPUT_ACTIVE_STATE;
		}
		else {
			outputs[index] = OUTPUT_INACTIVE_STATE;
		}
	}

	/* update states of oscillators */
	for (unsigned int index = 0; index < size(); index++) {
		pcnn_oscillator & oscillator = m_oscillators[index];

		oscillator.feeding = feeding[index];
		oscillator.linking = linking[index];
		oscillator.output = outputs[index];
		if(oscillator.output = OUTPUT_INACTIVE_STATE){
			oscillator.threshold -= m_params.step_value; // decrease threshold by step value
		}
		else{
			oscillator.threshold *= m_params.AT; // penalize if it has pulsed
		}
	}
	// std::cout<<"calculate_states ";
}

void pcnn::store_dynamic(const unsigned int step, pcnn_dynamic & dynamic) {
	pcnn_network_state & current_state = (pcnn_network_state &) dynamic[step];
	current_state.m_output.resize(size());

	current_state.m_time = step;
	for (size_t i = 0; i < m_oscillators.size(); i++) {
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