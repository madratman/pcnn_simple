// This main is just a demo to get familiar with the code base. 
// It's extract from the GTests of the ccore of pyclustering. (hence the static keyword)
// here: https://github.com/annoviko/pyclustering/blob/da8cdc231f615879bd6eb8c2b78a5f5cac2cd913/ccore/utcore/utest-pcnn.h

#include "pcnn_simple.h"
#include "network.h"
#include <unordered_set>
#include <iostream>

using namespace std;

vector< vector<double> > pcnn_dynamic_get_output(pcnn_dynamic dynamic, int num_osc) 
{
	vector< vector<double> > pcnn_data;
	pcnn_data.resize( dynamic.size() );

	for (int step = 0; step < dynamic.size(); step++)
	{
		pcnn_data[step].resize(num_osc);
		pcnn_network_state & current_state = dynamic[step];

		for (int i = 0; i < num_osc; i++) 
		{
			// pcnn_data[step][i] = current_state.m_output[i];
			pcnn_data[step][i] = dynamic.dynamic_oscillator_at(step, i);
			cout << "pcnn_dynamic_get_output(), " << "step " << step << ", osc_number "<<i <<", output " <<current_state.m_output[i]<<endl;
		}
	}

	return pcnn_data;
}

// generates, simulates PCNN. Returns a vector of outputs at a hardcoded timestep 
vector< vector<double> >  pcnn_ensemble_simulate(
		const unsigned int num_osc, 
		const unsigned int steps, 
		const conn_type type_conn, 
		const pcnn_stimulus & stimulus, 
		const pcnn_parameters * const params = nullptr) 
{
	// use default params if not specified
	pcnn_parameters parameters;
	if (params != nullptr) {
		parameters = *params;
	}

	// class pcnn has a vector of pcnn_oscillators, which are structs which have feeding, linking, output, threshold as members
	// class pcnn also has a member pcnn_params

	pcnn network(num_osc, type_conn, parameters);

	// pcnn_dynamic is a class that represents the dynamic output of pcnn. 
	// It's inherited from templated dynamic_data and uses "pcnn_network_state" (in place of the template)
	// dynamic_data::output_dynamics member stores a std::vector of pcnn_network_state 
	// pcnn_network_state is a struct with members std::vector<double> m_output and double m_time.

	pcnn_dynamic dynamic;

	// Performs static simulation of oscillatory network based on Hodgkin-Huxley neuron model.
	network.simulate(steps, stimulus, dynamic);

	// TODO : return a vector of vectors for all time steps
	vector< vector<double> >  pcnn_state; 

	// pcnn_state.resize( num_osc * dynamic.size() );
	pcnn_state = pcnn_dynamic_get_output(dynamic, num_osc);

	return pcnn_state;
}


static void template_dynamic_generation(
		const unsigned int num_osc, 
		const unsigned int steps, 
		const conn_type type_conn, 
		const pcnn_stimulus & stimulus) {

	pcnn_parameters parameters;
	pcnn network(num_osc, type_conn, parameters);

	pcnn_dynamic dynamic;
	network.simulate(steps, stimulus, dynamic);
	
	pcnn_time_signal time_signal;
	dynamic.allocate_time_signal(time_signal);
	pcnn_dynamic_get_output(dynamic, num_osc);
}

static void template_output_activity(
	const unsigned int num_osc,
	const unsigned int steps,
	const conn_type type_conn, 
	const pcnn_stimulus & stimulus,
	const bool activity_requirement,
	const pcnn_parameters * const params = nullptr) {

	pcnn_parameters parameters;
	if (params != nullptr) {
		parameters = *params;
	}

	pcnn network(num_osc, type_conn, parameters);

	pcnn_dynamic dynamic;
	network.simulate(steps, stimulus, dynamic);

	ensemble_data<pcnn_ensemble> sync_ensembles;
	ensemble_data<pcnn_ensemble> spike_ensembles;
	pcnn_time_signal time_signal;

	dynamic.allocate_sync_ensembles(sync_ensembles);
	dynamic.allocate_spike_ensembles(spike_ensembles);
	dynamic.allocate_time_signal(time_signal);

	/* check time signal for activity */
	bool output_activity = false;
	for (size_t i = 0; i < time_signal.size(); i++) {
		if (time_signal[i] > 0) {
			output_activity = true;
			break;
		}
	}

	/* if activity exists in time signal then at least one ensemble should be */}

static void template_ensemble_allocation(
	const unsigned int num_osc,
	const unsigned int steps,
	const conn_type type_conn, 
	const pcnn_stimulus & stimulus,
	const pcnn_parameters * const params = nullptr) {

	pcnn_parameters parameters;
	if (params != nullptr) {
		parameters = *params;
	}

	pcnn network(num_osc, type_conn, parameters);

	pcnn_dynamic dynamic;
	network.simulate(steps, stimulus, dynamic);

	ensemble_data<pcnn_ensemble> sync_ensembles;
	ensemble_data<pcnn_ensemble> spike_ensembles;
	pcnn_time_signal time_signal;

	dynamic.allocate_sync_ensembles(sync_ensembles);
	dynamic.allocate_spike_ensembles(spike_ensembles);
	dynamic.allocate_time_signal(time_signal);

	for (ensemble_data<pcnn_ensemble>::const_iterator iter = spike_ensembles.cbegin(); iter != spike_ensembles.cend(); iter++) {
		const pcnn_ensemble & ensemble = (*iter);
	}

	std::unordered_set<size_t> traverse_oscillators;

	for (ensemble_data<pcnn_ensemble>::const_iterator iter = sync_ensembles.cbegin(); iter != sync_ensembles.cend(); iter++) {
		const pcnn_ensemble & ensemble = (*iter);

		for (pcnn_ensemble::const_iterator iter_index = ensemble.cbegin(); iter_index != ensemble.cend(); iter_index++) {
			size_t index_oscillator = (*iter_index);
			traverse_oscillators.insert(index_oscillator);
		}
	}
}

int main()
{
	pcnn_stimulus stimulus { 10,1, 10,0}; 
	pcnn_ensemble_simulate(stimulus.size(), 10, conn_type::GRID_EIGHT, stimulus);
	// template_dynamic_generation(stimulus.size(), 20, conn_type::GRID_EIGHT, stimulus);
	// template_output_activity(stimulus.size(), 30, conn_type::GRID_EIGHT, stimulus, true);

	// pcnn_stimulus full_stimulus { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }; 
	// pcnn_stimulus partial_stimulus { 1, 0, 0, 1, 1, 1, 0, 0, 1, 1 };
	// pcnn_stimulus no_stimulus { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// pcnn_parameters params;
	// params.FAST_LINKING = true;

	// template_ensemble_allocation(stimulus.size(), 20, conn_type::ALL_TO_ALL, stimulus);
	// template_ensemble_allocation(full_stimulus.size(), 50, conn_type::ALL_TO_ALL, full_stimulus, &params);
}
