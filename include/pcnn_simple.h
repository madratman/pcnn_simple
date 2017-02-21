#ifndef _PCNN_H_
#define _PCNN_H_

#include "network.h"

#include <vector>

#define OUTPUT_ACTIVE_STATE			(double) 1.0
#define OUTPUT_INACTIVE_STATE		(double) 0.0

typedef struct pcnn_oscillator {
	double output;
	double feeding;
	double linking;
	double threshold;

	pcnn_oscillator(void) :
		output(0.0),
		feeding(0.0),
		linking(0.0),
		threshold(0.0) { }
} pcnn_oscillator;

// const int PCNN_NO_OF_STEPS = 10;
const int PCNN_NO_OF_STEPS = 5;

typedef struct pcnn_parameters {
    
    // Multiplier for the threshold at the current step.    
    double VT = 2;
    
    // Synaptic weight - neighbours influence on linking compartment
    std::vector<double> W = {sqrt(1/2), 1.0, sqrt(1/2),
							 1.0, 1.0, 1.0, 
						     sqrt(1/2), 1.0, sqrt(1/2)};

    // Linking strength in the network.
    double B = 0.2;

    // step_value to be subtracted from oscillator
    double step_value = 10; 
} pcnn_parameters;

typedef std::vector<unsigned int>		pcnn_ensemble;
typedef std::vector<double>				pcnn_stimulus;
typedef std::vector<unsigned int>		pcnn_time_signal;

typedef struct pcnn_network_state {
public:
	std::vector<double> m_output;

	double				m_time;

public:
	inline size_t size(void) const { return m_output.size(); }

	inline pcnn_network_state & operator=(const pcnn_network_state & other) {
		if (this != &other) {
			m_output.resize(other.size());
			std::copy(other.m_output.cbegin(), other.m_output.cend(), m_output.begin());

			m_time = other.m_time;
		}

		return *this;
	}
} pcnn_network_state;


class pcnn_dynamic : public dynamic_data<pcnn_network_state> {
public:
	pcnn_dynamic(void);

	pcnn_dynamic(const unsigned int number_oscillators, const unsigned int simulation_steps);

	~pcnn_dynamic(void);

public:
	void allocate_sync_ensembles(ensemble_data<pcnn_ensemble> & sync_ensembles) const;

	void allocate_spike_ensembles(ensemble_data<pcnn_ensemble> & spike_ensembles) const;

	void allocate_time_signal(pcnn_time_signal & time_signal) const;

public:
    /***********************************************************************************************
    *
    * @brief   Returns dynamic state of oscillator of the pulse-coupled neural network at the
    *          specified iteration step.
    *
    * @param[in] iteration (size_t): number of iteration at which oscillator state is required.
    * @param[in] index_oscillator (size_t): index of oscillator whose state is required.
    *
    * @return (double) dynamic state of the oscillator at the specified iteration.
    *
    ***********************************************************************************************/
    inline double dynamic_oscillator_at(const size_t iteration, const size_t index_oscillator) const { 
        return dynamic_at(iteration).m_output[index_oscillator]; 
    }
};


class pcnn : public network {
protected:
	std::vector<pcnn_oscillator> m_oscillators;

	pcnn_parameters m_params;

private:
	pcnn(void);

public:
	pcnn(const unsigned int num_osc, const conn_type connection_type, const pcnn_parameters & parameters);

    pcnn(const unsigned int num_osc, 
         const conn_type connection_type,
         const size_t height,
         const size_t width,
         const pcnn_parameters & parameters);

	virtual ~pcnn(void);

public:
	void simulate(const unsigned int steps, const pcnn_stimulus & stimulus, pcnn_dynamic & output_dynamic);


private:
	void calculate_states(const pcnn_stimulus & stimulus, const unsigned int current_step);

	void store_dynamic(const unsigned int step, pcnn_dynamic & dynamic);

	void fast_linking(const std::vector<double> & feeding, std::vector<double> & linking, std::vector<double> & output);
};

#endif
