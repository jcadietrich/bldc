/*
	Copyright 2016 - 2022 Benjamin Vedder	benjamin@vedder.se

	This file is part of the VESC firmware.

	The VESC firmware is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The VESC firmware is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "foc_math.h"
#include "utils_math.h"
#include <math.h>

// See http://cas.ensmp.fr/~praly/Telechargement/Journaux/2010-IEEE_TPEL-Lee-Hong-Nam-Ortega-Praly-Astolfi.pdf
void foc_observer_update(float v_alpha, float v_beta, float i_alpha, float i_beta,
		float dt, observer_state *state, float *phase, motor_all_state_t *motor) {

	mc_configuration *conf_now = motor->m_conf;

	float R = conf_now->foc_motor_r;
	float L = conf_now->foc_motor_l;
	float lambda = conf_now->foc_motor_flux_linkage;

	// Saturation compensation
	switch(conf_now->foc_sat_comp_mode) {
	case SAT_COMP_LAMBDA:
		// Here we assume that the inductance drops by the same amount as the flux linkage. I have
		// no idea if this is a valid or even a reasonable assumption.
		if (conf_now->foc_observer_type >= FOC_OBSERVER_ORTEGA_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXLEMMING_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			L = L * (state->lambda_est / lambda);
		}
		break;

	case SAT_COMP_FACTOR: {
		const float comp_fact = conf_now->foc_sat_comp * (motor->m_motor_state.i_abs_filter / conf_now->l_current_max);
		L -= L * comp_fact;
		lambda -= lambda * comp_fact;
	} break;

	case SAT_COMP_LAMBDA_AND_FACTOR: {
		if (conf_now->foc_observer_type >= FOC_OBSERVER_ORTEGA_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXLEMMING_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			L = L * (state->lambda_est / lambda);
		}
		const float comp_fact = conf_now->foc_sat_comp * (motor->m_motor_state.i_abs_filter / conf_now->l_current_max);
		L -= L * comp_fact;
	} break;

	default:
		break;
	}

	// Temperature compensation
	if (conf_now->foc_temp_comp) {
		R = motor->m_res_temp_comp;
	}

	float ld_lq_diff = conf_now->foc_motor_ld_lq_diff;
	float id = motor->m_motor_state.id;
	float iq = motor->m_motor_state.iq;

	// Adjust inductance for saliency.
	if (fabsf(id) > 0.1 || fabsf(iq) > 0.1) {
		L = L - ld_lq_diff / 2.0 + ld_lq_diff * SQ(iq) / (SQ(id) + SQ(iq));
	}

	float L_ia = L * i_alpha;
	float L_ib = L * i_beta;
	const float R_ia = R * i_alpha;
	const float R_ib = R * i_beta;
	const float gamma_half = motor->m_gamma_now * 0.5;

	switch (conf_now->foc_observer_type) {
	case FOC_OBSERVER_ORTEGA_ORIGINAL: {
		float err = SQ(lambda) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));

		// Forcing this term to stay negative helps convergence according to
		//
		// http://cas.ensmp.fr/Publications/Publications/Papers/ObserverPermanentMagnet.pdf
		// and
		// https://arxiv.org/pdf/1905.00833.pdf
		if (err > 0.0) {
			err = 0.0;
		}

		float x1_dot = v_alpha - R_ia + gamma_half * (state->x1 - L_ia) * err;
		float x2_dot = v_beta - R_ib + gamma_half * (state->x2 - L_ib) * err;

		state->x1 += x1_dot * dt;
		state->x2 += x2_dot * dt;
	} break;

	case FOC_OBSERVER_MXLEMMING:
	case FOC_OBSERVER_MXLEMMING_LAMBDA_COMP:
		// LICENCE NOTE:
		// This function deviates slightly from the BSD 3 clause licence.
		// The work here is entirely original to the MESC FOC project, and not based
		// on any appnotes, or borrowed from another project. This work is free to
		// use, as granted in BSD 3 clause, with the exception that this note must
		// be included in where this code is implemented/modified to use your
		// variable names, structures containing variables or other minor
		// rearrangements in place of the original names I have chosen, and credit
		// to David Molony as the original author must be noted.

		state->x1 += (v_alpha - R_ia) * dt - L * (i_alpha - state->i_alpha_last);
		state->x2 += (v_beta - R_ib) * dt - L * (i_beta - state->i_beta_last);

		if (conf_now->foc_observer_type == FOC_OBSERVER_MXLEMMING_LAMBDA_COMP) {
			float err = SQ(state->lambda_est) - (SQ(state->x1) + SQ(state->x2));
			state->lambda_est += 0.1 * gamma_half * state->lambda_est * -err * dt;
			utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

			utils_truncate_number_abs(&(state->x1), state->lambda_est);
			utils_truncate_number_abs(&(state->x2), state->lambda_est);
		} else {
			utils_truncate_number_abs(&(state->x1), lambda);
			utils_truncate_number_abs(&(state->x2), lambda);
		}

		// Set these to 0 to allow using the same atan2-code as for Ortega
		L_ia = 0.0;
		L_ib = 0.0;
		break;

	case FOC_OBSERVER_ORTEGA_LAMBDA_COMP: {
		float err = SQ(state->lambda_est) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));

		// FLux linkage observer. See:
		// https://cas.mines-paristech.fr/~praly/Telechargement/Conferences/2017_IFAC_Bernard-Praly.pdf
		state->lambda_est += 0.2 * gamma_half * state->lambda_est * -err * dt;

		// Clamp the observed flux linkage (not sure if this is needed)
		utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

		if (err > 0.0) {
			err = 0.0;
		}

		float x1_dot = v_alpha - R_ia + gamma_half * (state->x1 - L_ia) * err;
		float x2_dot = v_beta - R_ib + gamma_half * (state->x2 - L_ib) * err;

		state->x1 += x1_dot * dt;
		state->x2 += x2_dot * dt;
	} break;

	case FOC_OBSERVER_MXV:
	case FOC_OBSERVER_MXV_LAMBDA_COMP:
	case FOC_OBSERVER_MXV_LAMBDA_COMP_LIN:
		state->x1 += (v_alpha - R_ia) * dt;
		state->x2 += (v_beta - R_ib) * dt;

		if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
				float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
				UTILS_LP_FAST(state->lambda_est, mag, 0.1 * gamma_half * dt * SQ(state->lambda_est));
				utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

				if (mag > state->lambda_est) {
					state->x1 = (state->x1 / mag) * state->lambda_est;
					state->x2 = (state->x2 / mag) * state->lambda_est;
				}
			} else if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP) {
				float err = SQ(state->lambda_est) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));
				state->lambda_est += 0.2 * gamma_half * state->lambda_est * -err * dt;
				utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

				float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
				if (mag > state->lambda_est) {
					state->x1 = (state->x1 / mag) * state->lambda_est;
					state->x2 = (state->x2 / mag) * state->lambda_est;
				}
			}
		} else {
			float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
			if (mag > lambda) {
				state->x1 = (state->x1 / mag) * lambda;
				state->x2 = (state->x2 / mag) * lambda;
			}
		}
		break;

	default:
		break;
	}

	state->i_alpha_last = i_alpha;
	state->i_beta_last = i_beta;

	UTILS_NAN_ZERO(state->x1);
	UTILS_NAN_ZERO(state->x2);

	// Prevent the magnitude from getting too low, as that makes the angle very unstable.
	float mag = NORM2_f(state->x1, state->x2);
	if (mag < (lambda * 0.5)) {
		state->x1 *= 1.1;
		state->x2 *= 1.1;
	}

	if (phase) {
		*phase = utils_fast_atan2(state->x2 - L_ib, state->x1 - L_ia);
	}

	// Can we clamp the flux in dq with q flux = 0 and d flux is lambda
	// Then the state->x1 and state->x2 (which are the alpha and beta fluxes) are set as lambda*sin and lambda*cos
	// The d flux each time would have a residual after transform from ab to dq. This can be used as an input to the flux estimator
}

void foc_pll_run(float phase, float dt, float *phase_var,
					float *speed_var, mc_configuration *conf) {
	UTILS_NAN_ZERO(*phase_var);
	float delta_theta = phase - *phase_var;
	utils_norm_angle_rad(&delta_theta);
	UTILS_NAN_ZERO(*speed_var);
	*phase_var += (*speed_var + conf->foc_pll_kp * delta_theta) * dt;
	utils_norm_angle_rad((float*)phase_var);
	*speed_var += conf->foc_pll_ki * delta_theta * dt;
}

/**
 * @brief svm Space vector modulation. Magnitude must not be larger than sqrt(3)/2, or 0.866 to avoid overmodulation.
 *        See https://github.com/vedderb/bldc/pull/372#issuecomment-962499623 for a full description.
 * @param alpha voltage
 * @param beta Park transformed and normalized voltage
 * @param PWMFullDutyCycle is the peak value of the PWM counter.
 * @param tAout PWM duty cycle phase A (0 = off all of the time, PWMFullDutyCycle = on all of the time)
 * @param tBout PWM duty cycle phase B
 * @param tCout PWM duty cycle phase C
 */
void foc_svm(float alpha, float beta, uint32_t PWMFullDutyCycle,
				uint32_t* tAout, uint32_t* tBout, uint32_t* tCout, uint32_t *svm_sector) {
	uint32_t sector;

	if (beta >= 0.0f) {
		if (alpha >= 0.0f) {
			//quadrant I
			if (ONE_BY_SQRT3 * beta > alpha) {
				sector = 2;
			} else {
				sector = 1;
			}
		} else {
			//quadrant II
			if (-ONE_BY_SQRT3 * beta > alpha) {
				sector = 3;
			} else {
				sector = 2;
			}
		}
	} else {
		if (alpha >= 0.0f) {
			//quadrant IV5
			if (-ONE_BY_SQRT3 * beta > alpha) {
				sector = 5;
			} else {
				sector = 6;
			}
		} else {
			//quadrant III
			if (ONE_BY_SQRT3 * beta > alpha) {
				sector = 4;
			} else {
				sector = 5;
			}
		}
	}

	// PWM timings
	uint32_t tA, tB, tC;

	switch (sector) {

	// sector 1-2
	case 1: {
		// Vector on-times
		uint32_t t1 = (alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t2 = (TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tA = (PWMFullDutyCycle + t1 + t2) / 2;
		tB = tA - t1;
		tC = tB - t2;

		break;
	}

	// sector 2-3
	case 2: {
		// Vector on-times
		uint32_t t2 = (alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t3 = (-alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tB = (PWMFullDutyCycle + t2 + t3) / 2;
		tA = tB - t3;
		tC = tA - t2;

		break;
	}

	// sector 3-4
	case 3: {
		// Vector on-times
		uint32_t t3 = (TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t4 = (-alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tB = (PWMFullDutyCycle + t3 + t4) / 2;
		tC = tB - t3;
		tA = tC - t4;

		break;
	}

	// sector 4-5
	case 4: {
		// Vector on-times
		uint32_t t4 = (-alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t5 = (-TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tC = (PWMFullDutyCycle + t4 + t5) / 2;
		tB = tC - t5;
		tA = tB - t4;

		break;
	}

	// sector 5-6
	case 5: {
		// Vector on-times
		uint32_t t5 = (-alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t6 = (alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tC = (PWMFullDutyCycle + t5 + t6) / 2;
		tA = tC - t5;
		tB = tA - t6;

		break;
	}

	// sector 6-1
	case 6: {
		// Vector on-times
		uint32_t t6 = (-TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t1 = (alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tA = (PWMFullDutyCycle + t6 + t1) / 2;
		tC = tA - t1;
		tB = tC - t6;

		break;
	}
	}

	*tAout = tA;
	*tBout = tB;
	*tCout = tC;
	*svm_sector = sector;
}

void foc_run_pid_control_pos(bool index_found, float dt, motor_all_state_t *motor) {
	mc_configuration *conf_now = motor->m_conf;

	float angle_now = motor->m_pos_pid_now;
	float angle_set = motor->m_pos_pid_set;

	float p_term;
	float d_term;
	float d_term_proc;

	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_POS) {
		motor->m_pos_i_term = 0;
		motor->m_pos_prev_error = 0;
		motor->m_pos_prev_proc = angle_now;
		motor->m_pos_d_filter = 0.0;
		motor->m_pos_d_filter_proc = 0.0;
		return;
	}

	// Compute parameters
	float error = utils_angle_difference(angle_set, angle_now);
	float error_sign = 1.0;

	if (conf_now->m_sensor_port_mode != SENSOR_PORT_MODE_HALL) {
		if (conf_now->foc_encoder_inverted) {
			error_sign = -1.0;
		}
	}

	error *= error_sign;

	float kp = conf_now->p_pid_kp;
	float ki = conf_now->p_pid_ki;
	float kd = conf_now->p_pid_kd;
	float kd_proc = conf_now->p_pid_kd_proc;

	if (conf_now->p_pid_gain_dec_angle > 0.1) {
		float min_error = conf_now->p_pid_gain_dec_angle / conf_now->p_pid_ang_div;
		float error_abs = fabs(error);

		if (error_abs < min_error) {
			float scale = error_abs / min_error;
			kp *= scale;
			ki *= scale;
			kd *= scale;
			kd_proc *= scale;
		}
	}

	p_term = error * kp;
	motor->m_pos_i_term += error * (ki * dt);

	// Average DT for the D term when the error does not change. This likely
	// happens at low speed when the position resolution is low and several
	// control iterations run without position updates.
	// TODO: Are there problems with this approach?
	motor->m_pos_dt_int += dt;
	if (error == motor->m_pos_prev_error) {
		d_term = 0.0;
	} else {
		d_term = (error - motor->m_pos_prev_error) * (kd / motor->m_pos_dt_int);
		motor->m_pos_dt_int = 0.0;
	}

	// Filter D
	UTILS_LP_FAST(motor->m_pos_d_filter, d_term, conf_now->p_pid_kd_filter);
	d_term = motor->m_pos_d_filter;

	// Process D term
	motor->m_pos_dt_int_proc += dt;
	if (angle_now == motor->m_pos_prev_proc) {
		d_term_proc = 0.0;
	} else {
		d_term_proc = -utils_angle_difference(angle_now, motor->m_pos_prev_proc) * error_sign * (kd_proc / motor->m_pos_dt_int_proc);
		motor->m_pos_dt_int_proc = 0.0;
	}

	// Filter D process
	UTILS_LP_FAST(motor->m_pos_d_filter_proc, d_term_proc, conf_now->p_pid_kd_filter);
	d_term_proc = motor->m_pos_d_filter_proc;

	// I-term wind-up protection
	float p_tmp = p_term;
	utils_truncate_number_abs(&p_tmp, 1.0);
	utils_truncate_number_abs((float*)&motor->m_pos_i_term, 1.0 - fabsf(p_tmp));

	// Store previous error
	motor->m_pos_prev_error = error;
	motor->m_pos_prev_proc = angle_now;

	// Calculate output
	float output = p_term + motor->m_pos_i_term + d_term + d_term_proc;
	utils_truncate_number(&output, -1.0, 1.0);

	if (conf_now->m_sensor_port_mode != SENSOR_PORT_MODE_HALL) {
		if (index_found) {
			motor->m_iq_set = output * conf_now->l_current_max * conf_now->l_current_max_scale;;
		} else {
			// Rotate the motor with 40 % power until the encoder index is found.
			motor->m_iq_set = 0.4 * conf_now->l_current_max * conf_now->l_current_max_scale;;
		}
	} else {
		motor->m_iq_set = output * conf_now->l_current_max * conf_now->l_current_max_scale;;
	}
}

void foc_run_pid_control_speed(bool index_found, float dt, motor_all_state_t *motor) {
	mc_configuration *conf_now = motor->m_conf;
	float p_term;
	float d_term;

	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_SPEED) {
		motor->m_speed_i_term = 0.0;
		motor->m_speed_prev_error = 0.0;
		motor->m_speed_d_filter = 0.0;
		return;
	}

	if (conf_now->s_pid_ramp_erpms_s > 0.0) {
		utils_step_towards((float*)&motor->m_speed_pid_set_rpm, motor->m_speed_command_rpm, conf_now->s_pid_ramp_erpms_s * dt);
		if (!index_found) {
			utils_truncate_number_abs(&motor->m_speed_pid_set_rpm, conf_now->foc_openloop_rpm);
		}
		utils_truncate_number(&motor->m_speed_pid_set_rpm, conf_now->l_min_erpm, conf_now->l_max_erpm);
	}

	float rpm = 0.0;
	switch (conf_now->s_pid_speed_source) {
	case S_PID_SPEED_SRC_PLL:
		rpm = RADPS2RPM_f(motor->m_pll_speed);
		break;
	case S_PID_SPEED_SRC_FAST:
		rpm = RADPS2RPM_f(motor->m_speed_est_fast);
		break;
	case S_PID_SPEED_SRC_FASTER:
		rpm = RADPS2RPM_f(motor->m_speed_est_faster);
		break;
	}

	float error = motor->m_speed_pid_set_rpm - rpm;

	// Too low RPM set. Reset state, release motor and return.
	if (fabsf(motor->m_speed_pid_set_rpm) < conf_now->s_pid_min_erpm) {
		motor->m_speed_i_term = 0.0;
		motor->m_speed_prev_error = error;
		motor->m_iq_set = 0.0;
		return;
	}

	// Compute parameters
	p_term = error * conf_now->s_pid_kp * (1.0 / 20.0);
	d_term = (error - motor->m_speed_prev_error) * (conf_now->s_pid_kd / dt) * (1.0 / 20.0);

	// Filter D
	UTILS_LP_FAST(motor->m_speed_d_filter, d_term, conf_now->s_pid_kd_filter);
	d_term = motor->m_speed_d_filter;

	// Store previous error
	motor->m_speed_prev_error = error;

	// Calculate output
	float output = p_term + motor->m_speed_i_term + d_term;
	utils_truncate_number_abs(&output, 1.0);

	// Integrator windup protection
	motor->m_speed_i_term += error * conf_now->s_pid_ki * dt * (1.0 / 20.0);
	utils_truncate_number_abs(&motor->m_speed_i_term, 1.0);

	if (conf_now->s_pid_ki < 1e-9) {
		motor->m_speed_i_term = 0.0;
	}

	// Optionally disable braking
	if (!conf_now->s_pid_allow_braking) {
		if (rpm > 20.0 && output < 0.0) {
			output = 0.0;
		}

		if (rpm < -20.0 && output > 0.0) {
			output = 0.0;
		}
	}

	motor->m_iq_set = output * conf_now->lo_current_max * conf_now->l_current_max_scale;
}

float foc_correct_encoder(float obs_angle, float enc_angle, float speed,
							 float sl_erpm, motor_all_state_t *motor) {
	float rpm_abs = fabsf(RADPS2RPM_f(speed));

	// Hysteresis 5 % of total speed
	float hyst = sl_erpm * 0.05;
	if (motor->m_using_encoder) {
		if (rpm_abs > (sl_erpm + hyst)) {
			motor->m_using_encoder = false;
		}
	} else {
		if (rpm_abs < (sl_erpm- hyst)) {
			motor->m_using_encoder = true;
		}
	}

	return motor->m_using_encoder ? enc_angle : obs_angle;
}

float foc_correct_hall(float angle, float dt, motor_all_state_t *motor, int hall_val) {
	mc_configuration *conf_now = motor->m_conf;
	motor->m_hall_dt_diff_now += dt;

	float rpm_abs = fabsf(RADPS2RPM_f(motor->m_pll_speed));
	float rad_per_sec_hall = (M_PI / 3.0) / motor->m_hall_dt_diff_last;
	float rpm_abs_hall = fabsf(RADPS2RPM_f(rad_per_sec_hall));

	motor->m_using_hall = rpm_abs < conf_now->foc_sl_erpm;
	float angle_old = angle;

	int ang_hall_int = conf_now->foc_hall_table[hall_val];

	// Only override the observer if the hall sensor value is valid.
	if (ang_hall_int < 201) {
		// Scale to the circle and convert to radians
		float ang_hall_now = ((float)ang_hall_int / 200.0) * 2.0 * M_PI;

		if (motor->m_ang_hall_int_prev < 0) {
			// Previous angle not valid
			motor->m_ang_hall_int_prev = ang_hall_int;
			motor->m_ang_hall = ang_hall_now;
		} else if (ang_hall_int != motor->m_ang_hall_int_prev) {
			int diff = ang_hall_int - motor->m_ang_hall_int_prev;
			if (diff > 100) {
				diff -= 200;
			} else if (diff < -100) {
				diff += 200;
			}

			// This is only valid if the direction did not just change. If it did, we use the
			// last speed together with the sign right now.
			if (SIGN(diff) == SIGN(motor->m_hall_dt_diff_last)) {
				if (diff > 0) {
					motor->m_hall_dt_diff_last = motor->m_hall_dt_diff_now;
				} else {
					motor->m_hall_dt_diff_last = -motor->m_hall_dt_diff_now;
				}
			} else {
				motor->m_hall_dt_diff_last = -motor->m_hall_dt_diff_last;
			}

			motor->m_hall_dt_diff_now = 0.0;

			// A transition was just made. The angle is in the middle of the new and old angle.
			int ang_avg = motor->m_ang_hall_int_prev + diff / 2;
			ang_avg %= 200;

			// Scale to the circle and convert to radians
			motor->m_ang_hall = ((float)ang_avg / 200.0) * 2.0 * M_PI;
		}

		motor->m_ang_hall_int_prev = ang_hall_int;

		if (RADPS2RPM_f((M_PI / 3.0) / fmaxf(fabsf(motor->m_hall_dt_diff_now),
				fabsf(motor->m_hall_dt_diff_last))) < conf_now->foc_hall_interp_erpm) {
			// Don't interpolate on very low speed, just use the closest hall sensor. The reason is that we might
			// get stuck at 60 degrees off if a direction change happens between two steps.
			motor->m_ang_hall = ang_hall_now;
		} else {
			// Interpolate
			float diff = utils_angle_difference_rad(motor->m_ang_hall, ang_hall_now);
			if (fabsf(diff) < ((2.0 * M_PI) / 12.0) || SIGN(diff) != SIGN(rad_per_sec_hall)) {
				// Do interpolation
				motor->m_ang_hall += rad_per_sec_hall * dt;
			} else {
				// We are too far away with the interpolation
				motor->m_ang_hall -= diff * 0.01;
			}
		}

		utils_norm_angle_rad((float*)&motor->m_ang_hall);

		// Limit hall sensor rate of change. This will reduce current spikes in the current controllers when the angle estimation
		// changes fast.
		float angle_step = (fmaxf(rpm_abs_hall, conf_now->foc_hall_interp_erpm) / 60.0) * 2.0 * M_PI * dt * 1.5;
		float angle_diff = utils_angle_difference_rad(motor->m_ang_hall, motor->m_ang_hall_rate_limited);
		if (fabsf(angle_diff) < angle_step) {
			motor->m_ang_hall_rate_limited = motor->m_ang_hall;
		} else {
			motor->m_ang_hall_rate_limited += angle_step * SIGN(angle_diff);
		}

		utils_norm_angle_rad((float*)&motor->m_ang_hall_rate_limited);

		if (motor->m_using_hall) {
			angle = motor->m_ang_hall_rate_limited;
		}
	} else {
		// Invalid hall reading. Don't update angle.
		motor->m_ang_hall_int_prev = -1;

		// Also allow open loop in order to behave like normal sensorless
		// operation. Then the motor works even if the hall sensor cable
		// gets disconnected (when the sensor spacing is 120 degrees).
		if (motor->m_phase_observer_override && motor->m_state == MC_STATE_RUNNING) {
			angle = motor->m_phase_now_observer_override;
		}
	}

	// Map output angle between hall angle and observer angle in transition region to make
	// a smooth transition.
	if (angle_old != angle) {
		float weight_hall = utils_map(rpm_abs, conf_now->foc_sl_erpm_start, conf_now->foc_sl_erpm, 1.0, 0.0);
		utils_truncate_number(&weight_hall, 0.0, 1.0);
		angle = utils_interpolate_angles_rad(angle, angle_old, weight_hall);
	}

	return angle;
}

void foc_run_fw(motor_all_state_t *motor, float dt) {
	if (motor->m_conf->foc_fw_current_max < fmaxf(motor->m_conf->cc_min_current, 0.001)) {
		return;
	}

	// Field Weakening
	// FW is used in the current and speed control modes. If a different mode is used
	// this code also runs if field weakening was active before. This allows
	// changing control mode even while in field weakening.
	if (motor->m_state == MC_STATE_RUNNING &&
			(motor->m_control_mode == CONTROL_MODE_CURRENT ||
					motor->m_control_mode == CONTROL_MODE_CURRENT_BRAKE ||
					motor->m_control_mode == CONTROL_MODE_SPEED ||
					motor->m_i_fw_set > motor->m_conf->cc_min_current)) {
		float fw_current_now = 0.0;
		float duty_abs = motor->m_duty_abs_filtered;

		if (motor->m_conf->foc_fw_duty_start < 0.99 &&
			duty_abs > motor->m_conf->foc_fw_duty_start * motor->m_conf->l_max_duty
		) {
			//Modified Saluqi begin

			#define TABLEROWS 19 // number of rows in table minus 1
			#define TABLECOLUMNS 19 // number of columns in table minus 1
			const float MAXRPM = 15000.0; // Maximum motor speed boundary in the table
			const int RPMRATED = 500; // For motor speeds larger than this speed field weakening is activated
			const int REFVOLTAGE = 50; // DC bus reference voltage used to generate the iq and Id table

			//#define MAXRPM 15000 // Maximum motor speed boundary in the table
			//#define RPMRATED 500 // For motor speeds larger than this speed field weakening is activated
			//#define REFVOLTAGE 50 // DC bus reference voltage used to generate the iq and Id table


			//float idTable[TABLEROWS+1][TABLECOLUMNS+1] =
			//    	   {{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,-0.0089,-0.1724,-0.2945,-0.3806,-0.4482,-0.5032,-0.5489,-0.5874,-0.6204,-0.6490,-0.6739,-0.6960,-0.7155,-0.7330},
			//			{0.00000,-0.0018,-0.0038,-0.0061,-0.0087,-0.0125,-0.0294,-0.1805,-0.3027,-0.3946,-0.4609,-0.5220,-0.5584,-0.6006,-0.6340,-0.6655,-0.7020,-0.7034,-0.7287,-0.7711},
			//			{-0.0004,-0.0043,-0.0074,-0.0124,-0.0172,-0.0240,-0.0387,-0.1866,-0.3091,-0.4018,-0.4689,-0.5291,-0.5685,-0.6092,-0.6487,-0.6769,-0.7117,-0.7206,-0.7387,-0.7813},
			//			{-0.0028,-0.0067,-0.0076,-0.0137,-0.0180,-0.0234,-0.0413,-0.1941,-0.3162,-0.4056,-0.4754,-0.5323,-0.5813,-0.6169,-0.6638,-0.6870,-0.7152,-0.7455,-0.7505,-0.7845},
			//			{-0.0087,-0.0095,-0.0132,-0.0164,-0.0208,-0.0261,-0.0466,-0.2027,-0.3271,-0.4187,-0.4887,-0.5519,-0.6030,-0.6364,-0.6747,-0.7126,-0.7316,-0.7589,-0.7916,-0.8051},
			//			{-0.0149,-0.0140,-0.0168,-0.0204,-0.0258,-0.0324,-0.0522,-0.2130,-0.3405,-0.4344,-0.5059,-0.5698,-0.6231,-0.6602,-0.6949,-0.7412,-0.7653,-0.7914,-0.8280,-0.8405},
			//			{-0.0196,-0.0208,-0.0203,-0.0249,-0.0308,-0.0386,-0.0601,-0.2271,-0.3569,-0.4536,-0.5285,-0.5895,-0.6456,-0.6916,-0.7295,-0.7759,-0.8133,-0.8458,-0.8777,-0.9107},
			//			{-0.0230,-0.0286,-0.0281,-0.0332,-0.0345,-0.0398,-0.0701,-0.2459,-0.3770,-0.4774,-0.5565,-0.6169,-0.6765,-0.7349,-0.7818,-0.8237,-0.8680,-0.9191,-0.9519,-0.9579},
			//			{-0.0297,-0.0330,-0.0379,-0.0391,-0.0438,-0.0500,-0.0816,-0.2652,-0.4000,-0.5043,-0.5855,-0.6560,-0.7200,-0.7813,-0.8542,-0.9146,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0395,-0.0384,-0.0448,-0.0463,-0.0522,-0.0594,-0.0932,-0.2873,-0.4266,-0.5354,-0.6236,-0.7074,-0.7855,-0.8610,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0500,-0.0500,-0.0500,-0.0569,-0.0583,-0.0664,-0.1045,-0.3136,-0.4571,-0.5722,-0.6743,-0.7736,-0.8759,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0632,-0.0632,-0.0612,-0.0640,-0.0718,-0.0737,-0.1239,-0.3448,-0.4920,-0.6181,-0.7403,-0.8605,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0710,-0.0710,-0.0760,-0.0741,-0.0811,-0.0845,-0.1435,-0.3779,-0.5376,-0.6799,-0.8286,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0777,-0.0777,-0.0886,-0.0896,-0.0888,-0.0993,-0.1627,-0.4136,-0.5970,-0.7629,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.0934,-0.0934,-0.0963,-0.1050,-0.1063,-0.1081,-0.1867,-0.4603,-0.6888,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.1106,-0.1106,-0.1090,-0.1160,-0.1187,-0.1236,-0.2313,-0.5179,-0.7311,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.1295,-0.1295,-0.1291,-0.1278,-0.1310,-0.1403,-0.2931,-0.5919,-0.7311,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.1506,-0.1506,-0.1493,-0.1492,-0.1495,-0.1509,-0.3519,-0.6329,-0.7311,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.1620,-0.1620,-0.1739,-0.1739,-0.1736,-0.1723,-0.4212,-0.6329,-0.7311,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579},
			//			{-0.1867,-0.1867,-0.1867,-0.1867,-0.1867,-0.1867,-0.4429,-0.6329,-0.7311,-0.7898,-0.8358,-0.8607,-0.8890,-0.9062,-0.9190,-0.9302,-0.9394,-0.9455,-0.9519,-0.9579}};

			float idTable[TABLEROWS+1][TABLECOLUMNS+1] =
				{
					{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,-0.0664,-0.1864,-0.2675,-0.3320,-0.3853,-0.4308,-0.4701,-0.5045,-0.5347,-0.5616,-0.5857},
					{-0.0005,-0.0017,-0.0029,-0.0051,-0.0082,-0.0112,-0.0150,-0.0200,-0.0246,-0.0700,-0.1908,-0.2754,-0.3418,-0.3925,-0.4418,-0.4876,-0.5133,-0.5468,-0.5738,-0.6050},
					{-0.0021,-0.0042,-0.0065,-0.0115,-0.0177,-0.0241,-0.0314,-0.0411,-0.0494,-0.0721,-0.1955,-0.2812,-0.3482,-0.4004,-0.4491,-0.4965,-0.5228,-0.5575,-0.5839,-0.6134},
					{-0.0048,-0.0064,-0.0074,-0.0137,-0.0181,-0.0246,-0.0298,-0.0384,-0.0464,-0.0823,-0.2058,-0.2871,-0.3524,-0.4114,-0.4543,-0.5000,-0.5345,-0.5676,-0.5933,-0.6139},
					{-0.0085,-0.0091,-0.0119,-0.0167,-0.0215,-0.0278,-0.0333,-0.0400,-0.0478,-0.0916,-0.2142,-0.2970,-0.3635,-0.4219,-0.4708,-0.5142,-0.5474,-0.5845,-0.6095,-0.6322},
					{-0.0133,-0.0130,-0.0181,-0.0204,-0.0253,-0.0309,-0.0374,-0.0440,-0.0516,-0.1045,-0.2243,-0.3107,-0.3781,-0.4365,-0.4920,-0.5336,-0.5662,-0.6072,-0.6327,-0.6609},
					{-0.0192,-0.0191,-0.0245,-0.0246,-0.0296,-0.0340,-0.0418,-0.0501,-0.0575,-0.1206,-0.2376,-0.3282,-0.3947,-0.4580,-0.5139,-0.5555,-0.5947,-0.6361,-0.6643,-0.6952},
					{-0.0262,-0.0262,-0.0288,-0.0308,-0.0351,-0.0399,-0.0492,-0.0539,-0.0629,-0.1353,-0.2564,-0.3463,-0.4193,-0.4816,-0.5382,-0.5840,-0.6289,-0.6715,-0.7005,-0.7389},
					{-0.0346,-0.0344,-0.0335,-0.0399,-0.0419,-0.0471,-0.0567,-0.0596,-0.0687,-0.1529,-0.2772,-0.3681,-0.4474,-0.5108,-0.5699,-0.6220,-0.6710,-0.7194,-0.7567,-0.8045},
					{-0.0446,-0.0437,-0.0436,-0.0505,-0.0513,-0.0591,-0.0633,-0.0674,-0.0750,-0.1761,-0.2981,-0.3950,-0.4759,-0.5468,-0.6126,-0.6725,-0.7229,-0.7854,-0.8418,-0.9089},
					{-0.0502,-0.0545,-0.0543,-0.0578,-0.0635,-0.0668,-0.0729,-0.0753,-0.0836,-0.2011,-0.3294,-0.4262,-0.5124,-0.5899,-0.6627,-0.7347,-0.8054,-0.8871,-0.9107,-0.9147},
					{-0.0560,-0.0658,-0.0663,-0.0655,-0.0741,-0.0746,-0.0833,-0.0839,-0.0929,-0.2284,-0.3635,-0.4638,-0.5575,-0.6445,-0.7285,-0.8196,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.0688,-0.0765,-0.0805,-0.0791,-0.0815,-0.0901,-0.0916,-0.1006,-0.1083,-0.2579,-0.3996,-0.5078,-0.6117,-0.7169,-0.8298,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.0828,-0.0866,-0.0894,-0.0939,-0.0924,-0.0999,-0.1046,-0.1112,-0.1278,-0.2928,-0.4410,-0.5609,-0.6871,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.0979,-0.0975,-0.0982,-0.1100,-0.1100,-0.1091,-0.1209,-0.1209,-0.1531,-0.3344,-0.4920,-0.6317,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.1144,-0.1141,-0.1266,-0.1293,-0.1293,-0.1268,-0.1309,-0.1309,-0.1866,-0.3803,-0.5551,-0.7265,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.1326,-0.1326,-0.1413,-0.1413,-0.1413,-0.1465,-0.1446,-0.1447,-0.2194,-0.4365,-0.6365,-0.7265,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.1526,-0.1526,-0.1521,-0.1521,-0.1521,-0.1672,-0.1677,-0.1677,-0.2533,-0.5080,-0.6574,-0.7265,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.1739,-0.1739,-0.1728,-0.1728,-0.1765,-0.1929,-0.1930,-0.1931,-0.2923,-0.5536,-0.6574,-0.7265,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147},
					{-0.2013,-0.2013,-0.2013,-0.2013,-0.2013,-0.2013,-0.2013,-0.2015,-0.3594,-0.5536,-0.6574,-0.7265,-0.7753,-0.8123,-0.8425,-0.8661,-0.8789,-0.8984,-0.9107,-0.9147}
				};
				
			//Modified Saluqi end

			//			fw_current_now = utils_map(duty_abs,
			//					motor->m_conf->foc_fw_duty_start * motor->m_conf->l_max_duty,
			//					motor->m_conf->l_max_duty,
			//					0.0, motor->m_conf->foc_fw_current_max);

			// m_current_off_delay is used to not stop the modulation too soon after leaving FW. If axis decoupling
			// is not working properly an oscillation can occur on the modulation when changing the current
			// fast, which can make the estimated duty cycle drop below the FW threshold long enough to stop
			// modulation. When that happens the body diodes in the MOSFETs can see a lot of current and unexpected
			// braking happens. Therefore the modulation is left on for some time after leaving FW to give the
			// oscillation a chance to decay while the MOSFETs are still driven.

			motor->m_current_off_delay = 1.0;

			float rpm = fabsf(RADPS2RPM_f(motor->m_pll_speed));

			if (rpm > MAXRPM) {
				rpm = MAXRPM;
			}
			

			if (rpm > RPMRATED) {
			
				float currentNow = motor->m_motor_state.iq;
				float maxCurrent = motor->m_conf->foc_fw_current_max;

				if (currentNow > maxCurrent){
					currentNow = maxCurrent;
				}

				float currentStep = maxCurrent / (float)(TABLEROWS);

				float currentY = currentNow / currentStep;
				int currentY1 = floor(currentY);
				int currentY2 = ceil(currentY);
				
				// Check the indexes of currentY1 and currenY2.
				check_id_index(&currentY1, &currentY2);

				float rpmStep = MAXRPM / (float)(TABLECOLUMNS);

				float rpmX = rpm / rpmStep;
				int rpmX1 = floor(rpmX);
				int rpmX2 = ceil(rpmX);

				// Check the indexes of rpmX1 and rpmX2.
				check_id_index(&rpmX1, &rpmX2);

				float tableIdQ11 = idTable[currentY1][rpmX1];
				float tableIdQ21 = idTable[currentY1][rpmX2];
				float tableIdQ12 = idTable[currentY2][rpmX1];
				float tableIdQ22 = idTable[currentY2][rpmX2];

				float id_result = ((rpmX2 - rpmX) * (currentY2 - currentY)) / ((rpmX2 - rpmX1) * (currentY2 - currentY1)) * tableIdQ11 +
								  ((rpmX - rpmX1) * (currentY2 - currentY)) / ((rpmX2 - rpmX1) * (currentY2 - currentY1)) * tableIdQ21 +
								  ((rpmX2 - rpmX) * (currentY - currentY1)) / ((rpmX2 - rpmX1) * (currentY2 - currentY1)) * tableIdQ12 +
								  ((rpmX - rpmX1) * (currentY - currentY1)) / ((rpmX2 - rpmX1) * (currentY2 - currentY1)) * tableIdQ22;       

				fw_current_now = id_result * motor->m_conf->foc_fw_current_max;

				if (motor->m_conf->foc_fw_ramp_time < dt) {
					motor->m_i_fw_set = fw_current_now;
				} else {
					utils_step_towards((float*)&motor->m_i_fw_set, fw_current_now,
					(dt / motor->m_conf->foc_fw_ramp_time));
				}
			}
		}
	}
}

/*
Checks the indeces are not equal and if the index is in range of the interpolation table.
*/
void check_id_index(int *index1, int *index2) {
    if (*index1 == *index2) {
        if (*index2 < TABLEROWS){
            ++ *index2;
        }
        else {
            -- *index1;

            if (*index1 < 0) {
                *index1 = 0;
                *index2 = 1;
            }
        }
    }
}

void foc_hfi_adjust_angle(float ang_err, motor_all_state_t *motor, float dt) {
	mc_configuration *conf = motor->m_conf;
	utils_truncate_number_abs(&ang_err, conf->foc_hfi_max_err);

	// TODO: Check if ratio between these is sane or introduce separate gains
	const float gain_int = 4000.0 * conf->foc_hfi_gain;
	const float gain_int2 = 10.0 * conf->foc_hfi_gain;
	motor->m_hfi.double_integrator += ang_err * gain_int2;
	utils_truncate_number_abs(&motor->m_hfi.double_integrator, fabsf(motor->m_speed_est_fast));
	motor->m_hfi.angle -= dt * (gain_int * ang_err + motor->m_hfi.double_integrator);
	utils_norm_angle_rad((float*)&motor->m_hfi.angle);
	motor->m_hfi.ready = true;
}

void foc_precalc_values(motor_all_state_t *motor) {
	const mc_configuration *conf_now = motor->m_conf;
	motor->p_lq = conf_now->foc_motor_l + conf_now->foc_motor_ld_lq_diff * 0.5;
	motor->p_ld = conf_now->foc_motor_l - conf_now->foc_motor_ld_lq_diff * 0.5;
	motor->p_inv_ld_lq = (1.0 / motor->p_lq - 1.0 / motor->p_ld);
	motor->p_v2_v3_inv_avg_half = (0.5 / motor->p_lq + 0.5 / motor->p_ld) * 0.9; // With the 0.9 we undo the adjustment from the detection
	motor->m_observer_state.lambda_est = conf_now->foc_motor_flux_linkage;
}
