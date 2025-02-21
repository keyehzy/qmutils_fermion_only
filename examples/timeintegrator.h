#pragma once

class TimeIntegrator {
 public:
  TimeIntegrator(float initial_time, float final_time, size_t num_time_steps,
                 const arma::sp_cx_fmat& h_matrix)
      : initial_time_(initial_time),
        final_time_(final_time),
        num_time_steps_(num_time_steps),
        counter_(0),
        hamiltonian_matrix_(h_matrix) {
    time_step_ = (final_time_ - initial_time_) / num_time_steps_;

    identity_matrix_ = arma::speye<arma::sp_cx_fmat>(
        hamiltonian_matrix_.n_rows, hamiltonian_matrix_.n_cols);

    auto scaled_h =
        std::complex<float>(0, time_step_ * 0.5f) * hamiltonian_matrix_;

    crank_minus_ = identity_matrix_ - scaled_h;
    crank_plus_ = identity_matrix_ + scaled_h;

    if (!factoriser_.factorise(crank_minus_)) {
      std::cerr << "Matrix factorisation failed" << std::endl;
      std::abort();
    }
  }

  bool step(arma::cx_fvec& state_vector) {
    if (counter_ > num_time_steps_) {
      return false;
    }

    auto right_hand_side = crank_plus_ * state_vector;
    arma::cx_fvec result;

    if (!factoriser_.solve(result, right_hand_side)) {
      std::cerr << "Failed to solve linear system" << std::endl;
      std::abort();
    }

    state_vector = result;
    ++counter_;
    return true;
  }

  float time() const { return initial_time_ + counter_ * time_step_; }

 private:
  const float initial_time_;
  const float final_time_;
  const size_t num_time_steps_;
  float time_step_;
  size_t counter_;
  arma::sp_cx_fmat hamiltonian_matrix_;
  arma::sp_cx_fmat identity_matrix_;
  arma::sp_cx_fmat crank_minus_;
  arma::sp_cx_fmat crank_plus_;
  arma::spsolve_factoriser factoriser_;
};
