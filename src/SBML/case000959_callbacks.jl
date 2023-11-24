function get_callbacks_case000959(foo)

	function condition___parameter_ifelse2(u, t, integrator)
		t == 1
	end

	function affect___parameter_ifelse2!(integrator)
		integrator.p[4] = 1.0
	end

	function is_active_t0___parameter_ifelse2!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[4] = 0.0 # Default to being off
		if (t >= 1)
			p[4] = 1.0
		end
	end


	cb___parameter_ifelse2 = DiscreteCallback(condition___parameter_ifelse2, affect___parameter_ifelse2!, save_positions=(false, false))


	function condition___parameter_ifelse3(u, t, integrator)
		t == 0
	end

	function affect___parameter_ifelse3!(integrator)
		integrator.p[2] = 1.0
	end

	function is_active_t0___parameter_ifelse3!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[2] = 0.0 # Default to being off
		if !(t <= 0)
			p[2] = 1.0
		end
	end


	cb___parameter_ifelse3 = DiscreteCallback(condition___parameter_ifelse3, affect___parameter_ifelse3!, save_positions=(false, false))


	function condition___parameter_ifelse5(u, t, integrator)
		t == 1
	end

	function affect___parameter_ifelse5!(integrator)
		integrator.p[3] = 1.0
	end

	function is_active_t0___parameter_ifelse5!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[3] = 0.0 # Default to being off
		if !(t < 1)
			p[3] = 1.0
		end
	end


	cb___parameter_ifelse5 = DiscreteCallback(condition___parameter_ifelse5, affect___parameter_ifelse5!, save_positions=(false, false))


	function condition___parameter_ifelse1(u, t, integrator)
		t == 0
	end

	function affect___parameter_ifelse1!(integrator)
		integrator.p[5] = 1.0
	end

	function is_active_t0___parameter_ifelse1!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[5] = 0.0 # Default to being off
		if !(t <= 0)
			p[5] = 1.0
		end
	end


	cb___parameter_ifelse1 = DiscreteCallback(condition___parameter_ifelse1, affect___parameter_ifelse1!, save_positions=(false, false))


	function condition___parameter_ifelse4(u, t, integrator)
		t == 1
	end

	function affect___parameter_ifelse4!(integrator)
		integrator.p[1] = 1.0
	end

	function is_active_t0___parameter_ifelse4!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[1] = 0.0 # Default to being off
		if (t >= 1)
			p[1] = 1.0
		end
	end


	cb___parameter_ifelse4 = DiscreteCallback(condition___parameter_ifelse4, affect___parameter_ifelse4!, save_positions=(false, false))


	function condition___parameter_ifelse6(u, t, integrator)
		t == 1
	end

	function affect___parameter_ifelse6!(integrator)
		integrator.p[6] = 1.0
	end

	function is_active_t0___parameter_ifelse6!(u, p)
		t = 0.0 # Used to check conditions activated at t0=0
		p[6] = 0.0 # Default to being off
		if !(t < 1)
			p[6] = 1.0
		end
	end


	cb___parameter_ifelse6 = DiscreteCallback(condition___parameter_ifelse6, affect___parameter_ifelse6!, save_positions=(false, false))

	return CallbackSet(cb___parameter_ifelse2, cb___parameter_ifelse3, cb___parameter_ifelse5, cb___parameter_ifelse1, cb___parameter_ifelse4, cb___parameter_ifelse6), Function[is_active_t0___parameter_ifelse2!, is_active_t0___parameter_ifelse3!, is_active_t0___parameter_ifelse5!, is_active_t0___parameter_ifelse1!, is_active_t0___parameter_ifelse4!, is_active_t0___parameter_ifelse6!]
end


function compute_tstops(u::AbstractVector, p::AbstractVector)
	return Float64[1.0, -0.0, 1.0, -0.0, 1.0, 1.0]
end