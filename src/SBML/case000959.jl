function get_ODESystem_case000959(foo)
	# Model name: case000959

	ModelingToolkit.@variables t foo(t) P12(t) P18(t) P23(t) P16(t) P6(t) P7(t) P3(t) P22(t) P25(t) P14(t) P15(t) P21(t) P4(t) P17(t) P9(t) P19(t) P5(t) P24(t) P26(t) P8(t) P20(t) P2(t) P1(t) P11(t) P10(t) P13(t) 
	species = [foo, P12, P18, P23, P16, P6, P7, P3, P22, P25, P14, P15, P21, P4, P17, P9, P19, P5, P24, P26, P8, P20, P2, P1, P11, P10, P13]

	ModelingToolkit.@parameters __parameter_ifelse4 __parameter_ifelse3 __parameter_ifelse5 __parameter_ifelse2 __parameter_ifelse1 __parameter_ifelse6 
	parameters = [__parameter_ifelse4, __parameter_ifelse3, __parameter_ifelse5, __parameter_ifelse2, __parameter_ifelse1, __parameter_ifelse6]

    D = Differential(t)

	eqs = [
	D(foo) ~ 0.0,
	P12 ~ asec(-t-1),
	P18 ~ asinh(-t),
	P23 ~ acsch(t+1),
	P16 ~ acot(-t-0.001),
	P6 ~ cot(-t-0.001),
	P7 ~ sinh(t),
	P3 ~ csc(t+0.001),
	P22 ~ tanh(10 * (((1 - __parameter_ifelse1)*(1.0) + __parameter_ifelse1*(0.0))+((1 - __parameter_ifelse2)*(0.0) + __parameter_ifelse2*(1.0)))) * (0) + (1 - tanh(10 * (((1 - __parameter_ifelse1)*(1.0) + __parameter_ifelse1*(0.0))+((1 - __parameter_ifelse2)*(0.0) + __parameter_ifelse2*(1.0))))) * (asech(t)),
	P25 ~ acoth(t+1.001),
	P14 ~ acsc(-t-1),
	P15 ~ acot(t+0.001),
	P21 ~ ((1 - __parameter_ifelse5)*(atanh(-t)) + __parameter_ifelse5*(-10)),
	P4 ~ csc(-t-0.001),
	P17 ~ asinh(t),
	P9 ~ cosh(t),
	P19 ~ acosh(t+1),
	P5 ~ cot(t+0.001),
	P24 ~ acsch(-t-1),
	P26 ~ acoth(-t-1.001),
	P8 ~ sinh(-t),
	P20 ~ ((1 - __parameter_ifelse6)*(atanh(t)) + __parameter_ifelse6*(10)),
	P2 ~ sec(-t),
	P1 ~ sec(t),
	P11 ~ asec(t+1),
	P10 ~ cosh(-t),
	P13 ~ acsc(t+1),
	]

	@named sys = ODESystem(eqs, t, species, parameters)

	specie_map = [
	foo =>1.0,
	P12 => asec(-t-1),
	P18 => asinh(-t),
	P23 => acsch(t+1),
	P16 => acot(-t-0.001),
	P6 => cot(-t-0.001),
	P7 => sinh(t),
	P3 => csc(t+0.001),
	P22 => tanh(10 * (ifelse(t <= 0, 1.0, 0.0)+ifelse(t >= 1, 1.0, 0.0))) * (0) + (1 - tanh(10 * (ifelse(t <= 0, 1.0, 0.0)+ifelse(t >= 1, 1.0, 0.0)))) * (asech(t)),
	P25 => acoth(t+1.001),
	P14 => acsc(-t-1),
	P15 => acot(t+0.001),
	P21 => ifelse(t < 1, atanh(-t), -10),
	P4 => csc(-t-0.001),
	P17 => asinh(t),
	P9 => cosh(t),
	P19 => acosh(t+1),
	P5 => cot(t+0.001),
	P24 => acsch(-t-1),
	P26 => acoth(-t-1.001),
	P8 => sinh(-t),
	P20 => ifelse(t < 1, atanh(t), 10),
	P2 => sec(-t),
	P1 => sec(t),
	P11 => asec(t+1),
	P10 => cosh(-t),
	P13 => acsc(t+1),
	]

	# SBML file parameter values
	parameter_map = [
	__parameter_ifelse4 =>0.0,
	__parameter_ifelse3 =>0.0,
	__parameter_ifelse5 =>0.0,
	__parameter_ifelse2 =>0.0,
	__parameter_ifelse1 =>0.0,
	__parameter_ifelse6 =>0.0,
	]

    return sys, specie_map, parameter_map

end
