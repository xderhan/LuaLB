T = 0.5

rho_l = 2
rho_g = 1

--
--
--

findRoot = function(lhs, x0, x1, tol)
	local f0 = lhs(x0)
	local f1 = lhs(x1)

	assert(f0*f1 < 0)

	while math.abs(x0 - x1) > tol do
		local xc = (x0 + x1)/2
		local fc = lhs(xc)

		if fc*f0 < 0 then
			f1 = fc
			x1 = xc
		elseif fc*f1 < 0 then
			f0 = fc
			x0 = xc
		else
			break
		end
	end

	return 0.5*(x0 + x1)
end

auxLhs = function(b)
	local t1 = 1 - b*rho_l
	local t2 = 1 - b*rho_g

	return b - 2/(rho_l + rho_g)
	     + t1*t2*math.log((rho_l*t2)/(rho_g*t1))/(rho_l - rho_g)
end

--
--
--

b = findRoot(auxLhs, 0, 1/(rho_l + 0.1), 1.0e-8)
a = T/((rho_l + rho_g)*(1 - b*rho_l)*(1 - b*rho_g))

rho_c = 3*b
T_c = 8*a/(27*b)

io.write(string.format("for rho_g = %f, rho_l = %f and T = %f\n" ..
                       "  a = %f\n  b = %f\n  rho_c = %f\n  T_c = %f\n",
		       rho_g, rho_l, T, a, b, rho_c, T_c))

--
--  rho_l, rho_g
--
