

k = 14.38 #interaction constant to be fitted, in [eV/elem charge]
q = 0.08 #transition charge values, in [elem charge]
L = 7.74 #edge of the triangle within 10%, to be fitted, in [Angstrom]
ħ = 6.62607004*10^(-34)
c = 299792458
eV = 1.6*10^(-19)

function same_int(latt_const, mers, angle)
#  while angle > 360:
#    angle = angle - 360
#  while angle < 0:
#    angle = angle + 360
  angle = -angle/360*2*π #the minus is because I do the calculations with clockwise instead of anticlockwise angles
  z = latt_const*mers #distance along the chain
  tau = 2*π/3+angle #complementary angle 1
  theta = 2*π/3-angle #complementary angle 2



  V_same = 2*k*q^2/sqrt(2/3*(L^2)*(1-cos(angle))+z^2) - k*q^2/sqrt(2/3*L^2*(1-cos(tau))+z^2) - k*q^2/sqrt(2/3*L^2*(1-cos(theta))+z^2)
  return V_same

end

function cross_int(latt_const, mers, angle)
#  while angle > 360:
#    angle = angle - 360
#  while angle < 0:
#    angle = angle + 360
  angle = -angle/360*2*π #the minus is because I do the calculations with clockwise instead of anticlockwise angles
  z = latt_const*mers #distance along the chain
  tau = 2*π/3+angle #complementary angle 1
  theta = 2*π/3-angle #complementary angle 2



  V_cross = sqrt(3)*k*q^2/sqrt(2/3*L^2*(1-cos(tau))+z^2) - sqrt(3)*k*q^2/sqrt(2/3*L^2*(1-cos(theta))+z^2)
  return V_cross

end

function FC_overlap_0K(e_vibr, HR_fac)

  FCoverlap = 1/sqrt(factorial(e_vibr)) * exp(-HR_fac^2/2) * HR_fac^e_vibr * (-1)^(1*e_vibr)

  return FCoverlap
end

function gaussian(x, μ, σ)
   prob = 1/(sqrt(2*π)*σ)*exp( -((x - μ)^2/σ^2) /2)
   return prob
end

function OhmSD(E, wc)
    E = abs(E)
    SpecDens = E*exp(-(E/wc))
    return SpecDens
end

function Boltz_weight(E, kbT)
    if E > 0
        occup = exp(-E/(kbT))
    else
        occup = 1
    end
  return occup
end


function Bose_weight(E, kbT)
    if E > 0
        occup = 1/(exp((E+0.000001)/kbT)-1)
    else
        occup = 1+ 1/(exp((-E+0.000001)/kbT)-1)
    end
  return occup
end

function subdif_msd(t, alpha, A)
    msd = A*t^alpha

    return msd
end
